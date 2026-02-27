"""
Tier 3 Template: Functional State Annotation
=============================================
Code reference for phases/tier3.md rules.
Read tier3.md for rules/criteria FIRST, then use this code.
"""
import re
import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import json
import os
from pathlib import Path

# ── Step 1: Subset from Tier 2 ──────────────────────────────────────
def subset_by_dev_state(tier2_path, dev_state, tier2_col='tier2_annotation'):
    tier2_subset = sc.read_h5ad(tier2_path)
    subset = tier2_subset[tier2_subset.obs[tier2_col] == dev_state].copy()
    print(f"Subsetting → {dev_state}: {subset.n_obs} cells")
    return subset

# ── Tier 3 Input: Annotation + Cluster 결합 ─────────────────────────
def create_tier3_groups(subset, tier2_annot_col='tier2_annotation', tier2_cluster_col='tier2_cluster'):
    """
    Tier 3 입력: Tier 2 annotation + cluster 결합
    ex) Memory_B_1, Memory_B_2, Naive_B_3
    """
    subset.obs['tier3_group'] = (
        subset.obs[tier2_annot_col].astype(str) + '_' +
        subset.obs[tier2_cluster_col].astype(str)
    )
    print(f"Tier 3 groups: {subset.obs['tier3_group'].nunique()}")
    return subset

# ── Step 2: RE-COMPUTE DE (on tier3_group) ───────────────────────────
def compute_de(subset, groupby='tier3_group'):
    """DE on tier3_group (annotation + cluster combined labels)."""
    sc.tl.rank_genes_groups(subset, groupby=groupby, method='wilcoxon',
                            use_raw=True, pts=True)
    de_df = sc.get.rank_genes_groups_df(subset, group=None)
    print(f"⚠️ DE computed on tier3_group ({subset.n_obs} cells)")
    return de_df

# ── Step 3: Pathway Activity (MANDATORY) ─────────────────────────────
def compute_pathway_activity(subset, organism='human', top=500):
    """MANDATORY: Cannot proceed without pathway analysis."""
    print("🔬 Pathway Activity Analysis (MANDATORY)")
    net = dc.op.progeny(organism=organism, top=top)
    dc.run_mlm(mat=subset, net=net, source='source', target='target',
               weight='weight', verbose=True, use_raw=False)
    assert 'mlm_estimate' in subset.obsm, "❌ Pathway activity not computed!"
    print(f"✅ Pathway activity: {subset.obsm['mlm_estimate'].shape[1]} pathways")
    return subset

# ── Step 4.5: Pathway Outlier Detection ──────────────────────────────
def detect_pathway_outliers(subset, cluster_col='tier3_cluster'):
    """Detect pathway activity outliers using z-scores (data-driven)."""
    pw_key = 'score_mlm' if 'score_mlm' in subset.obsm else 'mlm_estimate'
    pval_key = 'padj_mlm' if 'padj_mlm' in subset.obsm else 'mlm_pvals'
    pw_scores = subset.obsm[pw_key]
    pw_pvals = subset.obsm[pval_key]

    outliers = []
    for cluster in subset.obs[cluster_col].unique():
        mask = subset.obs[cluster_col] == cluster
        cluster_pw_mean = pw_scores[mask].mean()
        cluster_pw_pval = pw_pvals[mask].mean()

        for pathway in pw_scores.columns:
            global_mean = pw_scores[pathway].mean()
            global_std = pw_scores[pathway].std()
            cluster_val = cluster_pw_mean[pathway]
            cluster_pval = cluster_pw_pval[pathway]

            if global_std > 1e-6:
                z_score = (cluster_val - global_mean) / global_std
            else:
                z_score = 0.0

            if abs(z_score) > 2.5 and abs(cluster_val) > 0.5 and cluster_pval < 0.05:
                outliers.append({
                    'cluster': str(cluster), 'pathway': pathway,
                    'activity': float(cluster_val), 'z_score': float(z_score),
                    'pval': float(cluster_pval),
                    'direction': 'HIGH' if z_score > 0 else 'LOW'
                })

    subset.uns['pathway_outliers'] = outliers
    print(f"   {'⚠️ ' + str(len(outliers)) + ' pathway outliers' if outliers else '✅ No outliers'}")
    return outliers

# ── Step 7: Filter Valid Markers (reasoning용, pct≥0.25) ─────────────
# (Steps 5.6, 8 are agent reasoning steps — see tier3.md)
def filter_valid_markers(de_df, cluster_id, top_n=50):
    """pct >= 25%, LFC >= 1, padj < 0.05, Top 50 (reasoning 단계용)"""
    cluster_df = de_df[de_df['group'] == cluster_id].copy()
    valid = cluster_df[
        (cluster_df['pct_nz_group'] >= 0.25) &
        (cluster_df['logfoldchanges'] >= 1.0) &
        (cluster_df['pvals_adj'] < 0.05)
    ]
    return valid.nlargest(top_n, 'logfoldchanges')

# ── HK removal + Dotplot-highest (DotPlot pool용) ─────────────────────
_HK_PREFIX_RE = re.compile(r'^(Rps|Rpl|Mrps|Mrpl|mt-|Hmgb|Gm\d)', re.IGNORECASE)
_HK_HISTONE_RE = re.compile(r'^H[1-4][a-zA-Z0-9\-]|^HIST[123]', re.IGNORECASE)
_HK_EXACT = frozenset({
    'Malat1', 'Tmsb4x', 'Tmsb10', 'Actb', 'Actg1',
    'Eef1a1', 'Eef2', 'Tpt1', 'Ppia', 'Fau', 'B2m',
    'Calm1', 'Cmss1', 'Psap', 'Ybx1', 'Ptma', 'Npm1',
    'MALAT1', 'TMSB4X', 'TMSB10', 'ACTB', 'ACTG1',
    'EEF1A1', 'EEF2', 'TPT1', 'PPIA', 'FAU', 'B2M',
    'CALM1', 'CMSS1', 'PSAP', 'YBX1', 'PTMA', 'NPM1',
})

def is_housekeeping(gene: str) -> bool:
    return (gene in _HK_EXACT or
            bool(_HK_PREFIX_RE.match(gene)) or
            bool(_HK_HISTONE_RE.match(gene)))

def remove_housekeeping_genes(de_df: pd.DataFrame) -> pd.DataFrame:
    mask = ~de_df['names'].apply(is_housekeeping)
    removed = (~mask).sum()
    if removed:
        print(f"   HK removed: {removed} genes "
              f"({', '.join(de_df.loc[~mask, 'names'].unique()[:5])}...)")
    return de_df[mask].copy()

def verify_dotplot_highest(de_df: pd.DataFrame, adata,
                           annotation_col: str) -> pd.DataFrame:
    """Keep genes where pct_nz is highest in the target annotation."""
    candidate_genes = de_df['names'].unique().tolist()
    var_names = list(adata.var_names)
    valid_genes = [g for g in candidate_genes if g in var_names]
    gene_idx = [var_names.index(g) for g in valid_genes]

    X_sub = adata.X[:, gene_idx]
    if hasattr(X_sub, 'toarray'):
        X_sub = X_sub.toarray()
    else:
        X_sub = np.asarray(X_sub)

    cell_type_labels = adata.obs[annotation_col].astype(str).values
    unique_cts = list(dict.fromkeys(cell_type_labels))

    pct_matrix = np.zeros((len(unique_cts), len(valid_genes)))
    for i, ct in enumerate(unique_cts):
        mask = cell_type_labels == ct
        pct_matrix[i, :] = (X_sub[mask, :] > 0).mean(axis=0)

    gene_to_highest = {
        gene: unique_cts[pct_matrix[:, j].argmax()]
        for j, gene in enumerate(valid_genes)
    }

    def _is_highest(row):
        return gene_to_highest.get(row['names'], '') == str(row['group'])

    mask = de_df.apply(_is_highest, axis=1)
    print(f"   Dotplot-highest: {mask.sum()} pass, {(~mask).sum()} removed")
    return de_df[mask].copy()

# ── Step 9: Per-Marker PMID Verification + Evidence ─────────────────
def build_per_marker_evidence(markers, cell_type, reasonings=None, species='human'):
    """Verify each marker and return per-marker evidence list for annotation_evidence.

    반환값은 annotation_evidence['markers'] 리스트로 직접 사용됩니다.
    report.md의 Marker Evidence 표에 Gene / pct / LFC / PMID / Title / Reasoning
    컬럼으로 출력됩니다.

    Args:
        markers: List of dicts [{'gene': ..., 'pct_in': ..., 'log2fc': ...}]
                 또는 gene 이름 list
        cell_type: Assigned cell type name (e.g., 'Cytotoxic')
        reasonings: dict {gene: reasoning_str} — agent가 작성한 선정 이유
        species: Species context

    Returns:
        list: [
          {
            'gene':      'Gzmb',
            'pct_in':    0.82,
            'log2fc':    3.1,
            'pmid':      '12345678',
            'title':     '...',
            'reasoning': '...'   # REQUIRED
          }, ...
        ]
    """
    from tools.dynamic_knowledge import pubmed_search
    import time

    reasonings = reasonings or {}

    if markers and isinstance(markers[0], str):
        markers = [{'gene': g, 'pct_in': 0, 'log2fc': 0} for g in markers]

    evidence = []
    for m in markers:
        gene = m.get('gene', m.get('names', ''))
        query = f"{gene} {cell_type} marker {species}"
        results = pubmed_search(query, max_results=3)
        time.sleep(0.15)

        top = results[0] if results else {}
        evidence.append({
            'gene':      gene,
            'pct_in':    m.get('pct_in', m.get('pct_nz_group', 0)),
            'log2fc':    m.get('log2fc', m.get('logfoldchanges', 0)),
            'pmid':      top.get('pmid', ''),
            'title':     top.get('title', ''),
            'reasoning': reasonings.get(gene, ''),
        })

    verified = sum(1 for e in evidence if e['pmid'])
    missing  = [e['gene'] for e in evidence if not e['reasoning']]
    print(f"   PMID verified: {verified}/{len(evidence)} markers")
    if missing:
        print(f"   ⚠️ reasoning missing for: {missing}")
    return evidence

# ── Step 6: Verify ───────────────────────────────────────────────────
def verify_functional_analysis(subset):
    errors = []
    if 'mlm_estimate' not in subset.obsm:
        errors.append("Pathway activity (mlm_estimate) not computed")
    if 'mlm_pvals' not in subset.obsm:
        errors.append("Pathway p-values (mlm_pvals) not computed")
    if errors:
        raise AssertionError("❌ INCOMPLETE:\n" + "\n".join(f"  - {e}" for e in errors))
    print("✅ Functional analysis verified")


# ── Step 10: 3-Iteration Reasoning ──────────────────────────────────
def run_3iteration_reasoning(subset, de_df, cluster_id,
                             group_col='tier3_group',
                             annotation_col='tier3_annotation'):
    """Execute 3-iteration reasoning for a single cluster/group.

    Reference: reasoning/integrated_format.md (Tier 3 section)

    Iteration 1 (Bioinformatician): Collect DE + Pathway evidence
    Iteration 2 (Computational Biologist): Evaluate candidates + literature
    Iteration 3 (PI): Integrated decision + confidence scoring

    Args:
        subset: AnnData with pathway activity computed
        de_df: DE results DataFrame
        cluster_id: Group/cluster to reason about
        group_col: Column with group labels (tier3_group)
        annotation_col: Column to store annotation

    Returns:
        dict: Reasoning result with annotation, confidence, evidence, candidates
    """
    # ── Iteration 1: Evidence Collection ─────────────────────────────
    # 1a. Top markers
    valid_markers = filter_valid_markers(de_df, cluster_id, top_n=50)
    top_markers = valid_markers[['names', 'pct_nz_group',
                                 'logfoldchanges', 'pvals_adj']].to_dict('records')

    # 1b. Pathway activity
    pw_key = 'mlm_estimate' if 'mlm_estimate' in subset.obsm else 'score_mlm'
    pval_key = 'mlm_pvals' if 'mlm_pvals' in subset.obsm else 'padj_mlm'
    mask = subset.obs[group_col] == cluster_id
    pw_mean = subset.obsm[pw_key][mask].mean()
    pw_pval = subset.obsm[pval_key][mask].mean()
    top_pathways = pw_mean.abs().nlargest(10)
    pathway_evidence = [
        {'pathway': pw, 'activity': float(pw_mean[pw]),
         'pval': float(pw_pval[pw]),
         'direction': 'HIGH' if pw_mean[pw] > 0 else 'LOW'}
        for pw in top_pathways.index
    ]

    # 1c. Outliers for this group
    pw_outliers = [o for o in subset.uns.get('pathway_outliers', [])
                   if str(o['cluster']) == str(cluster_id)]

    iteration1 = {
        'markers': top_markers,
        'pathway_activity': pathway_evidence,
        'pathway_outliers': pw_outliers,
        'n_valid_markers': len(valid_markers)
    }

    # ── Iteration 2: Candidate Reasoning ─────────────────────────────
    iteration2 = {
        'candidates': [],  # Agent populates
        'literature_queries': [],
        'instruction': ('Use markers + pathway activity to generate 2-3 candidates. '
                        'Search PubMed for each candidate. '
                        'Format: reasoning/integrated_format.md Tier 3 Iteration 2')
    }

    # ── Iteration 3: Decision Template ───────────────────────────────
    iteration3 = {
        'annotation': None,
        'confidence_score': 0,
        'confidence_level': 'INSUFFICIENT',
        'scoring': {
            'markers': 0,
            'references': 0,
            'pathway_consistency': 0,
            'literature': 0
        },
        'key_evidence': [],
        'references': [],
        'instruction': ('Integrate all evidence into final decision. '
                        'Score each criterion 0-3 pts. '
                        'Format: reasoning/integrated_format.md Tier 3 Iteration 3')
    }

    reasoning = {
        'cluster_id': str(cluster_id),
        'tier': 3,
        'iteration1': iteration1,
        'iteration2': iteration2,
        'iteration3': iteration3
    }

    print(f"   Group {cluster_id}: {len(top_markers)} markers, "
          f"{len(pathway_evidence)} pathways, {len(pw_outliers)} outliers")
    return reasoning

# ── Step 12: Assign Annotations ──────────────────────────────────────
def assign_annotations(subset, group_to_annotation, group_col='tier3_group',
                       annotation_col='tier3_annotation'):
    """Assign tier3 annotations based on group-to-annotation mapping.

    Args:
        subset: AnnData with tier3_group labels
        group_to_annotation: dict mapping group → annotation name
            e.g., {'Memory_B_1': 'Cytotoxic', 'Memory_B_2': 'Exhausted'}
        group_col: Column with group labels
        annotation_col: Column to write annotations to

    Returns:
        AnnData with annotation_col added
    """
    subset.obs[annotation_col] = (
        subset.obs[group_col].astype(str).map(group_to_annotation)
    )
    unmapped = subset.obs[annotation_col].isna().sum()
    if unmapped > 0:
        print(f"⚠️ {unmapped} cells have no annotation mapping")
        subset.obs[annotation_col] = subset.obs[annotation_col].fillna('Uncertain')
    subset.obs[annotation_col] = pd.Categorical(subset.obs[annotation_col])
    print(f"✅ Assigned {subset.obs[annotation_col].nunique()} annotations to '{annotation_col}'")
    return subset


# ── Step 12a: Annotation 기반 DE 재계산 ──────────────────────────────
def compute_de_by_annotation(subset, annotation_col='tier3_annotation'):
    """Annotation 단위 DE 재계산 — DotPlot pool 구성용.

    tier3_group 단위 DE(Step 2)와 별도로, 최종 annotation 기준 DE를 재계산.
    groupby=annotation_col, wilcoxon, use_raw=True, pts=True.

    Args:
        subset: tier3_annotation이 할당된 AnnData
        annotation_col: tier3 annotation 컬럼명

    Returns:
        de_df: annotation 단위 DE DataFrame (group = annotation name)
    """
    sc.tl.rank_genes_groups(subset, groupby=annotation_col, method='wilcoxon',
                            use_raw=True, pts=True)
    de_df = sc.get.rank_genes_groups_df(subset, group=None)
    n_groups = subset.obs[annotation_col].nunique()
    print(f"DE re-computed: {n_groups} annotations (groupby='{annotation_col}')")
    return de_df


# ── Step 12.5: DotPlot 후보 마커 Pool 구성 ───────────────────────────
def build_marker_pool(ann_de_df, adata,
                      annotation_col='tier3_annotation', pool_size=50):
    """Steps 12b-12d: annotation DE에서 DotPlot 후보 마커 pool 구성.

    Step 12b: pct≥0.25, LFC≥1, padj<0.05 strict filter → Top 50 per annotation
              (Tier 3는 dev-state subset → functional marker → pct 0.25 유지)
    Step 12c: HK removal — Rps/Rpl/mt-/Gm*/기타
    Step 12d: Dotplot-highest — argmax pct == target annotation

    ※ Step 12.5는 Agent가 수행:
       pool에서 PubMed 검색 → PMID 확인된 gene만 DotPlot 대상
       PMID 없는 gene은 DotPlot에 포함하지 않음
       3-5개 최종 선정 (PMID + reasoning 필수)

    Args:
        ann_de_df: Annotation-based DE DataFrame (from compute_de_by_annotation)
                   group = annotation name
        adata: AnnData with annotation_col and adata.X (log-normalized)
        annotation_col: Column with annotation labels
        pool_size: Candidate pool size per annotation (default 50)

    Returns:
        dict: {annotation_name: DataFrame (names, pct_nz_group, logfoldchanges, pvals_adj)}
    """
    # Step 12b: strict filter
    filtered = ann_de_df[
        (ann_de_df['pct_nz_group'] >= 0.25) &
        (ann_de_df['logfoldchanges'] >= 1.0) &
        (ann_de_df['pvals_adj'] < 0.05)
    ].copy()
    print(f"   Step 12b: {len(filtered)} genes pass pct≥0.25, LFC≥1, padj<0.05")

    # Step 12c: HK removal
    filtered = remove_housekeeping_genes(filtered)

    # Step 12d: Dotplot-highest
    filtered = verify_dotplot_highest(filtered, adata, annotation_col)

    # Per-annotation top N pool
    pool = {}
    for ann in adata.obs[annotation_col].cat.categories:
        ann_pool = (filtered[filtered['group'].astype(str) == str(ann)]
                    .sort_values('logfoldchanges', ascending=False)
                    .head(pool_size)
                    .reset_index(drop=True))
        if len(ann_pool) > 0:
            pool[ann] = ann_pool
            print(f"   [{ann}]: {len(ann_pool)} candidate markers in pool")
        else:
            print(f"   ⚠️ [{ann}]: no markers pass all filters")

    total = sum(len(v) for v in pool.values())
    print(f"   Pool ready: {total} total candidates across {len(pool)} annotations")
    print(f"   → Agent: PubMed 검색 후 PMID 확인된 gene만 DotPlot 대상 (3-5개/type)")
    return pool


def build_marker_dict_from_selections(agent_selections, annotation_order=None):
    """Agent 선정 완료 후 → DotPlot용 marker_dict 변환.

    Args:
        agent_selections: dict {annotation_name: [gene_list]}
                          Agent가 PMID 검증 후 선정한 3-5개 마커
        annotation_order: list — DotPlot 표시 순서 (없으면 삽입 순)

    Returns:
        dict: {annotation_name: [marker_gene_list]} (ordered)
    """
    if annotation_order:
        ordered = {ann: agent_selections[ann]
                   for ann in annotation_order if ann in agent_selections}
        for ann, markers in agent_selections.items():
            if ann not in ordered:
                ordered[ann] = markers
        return ordered
    return agent_selections


# ── Step 13: Save & Visualize ────────────────────────────────────────
def save_results(subset, major_type, dev_state, annotation_evidence, marker_dict,
                 annotation_col='tier3_annotation',
                 output_dir='annotation_output'):
    """Save subsets, evidence, AND generate visualizations with brackets.

    Args:
        subset: AnnData with tier3 annotations
        major_type: Major cell type name (e.g., 'B_cells')
        dev_state: Developmental state name (e.g., 'Naive')
        annotation_evidence: List of evidence dicts
                             markers[].pmid / .title / .reasoning 포함 필수
        marker_dict: {annotation: [genes]} — Agent가 선정한 3-5개 마커
                     (build_marker_dict_from_selections() 결과)
        annotation_col: Column with tier3 annotations
        output_dir: Output directory
    """
    safe_name = f"{major_type}_{dev_state}".replace(' ', '_')
    os.makedirs(f'{output_dir}/subsets/tier3', exist_ok=True)
    os.makedirs(f'{output_dir}/figures/tier3', exist_ok=True)
    os.makedirs(f'{output_dir}/references', exist_ok=True)

    # 1. Save h5ad
    subset.write(f'{output_dir}/subsets/tier3/{safe_name}.h5ad')

    # 2. Save evidence JSON
    evidence_file = f'{output_dir}/references/annotation_evidence.json'
    all_evidence = []
    if os.path.exists(evidence_file):
        with open(evidence_file, 'r') as f:
            all_evidence = json.load(f)
    all_evidence.extend(annotation_evidence)
    with open(evidence_file, 'w') as f:
        json.dump(all_evidence, f, indent=2, ensure_ascii=False)

    # 3. annotation_order from marker_dict
    annotation_order = list(marker_dict.keys())

    # 4. Generate visualizations (UMAP + bracketed DotPlot + report)
    from tools.visualization_template import save_all_visualizations
    save_all_visualizations(
        adata=subset,
        annotation_col=annotation_col,
        marker_dict=marker_dict,
        evidence_list=annotation_evidence,
        annotation_order=annotation_order,
        title_prefix=f'{safe_name} Tier 3',
        output_dir=f'{output_dir}/figures/tier3',
        report_kwargs={'lineage_name': f'{major_type} → {dev_state}'}
    )

    print(f"✅ Saved: {safe_name}.h5ad + evidence + figures")
