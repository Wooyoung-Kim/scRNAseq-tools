"""
Tier 2 Template: Developmental State Annotation
================================================
Code reference for phases/tier2.md rules.
Read tier2.md for rules/criteria FIRST, then use this code.

Usage: Adapt this template for each major cell type subset.
"""
import re
import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import json
import os
from pathlib import Path

# ── Step 1: Subset Data ──────────────────────────────────────────────
def subset_by_lineage(adata, major_type, tier1_col='tier1_annotation'):
    subset = adata[adata.obs[tier1_col] == major_type].copy()
    print(f"Subsetting {major_type}: {subset.n_obs} cells")
    return subset

# ── Step 2: Re-cluster ───────────────────────────────────────────────
def recluster(subset, resolution=0.1, use_rep='X_pca_harmony'):
    sc.pp.neighbors(subset, use_rep=use_rep, n_neighbors=15, n_pcs=50)
    sc.tl.leiden(subset, resolution=resolution, key_added='tier2_cluster')
    n = subset.obs['tier2_cluster'].nunique()
    print(f"Clusters: {n}")
    return subset


def check_cluster_count(subset, min_clusters=4, cluster_col='tier2_cluster'):
    """Re-cluster 직후 호출 — 클러스터 수 < 4이면 Tier 2/3 생략.

    Args:
        subset: re-clustered AnnData
        min_clusters: 최소 클러스터 수 (default 4)
        cluster_col: cluster column name

    Returns:
        bool: True → 진행, False → Tier 2/3 SKIP

    Usage:
        subset = recluster(subset)
        if not check_cluster_count(subset):
            return None  # SKIP — manifest에 기록
    """
    n = subset.obs[cluster_col].nunique()
    if n < min_clusters:
        print(f"   ⏭️ SKIP: {n} clusters < {min_clusters} → Tier 2/3 생략")
        subset.uns['tier2_skip'] = {
            'reason': f'n_clusters={n} < {min_clusters}',
            'action': 'tier1_annotation 유지, tier2/3 생략'
        }
        return False
    print(f"   ✅ {n} clusters >= {min_clusters} → Tier 2 진행")
    return True

# ── Step 3: RE-COMPUTE DE ────────────────────────────────────────────
def compute_de(subset, groupby='tier2_cluster'):
    sc.tl.rank_genes_groups(subset, groupby=groupby, method='wilcoxon',
                            use_raw=True, pts=True)
    de_df = sc.get.rank_genes_groups_df(subset, group=None)
    print(f"⚠️ DE computed for subset ({subset.n_obs} cells)")
    return de_df

# ── Step 4: TF Activity (MANDATORY) ─────────────────────────────────
def compute_tf_activity(subset, organism='human'):
    """MANDATORY: Cannot proceed without TF activity."""
    print("🔬 TF Activity Analysis (MANDATORY)")
    net = dc.op.collectri(organism=organism)
    dc.run_ulm(mat=subset, net=net, source='source', target='target',
               weight='weight', verbose=True, use_raw=False)
    assert 'ulm_estimate' in subset.obsm, "❌ TF activity not computed!"
    print(f"✅ TF activity: {subset.obsm['ulm_estimate'].shape[1]} TFs")
    return subset

# ── Step 5: Trajectory (MANDATORY if >= 2000 cells) ─────────────────
def compute_trajectory(subset, naive_markers=None):
    """MANDATORY if subset.n_obs >= 2000."""
    import palantir

    if subset.n_obs < 2000:
        print(f"   {subset.n_obs} cells < 2000, skipping trajectory")
        subset.obs['pseudotime'] = np.nan
        subset.obs['pseudotime_category'] = 'N/A'
        return subset

    print(f"   {subset.n_obs} cells >= 2000, computing trajectory...")
    palantir.utils.run_diffusion_maps(subset, n_components=10, pca_key='X_pca')
    palantir.utils.determine_multiscale_space(subset)

    # Find start cell
    if naive_markers is None:
        naive_markers = ['CCR7', 'TCF7', 'SELL', 'IL7R']
    naive_score = np.zeros(subset.n_obs)
    for marker in naive_markers:
        if marker in subset.var_names:
            expr = subset[:, marker].X.toarray().flatten() if hasattr(subset.X, 'toarray') else subset[:, marker].X.flatten()
            naive_score += (expr > 0).astype(float)
    start_cell = subset.obs_names[np.argmax(naive_score)]

    palantir.core.run_palantir(subset, start_cell, num_waypoints=500)
    subset.obs['pseudotime'] = subset.obs['palantir_pseudotime']
    subset.obs['pseudotime_category'] = pd.cut(
        subset.obs['pseudotime'], bins=[0, 0.33, 0.67, 1.0],
        labels=['Early', 'Mid', 'Late']
    )
    print(f"✅ Trajectory: pseudotime [{subset.obs['pseudotime'].min():.2f}, {subset.obs['pseudotime'].max():.2f}]")
    return subset

# ── Step 5.5: Statistical Outlier Detection ──────────────────────────
def detect_tf_outliers(subset, cluster_col='tier2_cluster'):
    """Detect TF activity outliers using z-scores (data-driven)."""
    tf_key = 'score_ulm' if 'score_ulm' in subset.obsm else 'ulm_estimate'
    pval_key = 'padj_ulm' if 'padj_ulm' in subset.obsm else 'ulm_pvals'
    tf_scores = subset.obsm[tf_key]
    tf_pvals = subset.obsm[pval_key]

    outliers = []
    for cluster in subset.obs[cluster_col].unique():
        mask = subset.obs[cluster_col] == cluster
        cluster_tf_mean = tf_scores[mask].mean()
        cluster_tf_pval = tf_pvals[mask].mean()

        for tf in tf_scores.columns:
            global_mean = tf_scores[tf].mean()
            global_std = tf_scores[tf].std()
            cluster_val = cluster_tf_mean[tf]
            cluster_pval = cluster_tf_pval[tf]

            if global_std > 1e-6:
                z_score = (cluster_val - global_mean) / global_std
            else:
                z_score = 0.0

            if abs(z_score) > 2.5 and cluster_val > 0.5 and cluster_pval < 0.05:
                outliers.append({
                    'cluster': str(cluster), 'tf': tf,
                    'activity': float(cluster_val), 'z_score': float(z_score),
                    'pval': float(cluster_pval),
                    'direction': 'HIGH' if z_score > 0 else 'LOW'
                })

    subset.uns['tf_outliers'] = outliers
    print(f"   {'⚠️ ' + str(len(outliers)) + ' TF outliers' if outliers else '✅ No outliers'}")
    return outliers

# ── Step 7: Filter Valid Markers (reasoning용, pct≥0.25) ─────────────
# (Steps 5.6, 8 are agent reasoning steps — see tier2.md)
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
    """Keep genes where pct_nz is highest in the target cell type."""
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
        cell_type: Assigned cell type name (e.g., 'Naive T cells')
        reasonings: dict {gene: reasoning_str} — agent가 작성한 선정 이유
        species: Species context

    Returns:
        list: [
          {
            'gene':      'Ccr7',
            'pct_in':    0.78,
            'log2fc':    2.1,
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

# ── Step 6: Verify Functional Analysis (run BEFORE annotation) ──────
def verify_functional_analysis(subset):
    errors = []
    if 'ulm_estimate' not in subset.obsm:
        errors.append("TF activity (ulm_estimate) not computed")
    if subset.n_obs >= 2000 and 'pseudotime' not in subset.obs.columns:
        errors.append("Trajectory not computed for >= 2000 cells")
    if errors:
        raise AssertionError("❌ INCOMPLETE:\n" + "\n".join(f"  - {e}" for e in errors))
    print("✅ Functional analysis verified")


# ── Step 10: 3-Iteration Reasoning ──────────────────────────────────
def run_3iteration_reasoning(subset, de_df, cluster_id,
                             cluster_col='tier2_cluster',
                             annotation_col='tier2_annotation'):
    """Execute 3-iteration reasoning for a single cluster.

    Reference: reasoning/integrated_format.md

    Iteration 1 (Bioinformatician): Collect DE + TF + Trajectory evidence
    Iteration 2 (Computational Biologist): Evaluate candidates + literature
    Iteration 3 (PI): Integrated decision + confidence scoring

    Args:
        subset: AnnData with TF activity and trajectory computed
        de_df: DE results DataFrame
        cluster_id: Cluster to reason about
        cluster_col: Column with cluster labels
        annotation_col: Column to store annotation

    Returns:
        dict: Reasoning result with annotation, confidence, evidence, candidates
    """
    # ── Iteration 1: Evidence Collection ─────────────────────────────
    # 1a. Top markers
    valid_markers = filter_valid_markers(de_df, cluster_id, top_n=50)
    top_markers = valid_markers[['names', 'pct_nz_group',
                                 'logfoldchanges', 'pvals_adj']].to_dict('records')

    # 1b. TF activity
    tf_key = 'ulm_estimate' if 'ulm_estimate' in subset.obsm else 'score_ulm'
    pval_key = 'ulm_pvals' if 'ulm_pvals' in subset.obsm else 'padj_ulm'
    mask = subset.obs[cluster_col] == cluster_id
    tf_mean = subset.obsm[tf_key][mask].mean()
    tf_pval = subset.obsm[pval_key][mask].mean()
    top_tfs = tf_mean.abs().nlargest(10)
    tf_evidence = [
        {'tf': tf, 'activity': float(tf_mean[tf]),
         'pval': float(tf_pval[tf]),
         'direction': 'HIGH' if tf_mean[tf] > 0 else 'LOW'}
        for tf in top_tfs.index
    ]

    # 1c. Trajectory
    trajectory = {'mean_pseudotime': None, 'category': 'N/A'}
    if 'pseudotime' in subset.obs.columns:
        pt = subset.obs.loc[mask, 'pseudotime']
        if pt.notna().any():
            mean_pt = float(pt.mean())
            if mean_pt <= 0.33:
                cat = 'Early'
            elif mean_pt <= 0.67:
                cat = 'Mid'
            else:
                cat = 'Late'
            trajectory = {'mean_pseudotime': mean_pt, 'category': cat}

    # 1d. Outliers for this cluster
    tf_outliers = [o for o in subset.uns.get('tf_outliers', [])
                   if str(o['cluster']) == str(cluster_id)]

    iteration1 = {
        'markers': top_markers,
        'tf_activity': tf_evidence,
        'trajectory': trajectory,
        'tf_outliers': tf_outliers,
        'n_valid_markers': len(valid_markers)
    }

    # ── Iteration 2: Candidate Reasoning ─────────────────────────────
    # Agent fills candidates based on literature search
    # Template structure for the agent to populate
    iteration2 = {
        'candidates': [],  # Agent populates: [{name, marker_support, tf_support, trajectory_support, confidence}]
        'literature_queries': [],  # Agent populates: [{query, pmids}]
        'instruction': ('Use markers + TF + trajectory to generate 2-3 candidates. '
                        'Search PubMed for each candidate. '
                        'Format: reasoning/integrated_format.md Iteration 2')
    }

    # ── Iteration 3: Decision Template ───────────────────────────────
    iteration3 = {
        'annotation': None,  # Agent fills
        'confidence_score': 0,  # 0-12
        'confidence_level': 'INSUFFICIENT',  # HIGH/MEDIUM/LOW/INSUFFICIENT
        'scoring': {
            'markers': 0,      # 0-3
            'references': 0,   # 0-3
            'tf_consistency': 0,  # 0-3
            'trajectory': 0    # 0-3
        },
        'key_evidence': [],
        'references': [],
        'instruction': ('Integrate all evidence into final decision. '
                        'Score each criterion 0-3 pts. '
                        'Format: reasoning/integrated_format.md Iteration 3')
    }

    reasoning = {
        'cluster_id': str(cluster_id),
        'tier': 2,
        'iteration1': iteration1,
        'iteration2': iteration2,
        'iteration3': iteration3
    }

    print(f"   Cluster {cluster_id}: {len(top_markers)} markers, "
          f"{len(tf_evidence)} TFs, trajectory={trajectory['category']}, "
          f"{len(tf_outliers)} outliers")
    return reasoning

# ── Step 11: Assign Annotations ──────────────────────────────────────
def assign_annotations(subset, cluster_to_annotation, cluster_col='tier2_cluster',
                       annotation_col='tier2_annotation'):
    """Assign tier2 annotations based on cluster-to-annotation mapping.

    Args:
        subset: AnnData with cluster labels
        cluster_to_annotation: dict mapping cluster_id → annotation name
            e.g., {'0': 'Naive_B', '1': 'Memory_B', '2': 'GC_B'}
        cluster_col: Column with cluster IDs
        annotation_col: Column to write annotations to

    Returns:
        AnnData with annotation_col added
    """
    subset.obs[annotation_col] = (
        subset.obs[cluster_col].astype(str).map(cluster_to_annotation)
    )
    unmapped = subset.obs[annotation_col].isna().sum()
    if unmapped > 0:
        print(f"⚠️ {unmapped} cells have no annotation mapping")
        subset.obs[annotation_col] = subset.obs[annotation_col].fillna('Uncertain')
    subset.obs[annotation_col] = pd.Categorical(subset.obs[annotation_col])
    print(f"✅ Assigned {subset.obs[annotation_col].nunique()} annotations to '{annotation_col}'")
    return subset


# ── Step 11a: Annotation 기반 DE 재계산 ──────────────────────────────
def compute_de_by_annotation(subset, annotation_col='tier2_annotation'):
    """Annotation 단위 DE 재계산 — DotPlot pool 구성용.

    cluster 단위 DE(Step 3)와 별도로, 최종 annotation 기준 DE를 재계산.
    groupby=annotation_col, wilcoxon, use_raw=True, pts=True.

    Args:
        subset: tier2_annotation이 할당된 AnnData
        annotation_col: tier2 annotation 컬럼명

    Returns:
        de_df: annotation 단위 DE DataFrame (group = annotation name)
    """
    sc.tl.rank_genes_groups(subset, groupby=annotation_col, method='wilcoxon',
                            use_raw=True, pts=True)
    de_df = sc.get.rank_genes_groups_df(subset, group=None)
    n_groups = subset.obs[annotation_col].nunique()
    print(f"DE re-computed: {n_groups} annotations (groupby='{annotation_col}')")
    return de_df


# ── Steps 11.5: DotPlot 후보 마커 Pool 구성 ──────────────────────────
def build_marker_pool(ann_de_df, adata,
                      annotation_col='tier2_annotation', pool_size=50):
    """Steps 11b-11d: annotation DE에서 DotPlot 후보 마커 pool 구성.

    Step 11b: pct≥0.25, LFC≥1, padj<0.05 strict filter → Top 50 per annotation
              (Tier 2는 lineage subset → developmental marker → pct 0.25 유지)
    Step 11c: HK removal — Rps/Rpl/mt-/Gm*/기타
    Step 11d: Dotplot-highest — argmax pct == target annotation

    ※ Step 11.5는 Agent가 수행:
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
    # Step 11b: strict filter
    filtered = ann_de_df[
        (ann_de_df['pct_nz_group'] >= 0.25) &
        (ann_de_df['logfoldchanges'] >= 1.0) &
        (ann_de_df['pvals_adj'] < 0.05)
    ].copy()
    print(f"   Step 11b: {len(filtered)} genes pass pct≥0.25, LFC≥1, padj<0.05")

    # Step 11c: HK removal
    filtered = remove_housekeeping_genes(filtered)

    # Step 11d: Dotplot-highest
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


# ── Step 12: Save & Visualize ────────────────────────────────────────
def save_results(subset, major_type, annotation_evidence, marker_dict,
                 annotation_col='tier2_annotation',
                 output_dir='annotation_output'):
    """Save subsets, evidence, AND generate visualizations with brackets.

    Args:
        subset: AnnData with tier2 annotations
        major_type: Major cell type name (e.g., 'B_cells')
        annotation_evidence: List of evidence dicts
                             markers[].pmid / .title / .reasoning 포함 필수
        marker_dict: {annotation: [genes]} — Agent가 선정한 3-5개 마커
                     (build_marker_dict_from_selections() 결과)
        annotation_col: Column with tier2 annotations
        output_dir: Output directory
    """
    safe_name = major_type.replace(" ", "_")
    os.makedirs(f'{output_dir}/subsets/tier2', exist_ok=True)
    os.makedirs(f'{output_dir}/figures/tier2', exist_ok=True)
    os.makedirs(f'{output_dir}/references', exist_ok=True)

    # 1. Save h5ad
    subset.write(f'{output_dir}/subsets/tier2/{safe_name}.h5ad')

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
        title_prefix=f'{safe_name} Tier 2',
        output_dir=f'{output_dir}/figures/tier2',
        report_kwargs={'lineage_name': major_type}
    )

    print(f"✅ Saved: {safe_name}.h5ad + evidence + figures")
