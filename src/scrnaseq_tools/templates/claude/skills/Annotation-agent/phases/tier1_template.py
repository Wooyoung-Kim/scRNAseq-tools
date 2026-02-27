"""
Tier 1 Template: Major Cell Type (Lineage) Annotation
=====================================================
Code reference for phases/tier1.md rules.
Read tier1.md for rules/criteria FIRST, then use this code.
"""
import re
import scanpy as sc
import pandas as pd
import numpy as np
import json
import os

# ── Step 1: Clustering ───────────────────────────────────────────────
def cluster(adata, resolution=0.8, use_rep=None):
    if use_rep is None:
        use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, n_pcs=50)
    sc.tl.leiden(adata, resolution=resolution, key_added='tier1_cluster')
    print(f"Clusters: {adata.obs['tier1_cluster'].nunique()}")
    return adata

# ── Step 2: Compute DE ───────────────────────────────────────────────
def compute_de(adata, groupby='tier1_cluster'):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method='wilcoxon',
                            use_raw=True, pts=True)
    de_df = sc.get.rank_genes_groups_df(adata, group=None)
    n_groups = adata.obs[groupby].nunique()
    print(f"DE computed: {n_groups} groups (groupby='{groupby}'), {adata.n_obs} cells")
    return de_df

# ── Step 7a: Annotation 기반 DE 재계산 ───────────────────────────────
def compute_de_by_annotation(adata, annotation_col='tier1_annotation'):
    """클러스터 단위 DE 대신 최종 cell type 단위로 DE 재계산.

    24개 클러스터 → 16개 cell type으로 병합 후, 실제 사용 레이블 기준의
    마커를 얻기 위해 groupby=annotation_col로 rank_genes_groups 재실행.

    Args:
        adata: tier1_annotation이 할당된 AnnData
        annotation_col: tier1 annotation 컬럼명

    Returns:
        de_df: annotation 단위 DE DataFrame (group = cell type name)
    """
    return compute_de(adata, groupby=annotation_col)

# ── Step 2.5: Statistical Outlier Detection ──────────────────────────
def detect_marker_outliers(adata, de_df, cluster_col='tier1_cluster'):
    """Detect marker expression outliers using z-scores (data-driven)."""
    marker_stats = {}
    for cluster in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == cluster
        cluster_cells = adata[mask]
        expr = cluster_cells.X.toarray() if hasattr(cluster_cells.X, 'toarray') else cluster_cells.X
        marker_stats[cluster] = pd.Series(expr.mean(axis=0), index=adata.var_names)

    global_means = pd.concat(marker_stats.values(), axis=1).mean(axis=1)
    global_stds = pd.concat(marker_stats.values(), axis=1).std(axis=1)

    outliers = []
    for cluster, cluster_expr in marker_stats.items():
        for gene in cluster_expr.index:
            if global_stds[gene] > 1e-6:
                z_score = (cluster_expr[gene] - global_means[gene]) / global_stds[gene]
                if z_score > 2.5 and cluster_expr[gene] > 0.1:
                    outliers.append({
                        'cluster': str(cluster), 'gene': gene,
                        'expression': float(cluster_expr[gene]),
                        'z_score': float(z_score), 'direction': 'HIGH'
                    })

    adata.uns['marker_outliers'] = outliers
    print(f"   {'⚠️ ' + str(len(outliers)) + ' marker outliers' if outliers else '✅ No outliers'}")
    return outliers

# ── Step 2.6: Cross-Lineage Contamination Detection ─────────────────
def detect_contamination(adata, de_df, cluster_col='tier1_cluster'):
    """Detect cross-lineage contamination (data-driven, no hardcoding)."""
    contamination_issues = []
    for cluster in adata.obs[cluster_col].unique():
        cluster_de = de_df[de_df['group'] == cluster]
        cluster_de = cluster_de[
            (cluster_de['pct_nz_group'] >= 0.25) &
            (cluster_de['logfoldchanges'] >= 1.0) &
            (cluster_de['pvals_adj'] < 0.05)
        ].nlargest(50, 'logfoldchanges')
        top_markers = cluster_de['names'].tolist()

        for other_cluster in adata.obs[cluster_col].unique():
            if other_cluster == cluster:
                continue
            other_cells = adata[adata.obs[cluster_col] == other_cluster]
            contamination_score = 0
            contaminating_markers = []

            for marker in top_markers:
                if marker in adata.var_names:
                    expr = other_cells[:, marker].X.toarray().flatten() if hasattr(other_cells.X, 'toarray') else other_cells[:, marker].X.flatten()
                    pct_expr = (expr > 0).mean()
                    if pct_expr > 0.10:
                        contamination_score += pct_expr
                        contaminating_markers.append(f"{marker} ({pct_expr:.1%})")

            if contamination_score > 0.5:
                contamination_issues.append({
                    'source_cluster': str(cluster),
                    'target_cluster': str(other_cluster),
                    'contamination_score': float(contamination_score),
                    'markers': contaminating_markers[:3],
                    'severity': 'HIGH' if contamination_score > 1.0 else 'MEDIUM'
                })

    adata.uns['contamination_issues'] = contamination_issues
    print(f"   {'⚠️ ' + str(len(contamination_issues)) + ' contamination issues' if contamination_issues else '✅ No contamination'}")
    return contamination_issues

# ── Step 3: Filter Valid Markers (initial clustering phase) ──────────
# (Steps 4, 5 are agent reasoning steps — see tier1.md)
def filter_valid_markers(de_df, cluster_id, top_n=50):
    """pct >= 25%, LFC >= 1, padj < 0.05, Top 50 (for initial cluster-level reasoning)"""
    cluster_df = de_df[de_df['group'] == cluster_id].copy()
    valid = cluster_df[
        (cluster_df['pct_nz_group'] >= 0.25) &
        (cluster_df['logfoldchanges'] >= 1.0) &
        (cluster_df['pvals_adj'] < 0.05)
    ]
    return valid.nlargest(top_n, 'logfoldchanges')

# ── Steps 7b-7d: Marker Selection Pipeline (post-annotation) ─────────

# Housekeeping gene patterns
_HK_PREFIX_RE = re.compile(
    r'^(Rps|Rpl|Mrps|Mrpl|mt-|Hmgb|Gm\d)',
    re.IGNORECASE
)
_HK_HISTONE_RE = re.compile(
    r'^H[1-4][a-zA-Z0-9\-]|^HIST[123]',
    re.IGNORECASE
)
_HK_EXACT = frozenset({
    'Malat1', 'Tmsb4x', 'Tmsb10', 'Actb', 'Actg1',
    'Eef1a1', 'Eef2', 'Tpt1', 'Ppia', 'Fau', 'B2m',
    'Calm1', 'Cmss1', 'Psap', 'Ybx1', 'Ptma', 'Npm1',
    # uppercase
    'MALAT1', 'TMSB4X', 'TMSB10', 'ACTB', 'ACTG1',
    'EEF1A1', 'EEF2', 'TPT1', 'PPIA', 'FAU', 'B2M',
    'CALM1', 'CMSS1', 'PSAP', 'YBX1', 'PTMA', 'NPM1',
})

def is_housekeeping(gene: str) -> bool:
    return (gene in _HK_EXACT or
            bool(_HK_PREFIX_RE.match(gene)) or
            bool(_HK_HISTONE_RE.match(gene)))

def remove_housekeeping_genes(de_df: pd.DataFrame) -> pd.DataFrame:
    """Step 7c: Remove ribosomal, mitochondrial, and other housekeeping genes."""
    mask = ~de_df['names'].apply(is_housekeeping)
    removed = (~mask).sum()
    if removed:
        print(f"   HK removed: {removed} genes "
              f"({', '.join(de_df.loc[~mask, 'names'].unique()[:5])}...)")
    return de_df[mask].copy()

def verify_dotplot_highest(de_df: pd.DataFrame, adata,
                           annotation_col: str) -> pd.DataFrame:
    """Step 7d: Keep genes where pct_nz is highest in the target cell type.

    For each candidate gene, computes pct_nz across all cell types and
    verifies argmax == target. Genes like Aif1 (highest in cDC1) or Ccr7
    (highest in T_cells) are automatically removed from wrong cell types.

    Uses adata.X (log-normalized) for pct computation.
    """
    candidate_genes = de_df['names'].unique().tolist()
    var_names = list(adata.var_names)

    # Restrict to genes present in adata
    valid_genes = [g for g in candidate_genes if g in var_names]
    gene_idx = [var_names.index(g) for g in valid_genes]

    # Extract candidate gene columns only (efficient)
    X_sub = adata.X[:, gene_idx]
    if hasattr(X_sub, 'toarray'):
        X_sub = X_sub.toarray()
    else:
        X_sub = np.asarray(X_sub)

    cell_type_labels = adata.obs[annotation_col].astype(str).values
    unique_cts = list(dict.fromkeys(cell_type_labels))  # ordered, unique

    # pct_matrix: shape (n_celltypes, n_genes)
    pct_matrix = np.zeros((len(unique_cts), len(valid_genes)))
    for i, ct in enumerate(unique_cts):
        mask = cell_type_labels == ct
        pct_matrix[i, :] = (X_sub[mask, :] > 0).mean(axis=0)

    # Gene → argmax cell type
    gene_to_highest = {
        gene: unique_cts[pct_matrix[:, j].argmax()]
        for j, gene in enumerate(valid_genes)
    }

    def _is_highest(row):
        return gene_to_highest.get(row['names'], '') == str(row['group'])

    mask = de_df.apply(_is_highest, axis=1)
    removed = (~mask).sum()
    print(f"   Dotplot-highest: {mask.sum()} pass, {removed} removed")
    return de_df[mask].copy()

# ── Step 8: Cross-Contamination Verification (post-annotation) ─────────
def verify_contamination_post_annotation(adata, annotation_col, de_df):
    """Verify cross-contamination using data-driven approach after annotation."""
    annotation_signatures = {}
    for ann in adata.obs[annotation_col].unique():
        ann_mask = adata.obs[annotation_col] == ann
        clusters_in_ann = adata.obs[ann_mask]['tier1_cluster'].unique()
        top_markers = []
        for cluster in clusters_in_ann:
            cluster_de = de_df[de_df['group'] == cluster]
            cluster_top = cluster_de[
                (cluster_de['pct_nz_group'] >= 0.50) &
                (cluster_de['logfoldchanges'] >= 2.0) &
                (cluster_de['pvals_adj'] < 0.001)
            ].nlargest(5, 'logfoldchanges')['names'].tolist()
            top_markers.extend(cluster_top)
        annotation_signatures[ann] = list(set(top_markers))[:10]

    issues = []
    for ann in adata.obs[annotation_col].unique():
        mask = adata.obs[annotation_col] == ann
        subset = adata[mask]
        for other_ann, other_markers in annotation_signatures.items():
            if other_ann == ann:
                continue
            contamination_count = 0
            contaminating_genes = []
            for marker in other_markers:
                if marker in subset.var_names:
                    expr = subset[:, marker].X.toarray().flatten() if hasattr(subset.X, 'toarray') else subset[:, marker].X.flatten()
                    pct = (expr > 0).mean()
                    if pct > 0.05:
                        contamination_count += 1
                        contaminating_genes.append(f"{marker} ({pct:.1%})")
            if contamination_count >= 2:
                issues.append({
                    'annotation': ann, 'contaminated_by': other_ann,
                    'n_markers': contamination_count,
                    'markers': contaminating_genes[:3],
                    'severity': 'HIGH' if contamination_count >= 3 else 'MEDIUM'
                })
    return issues

# ── Step 6 / Step 7e: Per-Marker PMID Verification + Evidence ────────
def build_per_marker_evidence(markers, cell_type, reasonings=None, species='human'):
    """Verify each marker and return per-marker evidence list for annotation_evidence.

    반환값은 annotation_evidence['markers'] 리스트로 직접 사용됩니다.
    report.md의 Marker Evidence 표에 Gene / pct / LFC / PMID / Title / Reasoning
    컬럼으로 출력됩니다.

    Args:
        markers: List of dicts with gene stats, e.g.:
                 [{'gene': 'Cd3e', 'pct_in': 0.85, 'log2fc': 2.5}, ...]
                 또는 gene 이름 list (pct/LFC는 0으로 설정됨)
        cell_type: Assigned cell type name (e.g., 'T cells')
        reasonings: dict {gene: reasoning_str} — agent가 작성한 선정 이유
                    (없으면 빈 문자열로 저장)
        species: Species context

    Returns:
        list: [
          {
            'gene': 'Cd3e',
            'pct_in': 0.85,      # fraction (0-1)
            'log2fc': 2.5,
            'pmid': '12345678',  # 최상위 PMID
            'title': '...',       # 논문 제목
            'reasoning': '...'    # agent 선정 이유 (REQUIRED)
          }, ...
        ]
    """
    from tools.dynamic_knowledge import pubmed_search
    import time

    reasonings = reasonings or {}

    # Normalize markers input
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

# ── Step 6.5: 3-Iteration Reasoning ───────────────────────────────
def run_3iteration_reasoning(adata, de_df, cluster_id,
                             cluster_col='tier1_cluster',
                             annotation_col='tier1_annotation'):
    """Execute 3-iteration reasoning for a single cluster.

    Reference: reasoning/integrated_format.md

    Tier 1 uses: DE markers + statistical outliers + contamination evidence.
    (No TF/pathway/trajectory at this tier.)

    Iteration 1 (Bioinformatician): Collect DE + outlier evidence
    Iteration 2 (Computational Biologist): Dynamic knowledge search + candidates
    Iteration 3 (PI): Integrated decision + confidence scoring

    Args:
        adata: AnnData with outlier/contamination detection done
        de_df: DE results DataFrame
        cluster_id: Cluster to reason about
        cluster_col: Column with cluster labels
        annotation_col: Column to store annotation

    Returns:
        dict: Reasoning result with annotation, confidence, evidence, candidates
    """
    # ── Iteration 1: Evidence Collection ─────────────────────────────
    # 1a. All valid markers (Top 50)
    valid_markers = filter_valid_markers(de_df, cluster_id, top_n=50)
    top_markers = valid_markers[['names', 'pct_nz_group',
                                 'logfoldchanges', 'pvals_adj']].to_dict('records')

    # 1b. Outliers for this cluster
    marker_outliers = [o for o in adata.uns.get('marker_outliers', [])
                       if str(o['cluster']) == str(cluster_id)]

    # 1c. Contamination for this cluster
    contamination = [c for c in adata.uns.get('contamination_issues', [])
                     if str(c['source_cluster']) == str(cluster_id)
                     or str(c['target_cluster']) == str(cluster_id)]

    iteration1 = {
        'markers': top_markers,
        'marker_outliers': marker_outliers,
        'contamination': contamination,
        'n_valid_markers': len(valid_markers)
    }

    # ── Iteration 2: Candidate Reasoning ─────────────────────────────
    iteration2 = {
        'candidates': [],  # Agent populates via dynamic_knowledge.py
        'literature_queries': [],
        'instruction': ('Use top DE markers to search dynamic_knowledge.py. '
                        'Generate 2-3 lineage candidates ranked by confidence. '
                        'NO hardcoded marker lists.')
    }

    # ── Iteration 3: Decision Template ───────────────────────────────
    iteration3 = {
        'annotation': None,
        'confidence_score': 0,
        'confidence_level': 'INSUFFICIENT',
        'scoring': {
            'markers': 0,
            'references': 0,
            'outlier_consistency': 0,
            'contamination_check': 0
        },
        'key_evidence': [],
        'references': [],
        'instruction': ('Integrate markers + literature + outliers into decision. '
                        'Score each criterion 0-3 pts. '
                        'Flag as NOVEL if no literature match.')
    }

    reasoning = {
        'cluster_id': str(cluster_id),
        'tier': 1,
        'iteration1': iteration1,
        'iteration2': iteration2,
        'iteration3': iteration3
    }

    print(f"   Cluster {cluster_id}: {len(top_markers)} markers, "
          f"{len(marker_outliers)} outliers, {len(contamination)} contamination")
    return reasoning


# ── Step 7: Assign Annotations ───────────────────────────────────────
def assign_annotations(adata, cluster_to_annotation, cluster_col='tier1_cluster',
                       annotation_col='tier1_annotation'):
    """Assign tier1 annotations based on cluster-to-annotation mapping.

    Args:
        adata: AnnData with cluster labels
        cluster_to_annotation: dict mapping cluster_id → annotation name
            e.g., {'0': 'T cells', '1': 'B lineage', '2': 'Myeloid'}
        cluster_col: Column with cluster IDs
        annotation_col: Column to write annotations to

    Returns:
        AnnData with annotation_col added
    """
    adata.obs[annotation_col] = (
        adata.obs[cluster_col].astype(str).map(cluster_to_annotation)
    )
    unmapped = adata.obs[annotation_col].isna().sum()
    if unmapped > 0:
        print(f"⚠️ {unmapped} cells have no annotation mapping")
        adata.obs[annotation_col] = adata.obs[annotation_col].fillna('Uncertain')
    adata.obs[annotation_col] = pd.Categorical(adata.obs[annotation_col])
    print(f"✅ Assigned {adata.obs[annotation_col].nunique()} annotations to '{annotation_col}'")
    return adata


# ── Steps 7b-7d: DotPlot 후보 마커 Pool 구성 (Code) ─────────────────
def build_marker_pool(ann_de_df, adata,
                      annotation_col='tier1_annotation', pool_size=50):
    """Steps 7b-7d: annotation DE에서 DotPlot 후보 마커 pool 구성.

    Step 7b: pct≥0.40, LFC≥1, padj<0.05 strict filter → Top 50 per cell type
             (annotation DE groupby=tier1_annotation 기준)
    Step 7c: HK removal — Rps/Rpl/mt-/Gm*/기타
    Step 7d: Dotplot-highest — argmax pct == target cell type

    ※ Step 7e는 Agent가 수행:
       pool에서 PubMed 검색 → PMID 확인된 gene만 DotPlot 대상
       PMID 없는 gene은 DotPlot에 포함하지 않음
       3-5개 최종 선정 (PMID + reasoning 필수)

    Args:
        ann_de_df: Annotation-based DE DataFrame (from compute_de_by_annotation)
                   group = cell type name
        adata: AnnData with annotation_col and adata.X (log-normalized)
        annotation_col: Column with cell type labels
        pool_size: Candidate pool size per cell type (default 50)

    Returns:
        dict: {annotation_name: DataFrame (names, pct_nz_group, logfoldchanges, pvals_adj)}
    """
    # ── Step 7b: strict filter ────────────────────────────────────────
    filtered = ann_de_df[
        (ann_de_df['pct_nz_group'] >= 0.40) &
        (ann_de_df['logfoldchanges'] >= 1.0) &
        (ann_de_df['pvals_adj'] < 0.05)
    ].copy()
    print(f"   Step 7b: {len(filtered)} genes pass pct≥0.40, LFC≥1, padj<0.05")

    # ── Step 7c: HK removal ───────────────────────────────────────────
    filtered = remove_housekeeping_genes(filtered)

    # ── Step 7d: Dotplot-highest ──────────────────────────────────────
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
    print(f"   Pool ready: {total} total candidates across {len(pool)} types")
    print(f"   → Agent: PubMed 검색 후 PMID 확인된 gene만 DotPlot 대상 (3-5개/type)")
    return pool


# ── Step 7e: Agent가 pool에서 최종 마커 선정 후 dict 구성 ──────────────
def build_marker_dict_from_selections(agent_selections, annotation_order=None):
    """Step 7e 완료 후 agent_selections → DotPlot용 marker_dict 변환.

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


def build_marker_dict(ann_markers, annotation_order=None):
    """Convert per-annotation markers to ordered dict for dotplot brackets."""
    if annotation_order:
        ordered = {}
        for ann in annotation_order:
            if ann in ann_markers:
                ordered[ann] = ann_markers[ann]
        for ann, markers in ann_markers.items():
            if ann not in ordered:
                ordered[ann] = markers
        return ordered
    return ann_markers


# ── Step 9: Save Results + Visualizations ─────────────────────────────
def save_results(adata, annotation_evidence, marker_dict,
                 annotation_col='tier1_annotation',
                 output_dir='annotation_output'):
    """Save annotated data, evidence, AND generate visualizations with brackets.

    Args:
        adata: AnnData with tier1 annotations
        annotation_evidence: List of evidence dicts
        marker_dict: {annotation: [genes]} — Agent가 Step 7e에서 선정한 3-5개 마커
                     (build_marker_dict_from_selections() 결과)
        annotation_col: Column with tier1 annotations
        output_dir: Output directory
    """
    os.makedirs(f'{output_dir}/subsets/tier1', exist_ok=True)
    os.makedirs(f'{output_dir}/figures/tier1', exist_ok=True)
    os.makedirs(f'{output_dir}/references', exist_ok=True)

    # 1. Save full annotated h5ad
    adata.write(f'{output_dir}/tier1_annotated.h5ad')

    # 2. Save per-type subsets
    for major_type in adata.obs[annotation_col].unique():
        subset = adata[adata.obs[annotation_col] == major_type].copy()
        safe_name = major_type.replace(' ', '_')
        subset.write(f'{output_dir}/subsets/tier1/{safe_name}.h5ad')
        print(f"Saved: {safe_name} ({subset.n_obs} cells)")

    # 3. Save evidence JSON
    evidence_file = f'{output_dir}/references/tier1_annotation_evidence.json'
    with open(evidence_file, 'w') as f:
        json.dump(annotation_evidence, f, indent=2, ensure_ascii=False)

    # 4. annotation_order from marker_dict
    annotation_order = list(marker_dict.keys())

    # 5. Generate visualizations (UMAP + bracketed DotPlot + report)
    from tools.visualization_template import save_all_visualizations
    save_all_visualizations(
        adata=adata,
        annotation_col=annotation_col,
        marker_dict=marker_dict,
        evidence_list=annotation_evidence,
        annotation_order=annotation_order,
        title_prefix='Tier 1',
        output_dir=f'{output_dir}/figures/tier1'
    )

    print(f"✅ Saved: tier1_annotated.h5ad + evidence + figures")
