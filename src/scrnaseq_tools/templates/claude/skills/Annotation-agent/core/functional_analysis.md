# Functional Analysis Module (v3 - DATA-DRIVEN)

**MANDATORY** functional analysis code for Tier 2 and Tier 3.

```
╔══════════════════════════════════════════════════════════════════════╗
║  🚨 v3 CHANGE: NO HARDCODED TF/PATHWAY ASSOCIATIONS                  ║
║  ALL interpretations driven by statistical analysis + LLM reasoning  ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Overview

```
╔══════════════════════════════════════════════════════════════════════╗
║  Tier 2: TF Activity (MANDATORY) + Trajectory (if >= 2000 cells)     ║
║  Tier 3: Pathway Activity (MANDATORY)                                ║
║                                                                      ║
║  Results stored in (decoupler v2):                                   ║
║  - adata.obsm['score_ulm'] / ['padj_ulm']    (TF)                    ║
║  - adata.obsm['score_mlm'] / ['padj_mlm']    (Pathway)               ║
║  - adata.obs['pseudotime']                    (Trajectory)            ║
║                                                                      ║
║  COMPATIBILITY: Create aliases for legacy code:                       ║
║  - adata.obsm['ulm_estimate'] = adata.obsm['score_ulm']              ║
║  - adata.obsm['ulm_pvals']    = adata.obsm['padj_ulm']               ║
║                                                                      ║
║  🆕 v3 ADDITIONS:                                                     ║
║  - compute_tf_statistics() → Global stats + z-scores                 ║
║  - compute_pathway_statistics() → Global stats + z-scores            ║
║  - interpret_tf_activity() → NO hardcoded associations               ║
║  - interpret_pathway_activity() → NO hardcoded associations          ║
║  - Returns reasoning_prompt to guide LLM literature search           ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Key Principles (v3)

**REMOVED:**
- ❌ `TF_ASSOCIATIONS` dictionary (hardcoded TF meanings)
- ❌ `PATHWAY_ASSOCIATIONS` dictionary (hardcoded pathway meanings)
- ❌ Automated state suggestions based on hardcoded rules

**ADDED:**
- ✅ Statistical analysis: z-scores, outlier detection
- ✅ Reasoning prompts: guide LLM to search literature
- ✅ Data-driven interpretation: let LLM + literature decide meaning
- ✅ Context-aware: considers parent cell type for specialized states

**Philosophy:**
```
OLD: "TBX21 high" → automatically interpret as "Th1/Effector"
NEW: "TBX21 z-score=3.5" → prompt LLM: "Search 'TBX21 [cell_lineage]'"
     → LLM finds: "TBX21 in B cells = age-associated B cells"
```

---

## 1. TF Activity Analysis (Tier 2 - MANDATORY)

### Run TF Activity

```python
import decoupler as dc
import pandas as pd
import numpy as np

def run_tf_activity(adata, organism='human', split_complexes=False):
    """
    Compute TF activity using decoupler + CollecTRI.

    MANDATORY for Tier 2 annotation.

    Parameters:
    -----------
    adata : AnnData
        Subset for one major type
    organism : str
        'human' or 'mouse'
    split_complexes : bool
        Split TF complexes into individual TFs

    Returns:
    --------
    adata : AnnData with TF activity in obsm
    """
    print("🔬 Computing TF activity (MANDATORY)...")

    # Get CollecTRI network (TF-target regulatory network)
    net = dc.op.collectri(organism=organism, split_complexes=split_complexes)
    print(f"   CollecTRI network: {len(net)} TF-target pairs")

    # Run ULM (Univariate Linear Model) — decoupler v2 API
    dc.mt.ulm(adata, net, verbose=True, raw=False)

    # decoupler v2 stores results as 'score_ulm'/'padj_ulm'
    # Create legacy aliases for backward compatibility
    if 'score_ulm' in adata.obsm:
        adata.obsm['ulm_estimate'] = adata.obsm['score_ulm']
        adata.obsm['ulm_pvals'] = adata.obsm['padj_ulm']

    # Verify results (check both v2 and legacy keys)
    assert 'score_ulm' in adata.obsm or 'ulm_estimate' in adata.obsm, \
        "❌ TF activity computation failed!"

    tf_key = 'score_ulm' if 'score_ulm' in adata.obsm else 'ulm_estimate'
    n_tfs = adata.obsm[tf_key].shape[1]
    print(f"✅ TF activity computed: {n_tfs} TFs")

    return adata


def get_tf_activity_per_cluster(adata, cluster_col, top_n=10):
    """
    Get top TF activities per cluster.

    Returns:
    --------
    dict : {cluster_id: DataFrame with TF, activity, pval}
    """
    tf_results = {}

    # Support both decoupler v2 (score_ulm) and legacy (ulm_estimate) key names
    tf_key = 'score_ulm' if 'score_ulm' in adata.obsm else 'ulm_estimate'
    pval_key = 'padj_ulm' if 'padj_ulm' in adata.obsm else 'ulm_pvals'
    tf_activities = adata.obsm[tf_key]
    tf_pvals = adata.obsm[pval_key]
    tf_names = tf_activities.columns.tolist()

    for cluster in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == cluster
        cluster_tf = tf_activities[mask].mean()
        cluster_pval = tf_pvals[mask].mean()

        # Create DataFrame
        df = pd.DataFrame({
            'tf': tf_names,
            'activity': cluster_tf.values,
            'pval': cluster_pval.values
        })

        # Filter significant and sort by activity
        df_sig = df[df['pval'] < 0.05].copy()
        df_sig = df_sig.sort_values('activity', ascending=False).head(top_n)

        tf_results[cluster] = df_sig

    return tf_results
```

### TF Activity Interpretation (DATA-DRIVEN, NO HARDCODING)

```python
def compute_tf_statistics(adata, cluster_col='tier2_cluster'):
    """
    Compute statistical measures for TF activity across clusters.
    NO hardcoded TF associations - purely data-driven.

    Returns:
    --------
    dict with global statistics and per-cluster z-scores
    """
    tf_key = 'score_ulm' if 'score_ulm' in adata.obsm else 'ulm_estimate'
    pval_key = 'padj_ulm' if 'padj_ulm' in adata.obsm else 'ulm_pvals'

    tf_scores = adata.obsm[tf_key]
    tf_pvals = adata.obsm[pval_key]

    # Global statistics
    global_stats = {
        'mean': tf_scores.mean(),
        'std': tf_scores.std(),
        'median': tf_scores.median()
    }

    # Per-cluster statistics with z-scores
    cluster_stats = {}
    for cluster in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == cluster
        cluster_tf = tf_scores[mask].mean()
        cluster_pval = tf_pvals[mask].mean()

        # Compute z-scores for each TF
        z_scores = {}
        for tf in tf_scores.columns:
            if global_stats['std'][tf] > 1e-6:
                z = (cluster_tf[tf] - global_stats['mean'][tf]) / global_stats['std'][tf]
            else:
                z = 0.0
            z_scores[tf] = z

        cluster_stats[cluster] = {
            'tf_activity': cluster_tf,
            'tf_pvals': cluster_pval,
            'z_scores': z_scores
        }

    return {
        'global_stats': global_stats,
        'cluster_stats': cluster_stats
    }


def interpret_tf_activity(tf_df, cluster_stats=None, cluster_id=None):
    """
    Interpret TF activity for reasoning - DATA-DRIVEN approach.
    NO hardcoded associations - returns statistical analysis for LLM to interpret.

    Parameters:
    -----------
    tf_df : DataFrame
        Top TFs with activity and p-values
    cluster_stats : dict (optional)
        Output from compute_tf_statistics()
    cluster_id : str (optional)
        Cluster ID for z-score lookup

    Returns:
    --------
    dict with statistical interpretation prompts for LLM reasoning
    """
    top_tfs = tf_df['tf'].tolist()[:10]

    # Collect TF data with statistical context
    tf_data = []
    for _, row in tf_df.iterrows():
        tf = row['tf']
        activity = row['activity']
        pval = row['pval']

        tf_info = {
            'tf': tf,
            'activity': activity,
            'pval': pval,
            'is_significant': pval < 0.05,
        }

        # Add z-score if available
        if cluster_stats and cluster_id and cluster_id in cluster_stats['cluster_stats']:
            z_score = cluster_stats['cluster_stats'][cluster_id]['z_scores'].get(tf, 0.0)
            tf_info['z_score'] = z_score
            tf_info['is_outlier'] = abs(z_score) > 2.5
        else:
            tf_info['z_score'] = None
            tf_info['is_outlier'] = False

        tf_data.append(tf_info)

    # Return structured data for LLM reasoning (NO interpretations)
    return {
        'top_tfs': top_tfs,
        'tf_data': tf_data,
        'n_significant': sum(1 for t in tf_data if t['is_significant']),
        'n_outliers': sum(1 for t in tf_data if t['is_outlier']),
        'outlier_tfs': [t['tf'] for t in tf_data if t['is_outlier']],
        'reasoning_prompt': generate_tf_reasoning_prompt(tf_data)
    }


def generate_tf_reasoning_prompt(tf_data):
    """
    Generate a reasoning prompt for LLM based on TF statistics.
    NO hardcoded interpretations - guides LLM to search literature.
    """
    prompt_parts = []

    # High activity TFs
    high_activity = [t for t in tf_data if t['activity'] > 1.0 and t['is_significant']]
    if high_activity:
        tf_names = ', '.join([t['tf'] for t in high_activity[:5]])
        prompt_parts.append(
            f"High TF activity: {tf_names}\n"
            f"→ Search literature: '[TF_name] [cell_lineage]' to identify cell state"
        )

    # Outlier TFs (statistically unusual)
    outliers = [t for t in tf_data if t['is_outlier']]
    if outliers:
        outlier_names = ', '.join([f"{t['tf']} (z={t['z_score']:.2f})" for t in outliers])
        prompt_parts.append(
            f"⚠️ STATISTICAL OUTLIERS: {outlier_names}\n"
            f"→ These TFs are unusually active compared to other clusters\n"
            f"→ PRIORITY: Search '[outlier_TF] [cell_lineage]' for specialized cell types"
        )

    # Multiple significant TFs
    if len([t for t in tf_data if t['is_significant']]) >= 5:
        prompt_parts.append(
            "Multiple TFs active - suggests complex transcriptional program\n"
            "→ Search combinations: '[TF1] [TF2] [cell_lineage]'"
        )

    return "\n\n".join(prompt_parts) if prompt_parts else "Standard TF pattern - proceed with marker-based annotation"
```

---

## 2. Trajectory Analysis (Tier 2 - MANDATORY if >= 2000 cells)

### Run Trajectory

```python
import palantir

def run_trajectory(adata, naive_markers, start_cell=None, n_waypoints=500, min_cells=2000):
    """
    Compute trajectory/pseudotime using Palantir.

    MANDATORY for Tier 2 if >= 2000 cells.

    Parameters:
    -----------
    adata : AnnData
        Subset for one major type
    naive_markers : list (REQUIRED)
        Lineage-specific early/progenitor markers for auto-detecting start cell.
        B lineage: ['DNTT', 'RAG1', 'ERG']
        T cells: ['CCR7', 'TCF7', 'SELL']
    start_cell : str
        Cell barcode for trajectory start (auto-detected if None)
    n_waypoints : int
        Number of waypoints for Palantir

    Returns:
    --------
    adata : AnnData with pseudotime in obs
    """
    if adata.n_obs < min_cells:
        print(f"⚠️ Skipping trajectory: {adata.n_obs} < {min_cells} cells")
        adata.obs['pseudotime'] = np.nan
        adata.obs['pseudotime_category'] = 'N/A'
        return adata

    print("🔬 Computing trajectory (MANDATORY)...")

    # Run diffusion maps
    palantir.utils.run_diffusion_maps(adata, n_components=10, pca_key='X_pca')

    # Determine multiscale space
    palantir.utils.determine_multiscale_space(adata)

    # Auto-detect start cell if not provided
    if start_cell is None:
        start_cell = find_earliest_cell(adata, naive_markers)
        print(f"   Auto-detected start cell: {start_cell}")

    # Run Palantir (results stored directly in AnnData)
    palantir.core.run_palantir(
        adata,
        start_cell,
        num_waypoints=n_waypoints
    )

    # Copy to standard key
    adata.obs['pseudotime'] = adata.obs['palantir_pseudotime']

    # Categorize pseudotime
    adata.obs['pseudotime_category'] = pd.cut(
        adata.obs['pseudotime'],
        bins=[0, 0.33, 0.67, 1.0],
        labels=['Early', 'Mid', 'Late']
    )

    # Verify
    assert 'palantir_pseudotime' in adata.obs.columns, "❌ Pseudotime computation failed!"
    print(f"✅ Trajectory computed: pseudotime range [{adata.obs['pseudotime'].min():.2f}, {adata.obs['pseudotime'].max():.2f}]")

    return adata


def find_earliest_cell(adata, naive_markers):
    """
    Auto-detect the earliest/most naive cell for trajectory start.

    Parameters:
    -----------
    naive_markers : list (REQUIRED, no default)
        Lineage-specific early/progenitor markers.
        B lineage: ['DNTT', 'RAG1', 'ERG']
        T cells: ['CCR7', 'TCF7', 'SELL']

    Alternative: Use palantir.utils.early_cell() when cell type info is available.
    """
    # Score cells by naive marker expression
    naive_score = np.zeros(adata.n_obs)

    for marker in naive_markers:
        if marker in adata.var_names:
            if hasattr(adata.X, 'toarray'):
                expr = adata[:, marker].X.toarray().flatten()
            else:
                expr = adata[:, marker].X.flatten()
            naive_score += (expr > 0).astype(float)

    # Return cell with highest naive score
    best_idx = np.argmax(naive_score)
    return adata.obs_names[best_idx]


def get_pseudotime_per_cluster(adata, cluster_col):
    """
    Get pseudotime statistics per cluster.

    Returns:
    --------
    dict : {cluster_id: {'mean': float, 'category': str, 'interpretation': str}}
    """
    results = {}

    for cluster in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == cluster
        pt = adata.obs.loc[mask, 'pseudotime']

        if pt.isna().all():
            results[cluster] = {
                'mean': np.nan,
                'std': np.nan,
                'category': 'N/A',
                'interpretation': 'Trajectory not computed'
            }
        else:
            mean_pt = pt.mean()
            std_pt = pt.std()

            if mean_pt < 0.33:
                category = 'Early'
                interp = 'Progenitor/Naive state (early differentiation)'
            elif mean_pt < 0.67:
                category = 'Mid'
                interp = 'Transitional state (active differentiation)'
            else:
                category = 'Late'
                interp = 'Terminal/Effector state (late differentiation)'

            results[cluster] = {
                'mean': mean_pt,
                'std': std_pt,
                'category': category,
                'interpretation': interp
            }

    return results
```

---

## 3. Pathway Activity Analysis (Tier 3 - MANDATORY)

### Run Pathway Activity

```python
def run_pathway_activity(adata, organism='human', top=500):
    """
    Compute pathway activity using decoupler + PROGENy.

    MANDATORY for Tier 3 annotation.

    Parameters:
    -----------
    adata : AnnData
        Subset for one developmental state
    organism : str
        'human' or 'mouse'
    top : int
        Number of top genes per pathway

    Returns:
    --------
    adata : AnnData with pathway activity in obsm
    """
    print("🔬 Computing Pathway activity (MANDATORY)...")

    # Get PROGENy network
    net = dc.op.progeny(organism=organism, top=top)
    print(f"   PROGENy network: {len(net['source'].unique())} pathways")

    # Run MLM (Multivariate Linear Model) — decoupler v2 API
    dc.mt.mlm(adata, net, verbose=True, raw=False)

    # decoupler v2 stores results as 'score_mlm'/'padj_mlm'
    # Create legacy aliases for backward compatibility
    if 'score_mlm' in adata.obsm:
        adata.obsm['mlm_estimate'] = adata.obsm['score_mlm']
        adata.obsm['mlm_pvals'] = adata.obsm['padj_mlm']

    # Verify results (check both v2 and legacy keys)
    assert 'score_mlm' in adata.obsm or 'mlm_estimate' in adata.obsm, \
        "❌ Pathway activity computation failed!"

    pw_key = 'score_mlm' if 'score_mlm' in adata.obsm else 'mlm_estimate'
    n_pathways = adata.obsm[pw_key].shape[1]
    print(f"✅ Pathway activity computed: {n_pathways} pathways")

    return adata


def get_pathway_activity_per_cluster(adata, cluster_col, top_n=10):
    """
    Get top pathway activities per cluster.

    Returns:
    --------
    dict : {cluster_id: DataFrame with pathway, activity, pval}
    """
    pathway_results = {}

    # Support both decoupler v2 (score_mlm) and legacy (mlm_estimate) key names
    pw_key = 'score_mlm' if 'score_mlm' in adata.obsm else 'mlm_estimate'
    pval_key = 'padj_mlm' if 'padj_mlm' in adata.obsm else 'mlm_pvals'
    pathway_activities = adata.obsm[pw_key]
    pathway_pvals = adata.obsm[pval_key]
    pathway_names = pathway_activities.columns.tolist()

    for cluster in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == cluster
        cluster_pathway = pathway_activities[mask].mean()
        cluster_pval = pathway_pvals[mask].mean()

        df = pd.DataFrame({
            'pathway': pathway_names,
            'activity': cluster_pathway.values,
            'pval': cluster_pval.values
        })

        df_sig = df[df['pval'] < 0.05].copy()
        df_sig = df_sig.sort_values('activity', ascending=False).head(top_n)

        pathway_results[cluster] = df_sig

    return pathway_results
```

### Pathway Interpretation (DATA-DRIVEN, NO HARDCODING)

```python
def compute_pathway_statistics(adata, cluster_col='tier3_cluster'):
    """
    Compute statistical measures for Pathway activity across clusters.
    NO hardcoded pathway associations - purely data-driven.

    Returns:
    --------
    dict with global statistics and per-cluster z-scores
    """
    pw_key = 'score_mlm' if 'score_mlm' in adata.obsm else 'mlm_estimate'
    pval_key = 'padj_mlm' if 'padj_mlm' in adata.obsm else 'mlm_pvals'

    pw_scores = adata.obsm[pw_key]
    pw_pvals = adata.obsm[pval_key]

    # Global statistics
    global_stats = {
        'mean': pw_scores.mean(),
        'std': pw_scores.std(),
        'median': pw_scores.median()
    }

    # Per-cluster statistics with z-scores
    cluster_stats = {}
    for cluster in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == cluster
        cluster_pw = pw_scores[mask].mean()
        cluster_pval = pw_pvals[mask].mean()

        # Compute z-scores for each pathway
        z_scores = {}
        for pw in pw_scores.columns:
            if global_stats['std'][pw] > 1e-6:
                z = (cluster_pw[pw] - global_stats['mean'][pw]) / global_stats['std'][pw]
            else:
                z = 0.0
            z_scores[pw] = z

        cluster_stats[cluster] = {
            'pathway_activity': cluster_pw,
            'pathway_pvals': cluster_pval,
            'z_scores': z_scores
        }

    return {
        'global_stats': global_stats,
        'cluster_stats': cluster_stats
    }


def interpret_pathway_activity(pathway_df, cluster_stats=None, cluster_id=None, parent_context=None):
    """
    Interpret pathway activity for reasoning - DATA-DRIVEN approach.
    NO hardcoded associations - returns statistical analysis for LLM to interpret.

    Parameters:
    -----------
    pathway_df : DataFrame
        Top pathways with activity and p-values
    cluster_stats : dict (optional)
        Output from compute_pathway_statistics()
    cluster_id : str (optional)
        Cluster ID for z-score lookup
    parent_context : str (optional)
        Parent cell type context (e.g., "T_Effector")

    Returns:
    --------
    dict with statistical interpretation prompts for LLM reasoning
    """
    top_pathways = pathway_df['pathway'].tolist()[:10]

    # Collect pathway data with statistical context
    pw_data = []
    for _, row in pathway_df.iterrows():
        pathway = row['pathway']
        activity = row['activity']
        pval = row['pval']

        pw_info = {
            'pathway': pathway,
            'activity': activity,
            'pval': pval,
            'is_significant': pval < 0.05,
            'direction': 'HIGH' if activity > 0 else 'LOW',
        }

        # Add z-score if available
        if cluster_stats and cluster_id and cluster_id in cluster_stats['cluster_stats']:
            z_score = cluster_stats['cluster_stats'][cluster_id]['z_scores'].get(pathway, 0.0)
            pw_info['z_score'] = z_score
            pw_info['is_outlier'] = abs(z_score) > 2.5
        else:
            pw_info['z_score'] = None
            pw_info['is_outlier'] = False

        pw_data.append(pw_info)

    # Return structured data for LLM reasoning (NO interpretations)
    return {
        'top_pathways': top_pathways,
        'pathway_data': pw_data,
        'n_significant': sum(1 for p in pw_data if p['is_significant']),
        'n_outliers': sum(1 for p in pw_data if p['is_outlier']),
        'outlier_pathways': [p['pathway'] for p in pw_data if p['is_outlier']],
        'reasoning_prompt': generate_pathway_reasoning_prompt(pw_data, parent_context)
    }


def generate_pathway_reasoning_prompt(pw_data, parent_context=None):
    """
    Generate a reasoning prompt for LLM based on Pathway statistics.
    NO hardcoded interpretations - guides LLM to search literature.
    """
    prompt_parts = []
    context_str = f" in {parent_context}" if parent_context else ""

    # High activity pathways
    high_activity = [p for p in pw_data if p['activity'] > 1.0 and p['is_significant']]
    if high_activity:
        pw_names = ', '.join([p['pathway'] for p in high_activity[:5]])
        prompt_parts.append(
            f"High pathway activity: {pw_names}\n"
            f"→ Search literature: '[pathway_name]{context_str}' to identify functional state"
        )

    # Outlier pathways (statistically unusual)
    outliers = [p for p in pw_data if p['is_outlier']]
    if outliers:
        outlier_high = [p for p in outliers if p['direction'] == 'HIGH']
        outlier_low = [p for p in outliers if p['direction'] == 'LOW']

        if outlier_high:
            high_names = ', '.join([f"{p['pathway']} (z={p['z_score']:.2f})" for p in outlier_high])
            prompt_parts.append(
                f"⚠️ UNUSUALLY HIGH PATHWAYS: {high_names}\n"
                f"→ These pathways are statistically higher than other clusters\n"
                f"→ PRIORITY: Search '[pathway]{context_str}' for specialized functional states"
            )

        if outlier_low:
            low_names = ', '.join([f"{p['pathway']} (z={p['z_score']:.2f})" for p in outlier_low])
            prompt_parts.append(
                f"⚠️ UNUSUALLY LOW PATHWAYS: {low_names}\n"
                f"→ Suppression/quiescence of these pathways is unusual\n"
                f"→ Consider: exhausted, anergic, or quiescent states"
            )

    # Multiple high pathways (activation signature)
    if len(high_activity) >= 3:
        prompt_parts.append(
            "Multiple pathways highly active - suggests complex functional program\n"
            f"→ Search combinations: '[pathway1] [pathway2]{context_str}'"
        )

    # Low overall activity (potential exhaustion/quiescence)
    significant_activities = [p['activity'] for p in pw_data if p['is_significant']]
    if significant_activities and max(significant_activities) < 0.5:
        prompt_parts.append(
            "Low overall pathway activity across all pathways\n"
            f"→ Consider: exhausted, anergic, quiescent, or senescent states{context_str}\n"
            "→ Search: 'exhaustion', 'anergy', 'quiescence'"
        )

    return "\n\n".join(prompt_parts) if prompt_parts else f"Standard pathway pattern{context_str} - proceed with marker-based annotation"
```

---

## 4. Combined Evidence Collector (UPDATED - DATA-DRIVEN)

```python
def collect_functional_evidence(adata, cluster_col, cluster_id, tier=2, parent_context=None):
    """
    Collect all functional evidence for a cluster - DATA-DRIVEN approach.
    NO hardcoded interpretations - returns statistical analysis + reasoning prompts.

    Parameters:
    -----------
    tier : int
        2 = TF + Trajectory
        3 = Pathway
    parent_context : str (optional)
        Parent cell type context for Tier 3 (e.g., "T_Effector")

    Returns:
    --------
    dict with all evidence for reasoning (statistical + prompts, NO interpretations)
    """
    evidence = {
        'cluster_id': cluster_id,
        'tier': tier,
        'n_cells': (adata.obs[cluster_col] == cluster_id).sum()
    }

    if tier == 2:
        # TF Activity (check both v2 and legacy keys)
        if 'score_ulm' in adata.obsm or 'ulm_estimate' in adata.obsm:
            # Compute statistics
            tf_stats = compute_tf_statistics(adata, cluster_col)

            # Get top TFs for this cluster
            tf_results = get_tf_activity_per_cluster(adata, cluster_col, top_n=10)
            cluster_tf_df = tf_results.get(cluster_id, pd.DataFrame())

            # Interpret (NO hardcoded associations)
            tf_interp = interpret_tf_activity(
                cluster_tf_df,
                cluster_stats=tf_stats,
                cluster_id=cluster_id
            )

            evidence['tf_activity'] = {
                'data': cluster_tf_df.to_dict('records'),
                'interpretation': tf_interp,
                'has_outliers': tf_interp['n_outliers'] > 0,
                'outlier_tfs': tf_interp['outlier_tfs'],
                'reasoning_prompt': tf_interp['reasoning_prompt']
            }
        else:
            evidence['tf_activity'] = {'error': 'TF activity not computed!'}

        # Trajectory
        if 'pseudotime' in adata.obs.columns:
            pt_results = get_pseudotime_per_cluster(adata, cluster_col)
            evidence['trajectory'] = pt_results.get(cluster_id, {})
        else:
            evidence['trajectory'] = {'error': 'Trajectory not computed'}

        # Outliers from uns (if computed in tier2.md workflow)
        if 'tf_outliers' in adata.uns:
            cluster_outliers = [o for o in adata.uns['tf_outliers'] if o['cluster'] == str(cluster_id)]
            evidence['statistical_outliers'] = cluster_outliers
        else:
            evidence['statistical_outliers'] = []

        # Conflicts from uns (if computed in tier2.md workflow)
        if 'evidence_conflicts' in adata.uns:
            cluster_conflicts = [c for c in adata.uns['evidence_conflicts'] if c['cluster'] == str(cluster_id)]
            evidence['conflicts'] = cluster_conflicts
        else:
            evidence['conflicts'] = []

    elif tier == 3:
        # Pathway Activity (check both v2 and legacy keys)
        if 'score_mlm' in adata.obsm or 'mlm_estimate' in adata.obsm:
            # Compute statistics
            pw_stats = compute_pathway_statistics(adata, cluster_col)

            # Get top pathways for this cluster
            pw_results = get_pathway_activity_per_cluster(adata, cluster_col, top_n=10)
            cluster_pw_df = pw_results.get(cluster_id, pd.DataFrame())

            # Interpret (NO hardcoded associations)
            pw_interp = interpret_pathway_activity(
                cluster_pw_df,
                cluster_stats=pw_stats,
                cluster_id=cluster_id,
                parent_context=parent_context
            )

            evidence['pathway_activity'] = {
                'data': cluster_pw_df.to_dict('records'),
                'interpretation': pw_interp,
                'has_outliers': pw_interp['n_outliers'] > 0,
                'outlier_pathways': pw_interp['outlier_pathways'],
                'reasoning_prompt': pw_interp['reasoning_prompt']
            }
        else:
            evidence['pathway_activity'] = {'error': 'Pathway activity not computed!'}

        # Outliers from uns (if computed in tier3.md workflow)
        if 'pathway_outliers' in adata.uns:
            cluster_outliers = [o for o in adata.uns['pathway_outliers'] if o['cluster'] == str(cluster_id)]
            evidence['statistical_outliers'] = cluster_outliers
        else:
            evidence['statistical_outliers'] = []

        # Conflicts from uns (if computed in tier3.md workflow)
        if 'evidence_conflicts' in adata.uns:
            cluster_conflicts = [c for c in adata.uns['evidence_conflicts'] if c['cluster'] == str(cluster_id)]
            evidence['conflicts'] = cluster_conflicts
        else:
            evidence['conflicts'] = []

    return evidence
```

---

## Verification Function

```python
def verify_functional_analysis(adata, tier):
    """
    Verify that required functional analysis has been computed.

    Raises AssertionError if missing.
    """
    errors = []

    if tier == 2:
        # Check both decoupler v2 (score_ulm) and legacy (ulm_estimate) key names
        has_tf = 'score_ulm' in adata.obsm or 'ulm_estimate' in adata.obsm
        has_tf_p = 'padj_ulm' in adata.obsm or 'ulm_pvals' in adata.obsm
        if not has_tf:
            errors.append("TF activity not computed (expected score_ulm or ulm_estimate)")
        if not has_tf_p:
            errors.append("TF p-values not computed (expected padj_ulm or ulm_pvals)")
        if adata.n_obs >= 2000 and 'pseudotime' not in adata.obs.columns:
            errors.append("Trajectory (pseudotime) not computed for >= 2000 cells")

    elif tier == 3:
        # Check both decoupler v2 (score_mlm) and legacy (mlm_estimate) key names
        has_pw = 'score_mlm' in adata.obsm or 'mlm_estimate' in adata.obsm
        has_pw_p = 'padj_mlm' in adata.obsm or 'mlm_pvals' in adata.obsm
        if not has_pw:
            errors.append("Pathway activity not computed (expected score_mlm or mlm_estimate)")
        if not has_pw_p:
            errors.append("Pathway p-values not computed (expected padj_mlm or mlm_pvals)")

    if errors:
        raise AssertionError(
            f"❌ FUNCTIONAL ANALYSIS INCOMPLETE (Tier {tier}):\n" +
            "\n".join(f"  - {e}" for e in errors)
        )

    print(f"✅ Functional analysis verified for Tier {tier}")
    return True
```

---

## Usage Examples (v3 - Data-Driven)

### Example 1: Tier 2 with TF Outlier Detection

```python
# After computing TF activity in Tier 2 workflow

# Step 1: Compute statistics (NO hardcoding)
tf_stats = compute_tf_statistics(subset, cluster_col='tier2_cluster')

# Step 2: Collect evidence for a cluster
evidence = collect_functional_evidence(
    subset,
    cluster_col='tier2_cluster',
    cluster_id='3',
    tier=2
)

# Step 3: Present to LLM
print("=" * 60)
print(f"Cluster 3 Evidence")
print("=" * 60)
print("\n## TF Activity (Top 10)")
for tf_data in evidence['tf_activity']['interpretation']['tf_data']:
    outlier_flag = "⚠️ OUTLIER" if tf_data['is_outlier'] else ""
    print(f"  {tf_data['tf']}: activity={tf_data['activity']:.2f}, "
          f"z={tf_data.get('z_score', 0):.2f} {outlier_flag}")

print("\n## Reasoning Instruction")
print(evidence['tf_activity']['reasoning_prompt'])

# Output example:
# ⚠️ STATISTICAL OUTLIERS: TBX21 (z=3.8)
# → These TFs are unusually active compared to other clusters
# → PRIORITY: Search 'TBX21 B lineage' for specialized cell types
```

### Example 2: Tier 3 with Pathway Outlier Detection

```python
# After computing Pathway activity in Tier 3 workflow

# Step 1: Compute statistics (NO hardcoding)
pw_stats = compute_pathway_statistics(subset, cluster_col='tier3_cluster')

# Step 2: Collect evidence for a cluster
evidence = collect_functional_evidence(
    subset,
    cluster_col='tier3_cluster',
    cluster_id='2',
    tier=3,
    parent_context='T_Effector'
)

# Step 3: Present to LLM
print("=" * 60)
print(f"Cluster 2 Evidence (T_Effector context)")
print("=" * 60)
print("\n## Pathway Activity (Top 10)")
for pw_data in evidence['pathway_activity']['interpretation']['pathway_data']:
    outlier_flag = "⚠️ OUTLIER" if pw_data['is_outlier'] else ""
    print(f"  {pw_data['pathway']}: activity={pw_data['activity']:.2f}, "
          f"z={pw_data.get('z_score', 0):.2f} {outlier_flag}")

print("\n## Reasoning Instruction")
print(evidence['pathway_activity']['reasoning_prompt'])

# Output example:
# ⚠️ UNUSUALLY HIGH PATHWAYS: WNT (z=3.2)
# → These pathways are statistically higher than other clusters
# → PRIORITY: Search 'WNT T_Effector' for specialized functional states
# → (LLM searches and finds: Precursor exhausted T cells (Tpex))
```

### Example 3: LLM Reasoning Workflow

```python
# This is how LLM uses the reasoning prompts

# 1. Receive evidence with reasoning_prompt
evidence = collect_functional_evidence(...)

# 2. LLM sees reasoning_prompt:
# "⚠️ STATISTICAL OUTLIERS: TBX21 (z=3.8)
#  → PRIORITY: Search 'TBX21 B lineage'"

# 3. LLM automatically performs MCP search:
mcp_result = search_pubmed("TBX21 B cell age-associated")
# Returns: PMID:25006127 - "Age-associated B cells express T-bet"

# 4. LLM annotates based on literature:
annotation = "ABC"  # Age-associated B cells
confidence = "HIGH"  # Outlier explained by literature

# 5. If LLM tried to use standard annotation:
annotation = "Memory_B"
confidence = "LOW"  # PENALTY: TBX21 outlier unexplained
max_score = 7  # Reduced from 12
```

---

## Comparison: Old vs New

### OLD (v2 - Hardcoded)

```python
# ❌ Hardcoded dictionary
TF_ASSOCIATIONS = {
    'TBX21': 'Th1/Effector, cytotoxic',
    ...
}

# ❌ Always interprets TBX21 as Th1/Effector
if tf == 'TBX21' and activity > 0.5:
    suggested_states.append('Effector')

# ❌ Misses ABC (age-associated B cells)
# TBX21+ in B cells → incorrectly suggests "Effector" (wrong lineage!)
```

### NEW (v3 - Data-Driven)

```python
# ✅ No hardcoding - statistical analysis
z_score = (cluster_tf - global_mean) / global_std

# ✅ Detect unusual pattern
if abs(z_score) > 2.5:
    prompt = f"⚠️ OUTLIER: {tf} (z={z_score:.2f})\n"
    prompt += f"→ Search '{tf} [cell_lineage]'"

# ✅ LLM searches and finds correct answer
# "TBX21 B cell" → PMID:25006127 → ABC discovered!
```
