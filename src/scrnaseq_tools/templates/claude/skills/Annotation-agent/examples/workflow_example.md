# Workflow Example (v2)

> **WARNING: B lineage-specific example only.** Marker genes (DNTT, RAG1, PAX5, etc.), TF associations, pathway names, and cell type labels below are illustrative for B cell lineage. Do NOT apply these markers or labels to other lineages without verifying enrichment in the current dataset's DE results. PMIDs in this file are placeholders.

Complete example with functional analysis integrated.

---

## Tier 2 Example: B Lineage Annotation

### Step 1: Setup and Load Data

```python
import scanpy as sc
import decoupler as dc
import palantir
import pandas as pd
import numpy as np

# Load Tier 1 annotated data
adata = sc.read_h5ad('annotation_output/tier1_annotated.h5ad')

# Subset to B lineage
major_type = 'B_lineage'
subset = adata[adata.obs['tier1_annotation'] == major_type].copy()
print(f"B lineage cells: {subset.n_obs}")
```

### Step 2: Re-cluster

```python
# Use existing harmony embedding
sc.pp.neighbors(subset, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=50)
sc.tl.leiden(subset, resolution=0.4, key_added='tier2_cluster')
print(f"Clusters: {subset.obs['tier2_cluster'].nunique()}")
```

### Step 3: Compute DE

```python
sc.tl.rank_genes_groups(subset, groupby='tier2_cluster', method='wilcoxon')
de_df = sc.get.rank_genes_groups_df(subset, group=None)
print(f"⚠️ DE computed for subset: {major_type}")
```

### Step 4: TF Activity (MANDATORY)

```python
print("🔬 Computing TF activity...")

# Get CollecTRI network (TF-target regulatory network)
net = dc.op.collectri(organism='human')

# Run ULM
dc.run_ulm(
    mat=subset,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
    use_raw=False
)

# Verify
assert 'ulm_estimate' in subset.obsm
print(f"✅ TF activity: {subset.obsm['ulm_estimate'].shape[1]} TFs")
```

### Step 5: Trajectory (MANDATORY if >= 2000 cells)

```python
print(f"🔬 Computing trajectory ({subset.n_obs} cells)...")

if subset.n_obs >= 2000:
    # Run diffusion maps
    palantir.utils.run_diffusion_maps(subset, n_components=10, pca_key='X_pca')

    # Determine multiscale space
    palantir.utils.determine_multiscale_space(subset)

    # Find start cell (Pro-B markers)
    pro_b_markers = ['DNTT', 'RAG1', 'ERG']
    pro_b_score = np.zeros(subset.n_obs)
    for marker in pro_b_markers:
        if marker in subset.var_names:
            expr = subset[:, marker].X.toarray().flatten() if hasattr(subset.X, 'toarray') else subset[:, marker].X.flatten()
            pro_b_score += (expr > 0).astype(float)
    start_cell = subset.obs_names[np.argmax(pro_b_score)]
    print(f"   Start cell: {start_cell}")

    # Run Palantir (results stored directly in AnnData)
    palantir.core.run_palantir(subset, start_cell, num_waypoints=500)
    subset.obs['pseudotime'] = subset.obs['palantir_pseudotime']

    # Categorize
    subset.obs['pseudotime_category'] = pd.cut(
        subset.obs['pseudotime'],
        bins=[0, 0.33, 0.67, 1.0],
        labels=['Early', 'Mid', 'Late']
    )
    print(f"✅ Trajectory computed")
else:
    subset.obs['pseudotime'] = np.nan
    print("⚠️ Skipped trajectory (< 2000 cells)")
```

### Step 6: Verify Functional Analysis

```python
# MUST pass before annotation
errors = []
if 'ulm_estimate' not in subset.obsm:
    errors.append("TF activity not computed")
if subset.n_obs >= 2000 and 'pseudotime' not in subset.obs.columns:
    errors.append("Trajectory not computed")

if errors:
    raise AssertionError(f"❌ INCOMPLETE: {errors}")

print("✅ Functional analysis verified")
```

### Step 7: Collect Evidence for Cluster 0

```python
cluster_id = '0'

# DE Markers
cluster_de = de_df[de_df['group'] == cluster_id].copy()
valid_markers = cluster_de[
    (cluster_de['pct_nz_group'] >= 0.25) &
    (cluster_de['logfoldchanges'] >= 1.0) &
    (cluster_de['pvals_adj'] < 0.05)
].nlargest(20, 'logfoldchanges')

print("## Top 20 DE Markers")
print(valid_markers[['names', 'pct_nz_group', 'logfoldchanges', 'pvals_adj']].to_string())

# TF Activity
tf_act = subset.obsm['ulm_estimate']
tf_pval = subset.obsm['ulm_pvals']
mask = subset.obs['tier2_cluster'] == cluster_id

cluster_tf = pd.DataFrame({
    'tf': tf_act.columns,
    'activity': tf_act[mask].mean().values,
    'pval': tf_pval[mask].mean().values
})
cluster_tf = cluster_tf[cluster_tf['pval'] < 0.05].sort_values('activity', ascending=False).head(10)

print("\n## Top 10 TF Activity")
print(cluster_tf.to_string())

# Trajectory
if 'pseudotime' in subset.obs.columns:
    pt = subset.obs.loc[mask, 'pseudotime']
    print(f"\n## Trajectory Position")
    print(f"Mean pseudotime: {pt.mean():.3f}")
    print(f"Category: {'Early' if pt.mean() < 0.33 else 'Mid' if pt.mean() < 0.67 else 'Late'}")
```

### Step 8: Integrated Reasoning

```markdown
## Cluster 0 Evidence Summary

### DE Markers
- DNTT (96.2%, LFC=12.1) - Pro-B marker
- RAG1 (91.8%, LFC=6.8) - V(D)J recombination
- ERG (88.1%, LFC=8.2) - Early B TF

### TF Activity
- PAX5 (1.8, p<0.001) - B cell identity
- EBF1 (1.5, p<0.001) - Early B cell factor
- E2A (1.2, p<0.01) - B cell development

### Trajectory Position
- Mean pseudotime: 0.15
- Category: EARLY
- Interpretation: Progenitor state

### Evidence Integration
| Type | Finding | Supports |
|------|---------|----------|
| Markers | DNTT+RAG1+ERG | Pro_B |
| TF | PAX5+EBF1+E2A | Pro_B/Pre_B |
| Trajectory | EARLY (0.15) | Pro_B |
| **Consensus** | | **Pro_B** |

### Final Decision
- **Annotation**: Pro_B
- **Confidence**: HIGH (12/12)
  - Markers: 3 pts (3+ markers)
  - Refs: 3 pts (PMID:22508441, PMID:39179932)
  - TF: 3 pts (PAX5+EBF1 strong match)
  - Trajectory: 3 pts (EARLY matches)
```

### Step 9: Save with Functional Analysis

```python
# Assign annotations (example mapping)
tier2_mapping = {
    '0': 'Pro_B',
    '1': 'Pre_B',
    '2': 'Transitional_B',
    # ... etc
}
subset.obs['tier2_annotation'] = subset.obs['tier2_cluster'].map(tier2_mapping)

# Save
subset.write('annotation_output/subsets/tier2/B_lineage.h5ad')

# Verify saved data
saved = sc.read_h5ad('annotation_output/subsets/tier2/B_lineage.h5ad')
assert 'ulm_estimate' in saved.obsm, "TF activity not saved!"
assert 'pseudotime' in saved.obs.columns, "Pseudotime not saved!"
print("✅ Saved with functional analysis")
```

---

## Tier 3 Example: Pro_B Functional States

### Step 1: Load Tier 2 Subset

```python
tier2_data = sc.read_h5ad('annotation_output/subsets/tier2/B_lineage.h5ad')
subset = tier2_data[tier2_data.obs['tier2_annotation'] == 'Pro_B'].copy()
print(f"Pro_B cells: {subset.n_obs}")
```

### Step 2: Re-cluster

```python
sc.pp.neighbors(subset, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=50)
sc.tl.leiden(subset, resolution=0.8, key_added='tier3_cluster')
```

### Step 3: Re-compute DE

```python
sc.tl.rank_genes_groups(subset, groupby='tier3_cluster', method='wilcoxon')
```

### Step 4: Pathway Activity (MANDATORY)

```python
print("🔬 Computing Pathway activity...")

net = dc.op.progeny(organism='human', top=500)
dc.run_mlm(
    mat=subset,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
    use_raw=False
)

assert 'mlm_estimate' in subset.obsm
print(f"✅ Pathway activity: {subset.obsm['mlm_estimate'].shape[1]} pathways")
```

### Step 5: Collect Evidence

```python
cluster_id = '0'

# Pathway Activity
pw_act = subset.obsm['mlm_estimate']
pw_pval = subset.obsm['mlm_pvals']
mask = subset.obs['tier3_cluster'] == cluster_id

cluster_pw = pd.DataFrame({
    'pathway': pw_act.columns,
    'activity': pw_act[mask].mean().values,
    'pval': pw_pval[mask].mean().values
})
cluster_pw = cluster_pw[cluster_pw['pval'] < 0.05].sort_values('activity', ascending=False)

print("## Pathway Activity")
print(cluster_pw.to_string())
```

### Step 6: Integrated Reasoning

```markdown
## Pro_B Cluster 0 Evidence Summary

### DE Markers
- PCLAF (92%, LFC=2.5) - Cell cycle
- MKI67 (85%, LFC=2.2) - Proliferation
- TOP2A (80%, LFC=2.0) - DNA replication

### Pathway Activity
- JAK-STAT (1.8, p<0.001) - Cytokine signaling
- MAPK (1.5, p<0.01) - Proliferation
- PI3K (1.2, p<0.01) - Survival/Growth

### Evidence Integration
| Type | Finding | Supports |
|------|---------|----------|
| Markers | MKI67+TOP2A+PCLAF | Proliferating |
| Pathway | JAK-STAT+MAPK+PI3K | Proliferating |
| **Consensus** | | **Proliferating** |

### Final Decision
- **Annotation**: Proliferating
- **Full Label**: B_lineage_Pro_B_Proliferating
- **Confidence**: HIGH (11/12)
```

---

## Output Verification

### Check Saved Data Contains Functional Analysis

```python
# Tier 2
tier2 = sc.read_h5ad('annotation_output/subsets/tier2/B_lineage.h5ad')
print("Tier 2 obsm keys:", list(tier2.obsm.keys()))
# Expected: ['ulm_estimate', 'ulm_pvals', 'X_pca', ...]

# Tier 3
tier3 = sc.read_h5ad('annotation_output/subsets/tier3/B_lineage_Pro_B.h5ad')
print("Tier 3 obsm keys:", list(tier3.obsm.keys()))
# Expected: ['mlm_estimate', 'mlm_pvals', 'X_pca', ...]
```

### Summary Table

| Tier | Required in obsm | Required in obs |
|------|------------------|-----------------|
| 2 | ulm_estimate, ulm_pvals | tier2_annotation, pseudotime (if >= 2000) |
| 3 | mlm_estimate, mlm_pvals | tier3_annotation, final_annotation |
