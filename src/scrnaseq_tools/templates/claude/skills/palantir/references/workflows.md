# Palantir Workflow Examples

## Basic Hematopoiesis Analysis

```python
import palantir
import scanpy as sc
import pandas as pd
import numpy as np

# ============================================
# 1. Load and Preprocess Data
# ============================================
adata = sc.read_h5ad('bone_marrow_data.h5ad')

sc.pp.normalize_per_cell(adata)
palantir.preprocess.log_transform(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor="cell_ranger")
sc.pp.pca(adata)

# ============================================
# 2. Diffusion Maps
# ============================================
palantir.utils.run_diffusion_maps(adata, n_components=5)
palantir.utils.determine_multiscale_space(adata)

# ============================================
# 3. Visualization
# ============================================
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color='cell_type')

# ============================================
# 4. Define Start and Terminal States
# ============================================
start_cell = palantir.utils.early_cell(adata, celltype='HSC', celltype_column='cell_type')
terminal_states = palantir.utils.find_terminal_states(
    adata,
    celltypes=['Erythrocyte', 'Monocyte', 'DC'],
    celltype_column='cell_type'
)

# ============================================
# 5. Run Palantir
# ============================================
palantir.core.run_palantir(adata, start_cell, terminal_states=terminal_states, num_waypoints=500)

# ============================================
# 6. Visualize Results
# ============================================
palantir.plot.plot_palantir_results(adata, s=3)
sc.pl.umap(adata, color=['palantir_pseudotime', 'palantir_entropy'], ncols=2)

# ============================================
# 7. Gene Trends
# ============================================
palantir.utils.run_magic_imputation(adata)
palantir.presults.select_branch_cells(adata, q=0.01, eps=0.01)
palantir.presults.compute_gene_trends(adata, expression_key='MAGIC_imputed_data')

genes = ['CD34', 'GATA1', 'MPO', 'IRF8']
palantir.plot.plot_gene_trends(adata, genes)
palantir.plot.plot_gene_trend_heatmaps(adata, genes)

# ============================================
# 8. Save
# ============================================
adata.write('hematopoiesis_palantir.h5ad')
```

---

## Manual Start/Terminal Selection

```python
import palantir
import pandas as pd
import numpy as np

# Find start cell from marker expression
cd34_expr = adata[:, 'CD34'].X.toarray().flatten()
start_candidates = adata.obs_names[cd34_expr > np.percentile(cd34_expr, 95)]

# Select one at diffusion map extreme
dc1 = adata.obsm['DM_EigenVectors_multiscaled'][:, 0]
dc1_series = pd.Series(dc1, index=adata.obs_names)
start_cell = dc1_series.loc[start_candidates].idxmax()

# Manual terminal selection based on markers
terminal_markers = {'Erythroid': 'GATA1', 'Myeloid': 'MPO', 'Lymphoid': 'CD79A'}

terminal_states = {}
for lineage, marker in terminal_markers.items():
    marker_expr = adata[:, marker].X.toarray().flatten()
    top_cells = adata.obs_names[marker_expr > np.percentile(marker_expr, 99)]
    dc_vals = dc1_series.loc[top_cells]
    terminal_states[dc_vals.idxmin()] = lineage

terminal_states = pd.Series(terminal_states)

# Run Palantir
palantir.core.run_palantir(adata, start_cell, terminal_states=terminal_states)
```

---

## Integration with Decoupler (TF Activity)

```python
import palantir
import decoupler as dc
from scipy import stats

# Run Palantir (assume done)
# ...

# Get TF activities
tf_net = dc.op.collectri(organism='human')
dc.run_ulm(adata, net=tf_net, source='source', target='target', weight='weight')
dc.pp.get_obsm(adata, obsm_key='ulm_estimate')

# Correlate TF activities with pseudotime
pseudotime = adata.obs['palantir_pseudotime']
tf_activities = adata.obsm['ulm_estimate']

correlations = {}
for tf in tf_activities.columns:
    corr, pval = stats.spearmanr(pseudotime, tf_activities[tf])
    correlations[tf] = {'correlation': corr, 'pvalue': pval}

corr_df = pd.DataFrame(correlations).T.sort_values('pvalue')
print("Early TFs (negative correlation):")
print(corr_df[corr_df['correlation'] < 0].head(10))
print("\nLate TFs (positive correlation):")
print(corr_df[corr_df['correlation'] > 0].head(10))
```

---

## Integration with CellRank

```python
import palantir
import cellrank as cr
from cellrank.kernels import PseudotimeKernel, ConnectivityKernel
from cellrank.estimators import GPCCA

# Run Palantir
palantir.core.run_palantir(adata, start_cell, num_waypoints=500)

# Use Palantir pseudotime in CellRank
pk = PseudotimeKernel(adata, time_key='palantir_pseudotime')
pk.compute_transition_matrix(threshold_scheme='soft')

ck = ConnectivityKernel(adata)
ck.compute_transition_matrix()

combined_kernel = 0.5 * ck + 0.5 * pk

# CellRank analysis
g = GPCCA(combined_kernel)
g.compute_schur()
g.compute_macrostates(n_states=5)
g.set_terminal_states_from_macrostates()
g.compute_fate_probabilities()

# Compare fate probabilities
palantir_fates = adata.obsm['palantir_fate_probabilities']
cellrank_fates = adata.obsm['lineages_fwd']
```

---

## Best Practices

### Start Cell Selection
1. Use biological knowledge (progenitor markers)
2. Verify with diffusion components (should be at extreme)
3. Test multiple starts for robustness

```python
# Verify start cell position
dc1 = adata.obsm['DM_EigenVectors_multiscaled'][:, 0]
start_dc1 = dc1[adata.obs_names.get_loc(start_cell)]
print(f"Start cell DC1 percentile: {stats.percentileofscore(dc1, start_dc1):.1f}%")
```

### Terminal State Selection
1. Use marker genes for terminal types
2. Let Palantir auto-detect if uncertain
3. Validate terminals correspond to known mature types

### Gene Trend Analysis
1. Always use MAGIC imputation (smooths noise)
2. Select branch cells before computing trends
3. Focus on variable genes for meaningful trends

---

## Common Pitfalls

### 1. Wrong Start Cell
```python
# WRONG
start_cell = adata.obs_names[0]  # Random cell

# CORRECT
start_cell = palantir.utils.early_cell(adata, 'progenitor', 'cell_type')
```

### 2. Missing Terminal States
```python
# WRONG - missing lineages
terminal_states = pd.Series(['Ery'], index=['cell1'])

# CORRECT - all lineages
terminal_states = pd.Series(['Ery', 'Mono', 'DC'], index=['c1', 'c2', 'c3'])
```

### 3. Using Raw Counts for Gene Trends
```python
# WRONG
gene_trends = palantir.presults.compute_gene_trends(adata)

# CORRECT
palantir.utils.run_magic_imputation(adata)
gene_trends = palantir.presults.compute_gene_trends(
    adata, expression_key='MAGIC_imputed_data'
)
```

### 4. Forgetting Branch Selection
```python
# WRONG
gene_trends = palantir.presults.compute_gene_trends(adata)

# CORRECT
palantir.presults.select_branch_cells(adata, q=0.01)
gene_trends = palantir.presults.compute_gene_trends(adata)
```

### 5. Too Few Components/Waypoints
```python
# WRONG
palantir.utils.run_diffusion_maps(adata, n_components=2)
palantir.core.run_palantir(adata, start_cell, num_waypoints=100)

# CORRECT
palantir.utils.run_diffusion_maps(adata, n_components=10)
palantir.core.run_palantir(adata, start_cell, num_waypoints=500)
```

---

## Parameter Tuning Guide

| Data Size | n_components | knn | num_waypoints |
|-----------|--------------|-----|---------------|
| < 5,000 | 5-10 | 30 | 500 |
| 5,000-20,000 | 10-15 | 30-50 | 1000 |
| > 20,000 | 10-20 | 50+ | 1200 |

### Checking Eigenvalue Gap
```python
import matplotlib.pyplot as plt

eigenvalues = adata.uns['DM_EigenValues']
plt.plot(eigenvalues.values, 'o-')
plt.xlabel('Component')
plt.ylabel('Eigenvalue')
plt.title('Eigenvalue Spectrum')
# Look for gap to determine n_components
```
