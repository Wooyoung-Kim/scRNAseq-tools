# Palantir Post-Results Functions (palantir.presults)

## select_branch_cells

**Select cells belonging to each branch/lineage.**

```python
masks = palantir.presults.select_branch_cells(
    ad,                      # AnnData
    pseudo_time_key='palantir_pseudotime',
    fate_prob_key='palantir_fate_probabilities',
    q=0.01,                  # Probability quantile threshold
    eps=0.01,                # Threshold adjustment
    masks_key='branch_masks',
    save_as_df=None
)

# Returns boolean mask array (cells x branches)
# Stores in adata.obsm['branch_masks']
```

**How it works:** Cells are selected based on:
1. Pseudotime ordering
2. Fate probability exceeding threshold

**Example:**
```python
palantir.presults.select_branch_cells(adata, q=0.01, eps=0.01)
palantir.plot.plot_branch_selection(adata)
```

---

## compute_gene_trends

**Compute gene expression trends along trajectories.**

```python
gene_trends = palantir.presults.compute_gene_trends(
    ad,                      # AnnData
    lineages=None,           # Subset of lineages (default: all)
    masks_key='branch_masks',
    expression_key=None,     # Expression layer (default: X or MAGIC)
    pseudo_time_key='palantir_pseudotime',
    gene_trend_key='gene_trends',
    save_as_df=None,
    **kwargs                 # mellon.FunctionEstimator args
)

# Returns DataFrame with gene trends per branch
# Stores in adata.uns['gene_trends']
```

**Example:**
```python
# Select branch cells first
palantir.presults.select_branch_cells(adata)

# Use MAGIC imputed expression
palantir.utils.run_magic_imputation(adata)

# Compute trends
gene_trends = palantir.presults.compute_gene_trends(
    adata,
    expression_key='MAGIC_imputed_data'
)

# Visualize
palantir.plot.plot_gene_trends(adata, genes=['CD34', 'GATA1', 'MPO'])
```

---

## compute_gene_trends_legacy

**Legacy GAM-based gene trend computation (requires pygam).**

```python
trends = palantir.presults.compute_gene_trends_legacy(
    data,                    # AnnData or PResults
    gene_exprs=None,         # Expression matrix (optional)
    lineages=None,           # Subset of lineages
    n_splines=4,             # Number of splines
    spline_order=2,          # Spline order
    n_jobs=-1,               # Parallel jobs
    expression_key='MAGIC_imputed_data',
    pseudo_time_key='palantir_pseudotime',
    fate_prob_key='palantir_fate_probabilities',
    gene_trend_key='palantir_gene_trends',
    save_as_df=None
)

# Returns dict: {branch: {'trends': DataFrame, 'std': DataFrame}}
```

---

## cluster_gene_trends

**Cluster genes by expression trend similarity.**

```python
clusters = palantir.presults.cluster_gene_trends(
    data,                    # AnnData or DataFrame
    branch_name,             # Target branch
    genes=None,              # Gene subset (default: all)
    gene_trend_key='gene_trends',
    n_neighbors=150,         # k-NN for clustering
    **kwargs                 # scanpy.tl.leiden args
)

# Returns pd.Series with cluster labels per gene
```

**Example:**
```python
# Cluster genes for erythroid lineage
clusters = palantir.presults.cluster_gene_trends(
    adata,
    branch_name='Erythroid',
    genes=adata.var_names[:500]
)

# Visualize
palantir.plot.plot_gene_trend_clusters(adata, 'Erythroid')
```

---

## gam_fit_predict

**Fit GAM for single gene (low-level, requires pygam).**

```python
predictions, std = palantir.presults.gam_fit_predict(
    x,                       # Pseudotime values
    y,                       # Gene expression values
    weights=None,            # Branch weights
    pred_x=None,             # Prediction points
    n_splines=4,
    spline_order=2
)

# Returns predictions and standard deviations
```

---

## Workflow for Gene Trends

```python
import palantir

# 1. Run Palantir (assume already done)
# palantir.core.run_palantir(adata, start_cell, ...)

# 2. MAGIC imputation for smooth expression
palantir.utils.run_magic_imputation(adata)

# 3. Select branch cells
palantir.presults.select_branch_cells(adata, q=0.01, eps=0.01)

# 4. Compute gene trends
palantir.presults.compute_gene_trends(
    adata,
    expression_key='MAGIC_imputed_data'
)

# 5. Visualize
genes = ['CD34', 'GATA1', 'MPO', 'IRF8']
palantir.plot.plot_gene_trends(adata, genes)
palantir.plot.plot_gene_trend_heatmaps(adata, genes)

# 6. Cluster trends (optional)
clusters = palantir.presults.cluster_gene_trends(adata, 'Erythroid')
palantir.plot.plot_gene_trend_clusters(adata, 'Erythroid')
```
