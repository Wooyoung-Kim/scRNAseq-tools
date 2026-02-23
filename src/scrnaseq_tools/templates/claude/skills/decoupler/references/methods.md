# Decoupler Enrichment Methods (dc.mt)

## Overview of Methods

| Method | Function | Best For | Speed |
|--------|----------|----------|-------|
| ULM | `dc.run_ulm()` | TF activity | Fast |
| MLM | `dc.run_mlm()` | Pathway activity | Fast |
| WSUM | `dc.run_wsum()` | General weighted sum | Fast |
| ORA | `dc.run_ora()` | Categorical enrichment | Fast |
| GSEA | `dc.run_gsea()` | Ranked enrichment | Medium |
| AUCELL | `dc.run_aucell()` | Single-cell scoring | Slow |
| GSVA | `dc.run_gsva()` | Gene set variation | Medium |
| VIPER | `dc.run_viper()` | Regulatory activity | Slow |
| Z-score | `dc.run_zscore()` | Simple scoring | Fast |

---

## ULM (Univariate Linear Model)

**Recommended for TF activity inference.**

```python
dc.run_ulm(
    mat,                    # AnnData, DataFrame, or array
    net,                    # Prior knowledge network
    source='source',        # Column with source nodes (TFs)
    target='target',        # Column with target nodes (genes)
    weight='weight',        # Column with edge weights
    batch_size=10000,       # Batch size for memory efficiency
    min_n=5,                # Minimum targets per source
    verbose=True,           # Print progress
    use_raw=False           # Use adata.raw
)

# Results in adata.obsm:
#   - 'ulm_estimate': Activity scores
#   - 'ulm_pvals': P-values
```

**How it works**: For each TF, fits a linear model predicting gene expression from the TF's target gene weights. The t-statistic estimates TF activity.

---

## MLM (Multivariate Linear Model)

**Recommended for pathway activity scoring.**

```python
dc.run_mlm(
    mat,                    # AnnData, DataFrame, or array
    net,                    # Prior knowledge network
    source='source',        # Column with source nodes (pathways)
    target='target',        # Column with target nodes (genes)
    weight='weight',        # Column with edge weights
    batch_size=10000,       # Batch size for memory efficiency
    min_n=5,                # Minimum targets per source
    verbose=True,           # Print progress
    use_raw=False           # Use adata.raw
)

# Results in adata.obsm:
#   - 'mlm_estimate': Activity scores
#   - 'mlm_pvals': P-values
```

**How it works**: Fits a multivariate linear model using all genes as predictors. More robust than ULM when pathways share many genes.

---

## WSUM (Weighted Sum)

**Simple weighted sum of target gene expression.**

```python
dc.run_wsum(
    mat,                    # AnnData, DataFrame, or array
    net,                    # Prior knowledge network
    source='source',        # Column with source nodes
    target='target',        # Column with target nodes
    weight='weight',        # Column with edge weights
    times=1000,             # Permutations for p-value
    batch_size=10000,       # Batch size
    min_n=5,                # Minimum targets per source
    seed=42,                # Random seed
    verbose=True            # Print progress
)

# Results in adata.obsm:
#   - 'wsum_estimate': Weighted sum scores
#   - 'wsum_pvals': Permutation p-values
#   - 'wsum_norm': Normalized scores (z-score)
```

---

## ORA (Over-Representation Analysis)

**For enrichment of gene lists.**

```python
dc.run_ora(
    mat,                    # AnnData, DataFrame, or list of genes
    net,                    # Gene set network
    source='source',        # Column with gene set names
    target='target',        # Column with gene names
    n_background=None,      # Background size (default: all genes)
    min_n=5,                # Minimum genes per set
    verbose=True            # Print progress
)

# Returns DataFrame with:
#   - source: Gene set name
#   - statistic: Fold enrichment
#   - p_value: Fisher's exact test p-value
#   - fdr: FDR-corrected p-value
```

**Example with ranked genes:**
```python
sc.tl.rank_genes_groups(adata, groupby='cell_type')
enr_results = dc.run_ora(mat=adata, net=net, source='source', target='target')
```

**Example with gene list:**
```python
gene_list = ['STAT1', 'IRF1', 'IRF7', 'ISG15', 'MX1']
enr_results = dc.run_ora(mat=gene_list, net=net, source='source', target='target')
```

---

## GSEA (Gene Set Enrichment Analysis)

**For ranked enrichment analysis.**

```python
dc.run_gsea(
    mat,                    # AnnData or DataFrame with ranking statistic
    net,                    # Gene set network
    source='source',        # Column with gene set names
    target='target',        # Column with gene names
    weight='weight',        # Column with edge weights (optional)
    times=1000,             # Permutations
    min_n=5,                # Minimum genes per set
    seed=42,                # Random seed
    verbose=True            # Print progress
)

# Returns DataFrame with:
#   - source: Gene set name
#   - statistic: Enrichment score (ES)
#   - p_value: Permutation p-value
#   - fdr: FDR-corrected p-value
#   - nes: Normalized enrichment score
```

**Example:**
```python
sc.tl.rank_genes_groups(adata, groupby='condition', method='wilcoxon')
dc.run_gsea(mat=adata, net=hallmark_net, source='source', target='target')
dc.pl.leading_edge(adata, net=hallmark_net, source='HALLMARK_INTERFERON_GAMMA_RESPONSE')
```

---

## AUCELL (Area Under Curve)

**For single-cell gene set scoring.**

```python
dc.run_aucell(
    mat,                    # AnnData or array
    net,                    # Gene set network
    source='source',        # Column with gene set names
    target='target',        # Column with gene names
    n_cells=None,           # Number of cells (for sampling)
    use_raw=False,          # Use raw counts
    verbose=True            # Print progress
)

# Results in adata.obsm:
#   - 'aucell_estimate': AUC scores per cell
```

---

## GSVA (Gene Set Variation Analysis)

**For gene set variation scoring.**

```python
dc.run_gsva(
    mat,                    # AnnData or array
    net,                    # Gene set network
    source='source',        # Column with gene set names
    target='target',        # Column with gene names
    kcdf='Gaussian',        # Kernel: 'Gaussian' or 'Poisson'
    mx_diff=True,           # Use max difference
    abs_ranking=False,      # Use absolute ranking
    min_n=5,                # Minimum genes per set
    verbose=True            # Print progress
)

# Results in adata.obsm:
#   - 'gsva_estimate': GSVA scores
```

---

## VIPER

**For regulatory activity inference.**

```python
dc.run_viper(
    mat,                    # AnnData or array
    net,                    # Regulatory network
    source='source',        # Column with regulators
    target='target',        # Column with targets
    weight='weight',        # Column with weights
    min_n=5,                # Minimum targets
    pleiotropy=True,        # Correct for pleiotropy
    eset_filter=False,      # Filter expression set
    verbose=True            # Print progress
)

# Results in adata.obsm:
#   - 'viper_estimate': Activity scores
#   - 'viper_pvals': P-values
```

---

## Z-score

**Simple z-score based scoring.**

```python
dc.run_zscore(
    mat,                    # AnnData or array
    net,                    # Gene set network
    source='source',        # Column with gene set names
    target='target',        # Column with gene names
    weight='weight',        # Column with weights (optional)
    min_n=5,                # Minimum genes per set
    verbose=True            # Print progress
)

# Results in adata.obsm:
#   - 'zscore_estimate': Z-scores
#   - 'zscore_pvals': P-values
```

---

## Decouple (Multiple Methods)

**Run multiple methods at once.**

```python
results = dc.decouple(
    mat,                    # AnnData or array
    net,                    # Prior knowledge network
    source='source',        # Column with source nodes
    target='target',        # Column with target nodes
    weight='weight',        # Column with edge weights
    methods=['ulm', 'mlm', 'wsum'],  # Methods to run
    consensus=True,         # Compute consensus score
    dense=True,             # Return dense arrays
    min_n=5,                # Minimum targets
    verbose=True            # Print progress
)

# Results stored in obsm with method name prefixes:
#   - 'ulm_estimate', 'ulm_pvals'
#   - 'mlm_estimate', 'mlm_pvals'
#   - 'consensus_estimate' (if consensus=True)
```

---

## Consensus

**Aggregate results from multiple methods.**

```python
dc.run_consensus(
    mat,                    # AnnData with method results in obsm
    obs_keys=None           # Keys to aggregate (default: all *_estimate)
)

# Results in adata.obsm:
#   - 'consensus_estimate': Aggregated scores (mean of normalized scores)
```

---

## Tool Functions (dc.tl)

### rankby_group

**Rank sources by activity differences between groups.**

```python
dc.tl.rankby_group(
    adata,                  # AnnData
    group_key,              # Column with groups
    obs_key,                # Key with activities
    reference='rest',       # Reference group or 'rest'
    method='t-test',        # Test method: 't-test', 'wilcoxon'
    corr_method='fdr_bh',   # Multiple testing correction
    ascending=False         # Sort order
)

# Results stored in adata.uns['rankby_group']
```

### rankby_order

**Rank sources along a continuous ordering (e.g., pseudotime).**

```python
dc.tl.rankby_order(
    adata,                  # AnnData
    order_key,              # Column with ordering (e.g., pseudotime)
    obs_key,                # Key with activities
    method='spearman',      # Correlation method: 'spearman', 'pearson'
    ascending=False         # Sort order
)

# Results stored in adata.uns['rankby_order']
```
