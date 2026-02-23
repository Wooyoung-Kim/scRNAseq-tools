# Decoupler Workflow Examples

## Single-Cell TF Activity Analysis

```python
import decoupler as dc
import scanpy as sc
import pandas as pd
import numpy as np

# ============================================
# 1. Load and Preprocess Data
# ============================================
adata = sc.read_h5ad('single_cell_data.h5ad')

# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Compute embeddings
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# ============================================
# 2. Load TF Network
# ============================================
tf_net = dc.op.collectri(organism='human')
print(f"Loaded {tf_net['source'].nunique()} TFs with {len(tf_net)} interactions")

# ============================================
# 3. Run TF Activity Inference
# ============================================
dc.run_ulm(
    mat=adata,
    net=tf_net,
    source='source',
    target='target',
    weight='weight',
    verbose=True
)

# ============================================
# 4. Extract and Visualize
# ============================================
dc.pp.get_obsm(adata, obsm_key='ulm_estimate')
sc.pl.umap(adata, color=['STAT1', 'NF-kB', 'MYC', 'leiden'], ncols=2)

# Dotplot by cluster
dc.pl.dotplot(
    adata,
    obs_keys='ulm_estimate',
    var_names=['STAT1', 'IRF1', 'NF-kB', 'MYC', 'E2F1', 'TP53', 'HIF1A'],
    groupby='leiden'
)

# ============================================
# 5. Find Cluster-Specific TFs
# ============================================
dc.tl.rankby_group(
    adata,
    group_key='leiden',
    obs_key='ulm_estimate',
    reference='rest'
)

results = adata.uns['rankby_group']
for cluster in adata.obs['leiden'].unique():
    print(f"\nCluster {cluster} top TFs:")
    top = results[results['group'] == cluster].nsmallest(5, 'pvalue')
    print(top[['source', 'logfoldchange', 'pvalue']].to_string(index=False))
```

---

## Pseudobulk Differential Activity Analysis

```python
import decoupler as dc
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
import numpy as np

# ============================================
# 1. Load Single-Cell Data
# ============================================
adata = sc.read_h5ad('single_cell_data.h5ad')

# Ensure counts are available
if 'counts' not in adata.layers:
    adata.layers['counts'] = adata.raw.X.copy()

# ============================================
# 2. Create Pseudobulk
# ============================================
pdata = dc.pp.pseudobulk(
    adata,
    sample_col='sample',
    groups_col='cell_type',
    obs_cols=['condition', 'batch'],
    layer='counts',
    min_cells=10,
    min_counts=1000
)

# ============================================
# 3. Filter Genes
# ============================================
genes = dc.pp.filter_by_expr(pdata, group='condition', min_count=10, min_total_count=15)
pdata = pdata[:, genes].copy()

# ============================================
# 4. Run DESeq2 per Cell Type
# ============================================
results_dict = {}

for cell_type in pdata.obs['cell_type'].unique():
    mask = pdata.obs['cell_type'] == cell_type
    pdata_ct = pdata[mask].copy()

    n_ctrl = (pdata_ct.obs['condition'] == 'control').sum()
    n_disease = (pdata_ct.obs['condition'] == 'disease').sum()

    if n_ctrl < 2 or n_disease < 2:
        continue

    counts = pd.DataFrame(pdata_ct.X.T, index=pdata_ct.var_names, columns=pdata_ct.obs_names)
    dds = DeseqDataSet(counts=counts, metadata=pdata_ct.obs, design_factors='condition')
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=['condition', 'disease', 'control'])
    stat_res.summary()
    results_dict[cell_type] = stat_res.results_df

# ============================================
# 5. TF Enrichment on DE Genes
# ============================================
tf_net = dc.op.collectri(organism='human')

for cell_type, de_results in results_dict.items():
    de_results['stat'] = -np.log10(de_results['pvalue']) * np.sign(de_results['log2FoldChange'])
    de_results = de_results.dropna(subset=['stat'])

    estimate, pvals = dc.run_ulm(
        mat=de_results[['stat']].T,
        net=tf_net,
        source='source',
        target='target',
        weight='weight',
        verbose=False
    )

    tf_res = pd.DataFrame({
        'tf': estimate.columns,
        'activity': estimate.values[0],
        'pvalue': pvals.values[0]
    }).sort_values('pvalue')

    print(f"\n{cell_type} top TFs: {', '.join(tf_res.head(5)['tf'].values)}")
```

---

## Pseudotime Trajectory Analysis

```python
import decoupler as dc
import scanpy as sc

# ============================================
# 1. Load Data with Trajectory
# ============================================
adata = sc.read_h5ad('trajectory_data.h5ad')

# Assume pseudotime computed (e.g., from Palantir)
# adata.obs['dpt_pseudotime'] exists

# ============================================
# 2. Run TF Activity Inference
# ============================================
tf_net = dc.op.collectri(organism='human')
dc.run_ulm(mat=adata, net=tf_net, source='source', target='target', weight='weight')
dc.pp.get_obsm(adata, obsm_key='ulm_estimate')

# ============================================
# 3. Find TFs Correlated with Pseudotime
# ============================================
dc.tl.rankby_order(
    adata,
    order_key='dpt_pseudotime',
    obs_key='ulm_estimate',
    method='spearman'
)

results = adata.uns['rankby_order']
print("TFs positively correlated (late):")
print(results[results['correlation'] > 0].nsmallest(10, 'pvalue'))

print("\nTFs negatively correlated (early):")
print(results[results['correlation'] < 0].nsmallest(10, 'pvalue'))
```

---

## Cell Type Scoring with Custom Markers

```python
import decoupler as dc
import scanpy as sc
import pandas as pd

# Custom markers
custom_markers = pd.DataFrame([
    {'source': 'T_cells', 'target': 'CD3D', 'weight': 1},
    {'source': 'T_cells', 'target': 'CD3E', 'weight': 1},
    {'source': 'B_cells', 'target': 'CD19', 'weight': 1},
    {'source': 'B_cells', 'target': 'MS4A1', 'weight': 1},
    {'source': 'Monocytes', 'target': 'CD14', 'weight': 1},
    {'source': 'Monocytes', 'target': 'LYZ', 'weight': 1},
])

# Score cell types
dc.run_aucell(mat=adata, net=custom_markers, source='source', target='target')
dc.pp.get_obsm(adata, obsm_key='aucell_estimate')

# Visualize
sc.pl.umap(adata, color=['T_cells', 'B_cells', 'Monocytes'])

# Assign cell types based on highest score
cell_types = adata.obsm['aucell_estimate'].idxmax(axis=1)
adata.obs['predicted_cell_type'] = cell_types
```

---

## Best Practices

### Data Preparation
1. **Use log-normalized data** for ULM, MLM, GSEA methods
2. **Keep raw counts** in a layer for pseudobulk analysis
3. **Check gene name format** - match your data to the network

### Method Selection
| Data Type | Recommended Methods |
|-----------|---------------------|
| Single-cell TF | ULM, VIPER |
| Single-cell pathway | MLM, AUCell |
| Bulk TF | ULM, VIPER |
| Bulk pathway | MLM, GSVA |
| Gene list enrichment | ORA |
| Ranked enrichment | GSEA |

### Statistical Considerations
1. **Multiple testing correction**: Use FDR-corrected p-values
2. **Effect sizes**: Report both p-values and activity scores
3. **Replication**: Validate with consensus scores across methods

---

## Common Pitfalls

### 1. Wrong Input Data Type
```python
# WRONG: Using counts for ULM
dc.run_ulm(mat=adata, net=net, ...)  # adata.X has raw counts

# CORRECT: Use log-normalized data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
dc.run_ulm(mat=adata, net=net, ...)
```

### 2. Gene Name Mismatch
```python
# Check overlap
data_genes = set(adata.var_names)
net_genes = set(net['target'])
overlap = data_genes & net_genes
print(f"Overlap: {len(overlap)} / {len(net_genes)} genes")
```

### 3. Insufficient Targets
```python
# Check targets per source
targets_per_source = net.groupby('source')['target'].nunique()
print(f"Sources with <5 targets: {(targets_per_source < 5).sum()}")

# Prune network
net = dc.pp.prune(net, min_n=5)
```

### 4. Not Filtering Before DESeq2
```python
# WRONG
dds = DeseqDataSet(counts=pdata.to_df().T, ...)

# CORRECT
genes = dc.pp.filter_by_expr(pdata, group='condition')
pdata = pdata[:, genes].copy()
dds = DeseqDataSet(counts=pdata.to_df().T, ...)
```

---

## Integration with Other Tools

### With Palantir (Trajectory Analysis)
```python
import palantir
import decoupler as dc

# Run Palantir
palantir.core.run_palantir(adata, start_cell, ...)

# Run decoupler
dc.run_ulm(adata, net=tf_net, ...)

# Correlate with pseudotime
dc.tl.rankby_order(adata, order_key='palantir_pseudotime', obs_key='ulm_estimate')
```

### With PyDESeq2
```python
# See Pseudobulk workflow above
```

### With scvi-tools
```python
import scvi

# scVI normalized expression
scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata)
model.train()
adata.layers['scvi_normalized'] = model.get_normalized_expression()

# Use for decoupler
dc.pp.swap_layer(adata, layer='scvi_normalized')
dc.run_ulm(mat=adata, net=tf_net, ...)
```
