# Phase 4: Normalization & HVG Selection

```
╔══════════════════════════════════════════════════════════════════════╗
║  GOAL: "표현량이 비교 가능하고 중요한 유전자가 선택되었는가?"             ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Workflow

```
Step 1: Normalize counts (CPM + log)
    ↓
Step 2: Select Highly Variable Genes
    ↓
Step 3: Scale (optional, for PCA)
    ↓
Step 4: PCA
```

---

## Step 1: Normalize Counts

### Standard Normalization (Recommended)

```python
import scanpy as sc

# Save raw counts for later use (DE analysis, etc.)
adata.raw = adata.copy()

# Total count normalization (CPM-like)
sc.pp.normalize_total(
    adata,
    target_sum=1e4,  # Normalize to 10,000 counts per cell
    inplace=True
)

# Log transformation
sc.pp.log1p(adata)

# Verify
print(f"Normalization complete:")
print(f"  X max: {adata.X.max():.2f}")
print(f"  X min: {adata.X.min():.2f}")
print(f"  X mean: {adata.X.mean():.2f}")
```

### Alternative: scran Normalization (for highly variable library sizes)

```python
# For R-based scran normalization via rpy2
# Generally not necessary for standard 10x data

# import scanpy.external as sce
# sce.pp.normalize_scran(adata)
# sc.pp.log1p(adata)
```

---

## Step 2: Select Highly Variable Genes

### Option A: Seurat v3 Method (Recommended)

```python
# Identify HVGs
sc.pp.highly_variable_genes(
    adata,
    flavor='seurat_v3',  # Uses raw counts internally
    n_top_genes=3000,    # Number of HVGs to select
    span=0.3,
    batch_key='sample' if 'sample' in adata.obs.columns else None  # Batch-aware
)

# Summary
n_hvg = adata.var['highly_variable'].sum()
print(f"HVGs selected: {n_hvg}")

# Visualize
sc.pl.highly_variable_genes(adata, save='_hvg.png')
```

### Option B: Cell Ranger Method

```python
sc.pp.highly_variable_genes(
    adata,
    flavor='cell_ranger',
    n_top_genes=3000
)
```

### Check HVG Selection

```python
# HVG summary
print(f"\nHVG statistics:")
print(f"  Total genes: {adata.n_vars}")
print(f"  HVGs: {adata.var['highly_variable'].sum()}")
print(f"  Non-HVGs: {(~adata.var['highly_variable']).sum()}")

# Top HVGs by mean expression
hvg_df = adata.var[adata.var['highly_variable']].copy()
hvg_df = hvg_df.sort_values('means', ascending=False)
print(f"\nTop 20 HVGs by expression:")
print(hvg_df.head(20)[['means', 'dispersions', 'dispersions_norm']])
```

---

## Step 3: Scale Data (for PCA)

```python
# Scale to unit variance (mean=0, std=1)
# Only scale for PCA computation
sc.pp.scale(adata, max_value=10)

# Verify scaling
print(f"After scaling:")
print(f"  X mean: {adata.X.mean():.4f}")
print(f"  X std: {adata.X[:, adata.var['highly_variable']].std():.4f}")
```

**Note**: Scaling is done on all genes, but PCA will use only HVGs.

---

## Step 4: PCA

```python
# Run PCA on HVGs
sc.tl.pca(
    adata,
    n_comps=50,           # Number of components
    use_highly_variable=True,
    svd_solver='arpack'
)

# Check variance explained
print(f"\nPCA complete:")
print(f"  Components: {adata.obsm['X_pca'].shape[1]}")
print(f"  Variance explained (PC1-10):")
for i in range(10):
    print(f"    PC{i+1}: {adata.uns['pca']['variance_ratio'][i]:.3f}")

# Cumulative variance
cumvar = adata.uns['pca']['variance_ratio'].cumsum()
n_pcs_90 = (cumvar < 0.9).sum() + 1
print(f"\n  PCs for 90% variance: {n_pcs_90}")

# Elbow plot
sc.pl.pca_variance_ratio(adata, n_pcs=50, save='_pca_variance.png')
```

---

## Summary

```python
# Log normalization parameters
adata.uns['normalization'] = {
    'method': 'log1p_10k',
    'target_sum': 1e4,
    'n_hvg': int(adata.var['highly_variable'].sum()),
    'hvg_flavor': 'seurat_v3',
    'n_pcs': adata.obsm['X_pca'].shape[1]
}

print("\n✅ Normalization and HVG selection complete")
print(f"   Ready for batch integration (if needed)")
```

---

## Verification

```python
# Verify key outputs
assert adata.X.max() < 20, "Data not log-normalized"
assert 'highly_variable' in adata.var.columns, "HVGs not selected"
assert 'X_pca' in adata.obsm, "PCA not computed"

hvg_count = adata.var['highly_variable'].sum()
assert 1000 < hvg_count < 5000, f"Unusual HVG count: {hvg_count}"

print("✅ All normalization checks passed")
```

---

## Next Step

→ Read: phases/batch_integration.md
