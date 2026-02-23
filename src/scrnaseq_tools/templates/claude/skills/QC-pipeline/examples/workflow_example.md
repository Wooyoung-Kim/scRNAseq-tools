# Complete QC Workflow Example

```
╔══════════════════════════════════════════════════════════════════════╗
║  Full QC pipeline from raw data to annotation-ready                  ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Example: Multi-Sample PBMC Data

```python
import scanpy as sc
import scrublet as scr
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Create output directories
os.makedirs('qc_output/figures', exist_ok=True)
os.makedirs('qc_output/reports', exist_ok=True)

print("=" * 60)
print("QC PIPELINE START")
print("=" * 60)
```

---

## Step 1: Load Data

```python
# Load multiple 10x samples
samples = ['sample1', 'sample2', 'sample3']
adatas = []

for sample in samples:
    path = f'raw_data/{sample}/filtered_feature_bc_matrix'
    adata_sample = sc.read_10x_mtx(path)
    adata_sample.obs['sample'] = sample
    adatas.append(adata_sample)
    print(f"Loaded {sample}: {adata_sample.n_obs} cells, {adata_sample.n_vars} genes")

# Concatenate
adata = adatas[0].concatenate(adatas[1:], batch_key='sample', batch_categories=samples)
print(f"\nCombined: {adata.n_obs} cells, {adata.n_vars} genes")
```

---

## Step 2: Calculate QC Metrics

```python
# Identify gene categories
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.str.match('^HB[^P]')

print(f"Gene categories:")
print(f"  Mitochondrial: {adata.var['mt'].sum()}")
print(f"  Ribosomal: {adata.var['ribo'].sum()}")
print(f"  Hemoglobin: {adata.var['hb'].sum()}")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt', 'ribo', 'hb'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# QC summary
print(f"\nQC metrics summary:")
print(f"  Genes/cell: [{adata.obs['n_genes_by_counts'].min()}, "
      f"{adata.obs['n_genes_by_counts'].max()}] "
      f"(median: {adata.obs['n_genes_by_counts'].median():.0f})")
print(f"  Counts/cell: [{adata.obs['total_counts'].min():.0f}, "
      f"{adata.obs['total_counts'].max():.0f}] "
      f"(median: {adata.obs['total_counts'].median():.0f})")
print(f"  Mito %: [{adata.obs['pct_counts_mt'].min():.1f}, "
      f"{adata.obs['pct_counts_mt'].max():.1f}] "
      f"(median: {adata.obs['pct_counts_mt'].median():.1f})")
```

---

## Step 3: Visualize QC

```python
# Violin plots
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_qc_raw.png'
)

# Per-sample QC
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='sample',
    rotation=45,
    save='_qc_by_sample.png'
)

# Scatter plots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
              color='pct_counts_mt', ax=axes[0], show=False)
axes[0].axhline(200, color='r', linestyle='--', alpha=0.5)
axes[0].axhline(6000, color='r', linestyle='--', alpha=0.5)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[1], show=False)
axes[1].axhline(20, color='r', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('qc_output/figures/qc_scatter_raw.png', dpi=150, bbox_inches='tight')
plt.show()
```

---

## Step 4: Cell Filtering

```python
n_cells_original = adata.n_obs

# Filter 1: Minimum genes
min_genes = 200
adata = adata[adata.obs['n_genes_by_counts'] > min_genes].copy()
print(f"After min_genes > {min_genes}: {adata.n_obs} cells")

# Filter 2: Maximum genes (potential doublets)
max_genes = 6000
adata = adata[adata.obs['n_genes_by_counts'] < max_genes].copy()
print(f"After max_genes < {max_genes}: {adata.n_obs} cells")

# Filter 3: Mitochondrial percentage
max_mito = 20
adata = adata[adata.obs['pct_counts_mt'] < max_mito].copy()
print(f"After mito < {max_mito}%: {adata.n_obs} cells")

print(f"\nCell filtering summary:")
print(f"  Original: {n_cells_original}")
print(f"  Filtered: {adata.n_obs}")
print(f"  Removed: {n_cells_original - adata.n_obs} ({(n_cells_original - adata.n_obs)/n_cells_original:.1%})")
```

---

## Step 5: Gene Filtering

```python
n_genes_original = adata.n_vars

# Filter genes expressed in < 3 cells
min_cells = 3
sc.pp.filter_genes(adata, min_cells=min_cells)

print(f"\nGene filtering:")
print(f"  Original: {n_genes_original}")
print(f"  Filtered: {adata.n_vars}")
print(f"  Removed: {n_genes_original - adata.n_vars}")
```

---

## Step 6: Doublet Detection

```python
print("\n" + "=" * 40)
print("DOUBLET DETECTION")
print("=" * 40)

adata.obs['doublet_score'] = np.nan
adata.obs['predicted_doublet'] = False

for sample in adata.obs['sample'].unique():
    print(f"\nProcessing {sample}...")
    mask = adata.obs['sample'] == sample
    subset = adata[mask]

    # Expected doublet rate
    n_cells = subset.n_obs
    expected_rate = min(0.08, n_cells / 1000 * 0.008)

    # Run Scrublet
    scrub = scr.Scrublet(subset.X, expected_doublet_rate=expected_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # Store
    adata.obs.loc[mask, 'doublet_score'] = doublet_scores
    adata.obs.loc[mask, 'predicted_doublet'] = predicted_doublets

    n_doublets = predicted_doublets.sum()
    print(f"  Cells: {n_cells}, Doublets: {n_doublets} ({n_doublets/n_cells:.1%})")

# Total doublets
total_doublets = adata.obs['predicted_doublet'].sum()
print(f"\nTotal doublets: {total_doublets} ({total_doublets/adata.n_obs:.1%})")

# Remove doublets
n_before = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"After doublet removal: {adata.n_obs} cells (-{n_before - adata.n_obs})")
```

---

## Step 7: Normalization & HVG

```python
print("\n" + "=" * 40)
print("NORMALIZATION")
print("=" * 40)

# Save raw counts
adata.raw = adata.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print(f"Normalization complete:")
print(f"  X range: [{adata.X.min():.2f}, {adata.X.max():.2f}]")

# HVG selection
sc.pp.highly_variable_genes(
    adata,
    flavor='seurat_v3',
    n_top_genes=3000,
    batch_key='sample'
)

n_hvg = adata.var['highly_variable'].sum()
print(f"\nHVGs selected: {n_hvg}")

# Visualize HVGs
sc.pl.highly_variable_genes(adata, save='_hvg.png')
```

---

## Step 8: PCA

```python
# Scale
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

print(f"\nPCA complete:")
print(f"  Shape: {adata.obsm['X_pca'].shape}")
print(f"  Variance explained (first 10 PCs):")
cumvar = 0
for i in range(10):
    var = adata.uns['pca']['variance_ratio'][i]
    cumvar += var
    print(f"    PC{i+1}: {var:.3f} (cumulative: {cumvar:.3f})")

# Elbow plot
sc.pl.pca_variance_ratio(adata, n_pcs=50, save='_pca_elbow.png')
```

---

## Step 9: Batch Integration (Harmony)

```python
print("\n" + "=" * 40)
print("BATCH INTEGRATION")
print("=" * 40)

# Pre-integration UMAP (for comparison)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep='X_pca')
sc.tl.umap(adata)

# Save pre-integration UMAP
adata.obsm['X_umap_uncorrected'] = adata.obsm['X_umap'].copy()

# Run Harmony
sce.pp.harmony_integrate(
    adata,
    key='sample',
    basis='X_pca',
    adjusted_basis='X_pca_harmony'
)

print(f"Harmony complete: X_pca_harmony shape = {adata.obsm['X_pca_harmony'].shape}")

# Post-integration neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep='X_pca_harmony')
sc.tl.umap(adata)
adata.obsm['X_umap_harmony'] = adata.obsm['X_umap'].copy()

# Visualize comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Pre-integration
ax1 = axes[0]
scatter1 = ax1.scatter(
    adata.obsm['X_umap_uncorrected'][:, 0],
    adata.obsm['X_umap_uncorrected'][:, 1],
    c=pd.factorize(adata.obs['sample'])[0],
    cmap='tab10',
    s=1,
    alpha=0.5
)
ax1.set_title('Pre-integration')
ax1.set_xlabel('UMAP1')
ax1.set_ylabel('UMAP2')

# Post-integration
sc.pl.umap(adata, color='sample', ax=axes[1], show=False, title='Post-integration (Harmony)')

plt.tight_layout()
plt.savefig('qc_output/figures/batch_integration.png', dpi=150, bbox_inches='tight')
plt.show()
```

---

## Step 10: Final Verification

```python
print("\n" + "=" * 60)
print("FINAL QC VERIFICATION")
print("=" * 60)

errors = []
warnings = []

# Check all required components
if 'n_genes_by_counts' not in adata.obs.columns:
    errors.append("Missing n_genes_by_counts")
if 'pct_counts_mt' not in adata.obs.columns:
    errors.append("Missing pct_counts_mt")
if 'doublet_score' not in adata.obs.columns:
    errors.append("Missing doublet_score")
if adata.obs['pct_counts_mt'].max() > 25:
    errors.append(f"High mito cells: {adata.obs['pct_counts_mt'].max():.1f}%")
if 'highly_variable' not in adata.var.columns:
    errors.append("HVGs not selected")
if 'X_pca_harmony' not in adata.obsm:
    errors.append("Harmony not computed")
if 'X_umap' not in adata.obsm:
    errors.append("UMAP not computed")
if adata.X.max() > 20:
    errors.append(f"Data not normalized: max={adata.X.max():.1f}")
if adata.raw is None:
    errors.append("Raw counts not preserved")

# Report
if errors:
    print("❌ ERRORS:")
    for e in errors:
        print(f"   - {e}")
else:
    print("✅ All checks passed!")
    print(f"\nFinal dataset:")
    print(f"   Cells: {adata.n_obs:,}")
    print(f"   Genes: {adata.n_vars:,}")
    print(f"   HVGs: {adata.var['highly_variable'].sum():,}")
    print(f"   Samples: {adata.obs['sample'].nunique()}")

# Save
if not errors:
    adata.write('qc_output/qc_complete.h5ad')
    print(f"\n✅ Saved: qc_output/qc_complete.h5ad")
    print("✅ Ready for Annotation-agent")
```

---

## QC Report Generation

```python
# Generate summary report
report = f"""# QC Report

## Dataset Summary
- **Cells (final)**: {adata.n_obs:,}
- **Genes (final)**: {adata.n_vars:,}
- **Samples**: {adata.obs['sample'].nunique()}

## Filtering Summary
| Step | Cells |
|------|-------|
| Original | {n_cells_original:,} |
| Cell filtering | ... |
| Doublet removal | {adata.n_obs:,} |

## QC Thresholds
| Metric | Threshold |
|--------|-----------|
| min_genes | > {min_genes} |
| max_genes | < {max_genes} |
| max_mito | < {max_mito}% |

## Per-Sample Summary
| Sample | Cells | Median Genes | Median Mito |
|--------|-------|--------------|-------------|
"""

for sample in adata.obs['sample'].unique():
    subset = adata[adata.obs['sample'] == sample]
    report += f"| {sample} | {subset.n_obs:,} | {subset.obs['n_genes_by_counts'].median():.0f} | {subset.obs['pct_counts_mt'].median():.1f}% |\n"

report += f"""
## Next Steps
1. Load `qc_complete.h5ad`
2. Run Annotation-agent (Tier 1 → Tier 2 → Tier 3)
3. Save annotated data
"""

with open('qc_output/reports/qc_summary.md', 'w') as f:
    f.write(report)

print("Report saved: qc_output/reports/qc_summary.md")
```
