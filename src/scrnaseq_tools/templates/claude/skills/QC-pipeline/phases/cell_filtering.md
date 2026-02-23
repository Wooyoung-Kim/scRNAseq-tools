# Phase 1: Cell Filtering

```
╔══════════════════════════════════════════════════════════════════════╗
║  GOAL: "이 세포가 viable하고 high-quality인가?"                        ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Workflow

```
Step 1: Calculate QC metrics
    ↓
Step 2: Visualize distributions
    ↓
Step 3: Determine thresholds
    ↓
Step 4: Apply filters
    ↓
Step 5: Verify and log
```

---

## Step 1: Calculate QC Metrics

```python
import scanpy as sc
import numpy as np

# Load data
adata = sc.read_10x_mtx('raw_data/')
# or
adata = sc.read_h5ad('raw_data.h5ad')

print(f"Raw data: {adata.n_obs} cells, {adata.n_vars} genes")

# Identify mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
print(f"Mitochondrial genes: {adata.var['mt'].sum()}")

# Identify ribosomal genes (optional)
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
print(f"Ribosomal genes: {adata.var['ribo'].sum()}")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt', 'ribo'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Check available metrics
print("QC metrics in obs:")
print([col for col in adata.obs.columns if 'counts' in col or 'genes' in col or 'pct' in col])
```

**Expected output in adata.obs:**
- `n_genes_by_counts`: Number of genes with positive counts
- `total_counts`: Total UMI counts
- `pct_counts_mt`: Percentage of mitochondrial counts
- `pct_counts_ribo`: Percentage of ribosomal counts

---

## Step 2: Visualize Distributions

```python
import matplotlib.pyplot as plt

# Violin plots
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_qc_violin.png'
)

# Scatter plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Counts vs Genes
sc.pl.scatter(
    adata, x='total_counts', y='n_genes_by_counts',
    color='pct_counts_mt', ax=axes[0], show=False
)
axes[0].axhline(200, color='r', linestyle='--', alpha=0.5)
axes[0].axhline(6000, color='r', linestyle='--', alpha=0.5)

# Counts vs Mito
sc.pl.scatter(
    adata, x='total_counts', y='pct_counts_mt',
    ax=axes[1], show=False
)
axes[1].axhline(20, color='r', linestyle='--', alpha=0.5)

# Genes vs Mito
sc.pl.scatter(
    adata, x='n_genes_by_counts', y='pct_counts_mt',
    ax=axes[2], show=False
)
axes[2].axhline(20, color='r', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('figures/qc_scatter.png', dpi=150, bbox_inches='tight')
plt.show()
```

---

## Step 3: Determine Thresholds

### Option A: Fixed Thresholds (Simple)

```python
# Default thresholds
min_genes = 200
max_genes = 6000
max_mito = 20
min_counts = 500
```

### Option B: MAD-based Thresholds (Recommended)

```python
def mad_threshold(series, nmads=5, direction='both'):
    """Calculate MAD-based threshold."""
    median = np.median(series)
    mad = np.median(np.abs(series - median))

    if direction == 'upper':
        return median + nmads * mad
    elif direction == 'lower':
        return median - nmads * mad
    else:
        return median - nmads * mad, median + nmads * mad

# Calculate thresholds
min_genes_mad, max_genes_mad = mad_threshold(adata.obs['n_genes_by_counts'])
max_mito_mad = mad_threshold(adata.obs['pct_counts_mt'], direction='upper')

print(f"MAD-based thresholds:")
print(f"  Genes: [{min_genes_mad:.0f}, {max_genes_mad:.0f}]")
print(f"  Mito: < {max_mito_mad:.1f}%")

# Use more conservative of fixed/MAD
min_genes = max(200, min_genes_mad)
max_genes = min(6000, max_genes_mad)
max_mito = min(20, max_mito_mad)
```

### Option C: Per-Sample Thresholds

```python
# For multi-sample data
if 'sample' in adata.obs.columns:
    for sample in adata.obs['sample'].unique():
        subset = adata[adata.obs['sample'] == sample]
        print(f"\n{sample}:")
        print(f"  Genes: [{subset.obs['n_genes_by_counts'].quantile(0.01):.0f}, "
              f"{subset.obs['n_genes_by_counts'].quantile(0.99):.0f}]")
        print(f"  Mito: median={subset.obs['pct_counts_mt'].median():.1f}%")
```

---

## Step 4: Apply Filters

```python
# Store original counts
n_cells_original = adata.n_obs

# Filter 1: Minimum genes
adata = adata[adata.obs['n_genes_by_counts'] > min_genes].copy()
n_after_min_genes = adata.n_obs
print(f"After min_genes > {min_genes}: {n_after_min_genes} cells "
      f"(-{n_cells_original - n_after_min_genes})")

# Filter 2: Maximum genes (potential doublets)
adata = adata[adata.obs['n_genes_by_counts'] < max_genes].copy()
n_after_max_genes = adata.n_obs
print(f"After max_genes < {max_genes}: {n_after_max_genes} cells "
      f"(-{n_after_min_genes - n_after_max_genes})")

# Filter 3: Mitochondrial percentage
adata = adata[adata.obs['pct_counts_mt'] < max_mito].copy()
n_after_mito = adata.n_obs
print(f"After mito < {max_mito}%: {n_after_mito} cells "
      f"(-{n_after_max_genes - n_after_mito})")

# Summary
print(f"\n=== Cell Filtering Summary ===")
print(f"Original cells: {n_cells_original}")
print(f"Final cells: {n_after_mito}")
print(f"Cells removed: {n_cells_original - n_after_mito} "
      f"({(n_cells_original - n_after_mito)/n_cells_original:.1%})")
```

---

## Step 5: Verify and Log

```python
# Post-filtering visualization
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_qc_violin_filtered.png'
)

# Log filtering parameters
adata.uns['qc_params'] = {
    'min_genes': min_genes,
    'max_genes': max_genes,
    'max_mito': max_mito,
    'cells_original': n_cells_original,
    'cells_filtered': adata.n_obs
}

# Check remaining distribution
print(f"\nPost-filtering QC metrics:")
print(f"  n_genes: [{adata.obs['n_genes_by_counts'].min():.0f}, "
      f"{adata.obs['n_genes_by_counts'].max():.0f}] "
      f"(median: {adata.obs['n_genes_by_counts'].median():.0f})")
print(f"  total_counts: [{adata.obs['total_counts'].min():.0f}, "
      f"{adata.obs['total_counts'].max():.0f}] "
      f"(median: {adata.obs['total_counts'].median():.0f})")
print(f"  pct_mito: [{adata.obs['pct_counts_mt'].min():.1f}, "
      f"{adata.obs['pct_counts_mt'].max():.1f}] "
      f"(median: {adata.obs['pct_counts_mt'].median():.1f}%)")
```

---

## Filtering Log Template

```markdown
# Cell Filtering Log

## Thresholds
| Metric | Threshold | Method |
|--------|-----------|--------|
| min_genes | {min_genes} | Fixed/MAD |
| max_genes | {max_genes} | Fixed/MAD |
| max_mito | {max_mito}% | Fixed/MAD |

## Results
| Step | Cells | Removed |
|------|-------|---------|
| Original | {n_original} | - |
| min_genes | {n_after_min} | {n_removed_min} |
| max_genes | {n_after_max} | {n_removed_max} |
| max_mito | {n_final} | {n_removed_mito} |
| **Total** | **{n_final}** | **{n_total_removed} ({pct}%)** |

## Per-Sample Breakdown
| Sample | Original | Final | % Removed |
|--------|----------|-------|-----------|
{per_sample_table}
```

---

## Next Step

→ Read: phases/gene_filtering.md
