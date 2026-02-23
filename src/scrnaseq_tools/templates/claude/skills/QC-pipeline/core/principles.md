# QC Principles and Thresholds

```
╔══════════════════════════════════════════════════════════════════════╗
║  데이터의 품질이 분석의 품질을 결정한다                                ║
║  "Garbage in, garbage out"                                           ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Core Principles

### 1. Context-Aware Thresholds

**절대적 기준 없음** - tissue/protocol에 따라 다름

| Tissue | Expected mito % | Expected genes/cell |
|--------|-----------------|---------------------|
| PBMC | 5-10% | 500-2000 |
| Solid tumor | 10-20% | 1000-5000 |
| Brain | 5-15% | 1000-4000 |
| Kidney | 15-25% | 1000-3000 |

### 2. Data-Driven Filtering

```python
# BAD: Arbitrary thresholds
adata = adata[adata.obs['pct_counts_mt'] < 10].copy()  # 왜 10%?

# GOOD: MAD-based thresholds
def mad_filter(adata, metric, nmads=5):
    """Filter cells based on median absolute deviation."""
    values = adata.obs[metric]
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    upper = median + nmads * mad
    lower = median - nmads * mad
    return (values > lower) & (values < upper)
```

### 3. Preserve Information

```
❌ 과도한 필터링: 세포 50% 손실
✅ 적절한 필터링: 세포 10-20% 손실
⚠️ 불충분한 필터링: 노이즈 포함
```

---

## Default Thresholds

### Cell Filtering

| Metric | Default | Flexible Range | Rationale |
|--------|---------|----------------|-----------|
| `n_genes_by_counts` | > 200 | 200-500 | Empty droplets |
| `n_genes_by_counts` | < 6000 | 4000-8000 | Doublets |
| `total_counts` | > 500 | 500-1000 | Low quality |
| `total_counts` | < 50000 | 30000-100000 | Doublets |
| `pct_counts_mt` | < 20% | 10-25% | Dying cells |

### Gene Filtering

| Metric | Default | Rationale |
|--------|---------|-----------|
| `min_cells` | 3 | Noise genes |
| `min_counts` | 3 | Low expression |

### Doublet Detection

| Metric | Default | Rationale |
|--------|---------|-----------|
| `expected_doublet_rate` | 0.06 | ~6% for 10k cells |
| `threshold` | 0.25 | Scrublet default |

---

## Threshold Selection Strategy

### Step 1: Visualize First

```python
# ALWAYS visualize before filtering
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
              color='pct_counts_mt')
```

### Step 2: Identify Outliers

```python
# Look for natural breakpoints
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Histogram of genes
axes[0].hist(adata.obs['n_genes_by_counts'], bins=50)
axes[0].axvline(200, color='r', linestyle='--', label='min_genes=200')
axes[0].set_xlabel('Genes per cell')

# Histogram of counts
axes[1].hist(adata.obs['total_counts'], bins=50)
axes[1].set_xlabel('Counts per cell')

# Histogram of mito
axes[2].hist(adata.obs['pct_counts_mt'], bins=50)
axes[2].axvline(20, color='r', linestyle='--', label='max_mito=20')
axes[2].set_xlabel('Mito %')

plt.tight_layout()
```

### Step 3: Apply Filters with Logging

```python
print(f"Cells before filtering: {adata.n_obs}")

# Apply filters
adata = adata[adata.obs['n_genes_by_counts'] > 200].copy()
print(f"  After min_genes > 200: {adata.n_obs}")

adata = adata[adata.obs['n_genes_by_counts'] < 6000].copy()
print(f"  After max_genes < 6000: {adata.n_obs}")

adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
print(f"  After mito < 20%: {adata.n_obs}")

print(f"Cells after filtering: {adata.n_obs}")
print(f"Cells removed: {original_cells - adata.n_obs} ({(original_cells - adata.n_obs)/original_cells:.1%})")
```

---

## Sample-Specific Considerations

### Multi-Sample Data

```python
# Check QC per sample
for sample in adata.obs['sample'].unique():
    subset = adata[adata.obs['sample'] == sample]
    print(f"\n{sample}:")
    print(f"  Cells: {subset.n_obs}")
    print(f"  Median genes: {subset.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Median mito: {subset.obs['pct_counts_mt'].median():.1f}%")
```

### Batch-Aware Filtering

```python
# If one sample has much higher mito, consider:
# 1. Different tissue handling → adjust threshold for that sample
# 2. Technical issue → flag but don't adjust
# 3. Biology (tumor vs normal) → document and proceed
```

---

## QC Report Template

```markdown
# QC Report: {dataset_name}

## Summary
- Total cells loaded: {n_cells_raw}
- Cells after filtering: {n_cells_filtered}
- Cells removed: {n_cells_removed} ({pct_removed}%)

## Filtering Thresholds
| Metric | Threshold | Cells Removed |
|--------|-----------|---------------|
| min_genes | > 200 | {n_low_genes} |
| max_genes | < 6000 | {n_high_genes} |
| max_mito | < 20% | {n_high_mito} |
| doublets | score > 0.25 | {n_doublets} |

## Per-Sample Summary
| Sample | Cells In | Cells Out | Mito Median | Genes Median |
|--------|----------|-----------|-------------|--------------|
{sample_table}

## QC Status
- [x] Cell filtering complete
- [x] Gene filtering complete
- [x] Doublet detection complete
- [x] Ready for annotation
```

---

## Warning Signs

### Over-filtering Indicators

```
⚠️ > 30% cells removed
⚠️ One sample has > 50% removal
⚠️ Remaining cells < 1000 per sample
```

### Under-filtering Indicators

```
⚠️ Max mito > 30% still present
⚠️ Doublet score distribution bimodal but not filtered
⚠️ Clear outliers in scatter plots
```

### Action Items

| Warning | Action |
|---------|--------|
| Over-filtering | Relax thresholds, document rationale |
| Under-filtering | Tighten thresholds, check visualizations |
| Sample imbalance | Consider sample-specific thresholds |
| Unexpected loss | Check data loading, gene names |
