# Phase 3: Doublet Detection

```
╔══════════════════════════════════════════════════════════════════════╗
║  GOAL: "이 세포가 실제 single cell인가?"                               ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Why Doublet Detection?

- **Doublets**: Two or more cells captured in one droplet
- **Rate**: ~6% for 10,000 cells (increases with loading)
- **Impact**: Can appear as "novel" or "intermediate" cell types
- **Detection**: Scrublet simulates doublets and scores real cells

---

## Workflow

```
Step 1: Run Scrublet (per sample)
    ↓
Step 2: Visualize doublet scores
    ↓
Step 3: Apply threshold
    ↓
Step 4: Decision - Remove or Flag
```

---

## Step 1: Run Scrublet

```python
import scrublet as scr
import numpy as np
import pandas as pd

# IMPORTANT: Run per sample to avoid cross-sample artifacts
sample_col = 'sample' if 'sample' in adata.obs.columns else 'batch'

# Initialize columns
adata.obs['doublet_score'] = np.nan
adata.obs['predicted_doublet'] = False

if sample_col in adata.obs.columns:
    # Multi-sample: run separately
    print("Running Scrublet per sample...")

    for sample in adata.obs[sample_col].unique():
        print(f"\n  Processing {sample}...")

        # Get sample mask
        mask = adata.obs[sample_col] == sample
        subset = adata[mask]

        # Expected doublet rate (adjust based on loading)
        # ~0.8% per 1000 cells loaded
        n_cells = subset.n_obs
        expected_doublet_rate = min(0.08, n_cells / 1000 * 0.008)

        print(f"    Cells: {n_cells}")
        print(f"    Expected doublet rate: {expected_doublet_rate:.1%}")

        # Run Scrublet
        scrub = scr.Scrublet(
            subset.X,
            expected_doublet_rate=expected_doublet_rate
        )

        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30
        )

        # Store results
        adata.obs.loc[mask, 'doublet_score'] = doublet_scores
        adata.obs.loc[mask, 'predicted_doublet'] = predicted_doublets

        # Summary
        n_doublets = predicted_doublets.sum()
        print(f"    Detected doublets: {n_doublets} ({n_doublets/n_cells:.1%})")

else:
    # Single sample
    print("Running Scrublet on full dataset...")

    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets

print("\n✅ Scrublet complete")
```

---

## Step 2: Visualize Doublet Scores

```python
import matplotlib.pyplot as plt

# Histogram of doublet scores
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Overall distribution
axes[0].hist(adata.obs['doublet_score'], bins=50, edgecolor='black')
axes[0].axvline(0.25, color='r', linestyle='--', label='threshold=0.25')
axes[0].set_xlabel('Doublet Score')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Doublet Score Distribution')
axes[0].legend()

# Per-sample distribution
if sample_col in adata.obs.columns:
    for sample in adata.obs[sample_col].unique():
        subset = adata[adata.obs[sample_col] == sample]
        axes[1].hist(subset.obs['doublet_score'], bins=30, alpha=0.5, label=sample)
    axes[1].set_xlabel('Doublet Score')
    axes[1].set_ylabel('Number of Cells')
    axes[1].set_title('Per-Sample Distribution')
    axes[1].legend()

plt.tight_layout()
plt.savefig('figures/doublet_scores.png', dpi=150, bbox_inches='tight')
plt.show()

# If UMAP exists, visualize doublets
if 'X_umap' in adata.obsm:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sc.pl.umap(adata, color='doublet_score', ax=axes[0], show=False,
               title='Doublet Score')
    sc.pl.umap(adata, color='predicted_doublet', ax=axes[1], show=False,
               title='Predicted Doublets')

    plt.tight_layout()
    plt.savefig('figures/doublet_umap.png', dpi=150, bbox_inches='tight')
    plt.show()
```

---

## Step 3: Apply Threshold

```python
# Check Scrublet's automatic threshold
print("Doublet detection summary:")
print(f"  Total cells: {adata.n_obs}")
print(f"  Predicted doublets: {adata.obs['predicted_doublet'].sum()} "
      f"({adata.obs['predicted_doublet'].mean():.1%})")

# Optional: Adjust threshold manually
# If bimodal distribution is clear, you can adjust
custom_threshold = 0.25  # Scrublet default

adata.obs['predicted_doublet_manual'] = adata.obs['doublet_score'] > custom_threshold
print(f"  Manual threshold ({custom_threshold}): "
      f"{adata.obs['predicted_doublet_manual'].sum()} "
      f"({adata.obs['predicted_doublet_manual'].mean():.1%})")
```

---

## Step 4: Decision - Remove or Flag

### Option A: Remove Doublets (RECOMMENDED for most analyses)

```python
# Remove predicted doublets
n_before = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
n_after = adata.n_obs

print(f"\nDoublet removal:")
print(f"  Before: {n_before}")
print(f"  After: {n_after}")
print(f"  Removed: {n_before - n_after} ({(n_before - n_after)/n_before:.1%})")
```

### Option B: Flag Only (for later inspection)

```python
# Keep doublets but flag them
# Useful when you want to inspect what doublets look like

# They will be visible in annotation as potentially mixed cell types
print("Doublets flagged but NOT removed")
print("  predicted_doublet column can be used to filter later")
print("  Consider removing after annotation if clusters look suspicious")
```

---

## Doublet Summary

```python
# Store summary
adata.uns['doublet_detection'] = {
    'method': 'scrublet',
    'n_doublets_detected': int(adata.obs['predicted_doublet'].sum()) if 'predicted_doublet' in adata.obs.columns else 0,
    'doublet_rate': float(adata.obs['predicted_doublet'].mean()) if 'predicted_doublet' in adata.obs.columns else 0,
    'threshold': 0.25,
    'removed': True  # or False if only flagged
}
```

---

## Troubleshooting

### Issue: No clear bimodal distribution

```python
# Check if there's a clear separation
scores = adata.obs['doublet_score']
print(f"Score range: [{scores.min():.3f}, {scores.max():.3f}]")
print(f"Median: {scores.median():.3f}")
print(f"95th percentile: {scores.quantile(0.95):.3f}")

# If no clear separation, consider:
# 1. Adjust expected_doublet_rate
# 2. Use a fixed percentile threshold (e.g., top 5%)
```

### Issue: Very high doublet rate (>15%)

```python
# This might indicate:
# 1. Over-loading during library prep
# 2. Cell clumping
# 3. Technical issues

# Consider:
# - Checking per-sample rates
# - Using max_genes filtering as additional doublet filter
```

---

## Next Step

→ Read: phases/normalization.md
