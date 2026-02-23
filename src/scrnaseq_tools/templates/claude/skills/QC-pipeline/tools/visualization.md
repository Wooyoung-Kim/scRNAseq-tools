# QC Visualization Tools

---

## Standard QC Plots

### 1. Violin Plots

```python
import scanpy as sc
import matplotlib.pyplot as plt

# Basic QC violin
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_qc_violin.png'
)

# Per-sample violin
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='sample',
    rotation=45,
    save='_qc_by_sample.png'
)
```

### 2. Scatter Plots

```python
# Counts vs Genes colored by mito
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
              color='pct_counts_mt', ax=axes[0], show=False)
axes[0].axhline(200, color='r', linestyle='--', alpha=0.5, label='min_genes')
axes[0].axhline(6000, color='r', linestyle='--', alpha=0.5, label='max_genes')
axes[0].legend()

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',
              ax=axes[1], show=False)
axes[1].axhline(20, color='r', linestyle='--', alpha=0.5, label='max_mito')
axes[1].legend()

sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt',
              ax=axes[2], show=False)
axes[2].axhline(20, color='r', linestyle='--', alpha=0.5, label='max_mito')
axes[2].legend()

plt.tight_layout()
plt.savefig('figures/qc_scatter.png', dpi=150, bbox_inches='tight')
```

### 3. Histograms with Thresholds

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Genes histogram
axes[0].hist(adata.obs['n_genes_by_counts'], bins=100)
axes[0].axvline(200, color='r', linestyle='--', label='min=200')
axes[0].axvline(6000, color='r', linestyle='--', label='max=6000')
axes[0].set_xlabel('Genes per cell')
axes[0].set_ylabel('Count')
axes[0].legend()

# Counts histogram
axes[1].hist(adata.obs['total_counts'], bins=100)
axes[1].set_xlabel('UMIs per cell')
axes[1].set_ylabel('Count')

# Mito histogram
axes[2].hist(adata.obs['pct_counts_mt'], bins=100)
axes[2].axvline(20, color='r', linestyle='--', label='max=20%')
axes[2].set_xlabel('Mitochondrial %')
axes[2].set_ylabel('Count')
axes[2].legend()

plt.tight_layout()
plt.savefig('figures/qc_histograms.png', dpi=150, bbox_inches='tight')
```

---

## Doublet Visualization

```python
# Doublet score histogram
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].hist(adata.obs['doublet_score'], bins=50, edgecolor='black')
axes[0].axvline(0.25, color='r', linestyle='--', label='threshold')
axes[0].set_xlabel('Doublet Score')
axes[0].set_ylabel('Count')
axes[0].set_title('Doublet Score Distribution')
axes[0].legend()

# Per-sample
for sample in adata.obs['sample'].unique():
    subset = adata[adata.obs['sample'] == sample]
    axes[1].hist(subset.obs['doublet_score'], bins=30, alpha=0.5, label=sample)
axes[1].set_xlabel('Doublet Score')
axes[1].set_ylabel('Count')
axes[1].set_title('Per-Sample Distribution')
axes[1].legend()

plt.tight_layout()
plt.savefig('figures/doublet_scores.png', dpi=150, bbox_inches='tight')

# UMAP with doublet scores (if UMAP exists)
if 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color=['doublet_score', 'predicted_doublet'],
               ncols=2, save='_doublet_umap.png')
```

---

## Batch Integration Visualization

```python
# Pre vs Post integration comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Pre-integration (if stored)
if 'X_umap_uncorrected' in adata.obsm:
    for i, sample in enumerate(adata.obs['sample'].unique()):
        mask = adata.obs['sample'] == sample
        axes[0].scatter(
            adata.obsm['X_umap_uncorrected'][mask, 0],
            adata.obsm['X_umap_uncorrected'][mask, 1],
            s=1, alpha=0.5, label=sample
        )
    axes[0].set_title('Pre-integration')
    axes[0].legend(markerscale=5)

# Post-integration
sc.pl.umap(adata, color='sample', ax=axes[1], show=False,
           title='Post-integration (Harmony)')

plt.tight_layout()
plt.savefig('figures/batch_comparison.png', dpi=150, bbox_inches='tight')
```

---

## Combined QC Summary Figure

```python
def create_qc_summary_figure(adata, save_path='figures/qc_summary.png'):
    """Create comprehensive QC summary figure."""
    fig = plt.figure(figsize=(16, 12))

    # Layout: 3x3 grid
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Row 1: Basic QC
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(adata.obs['n_genes_by_counts'], bins=50)
    ax1.axvline(200, color='r', linestyle='--')
    ax1.set_xlabel('Genes/cell')
    ax1.set_title('Gene Distribution')

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(adata.obs['total_counts'], bins=50)
    ax2.set_xlabel('UMIs/cell')
    ax2.set_title('UMI Distribution')

    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(adata.obs['pct_counts_mt'], bins=50)
    ax3.axvline(20, color='r', linestyle='--')
    ax3.set_xlabel('Mito %')
    ax3.set_title('Mito Distribution')

    # Row 2: Scatter + Doublet
    ax4 = fig.add_subplot(gs[1, 0])
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                  color='pct_counts_mt', ax=ax4, show=False)
    ax4.set_title('Counts vs Genes')

    ax5 = fig.add_subplot(gs[1, 1])
    ax5.hist(adata.obs['doublet_score'], bins=50)
    ax5.axvline(0.25, color='r', linestyle='--')
    ax5.set_xlabel('Doublet Score')
    ax5.set_title('Doublet Scores')

    ax6 = fig.add_subplot(gs[1, 2])
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='doublet_score', ax=ax6, show=False)
    ax6.set_title('UMAP: Doublet Score')

    # Row 3: UMAP by sample and QC
    ax7 = fig.add_subplot(gs[2, 0])
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='sample', ax=ax7, show=False)
    ax7.set_title('UMAP: Sample')

    ax8 = fig.add_subplot(gs[2, 1])
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='n_genes_by_counts', ax=ax8, show=False)
    ax8.set_title('UMAP: Genes')

    ax9 = fig.add_subplot(gs[2, 2])
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='pct_counts_mt', ax=ax9, show=False)
    ax9.set_title('UMAP: Mito %')

    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Saved: {save_path}")


# Usage
create_qc_summary_figure(adata)
```
