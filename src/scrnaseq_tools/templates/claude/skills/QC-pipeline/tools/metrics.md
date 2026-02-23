# QC Metrics Reference

---

## Cell-Level Metrics (adata.obs)

| Metric | Description | Good Range | Calculated By |
|--------|-------------|------------|---------------|
| `n_genes_by_counts` | Genes with > 0 counts | 500-5000 | `sc.pp.calculate_qc_metrics()` |
| `total_counts` | Total UMI counts | 1000-50000 | `sc.pp.calculate_qc_metrics()` |
| `pct_counts_mt` | % mitochondrial | < 20% | `sc.pp.calculate_qc_metrics()` |
| `pct_counts_ribo` | % ribosomal | 5-60% | `sc.pp.calculate_qc_metrics()` |
| `pct_counts_hb` | % hemoglobin | < 1% (except blood) | `sc.pp.calculate_qc_metrics()` |
| `doublet_score` | Scrublet doublet score | < 0.25 | `scrublet` |
| `predicted_doublet` | Boolean doublet flag | False | `scrublet` |

---

## Gene-Level Metrics (adata.var)

| Metric | Description | Filter Criteria | Calculated By |
|--------|-------------|-----------------|---------------|
| `n_cells_by_counts` | Cells expressing gene | > 3 | `sc.pp.calculate_qc_metrics()` |
| `total_counts` | Total expression | - | `sc.pp.calculate_qc_metrics()` |
| `mean_counts` | Mean expression | - | `sc.pp.calculate_qc_metrics()` |
| `highly_variable` | HVG flag | Top 2000-4000 | `sc.pp.highly_variable_genes()` |
| `mt` | Mitochondrial gene | - | Manual assignment |
| `ribo` | Ribosomal gene | - | Manual assignment |

---

## Tissue-Specific Guidelines

### PBMC

```python
# Standard PBMC thresholds
min_genes = 200
max_genes = 4000
max_mito = 10  # PBMC should have low mito
```

### Solid Tissue (Tumor, Kidney, Heart)

```python
# Higher mito tolerance
min_genes = 500
max_genes = 6000
max_mito = 25  # Tissue damage during dissociation
```

### Brain

```python
# Brain-specific
min_genes = 500
max_genes = 8000  # Neurons can have high gene counts
max_mito = 15
```

---

## QC Summary Function

```python
def print_qc_summary(adata):
    """Print comprehensive QC summary."""
    print("=" * 60)
    print("QC METRICS SUMMARY")
    print("=" * 60)

    # Dataset size
    print(f"\n📊 Dataset Size:")
    print(f"   Cells: {adata.n_obs:,}")
    print(f"   Genes: {adata.n_vars:,}")

    # Cell metrics
    print(f"\n📈 Cell Metrics:")
    for metric in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']:
        if metric in adata.obs.columns:
            values = adata.obs[metric]
            print(f"   {metric}:")
            print(f"      Min: {values.min():.1f}")
            print(f"      Max: {values.max():.1f}")
            print(f"      Median: {values.median():.1f}")
            print(f"      Mean: {values.mean():.1f}")

    # Doublet info
    if 'doublet_score' in adata.obs.columns:
        print(f"\n🔍 Doublet Detection:")
        print(f"   Score range: [{adata.obs['doublet_score'].min():.3f}, "
              f"{adata.obs['doublet_score'].max():.3f}]")
        if 'predicted_doublet' in adata.obs.columns:
            n_doublets = adata.obs['predicted_doublet'].sum()
            print(f"   Predicted doublets: {n_doublets} ({n_doublets/adata.n_obs:.1%})")

    # Gene categories
    print(f"\n🧬 Gene Categories:")
    for cat in ['mt', 'ribo', 'hb', 'ig', 'tcr']:
        if cat in adata.var.columns:
            print(f"   {cat}: {adata.var[cat].sum()}")

    # HVGs
    if 'highly_variable' in adata.var.columns:
        print(f"\n⭐ Highly Variable Genes:")
        print(f"   Count: {adata.var['highly_variable'].sum()}")

    # Embeddings
    print(f"\n🗺️ Embeddings:")
    for key in ['X_pca', 'X_pca_harmony', 'X_umap']:
        if key in adata.obsm:
            print(f"   {key}: {adata.obsm[key].shape}")

    # Samples
    for col in ['sample', 'batch']:
        if col in adata.obs.columns:
            print(f"\n📦 {col.capitalize()}s:")
            for val in adata.obs[col].unique():
                n = (adata.obs[col] == val).sum()
                print(f"   {val}: {n:,} cells")

    print("=" * 60)


# Usage
print_qc_summary(adata)
```

---

## Export QC Metrics

```python
def export_qc_metrics(adata, output_dir='qc_output'):
    """Export QC metrics to CSV files."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Cell metrics
    cell_metrics = adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].copy()
    if 'doublet_score' in adata.obs.columns:
        cell_metrics['doublet_score'] = adata.obs['doublet_score']
    if 'sample' in adata.obs.columns:
        cell_metrics['sample'] = adata.obs['sample']

    cell_metrics.to_csv(f'{output_dir}/cell_qc_metrics.csv')
    print(f"Saved: {output_dir}/cell_qc_metrics.csv")

    # Gene metrics
    gene_metrics = adata.var[['n_cells_by_counts']].copy()
    if 'highly_variable' in adata.var.columns:
        gene_metrics['highly_variable'] = adata.var['highly_variable']

    gene_metrics.to_csv(f'{output_dir}/gene_qc_metrics.csv')
    print(f"Saved: {output_dir}/gene_qc_metrics.csv")

    # Sample summary
    if 'sample' in adata.obs.columns:
        sample_summary = adata.obs.groupby('sample').agg({
            'n_genes_by_counts': ['median', 'mean', 'std'],
            'total_counts': ['median', 'mean', 'std'],
            'pct_counts_mt': ['median', 'mean', 'std']
        })
        sample_summary.to_csv(f'{output_dir}/sample_summary.csv')
        print(f"Saved: {output_dir}/sample_summary.csv")


# Usage
export_qc_metrics(adata)
```
