# QC Verification Checklist

```
╔══════════════════════════════════════════════════════════════════════╗
║  ⚠️ MUST pass verification before proceeding to annotation           ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Verification Function

```python
def verify_qc_complete(adata, verbose=True):
    """
    Verify QC pipeline completion before annotation.
    Raises AssertionError if any check fails.
    """
    errors = []
    warnings = []

    # === REQUIRED CHECKS ===

    # 1. Cell filtering metrics exist
    required_obs = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    for col in required_obs:
        if col not in adata.obs.columns:
            errors.append(f"Missing obs column: {col}")

    # 2. Mito threshold check
    if 'pct_counts_mt' in adata.obs.columns:
        max_mito = adata.obs['pct_counts_mt'].max()
        if max_mito > 25:
            errors.append(f"High mito cells present: max={max_mito:.1f}%")
        elif max_mito > 20:
            warnings.append(f"Borderline mito: max={max_mito:.1f}%")

    # 3. Gene filtering check
    if adata.obs['n_genes_by_counts'].min() < 100:
        errors.append(f"Very low gene cells present: min={adata.obs['n_genes_by_counts'].min()}")
    elif adata.obs['n_genes_by_counts'].min() < 200:
        warnings.append(f"Low gene cells present: min={adata.obs['n_genes_by_counts'].min()}")

    # 4. Doublet detection
    if 'doublet_score' not in adata.obs.columns:
        errors.append("Doublet detection not performed (missing 'doublet_score')")
    else:
        if 'predicted_doublet' not in adata.obs.columns:
            warnings.append("Doublet scores present but no predictions")

    # 5. Normalization check
    if adata.X.max() > 100:
        errors.append(f"Data not normalized (max value: {adata.X.max():.1f})")

    # 6. HVG selection
    if 'highly_variable' not in adata.var.columns:
        errors.append("HVG selection not performed")
    else:
        n_hvg = adata.var['highly_variable'].sum()
        if n_hvg < 1000:
            warnings.append(f"Few HVGs selected: {n_hvg}")
        elif n_hvg > 5000:
            warnings.append(f"Many HVGs selected: {n_hvg}")

    # 7. PCA
    if 'X_pca' not in adata.obsm:
        errors.append("PCA not computed")
    else:
        n_pcs = adata.obsm['X_pca'].shape[1]
        if n_pcs < 30:
            warnings.append(f"Few PCs: {n_pcs}")

    # 8. Batch integration (if multi-batch)
    if 'batch' in adata.obs.columns or 'sample' in adata.obs.columns:
        if 'X_pca_harmony' not in adata.obsm:
            errors.append("Multi-sample data but Harmony not run")

    # 9. UMAP
    if 'X_umap' not in adata.obsm:
        errors.append("UMAP not computed")

    # === REPORT ===
    if verbose:
        print("=" * 60)
        print("QC VERIFICATION REPORT")
        print("=" * 60)

        # Summary
        print(f"\n📊 Dataset Summary:")
        print(f"   Cells: {adata.n_obs:,}")
        print(f"   Genes: {adata.n_vars:,}")
        if 'highly_variable' in adata.var.columns:
            print(f"   HVGs: {adata.var['highly_variable'].sum():,}")

        # Errors
        if errors:
            print(f"\n❌ ERRORS ({len(errors)}):")
            for e in errors:
                print(f"   - {e}")

        # Warnings
        if warnings:
            print(f"\n⚠️ WARNINGS ({len(warnings)}):")
            for w in warnings:
                print(f"   - {w}")

        # Status
        if errors:
            print(f"\n❌ QC INCOMPLETE - Cannot proceed to annotation")
        else:
            print(f"\n✅ QC COMPLETE - Ready for annotation")

        print("=" * 60)

    # Raise if errors
    if errors:
        raise AssertionError(
            "QC verification failed:\n" +
            "\n".join(f"  - {e}" for e in errors)
        )

    return True


def print_qc_summary(adata):
    """Print detailed QC summary."""
    print("\n" + "=" * 60)
    print("QC SUMMARY")
    print("=" * 60)

    # Cell metrics
    print("\n📊 Cell Metrics:")
    print(f"   n_genes_by_counts: {adata.obs['n_genes_by_counts'].median():.0f} (median)")
    print(f"   total_counts: {adata.obs['total_counts'].median():.0f} (median)")
    print(f"   pct_counts_mt: {adata.obs['pct_counts_mt'].median():.1f}% (median)")

    # Doublet info
    if 'doublet_score' in adata.obs.columns:
        print("\n🔍 Doublet Detection:")
        print(f"   Score range: [{adata.obs['doublet_score'].min():.3f}, {adata.obs['doublet_score'].max():.3f}]")
        if 'predicted_doublet' in adata.obs.columns:
            n_doublets = adata.obs['predicted_doublet'].sum()
            print(f"   Predicted doublets: {n_doublets} ({n_doublets/adata.n_obs:.1%})")

    # Normalization
    print("\n📈 Normalization:")
    print(f"   X max value: {adata.X.max():.2f}")
    print(f"   X min value: {adata.X.min():.2f}")

    # Embeddings
    print("\n🗺️ Embeddings:")
    for key in ['X_pca', 'X_pca_harmony', 'X_umap']:
        if key in adata.obsm:
            print(f"   {key}: {adata.obsm[key].shape}")

    # Sample info
    sample_col = 'batch' if 'batch' in adata.obs.columns else 'sample' if 'sample' in adata.obs.columns else None
    if sample_col:
        print(f"\n📦 Samples ({sample_col}):")
        for sample in adata.obs[sample_col].unique():
            n = (adata.obs[sample_col] == sample).sum()
            print(f"   {sample}: {n:,} cells")

    print("=" * 60)
```

---

## Quick Verification

```python
# One-liner verification before annotation
verify_qc_complete(adata)

# If passes, proceed to annotation
# -> Load Annotation-agent skill
```

---

## Checklist (Manual)

```
CELL FILTERING:
- [ ] n_genes_by_counts > 200 for all cells
- [ ] pct_counts_mt < 20% for all cells
- [ ] No extreme outliers in scatter plot

GENE FILTERING:
- [ ] Genes expressed in >= 3 cells
- [ ] Mitochondrial genes identified

DOUBLET DETECTION:
- [ ] Scrublet run successfully
- [ ] doublet_score in adata.obs
- [ ] Doublets flagged or removed

NORMALIZATION:
- [ ] Total count normalized
- [ ] Log-transformed
- [ ] HVGs selected (2000-4000)

BATCH INTEGRATION:
- [ ] Harmony run (if multi-sample)
- [ ] X_pca_harmony in adata.obsm
- [ ] Batch mixing verified in UMAP

EMBEDDINGS:
- [ ] PCA computed (50 components)
- [ ] UMAP computed
- [ ] Neighbors computed

READY FOR ANNOTATION:
- [ ] All checks passed
- [ ] verify_qc_complete() returns True
```

---

## Common Issues and Fixes

### Issue: "High mito cells present"

```python
# Check distribution
sc.pl.violin(adata, 'pct_counts_mt')

# Option 1: Filter more strictly
adata = adata[adata.obs['pct_counts_mt'] < 15].copy()

# Option 2: If tissue-specific (kidney, heart)
# Document and proceed with higher threshold
```

### Issue: "Doublet detection not performed"

```python
import scrublet as scr

# Run scrublet per sample
for sample in adata.obs['sample'].unique():
    mask = adata.obs['sample'] == sample
    subset = adata[mask]

    scrub = scr.Scrublet(subset.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs.loc[mask, 'doublet_score'] = doublet_scores
    adata.obs.loc[mask, 'predicted_doublet'] = predicted_doublets
```

### Issue: "Data not normalized"

```python
# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Verify
print(f"Max value after normalization: {adata.X.max()}")
```

### Issue: "Harmony not run"

```python
import scanpy.external as sce

# Run Harmony
sce.pp.harmony_integrate(adata, key='batch', basis='X_pca', adjusted_basis='X_pca_harmony')

# Verify
print('X_pca_harmony' in adata.obsm)
```
