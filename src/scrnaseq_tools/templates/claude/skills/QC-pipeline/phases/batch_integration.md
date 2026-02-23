# Phase 5: Batch Integration

```
╔══════════════════════════════════════════════════════════════════════╗
║  GOAL: "batch effect가 제거되고 생물학적 variation이 보존되었는가?"      ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## When to Run Batch Integration

| Condition | Action |
|-----------|--------|
| Single sample | Skip (no batch effect) |
| Multiple samples, same protocol | Run Harmony |
| Multiple protocols/technologies | Run Harmony + careful QC |
| Single-cell + bulk | Different pipeline needed |

---

## Workflow

```
Step 1: Check batch structure
    ↓
Step 2: Pre-integration UMAP (for comparison)
    ↓
Step 3: Run Harmony
    ↓
Step 4: Post-integration UMAP
    ↓
Step 5: Verify batch mixing
```

---

## Step 1: Check Batch Structure

```python
import scanpy as sc
import pandas as pd

# Identify batch column
batch_col = None
for col in ['batch', 'sample', 'donor', 'patient', 'library']:
    if col in adata.obs.columns:
        batch_col = col
        break

if batch_col is None:
    print("⚠️ No batch column found - skipping batch integration")
    # Create dummy batch if needed
    # adata.obs['batch'] = 'batch1'
else:
    print(f"Batch column: '{batch_col}'")
    print(f"\nBatch distribution:")
    print(adata.obs[batch_col].value_counts())

    # Check balance
    batch_counts = adata.obs[batch_col].value_counts()
    min_batch = batch_counts.min()
    max_batch = batch_counts.max()
    print(f"\nBatch balance:")
    print(f"  Min batch size: {min_batch}")
    print(f"  Max batch size: {max_batch}")
    print(f"  Ratio: {max_batch/min_batch:.1f}x")
```

---

## Step 2: Pre-Integration UMAP

```python
# Compute neighbors on uncorrected PCA
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep='X_pca')
sc.tl.umap(adata)

# Visualize batch effect
sc.pl.umap(
    adata,
    color=batch_col,
    title='Pre-integration (uncorrected)',
    save='_pre_integration.png'
)
```

---

## Step 3: Run Harmony

### Option A: scanpy.external (Recommended)

```python
import scanpy.external as sce

print("Running Harmony integration...")

# Run Harmony on PCA
sce.pp.harmony_integrate(
    adata,
    key=batch_col,
    basis='X_pca',
    adjusted_basis='X_pca_harmony',
    max_iter_harmony=20
)

# Verify output
assert 'X_pca_harmony' in adata.obsm
print(f"✅ Harmony complete: X_pca_harmony shape = {adata.obsm['X_pca_harmony'].shape}")
```

### Option B: harmonypy directly

```python
import harmonypy as hm
import pandas as pd

# Run Harmony
ho = hm.run_harmony(
    adata.obsm['X_pca'],
    adata.obs,
    batch_col,
    max_iter_harmony=20
)

# Store result
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
```

---

## Step 4: Post-Integration UMAP

```python
# Compute neighbors on CORRECTED PCA
sc.pp.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=50,
    use_rep='X_pca_harmony'
)
sc.tl.umap(adata)

# Visualize
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Pre-integration (store old UMAP if needed)
# axes[0].scatter(adata.obsm['X_umap_uncorrected'][:, 0], ...)

# Post-integration
sc.pl.umap(
    adata,
    color=batch_col,
    ax=axes[1],
    show=False,
    title='Post-integration (Harmony)'
)

plt.tight_layout()
plt.savefig('figures/batch_integration_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
```

---

## Step 5: Verify Batch Mixing

### Visual Check

```python
# Side-by-side comparison
sc.pl.umap(
    adata,
    color=[batch_col, 'n_genes_by_counts', 'pct_counts_mt'],
    ncols=3,
    save='_batch_qc.png'
)
```

### Quantitative Check: kBET or LISI

```python
# Simple batch mixing score (cells per neighborhood)
def calculate_batch_mixing(adata, batch_key, n_neighbors=30):
    """Calculate batch mixing score."""
    from sklearn.neighbors import NearestNeighbors

    # Get neighbors
    X = adata.obsm['X_pca_harmony']
    nn = NearestNeighbors(n_neighbors=n_neighbors)
    nn.fit(X)
    _, indices = nn.kneighbors(X)

    # Calculate entropy per cell
    batches = adata.obs[batch_key].values
    n_batches = len(np.unique(batches))

    scores = []
    for i in range(len(adata)):
        neighbor_batches = batches[indices[i]]
        batch_counts = pd.Series(neighbor_batches).value_counts(normalize=True)
        # Entropy
        entropy = -(batch_counts * np.log(batch_counts + 1e-10)).sum()
        scores.append(entropy)

    max_entropy = np.log(n_batches)
    mixing_score = np.mean(scores) / max_entropy

    return mixing_score

mixing = calculate_batch_mixing(adata, batch_col)
print(f"Batch mixing score: {mixing:.3f}")
print(f"  (0 = no mixing, 1 = perfect mixing)")
```

---

## Integration QC Summary

```python
# Store integration info
adata.uns['batch_integration'] = {
    'method': 'harmony',
    'batch_key': batch_col,
    'n_batches': adata.obs[batch_col].nunique(),
    'n_iterations': 20
}

# Final summary
print("\n" + "=" * 50)
print("BATCH INTEGRATION SUMMARY")
print("=" * 50)
print(f"Method: Harmony")
print(f"Batch key: {batch_col}")
print(f"Number of batches: {adata.obs[batch_col].nunique()}")
print(f"Output: X_pca_harmony ({adata.obsm['X_pca_harmony'].shape})")
print("=" * 50)
```

---

## Verification

```python
# Final QC checks
assert 'X_pca_harmony' in adata.obsm, "Harmony not computed!"
assert 'X_umap' in adata.obsm, "UMAP not computed!"

print("✅ Batch integration complete")
print("✅ Ready for Annotation-agent")
```

---

## When Harmony Isn't Enough

```python
# If batch effects remain visible:

# Option 1: Increase iterations
sce.pp.harmony_integrate(adata, key=batch_col, max_iter_harmony=50)

# Option 2: Try scVI (for complex batch effects)
# import scvi
# scvi.model.SCVI.setup_anndata(adata, batch_key=batch_col)
# model = scvi.model.SCVI(adata)
# model.train()
# adata.obsm['X_scvi'] = model.get_latent_representation()

# Option 3: Check if batches have biological differences
# (e.g., different conditions, tissues)
```

---

## Next Step

→ Read: connection/to_annotation.md
