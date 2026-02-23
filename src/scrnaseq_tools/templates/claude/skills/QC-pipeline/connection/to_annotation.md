# Connection: QC-pipeline → Annotation-agent

```
╔══════════════════════════════════════════════════════════════════════╗
║  QC-pipeline output is the INPUT for Annotation-agent                ║
║  This document defines the handoff                                   ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Pipeline Flow

```
Raw Data (10x/h5ad)
        ↓
┌───────────────────────┐
│    QC-pipeline        │
├───────────────────────┤
│ 1. Cell filtering     │
│ 2. Gene filtering     │
│ 3. Doublet detection  │
│ 4. Normalization      │
│ 5. Batch integration  │
└───────────────────────┘
        ↓
    qc_complete.h5ad
        ↓
┌───────────────────────┐
│  Annotation-agent     │
├───────────────────────┤
│ Tier 1: Major types   │ → TF Activity (decoupler)
│ Tier 2: Dev states    │ → + Trajectory (palantir)
│ Tier 3: Func states   │ → + Pathway (decoupler)
└───────────────────────┘
        ↓
    annotated.h5ad
```

---

## Required Output from QC-pipeline

### adata.obs (cell metadata)

| Column | Type | Description | Used by |
|--------|------|-------------|---------|
| `n_genes_by_counts` | int | Genes per cell | QC verification |
| `total_counts` | int | UMIs per cell | QC verification |
| `pct_counts_mt` | float | Mito percentage | QC verification |
| `doublet_score` | float | Scrublet score | QC verification |
| `sample` or `batch` | str | Sample ID | Batch integration |

### adata.var (gene metadata)

| Column | Type | Description | Used by |
|--------|------|-------------|---------|
| `highly_variable` | bool | HVG flag | PCA, clustering |
| `mt` | bool | Mitochondrial gene | QC metrics |
| `n_cells_by_counts` | int | Cells expressing | Gene filtering |

### adata.obsm (embeddings)

| Key | Shape | Description | Used by |
|-----|-------|-------------|---------|
| `X_pca` | (n_cells, 50) | PCA | Backup if no harmony |
| `X_pca_harmony` | (n_cells, 50) | Harmony-corrected | **Clustering, neighbors** |
| `X_umap` | (n_cells, 2) | UMAP | Visualization |

### adata.X (expression matrix)

- **Log-normalized** (max ~10-15)
- **Scaled** for PCA

### adata.raw

- **Raw counts** preserved for DE analysis

---

## Transition Checkpoint

```python
def verify_qc_for_annotation(adata):
    """
    Verify QC output is ready for Annotation-agent.
    Run this before starting annotation.
    """
    errors = []

    # === obs checks ===
    required_obs = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    for col in required_obs:
        if col not in adata.obs.columns:
            errors.append(f"Missing obs: {col}")

    if 'doublet_score' not in adata.obs.columns:
        errors.append("Doublet detection not run")

    if adata.obs['pct_counts_mt'].max() > 25:
        errors.append(f"High mito cells present: {adata.obs['pct_counts_mt'].max():.1f}%")

    # === var checks ===
    if 'highly_variable' not in adata.var.columns:
        errors.append("HVGs not selected")

    # === obsm checks ===
    if 'X_pca' not in adata.obsm:
        errors.append("PCA not computed")

    # Prefer harmony-corrected if multi-sample
    if 'sample' in adata.obs.columns or 'batch' in adata.obs.columns:
        if 'X_pca_harmony' not in adata.obsm:
            errors.append("Multi-sample data but Harmony not run")

    if 'X_umap' not in adata.obsm:
        errors.append("UMAP not computed")

    # === X checks ===
    if adata.X.max() > 50:
        errors.append(f"Data not log-normalized (max={adata.X.max():.1f})")

    # === raw checks ===
    if adata.raw is None:
        errors.append("Raw counts not preserved (adata.raw is None)")

    # === Report ===
    if errors:
        print("❌ QC INCOMPLETE - Cannot proceed to annotation:")
        for e in errors:
            print(f"   - {e}")
        raise AssertionError("QC verification failed")
    else:
        print("✅ QC COMPLETE - Ready for Annotation-agent")
        print(f"   Cells: {adata.n_obs:,}")
        print(f"   Genes: {adata.n_vars:,}")
        print(f"   HVGs: {adata.var['highly_variable'].sum():,}")
        if 'X_pca_harmony' in adata.obsm:
            print(f"   Using: X_pca_harmony (batch-corrected)")
        else:
            print(f"   Using: X_pca (no batch correction)")

    return True


# Usage
# verify_qc_for_annotation(adata)
# -> If passes, proceed to Annotation-agent
```

---

## Starting Annotation-agent

After QC verification passes:

```python
# Save QC-complete data
adata.write('qc_output/qc_complete.h5ad')

# === TRANSITION TO ANNOTATION ===

# Load data for annotation
adata = sc.read_h5ad('qc_output/qc_complete.h5ad')

# Verify (redundant but safe)
verify_qc_for_annotation(adata)

# Now follow Annotation-agent workflow
# -> Read: Annotation-agent/SKILL.md
# -> Read: Annotation-agent/phases/tier1.md

# Key settings for annotation
use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'

# Tier 1: Major cell types
sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, n_pcs=50)
sc.tl.leiden(adata, resolution=0.8, key_added='tier1_cluster')
sc.tl.rank_genes_groups(adata, groupby='tier1_cluster', method='wilcoxon')

# ... continue with annotation workflow
```

---

## Troubleshooting Transition Issues

### Issue: "High mito cells present"

```python
# Option 1: Filter more strictly in QC
adata = adata[adata.obs['pct_counts_mt'] < 15].copy()

# Option 2: Document and proceed (if biologically justified)
print("⚠️ Proceeding with high-mito cells - documented as tissue artifact")
```

### Issue: "Harmony not run"

```python
# Run Harmony
import scanpy.external as sce
sce.pp.harmony_integrate(adata, key='batch', basis='X_pca', adjusted_basis='X_pca_harmony')

# Re-compute UMAP
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
```

### Issue: "Raw counts not preserved"

```python
# If adata.raw is None, you need to go back to before normalization
# Or reload the QC-filtered but pre-normalized data

# Best practice: always save raw before normalization
# adata.raw = adata.copy()  # BEFORE sc.pp.normalize_total()
```

---

## Summary: QC → Annotation Handoff

```
QC-pipeline OUTPUTS:
├── adata.obs
│   ├── n_genes_by_counts ✓
│   ├── total_counts ✓
│   ├── pct_counts_mt ✓ (< 20%)
│   ├── doublet_score ✓
│   └── sample/batch ✓
├── adata.var
│   └── highly_variable ✓
├── adata.obsm
│   ├── X_pca ✓
│   ├── X_pca_harmony ✓ (if multi-sample)
│   └── X_umap ✓
├── adata.X
│   └── log-normalized ✓
└── adata.raw
    └── raw counts ✓

→ verify_qc_for_annotation(adata)

Annotation-agent INPUTS:
├── Clean cells (no doublets, low-mito)
├── Informative genes (HVGs)
├── Normalized expression (for DE)
├── Batch-corrected embedding (for clustering)
└── Raw counts (for future analyses)
```
