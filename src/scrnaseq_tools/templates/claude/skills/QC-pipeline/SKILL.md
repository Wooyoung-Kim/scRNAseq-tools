---
name: QC-pipeline
description: "Single-cell RNA-seq QC: cell/gene filtering, doublet detection, normalization, batch integration. Prerequisite for Annotation-agent."
---

# QC-pipeline Skill

```
╔══════════════════════════════════════════════════════════════════════╗
║  Single-cell RNA-seq Quality Control Pipeline                        ║
║  Prerequisite for Annotation-agent                                   ║
╚══════════════════════════════════════════════════════════════════════╝
```

## Overview

This skill provides comprehensive QC workflow for scRNA-seq data before annotation.
**MUST be completed before running Annotation-agent.**

## Goals

| QC Phase | Question | Output |
|----------|----------|--------|
| **Cell Filtering** | "이 세포가 viable한가?" | Filtered cells |
| **Gene Filtering** | "이 유전자가 informative한가?" | Filtered genes |
| **Doublet Detection** | "이 세포가 실제 single cell인가?" | Doublet scores |
| **Normalization** | "표현량이 비교 가능한가?" | Normalized matrix |
| **Batch Integration** | "batch effect가 제거되었는가?" | Harmony embedding |

---

## Quick Start

```python
import scanpy as sc
import scrublet as scr
import numpy as np

# 1. Load raw data
adata = sc.read_10x_mtx('raw_data/')

# 2. Run QC pipeline
# -> Read: phases/cell_filtering.md
# -> Read: phases/gene_filtering.md
# -> Read: phases/doublet_detection.md
# -> Read: phases/normalization.md
# -> Read: phases/batch_integration.md

# 3. Verify QC completion
# -> Read: core/qc_verification.md

# 4. Save and proceed to annotation
adata.write('qc_output/qc_complete.h5ad')
# -> Now ready for Annotation-agent
```

---

## Connection to Annotation

This skill is a **prerequisite** for Annotation-agent.

```
QC-pipeline                    Annotation-agent
├── Cell filtering        →    Uses clean cells
├── Gene filtering        →    Uses informative genes
├── Doublet detection     →    Removes doublet contamination
├── Normalization         →    Normalized expression for DE
└── Batch integration     →    X_pca_harmony for clustering
```

**Before annotation, verify:**
```python
# QC completion check
assert 'X_pca_harmony' in adata.obsm, "Run batch integration first!"
assert 'doublet_score' in adata.obs.columns, "Run doublet detection first!"
assert adata.obs['pct_counts_mt'].max() < 20, "High mito cells detected!"
print("✅ QC complete - ready for annotation")
```

---

## Skill Structure

```
QC-pipeline/
├── SKILL.md                    # This file
├── core/
│   ├── principles.md           # QC principles and thresholds
│   └── qc_verification.md      # Verification checklist
├── phases/
│   ├── cell_filtering.md       # Cell QC (mito, genes, counts)
│   ├── gene_filtering.md       # Gene filtering
│   ├── doublet_detection.md    # Scrublet workflow
│   ├── normalization.md        # Log-normalization + HVG
│   └── batch_integration.md    # Harmony integration
├── tools/
│   ├── visualization.md        # QC plots
│   └── metrics.md              # QC metrics table
├── examples/
│   └── workflow_example.md     # Complete example
└── connection/
    └── to_annotation.md        # How to proceed to annotation
```

---

## Required Packages

```python
# Core
scanpy>=1.9.0
anndata>=0.8.0

# Doublet detection
scrublet>=0.2.3

# Batch integration
harmonypy>=0.0.9
# or scanpy-harmony

# Optional
scipy>=1.9.0
pandas>=1.5.0
numpy>=1.21.0
```

---

## Workflow Overview

```
Step 1: Load Data
    ↓
Step 2: Calculate QC Metrics
    ↓
Step 3: Cell Filtering (mito, genes, counts)
    ↓
Step 4: Gene Filtering (min_cells)
    ↓
Step 5: Doublet Detection (scrublet)
    ↓
Step 6: Normalization + HVG Selection
    ↓
Step 7: PCA + Batch Integration (Harmony)
    ↓
Step 8: UMAP + Verification
    ↓
Step 9: Save → Annotation-agent
```

---

## Success Criteria

| Metric | Pass Criteria | Check |
|--------|---------------|-------|
| Mito % | All cells < 20% | `adata.obs['pct_counts_mt'].max() < 20` |
| Genes/cell | All cells > 200 | `adata.obs['n_genes_by_counts'].min() > 200` |
| Doublets | Removed or flagged | `'doublet_score' in adata.obs` |
| Normalization | Log-normalized | `adata.X.max() < 20` (log scale) |
| Harmony | Batch integrated | `'X_pca_harmony' in adata.obsm` |
| UMAP | Computed | `'X_umap' in adata.obsm` |

---

## Output Files

```
qc_output/
├── qc_complete.h5ad          # Main output for annotation
├── qc_metrics.csv            # Per-cell QC metrics
├── figures/
│   ├── violin_qc.png         # QC violin plots
│   ├── scatter_qc.png        # QC scatter plots
│   ├── doublet_histogram.png # Doublet score distribution
│   └── umap_batch.png        # UMAP colored by batch
└── reports/
    └── qc_summary.md         # QC summary report
```
