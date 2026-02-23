---
name: palantir
description: "Trajectory inference and pseudotime analysis for single-cell RNA-seq. Use for modeling differentiation trajectories, computing pseudotime ordering, inferring cell fate probabilities, gene expression trends along trajectories, and identifying terminal states."
---

# Palantir: Trajectory Inference for Single-Cell Data

## Overview

Palantir is a Python package for trajectory inference in single-cell data. It models cellular differentiation as a stochastic process, computing pseudotime ordering, cell fate probabilities, and gene expression dynamics along trajectories.

## When to Use This Skill

Use Palantir when you need to:
- Infer differentiation trajectories from single-cell data
- Compute pseudotime ordering of cells
- Determine cell fate probabilities for multiple lineages
- Identify terminal/progenitor states in differentiation
- Analyze gene expression dynamics along trajectories
- Quantify differentiation potential (entropy) of cells

## Core Capabilities

### 1. Core Functions (palantir.core)
Main trajectory inference algorithms. See `references/core.md` for:
- **run_palantir**: Main trajectory inference with pseudotime and fate probabilities
- **identify_terminal_states**: Automatic terminal state detection

### 2. Utility Functions (palantir.utils)
Preprocessing and helper functions. See `references/utils.md` for:
- **run_diffusion_maps**: Compute diffusion maps (core of Palantir)
- **determine_multiscale_space**: Automatic component selection
- **run_magic_imputation**: MAGIC imputation for smooth expression
- **early_cell/find_terminal_states**: Identify start and terminal cells

### 3. Post-Results Functions (palantir.presults)
Gene trend analysis. See `references/presults.md` for:
- **select_branch_cells**: Select cells per lineage
- **compute_gene_trends**: Gene expression along trajectories
- **cluster_gene_trends**: Cluster genes by trend similarity

### 4. Visualization (palantir.plot)
Plotting functions. See `references/plotting.md` for:
- **plot_palantir_results**: Main results visualization
- **plot_trajectory/plot_trajectories**: Trajectory arrows on embedding
- **plot_gene_trends**: Gene expression trends

## Quick Start

```python
import palantir
import scanpy as sc
import pandas as pd

# Load and preprocess data
adata = sc.read_h5ad("data.h5ad")
sc.pp.normalize_per_cell(adata)
palantir.preprocess.log_transform(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor="cell_ranger")
sc.pp.pca(adata)

# Diffusion maps
palantir.utils.run_diffusion_maps(adata, n_components=5)
palantir.utils.determine_multiscale_space(adata)

# Visualization
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Define start and terminal states
start_cell = palantir.utils.early_cell(adata, 'HSC', 'cell_type')
terminal_states = pd.Series(['Ery', 'Mono', 'DC'], index=['cell1', 'cell2', 'cell3'])

# Run Palantir
palantir.core.run_palantir(adata, start_cell, terminal_states=terminal_states, num_waypoints=500)

# Visualize
palantir.plot.plot_palantir_results(adata)
```

## Typical Workflow

```python
import palantir
import scanpy as sc

# 1. Preprocessing
adata = sc.read_h5ad("data.h5ad")
sc.pp.normalize_per_cell(adata)
palantir.preprocess.log_transform(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500)
sc.pp.pca(adata)

# 2. Diffusion maps (core of Palantir)
palantir.utils.run_diffusion_maps(adata, n_components=10, knn=30)
palantir.utils.determine_multiscale_space(adata)

# 3. Visualization setup
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# 4. Identify start and terminal cells
start_cell = palantir.utils.early_cell(adata, 'progenitor', 'cell_type')
terminal_states = palantir.utils.find_terminal_states(
    adata, celltypes=['Ery', 'Mono', 'DC'], celltype_column='cell_type'
)

# 5. Run Palantir
palantir.core.run_palantir(adata, start_cell, terminal_states=terminal_states)

# 6. Results stored in adata:
#   - adata.obs['palantir_pseudotime']
#   - adata.obs['palantir_entropy']
#   - adata.obsm['palantir_fate_probabilities']

# 7. Gene expression trends
palantir.utils.run_magic_imputation(adata)
palantir.presults.select_branch_cells(adata, q=0.01)
palantir.presults.compute_gene_trends(adata, expression_key='MAGIC_imputed_data')
palantir.plot.plot_gene_trends(adata, genes=['CD34', 'GATA1', 'MPO'])
```

**Key Design Principles:**
- **Diffusion maps**: Core representation for trajectory inference
- **Start cell required**: Must specify early/progenitor cell
- **Terminal states optional**: Can be auto-detected or manually specified
- **AnnData integration**: Results stored in obs, obsm, uns

## Results Interpretation

| Output | Location | Description |
|--------|----------|-------------|
| Pseudotime | `adata.obs['palantir_pseudotime']` | Ordering along differentiation |
| Entropy | `adata.obs['palantir_entropy']` | Differentiation potential (high=multipotent) |
| Fate probabilities | `adata.obsm['palantir_fate_probabilities']` | Probability of reaching each terminal state |
| Waypoints | `adata.uns['palantir_waypoints']` | Representative cells for computation |

## Parameter Guide

| Data Size | n_components | knn | num_waypoints |
|-----------|--------------|-----|---------------|
| < 5,000 cells | 5-10 | 30 | 500 |
| 5,000-20,000 | 10-15 | 30-50 | 1000 |
| > 20,000 cells | 10-20 | 50+ | 1200 |

## Best Practices

1. **Start cell selection**: Use biological knowledge (progenitor markers) or diffusion map extremes
2. **Terminal states**: Should correspond to known mature cell types
3. **Use MAGIC imputation** for gene trends (smooths noisy expression)
4. **Check entropy distribution**: High entropy = multipotent, Low = committed
5. **Run with different starts** to verify trajectory robustness

## Additional Resources

- **Core Functions**: `references/core.md` - run_palantir parameters
- **Utility Functions**: `references/utils.md` - Diffusion maps, MAGIC, cell finding
- **Post-Results**: `references/presults.md` - Gene trends and clustering
- **Visualization**: `references/plotting.md` - All plotting functions
- **Workflows**: `references/workflows.md` - Complete examples
- **Official Documentation**: https://palantir.readthedocs.io/

## Installation

```bash
pip install palantir
pip install palantir[all]  # With all dependencies
pip install palantir pygam  # For legacy GAM-based gene trends
```
