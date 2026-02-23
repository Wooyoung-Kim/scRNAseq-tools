---
name: decoupler
description: "Pathway and transcription factor activity inference from omics data. Use for TF activity scoring, pathway enrichment, functional analysis of single-cell and bulk RNA-seq data, pseudobulk analysis, and trajectory-based activity inference."
---

# Decoupler: Activity Inference for Omics Data

## Overview

Decoupler is a Python package for inferring biological activities from omics data using prior knowledge. It provides a unified framework for transcription factor (TF) activity inference, pathway activity scoring, and functional enrichment analysis.

## When to Use This Skill

Use decoupler when you need to:
- Infer TF activities from gene expression data
- Score pathway activities across cells or samples
- Perform enrichment analysis (ORA, GSEA) on gene lists
- Aggregate single-cell data to pseudobulk for differential analysis
- Analyze activities along trajectories (pseudotime)
- Compare activities between conditions or cell types

## Core Capabilities

### 1. Enrichment Methods (dc.mt)
Multiple statistical methods for activity inference. See `references/methods.md` for:
- **ULM/MLM**: Linear models for TF and pathway activity (recommended)
- **ORA/GSEA**: Over-representation and gene set enrichment analysis
- **AUCell/GSVA/VIPER**: Alternative scoring methods
- **Consensus**: Aggregate results from multiple methods

### 2. Preprocessing (dc.pp)
Data preparation and aggregation functions. See `references/preprocessing.md` for:
- **pseudobulk**: Aggregate single-cell to pseudobulk
- **get_obsm**: Extract activities to obs for plotting
- **filter_by_expr**: Filter genes for DESeq2 compatibility

### 3. Visualization (dc.pl)
Plotting functions for activity analysis. See `references/plotting.md` for:
- **barplot/dotplot**: Activity score visualization
- **volcano**: Differential activity plots
- **network**: TF-target network visualization

### 4. Prior Knowledge Databases (dc.op)
Built-in biological networks. See `references/databases.md` for:
- **CollecTRI**: TF-target regulatory network
- **PROGENy**: Pathway responsive genes
- **Hallmark**: MSigDB hallmark gene sets

## Quick Start

### TF Activity Scoring

```python
import decoupler as dc
import scanpy as sc

# Load data (log-normalized)
adata = sc.read_h5ad("data.h5ad")

# Get TF network
net = dc.op.collectri(organism='human')

# Run ULM (recommended for TF activity)
dc.run_ulm(
    mat=adata,
    net=net,
    source='source',
    target='target',
    weight='weight'
)

# Extract to obs and visualize
dc.pp.get_obsm(adata, obsm_key='ulm_estimate')
sc.pl.umap(adata, color=['STAT1', 'NF-kB'])
```

### Pathway Activity Scoring

```python
# Get pathway network
model = dc.op.progeny(organism='human', top=500)

# Run MLM (recommended for pathways)
dc.run_mlm(
    mat=adata,
    net=model,
    source='source',
    target='target',
    weight='weight'
)

# Visualize
dc.pp.get_obsm(adata, obsm_key='mlm_estimate')
sc.pl.umap(adata, color=['JAK-STAT', 'NFkB', 'Hypoxia'])
```

## Typical Workflow

```python
import decoupler as dc
import scanpy as sc

# 1. Load and ensure log-normalized data
adata = sc.read_h5ad("data.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 2. Load prior knowledge network
tf_net = dc.op.collectri(organism='human')

# 3. Run enrichment method
dc.run_ulm(mat=adata, net=tf_net, source='source', target='target', weight='weight')

# 4. Extract and visualize
dc.pp.get_obsm(adata, obsm_key='ulm_estimate')
sc.pl.umap(adata, color=['STAT1', 'MYC', 'leiden'])

# 5. Statistical comparison between groups
dc.tl.rankby_group(adata, group_key='condition', obs_key='ulm_estimate', reference='control')
```

**Key Design Principles:**
- **Log-normalized data**: Most methods expect log-normalized expression
- **Prior knowledge networks**: DataFrame with source, target, weight columns
- **AnnData integration**: Results stored in adata.obsm (estimates, p-values)
- **Scanpy compatibility**: Use scanpy plotting after extracting to obs

## Method Selection Guide

| Analysis Type | Recommended Method | Network |
|---------------|-------------------|---------|
| TF activity (single-cell) | `run_ulm()` | CollecTRI |
| TF activity (bulk) | `run_ulm()` | CollecTRI |
| Pathway activity | `run_mlm()` | PROGENy |
| Gene list enrichment | `run_ora()` | Hallmark, GO |
| Ranked enrichment | `run_gsea()` | Hallmark, GO |
| Cell type scoring | `run_aucell()` | Custom markers |

## Pseudobulk Analysis

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Aggregate to pseudobulk
pdata = dc.pp.pseudobulk(
    adata,
    sample_col='sample',
    groups_col='cell_type',
    layer='counts',
    min_cells=10
)

# Filter genes and run DESeq2
genes = dc.pp.filter_by_expr(pdata, group='condition')
pdata = pdata[:, genes].copy()

# ... DESeq2 analysis ...

# Enrich on DE results
dc.run_ulm(mat=de_results[['stat']].T, net=tf_net, ...)
```

## Best Practices

1. **Use log-normalized data** for ULM, MLM, GSEA methods
2. **Keep raw counts** in a layer for pseudobulk analysis
3. **Check gene name overlap** between data and network
4. **Filter lowly expressed genes** before DESeq2
5. **Use consensus scores** for robust inference across methods

## Additional Resources

- **Methods Reference**: `references/methods.md` - All enrichment methods with parameters
- **Preprocessing**: `references/preprocessing.md` - Data preparation functions
- **Visualization**: `references/plotting.md` - Plotting functions
- **Databases**: `references/databases.md` - Prior knowledge networks
- **Workflows**: `references/workflows.md` - Complete examples and best practices
- **Official Documentation**: https://decoupler-py.readthedocs.io/

## Installation

```bash
pip install decoupler
pip install decoupler[all]  # With all dependencies
```
