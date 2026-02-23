# Phase 2: Gene Filtering

```
╔══════════════════════════════════════════════════════════════════════╗
║  GOAL: "이 유전자가 informative한가?"                                  ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Workflow

```
Step 1: Filter lowly-expressed genes
    ↓
Step 2: Identify gene categories
    ↓
Step 3: Log filtering
```

---

## Step 1: Filter Lowly-Expressed Genes

```python
import scanpy as sc
import numpy as np

# Store original gene count
n_genes_original = adata.n_vars

# Calculate genes per cell (already in var if calculate_qc_metrics was run)
if 'n_cells_by_counts' not in adata.var.columns:
    sc.pp.filter_genes(adata, min_cells=0)  # Just calculates metrics

# Check distribution
print(f"Gene expression distribution:")
print(f"  Min cells: {adata.var['n_cells_by_counts'].min()}")
print(f"  Max cells: {adata.var['n_cells_by_counts'].max()}")
print(f"  Median cells: {adata.var['n_cells_by_counts'].median():.0f}")

# Visualize
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
ax.hist(adata.var['n_cells_by_counts'], bins=100, log=True)
ax.axvline(3, color='r', linestyle='--', label='min_cells=3')
ax.set_xlabel('Number of cells expressing gene')
ax.set_ylabel('Number of genes (log scale)')
ax.legend()
plt.savefig('figures/gene_filtering.png', dpi=150, bbox_inches='tight')
plt.show()

# Filter genes expressed in fewer than min_cells
min_cells = 3  # Default: gene must be in at least 3 cells
sc.pp.filter_genes(adata, min_cells=min_cells)

n_genes_filtered = adata.n_vars
print(f"\nGene filtering:")
print(f"  Original genes: {n_genes_original}")
print(f"  Filtered genes: {n_genes_filtered}")
print(f"  Removed: {n_genes_original - n_genes_filtered} "
      f"({(n_genes_original - n_genes_filtered)/n_genes_original:.1%})")
```

---

## Step 2: Identify Gene Categories

```python
# Mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

# Hemoglobin genes (for blood samples)
adata.var['hb'] = adata.var_names.str.match('^HB[^P]')

# Immunoglobulin genes (for immune samples)
adata.var['ig'] = adata.var_names.str.match('^IG[HKL][VDJC]')

# TCR genes
adata.var['tcr'] = adata.var_names.str.match('^TR[AB][VDJC]')

# Summary
print("Gene categories:")
print(f"  Mitochondrial: {adata.var['mt'].sum()}")
print(f"  Ribosomal: {adata.var['ribo'].sum()}")
print(f"  Hemoglobin: {adata.var['hb'].sum()}")
print(f"  Immunoglobulin: {adata.var['ig'].sum()}")
print(f"  TCR: {adata.var['tcr'].sum()}")
```

---

## Step 3: Optional - Remove Specific Gene Categories

**Generally NOT recommended** - these genes can be informative

```python
# Option A: Keep all genes (RECOMMENDED)
# Mitochondrial/ribosomal can indicate cell state
# IG/TCR are critical for immune cell annotation

# Option B: Remove specific categories (if justified)
# Only do this if you have a specific reason
if remove_mt:
    adata = adata[:, ~adata.var['mt']].copy()
    print(f"Removed MT genes: {adata.n_vars} remaining")

if remove_hb:
    adata = adata[:, ~adata.var['hb']].copy()
    print(f"Removed HB genes: {adata.n_vars} remaining")
```

---

## Log Gene Filtering

```python
# Store in uns
adata.uns['gene_filtering'] = {
    'n_genes_original': n_genes_original,
    'n_genes_filtered': adata.n_vars,
    'min_cells': min_cells,
    'mt_genes': int(adata.var['mt'].sum()),
    'ribo_genes': int(adata.var['ribo'].sum())
}

print("\nGene filtering complete")
print(f"  Final genes: {adata.n_vars}")
```

---

## Gene Filtering Summary

```markdown
# Gene Filtering Log

## Parameters
- min_cells: {min_cells}
- Remove MT: No
- Remove Ribo: No

## Results
| Category | Count |
|----------|-------|
| Original | {n_original} |
| Filtered | {n_filtered} |
| Removed | {n_removed} ({pct}%) |

## Gene Categories
| Category | Count | Percentage |
|----------|-------|------------|
| Mitochondrial | {n_mt} | {pct_mt}% |
| Ribosomal | {n_ribo} | {pct_ribo}% |
| Hemoglobin | {n_hb} | {pct_hb}% |
| Immunoglobulin | {n_ig} | {pct_ig}% |
| TCR | {n_tcr} | {pct_tcr}% |
```

---

## Next Step

→ Read: phases/doublet_detection.md
