# Compact Working Examples

> **WARNING: T cell-specific example only.** All marker genes (CD3D, GZMB, PDCD1, etc.), thresholds, and cell type labels below are illustrative for T cell lineage. Do NOT apply these markers or labels to other lineages (B cell, myeloid, etc.). Always derive markers from DE results of the current dataset. PMIDs in this file are placeholders.

Minimal examples demonstrating the hierarchical annotation workflow.

---

## Example 1: Tier 1 Annotation (T cells)

### Input
```
Cluster 0, Top 10 markers:
CD3D (pct=0.95, LFC=3.2, padj=1e-100)
CD3E (pct=0.92, LFC=2.8, padj=1e-95)
TRAC (pct=0.88, LFC=2.5, padj=1e-90)
CD2 (pct=0.85, LFC=2.1, padj=1e-80)
IL7R (pct=0.45, LFC=1.5, padj=1e-50)
...
```

### Reasoning (Compact)
```
## Cluster 0 Annotation

### Evidence
- Top markers: CD3D, CD3E, TRAC, CD2, IL7R
- All T cell lineage markers

### Candidates
1. T cells: CD3D + CD3E + TRAC (classic T cell signature)

### Literature
- Query: "CD3D CD3E TRAC T cell" -> PMID:XXXXXXXX (placeholder) (T cell markers)
- Verified: DOUBLE_VERIFIED

### Decision
- **Annotation**: T cells
- **Markers**: CD3D + CD3E + TRAC
- **Confidence**: High
- **References**: PMID:XXXXXXXX (placeholder)
```

---

## Example 2: Tier 2 Annotation (T cell -> Effector)

### Input (After subsetting to T cells)
```
Cluster 2, Top 10 markers (WITHIN T cells):
GZMB (pct=0.85, LFC=2.5, padj=1e-50)
PRF1 (pct=0.78, LFC=2.1, padj=1e-45)
NKG7 (pct=0.72, LFC=1.9, padj=1e-40)
GNLY (pct=0.68, LFC=1.8, padj=1e-38)
IFNG (pct=0.55, LFC=1.5, padj=1e-30)

TF Activity: TBX21 (1.8), EOMES (1.2)
Pseudotime: 0.75 (late in trajectory)
```

### Reasoning (Compact)
```
## T cells -> Cluster 2 Annotation

### Evidence
- Top markers: GZMB, PRF1, NKG7, GNLY, IFNG (cytotoxic signature)
- TF: TBX21 (effector TF), EOMES (effector TF)
- Trajectory: Late pseudotime (differentiated state)

### Candidates
1. Effector: GZMB + PRF1 + NKG7 (cytotoxic effector)
2. Activated: IFNG + NKG7 (early activation)

### Literature
- Query: "GZMB PRF1 NKG7 effector T cell" -> PMID:XXXXXXXX (placeholder) (effector signature)
- Query: "TBX21 EOMES cytotoxic" -> PMID:XXXXXXXX (placeholder) (TF confirms)

### Reasoning
- Effector: [+] High cytotoxic markers, [+] TBX21/EOMES high, [+] Late pseudotime
- Activated: [-] Too many cytotoxic markers for early activation

### Decision
- **Annotation**: Effector
- **Markers**: GZMB + PRF1 + NKG7
- **Confidence**: High
- **References**: PMID:XXXXXXXX, PMID:XXXXXXXX (placeholders)
```

---

## Example 3: Tier 3 Annotation (Effector -> Exhausted)

### Input (After subsetting to Effector)
```
Cluster 1, Top 10 markers (WITHIN Effector):
PDCD1 (pct=0.65, LFC=1.8, padj=1e-35)
LAG3 (pct=0.55, LFC=1.5, padj=1e-28)
HAVCR2 (pct=0.48, LFC=1.3, padj=1e-22)
TOX (pct=0.45, LFC=1.4, padj=1e-20)
TIGIT (pct=0.42, LFC=1.2, padj=1e-18)

Pathway Activity: Apoptosis (1.5), Hypoxia (1.2)
```

### Reasoning (Compact)
```
## T cells -> Effector -> Cluster 1 Annotation

### Evidence
- Top markers: PDCD1, LAG3, HAVCR2, TOX, TIGIT (exhaustion markers)
- Pathway: Apoptosis elevated (stress signature)
- Note: Still has GZMB (from parent) but NOT differentially expressed HERE

### Candidates
1. Exhausted: PDCD1 + LAG3 + TOX (terminal exhaustion)
2. Tpex: PDCD1 + TCF7 + TOX (progenitor exhausted)

### Literature
- Query: "PDCD1 LAG3 TOX exhausted T cell" -> PMID:XXXXXXXX (placeholder)
- Query: "HAVCR2 TIGIT dysfunction" -> PMID:XXXXXXXX (placeholder)

### Reasoning
- Exhausted: [+] All exhaustion markers high, [+] No TCF7 (not progenitor)
- Tpex: [-] TCF7 not in top markers (not stem-like)

### Decision
- **Annotation**: Exhausted
- **Markers**: PDCD1 + LAG3 + TOX + HAVCR2
- **Confidence**: High
- **References**: PMID:XXXXXXXX, PMID:XXXXXXXX (placeholders)
- **Final Label**: T cells_Effector_Exhausted
```

---

## Example 4: Novel Population Detection

### Input (After subsetting to Effector)
```
Cluster 3, Top 10 markers (WITHIN Effector):
GZMB (pct=0.82, LFC=2.0, padj=1e-40)
TCF7 (pct=0.55, LFC=1.6, padj=1e-30)
IL7R (pct=0.48, LFC=1.4, padj=1e-25)
PDCD1 (pct=0.35, LFC=1.1, padj=1e-15)
MKI67 (pct=0.30, LFC=1.2, padj=1e-12)

TF Activity: TBX21 (1.5), TCF7 (1.8), TOX (0.5)
Pseudotime: 0.45 (intermediate)
```

### Reasoning (Compact)
```
## T cells -> Effector -> Cluster 3 Annotation

### Evidence
- Top markers: GZMB, TCF7, IL7R, PDCD1, MKI67
- Unusual: Cytotoxic (GZMB) + Stem-like (TCF7, IL7R) + Exhaustion (PDCD1)
- Trajectory: Intermediate pseudotime

### Candidates
1. Cytotoxic: GZMB + ... -> but TCF7/IL7R unusual
2. Stem-like: TCF7 + IL7R -> but GZMB too high
3. Tpex: PDCD1 + TCF7 -> but GZMB very high, TOX low

### Literature (PRE)
- Query: "GZMB TCF7 IL7R T cell" -> 0 relevant hits
- Query: "cytotoxic stem-like T cell" -> PMID:XXXXXXXX (placeholder) (not exact match)
- Query: "GZMB TCF7 PDCD1 single cell" -> 0 relevant hits

### Novel Check
- [x] No literature match (3 queries failed)
- [x] Opposing markers: Cytotoxic (GZMB) + Stem (TCF7, IL7R)
- [x] Trajectory: Intermediate position
- [x] Unique TF: TBX21 + TCF7 high together (unusual)

### Novel Requirements
- [x] >= 3 distinguishing markers: GZMB, TCF7, IL7R
- [x] Unique signature: TBX21+TCF7 combination
- [x] >= 50 cells: 127 cells in cluster
- [x] 3 failed queries documented

### Decision
- **Annotation**: Novel_GZMB+TCF7+IL7R
- **Markers**: GZMB + TCF7 + IL7R
- **Confidence**: Medium (novel, needs validation)
- **Notes**: Potentially transitional population with cytotoxic capacity
  and stem-like properties. May represent tumor-specific population.
- **Final Label**: T cells_Effector_Novel_GZMB+TCF7+IL7R
```

---

## Example 5: Low Confidence Case

### Input
```
Cluster 5, Top 10 markers:
RPS12 (pct=0.95, LFC=1.1, padj=1e-20)
RPL13 (pct=0.93, LFC=1.0, padj=1e-18)
MALAT1 (pct=0.90, LFC=0.9, padj=1e-15)
MT-CO1 (pct=0.88, LFC=0.8, padj=0.01)
...

Cell count: 45 cells
Doublet score: 0.35 (elevated)
```

### Reasoning (Compact)
```
## Cluster 5 Annotation

### Evidence
- Top markers: Ribosomal + Mitochondrial genes
- Only 45 cells (< 50 threshold)
- Elevated doublet score

### Issues
- No lineage-specific markers meeting criteria
- Technical artifacts dominate
- Small cluster size
- High doublet score

### Decision
- **Annotation**: Low_Quality_Unassigned
- **Confidence**: Low
- **Action**: Flag for review, consider removing from analysis
- **Notes**: Likely doublets or stressed cells
```

---

## Code Example: Full Tier Workflow

```python
import scanpy as sc

def annotate_tier1(adata):
    """Example Tier 1 workflow."""
    # Cluster
    sc.tl.leiden(adata, resolution=0.4, key_added='tier1_cluster')

    # DE
    sc.tl.rank_genes_groups(adata, groupby='tier1_cluster', method='wilcoxon')

    # For each cluster
    results = {}
    for cluster in adata.obs['tier1_cluster'].unique():
        # Get markers
        de_df = sc.get.rank_genes_groups_df(adata, group=cluster)
        valid = de_df[
            (de_df['pct_nz_group'] >= 0.25) &
            (de_df['logfoldchanges'] >= 1.0) &
            (de_df['pvals_adj'] < 0.05)
        ].head(50)

        # [Literature search + Reasoning would happen here]

        # Store result
        results[cluster] = {
            'annotation': 'T cells',  # Example
            'markers': ['CD3D', 'CD3E', 'TRAC'],
            'confidence': 'High',
            'references': ['PMID:XXXXXXXX']  # placeholder
        }

    # Apply annotations
    mapping = {k: v['annotation'] for k, v in results.items()}
    adata.obs['tier1_annotation'] = adata.obs['tier1_cluster'].map(mapping)

    return adata, results
```
