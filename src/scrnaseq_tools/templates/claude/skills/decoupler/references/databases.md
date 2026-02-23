# Decoupler Prior Knowledge Databases (dc.op)

## CollecTRI (TF-Target Network)

**Curated TF-target gene regulatory network.**

```python
net = dc.op.collectri(
    organism='human',       # 'human' or 'mouse'
    split_complexes=False   # Split TF complexes into individual TFs
)

# Returns DataFrame with columns:
#   - source: TF name
#   - target: Target gene name
#   - weight: Interaction weight (-1 or 1)
#   - PMID: PubMed ID reference
```

**Example:**
```python
tf_net = dc.op.collectri(organism='human')
print(f"TFs: {tf_net['source'].nunique()}")
print(f"Targets: {tf_net['target'].nunique()}")
print(f"Interactions: {len(tf_net)}")
```

**Best for:** TF activity inference with ULM or VIPER

---

## PROGENy (Pathway Signatures)

**Pathway responsive genes (14 cancer pathways).**

```python
model = dc.op.progeny(
    organism='human',       # 'human' or 'mouse'
    top=500                 # Top n genes per pathway (100, 200, 300, 500, 1000)
)

# Returns DataFrame with columns:
#   - source: Pathway name
#   - target: Gene name
#   - weight: Gene weight
#   - p_value: Significance
```

**Pathways included:**
- Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB
- p53, PI3K, TGFb, TNFa, Trail, VEGF, WNT

**Example:**
```python
pathway_net = dc.op.progeny(organism='human', top=500)
print(f"Pathways: {pathway_net['source'].unique()}")
```

**Best for:** Pathway activity scoring with MLM

---

## MSigDB Hallmark

**Hallmark gene sets from MSigDB.**

```python
hallmark = dc.op.hallmark(
    organism='human'        # 'human' or 'mouse'
)

# Returns DataFrame with columns:
#   - source: Gene set name (e.g., 'HALLMARK_INTERFERON_GAMMA_RESPONSE')
#   - target: Gene name
#   - weight: 1 (uniform)
```

**Example:**
```python
hallmark_net = dc.op.hallmark(organism='human')
print(f"Gene sets: {hallmark_net['source'].nunique()}")
```

**Best for:** ORA and GSEA enrichment analysis

---

## PanglaoDB (Cell Type Markers)

**Cell type marker genes.**

```python
markers = dc.op.PanglaoDB()

# Returns DataFrame with columns:
#   - source: Cell type name
#   - target: Marker gene name
#   - weight: 1 (uniform)
#   - organism: 'Hs' (human) or 'Mm' (mouse)
```

**Example:**
```python
markers = dc.op.PanglaoDB()
markers_human = markers[markers['organism'] == 'Hs']
```

**Best for:** Cell type scoring with AUCell

---

## Show Resources

**List all available resources.**

```python
dc.op.show_resources()

# Returns DataFrame with available resources and descriptions
```

---

## Generic Resource Access

**Access any OmniPath resource.**

```python
net = dc.op.resource(
    name,                   # Resource name
    organism='human',       # Organism
    **kwargs                # Additional parameters
)
```

**Examples:**
```python
dorothea = dc.op.resource('dorothea', organism='human')
tfnet = dc.op.resource('collectri', organism='mouse')
```

---

## Custom Networks

### From GMT File

```python
custom_net = dc.pp.read_gmt('my_genesets.gmt')
```

### From CSV/TSV

```python
import pandas as pd
custom_net = pd.read_csv('my_network.csv')
# Ensure columns: source, target, weight
```

### From Dictionary

```python
import pandas as pd
genesets = {
    'pathway1': ['gene1', 'gene2', 'gene3'],
    'pathway2': ['gene4', 'gene5', 'gene6']
}
custom_net = pd.DataFrame([
    {'source': pathway, 'target': gene, 'weight': 1}
    for pathway, genes in genesets.items()
    for gene in genes
])
```

---

## Network Selection Guide

| Analysis | Network | Function |
|----------|---------|----------|
| TF activity | CollecTRI | `dc.op.collectri()` |
| Pathway activity | PROGENy | `dc.op.progeny()` |
| Gene set enrichment | Hallmark, GO, KEGG | `dc.op.hallmark()` |
| Cell type scoring | PanglaoDB, custom | `dc.op.PanglaoDB()` |

---

## Checking Network Overlap

```python
# Check overlap between data and network
data_genes = set(adata.var_names)
net_genes = set(net['target'])
overlap = data_genes & net_genes
print(f"Overlap: {len(overlap)} / {len(net_genes)} network genes ({100*len(overlap)/len(net_genes):.1f}%)")

# If low overlap, check gene name format
print(f"Data genes example: {list(adata.var_names[:5])}")
print(f"Network genes example: {list(net['target'].head())}")
```
