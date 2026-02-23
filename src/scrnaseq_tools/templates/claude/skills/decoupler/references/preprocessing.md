# Decoupler Preprocessing Functions (dc.pp)

## Pseudobulk Aggregation

**Aggregate single-cell data to pseudobulk samples.**

```python
pdata = dc.pp.pseudobulk(
    adata,                  # AnnData object
    sample_col,             # Column for sample/donor ID
    groups_col=None,        # Column for cell type/group
    obs_cols=None,          # Additional obs columns to keep
    layer=None,             # Layer with counts (default: X)
    use_raw=False,          # Use adata.raw
    mode='sum',             # Aggregation: 'sum', 'mean', 'median'
    min_cells=10,           # Minimum cells per pseudobulk
    min_counts=1000,        # Minimum total counts per pseudobulk
    dtype=np.float32,       # Output dtype
    skip_checks=False       # Skip input checks
)

# Returns AnnData with:
#   - obs: sample_col, groups_col, obs_cols, psbulk_n_cells, psbulk_counts
#   - X: Aggregated counts
#   - var: Gene information
```

**Examples:**
```python
# Simple pseudobulk by sample
pdata = dc.pp.pseudobulk(adata, sample_col='patient_id', layer='counts')

# Pseudobulk by sample and cell type
pdata = dc.pp.pseudobulk(
    adata,
    sample_col='patient_id',
    groups_col='cell_type',
    layer='counts',
    min_cells=10,
    min_counts=1000
)

# Keep additional metadata
pdata = dc.pp.pseudobulk(
    adata,
    sample_col='patient_id',
    groups_col='cell_type',
    obs_cols=['condition', 'batch', 'age'],
    layer='counts'
)
```

---

## Get Obsm

**Extract obsm values to obs for plotting.**

```python
dc.pp.get_obsm(
    adata,                  # AnnData object
    obsm_key,               # Key in adata.obsm
    obs_key=None,           # Output key in adata.obs (default: obsm_key)
    var_names=None          # Specific variables to extract (default: all)
)
```

**Examples:**
```python
# Extract all TF activities to obs
dc.pp.get_obsm(adata, obsm_key='ulm_estimate')
# Now adata.obs contains columns for each TF activity

# Extract specific TFs only
dc.pp.get_obsm(adata, obsm_key='ulm_estimate', var_names=['STAT1', 'NF-kB'])
```

---

## Filter by Expression

**Filter genes based on expression levels (for DESeq2 input).**

```python
genes = dc.pp.filter_by_expr(
    adata,                  # AnnData with counts
    group=None,             # Column for grouping
    lib_size=None,          # Library sizes (default: sum of counts)
    min_count=10,           # Minimum count per gene
    min_total_count=15,     # Minimum total count across samples
    large_n=10,             # Large n for proportion filtering
    min_prop=0.7            # Minimum proportion of samples
)

# Returns list of gene names passing filter
# Apply filter:
pdata = pdata[:, genes].copy()
```

---

## Filter by Proportion

**Filter genes by proportion of cells expressing.**

```python
genes = dc.pp.filter_by_prop(
    adata,                  # AnnData object
    min_prop=0.1,           # Minimum proportion of cells
    min_smpls=2,            # Minimum samples meeting threshold
    group=None              # Column for grouping
)

# Returns list of gene names passing filter
```

---

## Filter Samples

**Filter samples based on criteria.**

```python
adata = dc.pp.filter_samples(
    adata,                  # AnnData object
    min_cells=10,           # Minimum cells per sample
    min_counts=1000,        # Minimum counts per sample
    sample_col='sample'     # Column with sample IDs
)
```

---

## Swap Layer

**Swap adata.X with a layer.**

```python
dc.pp.swap_layer(
    adata,                  # AnnData object
    layer='counts',         # Layer to swap with X
    X_layer_key=None        # Key to store current X
)

# Example: Use counts for pseudobulk
dc.pp.swap_layer(adata, layer='counts', X_layer_key='normalized')
```

---

## Read GMT

**Read Gene Matrix Transposed (GMT) files.**

```python
net = dc.pp.read_gmt(
    path,                   # Path to GMT file
    sep='\t'                # Separator
)

# Returns DataFrame with:
#   - source: Gene set name
#   - target: Gene name
#   - weight: 1 (uniform)

# Example
hallmark = dc.pp.read_gmt('h.all.v2023.2.Hs.symbols.gmt')
```

---

## Network Operations

### Prune Network

```python
net = dc.pp.prune(
    net,                    # Network DataFrame
    min_n=5,                # Minimum targets per source
    max_n=500               # Maximum targets per source (optional)
)
```

### Adjacency Matrix

```python
adjmat = dc.pp.adjmat(
    net,                    # Network DataFrame
    source='source',
    target='target',
    weight='weight'
)
```

### Shuffle Network

```python
shuffled = dc.pp.shuffle_net(
    net,                    # Network DataFrame
    source='source',
    target='target',
    weight='weight',
    seed=42
)
```

### Network Correlation

```python
corr = dc.pp.net_corr(
    net1,                   # First network
    net2,                   # Second network
    source='source',
    target='target',
    weight='weight'
)
```

---

## Bin Order

**Bin cells by a continuous variable (e.g., pseudotime).**

```python
dc.pp.bin_order(
    adata,                  # AnnData object
    order_key,              # Column with ordering (e.g., 'dpt_pseudotime')
    bins=50,                # Number of bins
    layer=None,             # Layer to use
    use_raw=False           # Use raw
)

# Returns binned expression matrix
```

---

## KNN Smoothing

**Smooth expression using k-nearest neighbors.**

```python
dc.pp.knn(
    adata,                  # AnnData object
    k=15,                   # Number of neighbors
    use_rep='X_pca',        # Representation for neighbors
    layer=None,             # Layer to smooth
    use_raw=False           # Use raw
)

# Smooths adata.X using KNN averaging
```

---

## Extract

**Extract subset of AnnData.**

```python
sub_adata = dc.pp.extract(
    adata,                  # AnnData object
    source=None,            # Sources to extract
    obs_key=None,           # Key in obs
    obsm_key=None           # Key in obsm
)
```
