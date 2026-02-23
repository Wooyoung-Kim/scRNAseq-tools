# Decoupler Visualization Functions (dc.pl)

## Barplot

**Bar plot of activity scores.**

```python
dc.pl.barplot(
    adata,                  # AnnData object
    obs_keys,               # Key in adata.obs or list of sources
    var_names=None,         # Variables to plot
    groupby=None,           # Group by column
    top_n=10,               # Top n sources to show
    sort_by='mean',         # Sort by: 'mean', 'var', 'std'
    ascending=False,        # Sort order
    figsize=(10, 6),        # Figure size
    return_fig=False        # Return figure
)
```

**Examples:**
```python
# Plot top TF activities
dc.pl.barplot(adata, obs_keys='ulm_estimate', top_n=15, sort_by='mean')

# Plot specific TFs by cell type
dc.pl.barplot(
    adata,
    obs_keys='ulm_estimate',
    var_names=['STAT1', 'NF-kB', 'MYC'],
    groupby='cell_type'
)
```

---

## Dotplot

**Dot plot of activity scores.**

```python
dc.pl.dotplot(
    adata,                  # AnnData object
    obs_keys,               # Key in adata.obs or obsm
    var_names=None,         # Variables to plot
    groupby=None,           # Group by column
    cmap='RdBu_r',          # Colormap
    vmin=None,              # Min value for color
    vmax=None,              # Max value for color
    figsize=None,           # Figure size
    return_fig=False        # Return figure
)
```

**Example:**
```python
dc.pp.get_obsm(adata, obsm_key='ulm_estimate')
dc.pl.dotplot(
    adata,
    obs_keys='ulm_estimate',
    var_names=['STAT1', 'IRF1', 'NF-kB', 'MYC', 'E2F1'],
    groupby='cell_type',
    cmap='RdBu_r'
)
```

---

## Volcano Plot

**Volcano plot for differential activity.**

```python
dc.pl.volcano(
    data,                   # DataFrame with results
    x='log2FoldChange',     # Column for x-axis
    y='pvalue',             # Column for y-axis
    top_n=10,               # Label top n points
    sign_thr=0.05,          # Significance threshold
    lfc_thr=1.0,            # Log fold change threshold
    figsize=(8, 6),         # Figure size
    return_fig=False        # Return figure
)
```

**Example:**
```python
dc.pl.volcano(
    results_df,
    x='log2FoldChange',
    y='pvalue',
    top_n=10,
    sign_thr=0.05,
    lfc_thr=0.5
)
```

---

## Network Plot

**Visualize TF-target network.**

```python
dc.pl.network(
    net,                    # Network DataFrame
    source='source',        # Source column
    target='target',        # Target column
    weight='weight',        # Weight column
    source_name=None,       # Highlight specific source
    top_n=25,               # Top targets to show
    figsize=(8, 8),         # Figure size
    return_fig=False        # Return figure
)
```

**Example:**
```python
# Visualize STAT1 regulon
dc.pl.network(
    tf_net,
    source='source',
    target='target',
    weight='weight',
    source_name='STAT1',
    top_n=30
)
```

---

## Leading Edge Plot

**GSEA leading edge plot.**

```python
dc.pl.leading_edge(
    adata,                  # AnnData with GSEA results
    net,                    # Network used for GSEA
    source=None,            # Gene set name
    top_n=20,               # Top genes to show
    figsize=(10, 6),        # Figure size
    return_fig=False        # Return figure
)
```

**Example:**
```python
dc.run_gsea(adata, net=hallmark_net, source='source', target='target')
dc.pl.leading_edge(adata, net=hallmark_net, source='HALLMARK_INTERFERON_GAMMA_RESPONSE')
```

---

## Filter Plots

**Visualize filtering thresholds.**

```python
# Expression filter visualization
dc.pl.filter_by_expr(
    adata,                  # AnnData
    group=None,             # Group column
    min_count=10,           # Threshold
    figsize=(8, 4)
)

# Proportion filter visualization
dc.pl.filter_by_prop(
    adata,                  # AnnData
    min_prop=0.1,           # Threshold
    group=None,             # Group column
    figsize=(8, 4)
)

# Sample filter visualization
dc.pl.filter_samples(
    adata,                  # AnnData
    sample_col='sample',    # Sample column
    min_cells=10,           # Threshold
    figsize=(8, 4)
)
```

---

## Obsbar Plot

**Bar plot for obs annotations.**

```python
dc.pl.obsbar(
    adata,                  # AnnData
    obs_key,                # Key in adata.obs
    groupby=None,           # Group by column
    figsize=(8, 4),         # Figure size
    return_fig=False        # Return figure
)
```

---

## Obsm Plot

**Plot obsm values.**

```python
dc.pl.obsm(
    adata,                  # AnnData
    obsm_key,               # Key in adata.obsm
    var_names=None,         # Variables to plot
    groupby=None,           # Group by column
    cmap='viridis',         # Colormap
    figsize=None,           # Figure size
    return_fig=False        # Return figure
)
```

---

## Order Plots

**Plots for ordered/trajectory data.**

```python
# Plot activities along order
dc.pl.order(
    adata,                  # AnnData
    order_key,              # Ordering column (e.g., pseudotime)
    obs_keys,               # Activity keys
    var_names=None,         # Variables to plot
    smooth=True,            # Apply smoothing
    figsize=(10, 4),        # Figure size
    return_fig=False        # Return figure
)

# Plot target activities along order
dc.pl.order_targets(
    adata,                  # AnnData
    net,                    # Network
    order_key,              # Ordering column
    source_name,            # Source to plot targets for
    top_n=10,               # Top targets
    figsize=(10, 6),        # Figure size
    return_fig=False        # Return figure
)
```

---

## Source Targets Plot

**Visualize source and its targets.**

```python
dc.pl.source_targets(
    adata,                  # AnnData
    net,                    # Network
    source_name,            # Source (e.g., TF name)
    obs_key=None,           # Key with activities
    top_n=10,               # Top targets to show
    groupby=None,           # Group by column
    figsize=(10, 6),        # Figure size
    return_fig=False        # Return figure
)
```
