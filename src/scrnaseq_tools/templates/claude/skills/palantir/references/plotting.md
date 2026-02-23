# Palantir Plotting Functions (palantir.plot)

## Main Results

### plot_palantir_results

**Main visualization of Palantir results.**

```python
fig = palantir.plot.plot_palantir_results(
    data,                    # AnnData or DataFrame
    pr_res=None,             # PResults (optional if AnnData)
    embedding_basis='X_umap',
    pseudo_time_key='palantir_pseudotime',
    entropy_key='palantir_entropy',
    fate_prob_key='palantir_fate_probabilities',
    **kwargs
)

# Creates multi-panel figure with:
#   - Pseudotime on embedding
#   - Entropy/differentiation potential
#   - Fate probabilities for each terminal state
```

---

## Trajectory Visualization

### plot_trajectory

**Plot single trajectory with arrows.**

```python
ax = palantir.plot.plot_trajectory(
    ad,                      # AnnData
    branch,                  # Branch name
    ax=None,                 # Matplotlib axes
    pseudo_time_key='palantir_pseudotime',
    masks_key='branch_masks',
    embedding_basis='X_umap',
    cell_color='branch_selection',
    smoothness=1.0,          # Arrow smoothness
    pseudotime_interval=None,  # (start, end) range
    n_arrows=5,              # Number of direction arrows
    arrowprops={},           # Arrow properties
    scanpy_kwargs={},        # scanpy plotting args
    figsize=(5, 5),
    **kwargs
)
```

### plot_trajectories

**Plot multiple trajectories simultaneously.**

```python
ax = palantir.plot.plot_trajectories(
    ad,                      # AnnData
    groups=None,             # Branch names (default: all)
    pseudo_time_key='palantir_pseudotime',
    masks_key='branch_masks',
    embedding_basis='X_umap',
    cell_color='palantir_pseudotime',
    palantir_fates_colors=None,  # Color dict/list
    smoothness=1.0,
    pseudotime_interval=None,
    n_arrows=5,
    ax=None,
    figsize=(5, 5),
    show_legend=True,
    **kwargs
)
```

**Example:**
```python
palantir.plot.plot_trajectories(
    adata,
    pseudotime_interval=(0, 0.9),
    n_arrows=5
)
```

---

## Gene Trends

### plot_gene_trends

**Plot gene expression trends across branches.**

```python
fig = palantir.plot.plot_gene_trends(
    data,                    # Dict or AnnData
    genes=None,              # Gene list
    gene_trend_key='gene_trends',
    branch_names='branch_masks'
)

# Creates grid: rows=genes, columns=branches
```

### plot_trend

**Plot single gene trend with scatter overlay.**

```python
fig, ax = palantir.plot.plot_trend(
    ad,                      # AnnData
    branch_name,             # Target branch
    gene,                    # Gene name
    color=None,              # Color by variable
    masks_key='branch_masks',
    gene_trend_key='gene_trends',
    pseudo_time_key='palantir_pseudotime',
    figsize=(12, 4),
    **kwargs
)
```

### plot_gene_trend_heatmaps

**Heatmap of gene trends.**

```python
fig = palantir.plot.plot_gene_trend_heatmaps(
    data,                    # AnnData or Dict
    genes=None,              # Gene list
    gene_trend_key='gene_trends',
    branch_names='branch_masks',
    scaling='z-score',       # 'none', 'z-score', 'quantile', 'percent'
    basefigsize=(7, 0.7),
    **kwargs
)
```

### plot_gene_trend_clusters

**Plot clustered gene trends.**

```python
fig = palantir.plot.plot_gene_trend_clusters(
    data,                    # AnnData or DataFrame
    branch_name='',          # Target branch
    clusters=None,           # Cluster labels
    gene_trend_key='gene_trends'
)
```

---

## Branch and Cell Selection

### plot_branch_selection

**Visualize branch cell selection.**

```python
fig = palantir.plot.plot_branch_selection(
    ad,                      # AnnData
    pseudo_time_key='palantir_pseudotime',
    fate_prob_key='palantir_fate_probabilities',
    masks_key='branch_masks',
    fates=None,              # Subset of fates
    embedding_basis='X_umap',
    figsize=(15, 5),
    **kwargs
)
```

### plot_branch

**Scatter plot of cells along pseudotime.**

```python
palantir.plot.plot_branch(
    ad,                      # AnnData
    branch_name,             # Target branch
    position,                # Y-axis variable
    color=None,              # Color variable
    masks_key='branch_masks',
    pseudo_time_key='palantir_pseudotime',
    figsize=(12, 4),
    **kwargs
)

# X-axis: pseudotime, Y-axis: position variable
```

---

## Other Plots

### plot_terminal_state_probs

**Plot fate probability bars for specific cells.**

```python
fig = palantir.plot.plot_terminal_state_probs(
    data,                    # AnnData or DataFrame
    cells,                   # Cell barcodes
    **kwargs
)
```

### plot_diffusion_components

**Visualize diffusion components on embedding.**

```python
fig, ax = palantir.plot.plot_diffusion_components(
    data,                    # AnnData or DataFrame
    dm_res='DM_EigenVectors',
    embedding_basis='X_umap',
    **kwargs
)
```

### highlight_cells_on_umap

**Highlight specific cells on UMAP.**

```python
palantir.plot.highlight_cells_on_umap(
    data,                    # AnnData or DataFrame
    cells,                   # Cell barcodes to highlight
    **kwargs
)
```

### plot_stats

**Generic scatter plot function.**

```python
fig, ax = palantir.plot.plot_stats(
    ad,                      # AnnData
    x,                       # X-axis variable
    y,                       # Y-axis variable
    color=None,              # Color variable
    branch_name=None,
    masks_key='branch_masks',
    figsize=(6, 6),
    **kwargs
)
```

---

## FigureGrid Helper

**Helper class for subplot grids.**

```python
fg = palantir.plot.FigureGrid(n_rows, n_cols, figsize=None)

for ax in fg:
    ax.plot(...)

fig = fg.figure
```
