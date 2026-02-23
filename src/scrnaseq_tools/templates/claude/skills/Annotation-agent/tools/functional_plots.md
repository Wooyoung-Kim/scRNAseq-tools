# Functional Analysis Visualization

Publication-ready plots for TF activity, Pathway activity, and Trajectory.

---

## 1. TF Activity Visualization

### Heatmap per Cluster

```python
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def plot_tf_activity_heatmap(adata, cluster_col, top_n=15, figsize=(8, 6), save_path=None):
    """
    Plot TF activity heatmap per cluster.

    Parameters:
    -----------
    adata : AnnData with ulm_estimate in obsm
    cluster_col : str, cluster column name
    top_n : int, top TFs to show
    """
    if 'ulm_estimate' not in adata.obsm:
        raise ValueError("TF activity not computed! Run decoupler first.")

    tf_activities = adata.obsm['ulm_estimate']
    tf_pvals = adata.obsm['ulm_pvals']

    # Compute mean TF activity per cluster
    clusters = adata.obs[cluster_col].unique()
    mean_activity = pd.DataFrame(index=clusters, columns=tf_activities.columns)

    for cluster in clusters:
        mask = adata.obs[cluster_col] == cluster
        mean_activity.loc[cluster] = tf_activities[mask].mean()

    # Select top variable TFs
    tf_var = mean_activity.var()
    top_tfs = tf_var.nlargest(top_n).index.tolist()

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        mean_activity[top_tfs].astype(float),
        cmap='RdBu_r',
        center=0,
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'TF Activity Score'}
    )
    ax.set_xlabel('Transcription Factor')
    ax.set_ylabel('Cluster')
    ax.set_title('TF Activity per Cluster')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if save_path:
        fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}.svg', bbox_inches='tight')
        print(f"Saved: {save_path}.png/.svg")

    plt.close(fig)
    return mean_activity[top_tfs]
```

### TF Activity Dotplot

```python
def plot_tf_dotplot(adata, cluster_col, tfs=None, top_n=10, figsize=(10, 6), save_path=None):
    """
    Plot TF activity as dotplot (like gene dotplot).
    Size = significance (-log10 pval)
    Color = activity score
    """
    if 'ulm_estimate' not in adata.obsm:
        raise ValueError("TF activity not computed!")

    tf_act = adata.obsm['ulm_estimate']
    tf_pval = adata.obsm['ulm_pvals']

    clusters = adata.obs[cluster_col].unique()

    # Select TFs
    if tfs is None:
        # Top variable TFs
        mean_act = pd.DataFrame(index=clusters)
        for cluster in clusters:
            mask = adata.obs[cluster_col] == cluster
            mean_act.loc[cluster] = tf_act[mask].mean()
        tf_var = mean_act.var()
        tfs = tf_var.nlargest(top_n).index.tolist()

    # Prepare data
    data = []
    for cluster in clusters:
        mask = adata.obs[cluster_col] == cluster
        for tf in tfs:
            activity = tf_act.loc[mask, tf].mean()
            pval = tf_pval.loc[mask, tf].mean()
            data.append({
                'Cluster': cluster,
                'TF': tf,
                'Activity': activity,
                'pval': pval,
                '-log10(pval)': -np.log10(pval + 1e-10)
            })

    df = pd.DataFrame(data)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        df['TF'].astype('category').cat.codes,
        df['Cluster'].astype('category').cat.codes,
        s=df['-log10(pval)'] * 30,  # Size = significance
        c=df['Activity'],  # Color = activity
        cmap='RdBu_r',
        vmin=-df['Activity'].abs().max(),
        vmax=df['Activity'].abs().max(),
        alpha=0.8
    )

    ax.set_xticks(range(len(tfs)))
    ax.set_xticklabels(tfs, rotation=45, ha='right')
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters)

    ax.set_xlabel('Transcription Factor')
    ax.set_ylabel('Cluster')
    ax.set_title('TF Activity Dotplot')

    cbar = plt.colorbar(scatter, ax=ax, label='Activity Score')

    # Size legend
    sizes = [1, 2, 3, 4]
    size_labels = ['1', '2', '3', '4']
    for i, (s, label) in enumerate(zip(sizes, size_labels)):
        ax.scatter([], [], s=s*30, c='gray', alpha=0.5, label=f'-log10(p)={label}')
    ax.legend(loc='upper left', bbox_to_anchor=(1.15, 1), title='Significance')

    plt.tight_layout()

    if save_path:
        fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}.svg', bbox_inches='tight')

    plt.close(fig)
```

---

## 2. Pathway Activity Visualization

### Pathway Heatmap

```python
def plot_pathway_activity_heatmap(adata, cluster_col, figsize=(8, 6), save_path=None):
    """
    Plot pathway activity heatmap per cluster.
    """
    if 'mlm_estimate' not in adata.obsm:
        raise ValueError("Pathway activity not computed!")

    pathway_act = adata.obsm['mlm_estimate']
    clusters = adata.obs[cluster_col].unique()

    # Mean pathway activity per cluster
    mean_activity = pd.DataFrame(index=clusters, columns=pathway_act.columns)
    for cluster in clusters:
        mask = adata.obs[cluster_col] == cluster
        mean_activity.loc[cluster] = pathway_act[mask].mean()

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        mean_activity.astype(float),
        cmap='RdBu_r',
        center=0,
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'Pathway Activity Score'},
        annot=True,
        fmt='.1f',
        annot_kws={'size': 8}
    )
    ax.set_xlabel('Pathway')
    ax.set_ylabel('Cluster')
    ax.set_title('Pathway Activity per Cluster')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if save_path:
        fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}.svg', bbox_inches='tight')

    plt.close(fig)
    return mean_activity
```

### Pathway Activity Barplot per Cluster

```python
def plot_pathway_bars(adata, cluster_col, cluster_id, top_n=10, figsize=(6, 4), save_path=None):
    """
    Bar plot of pathway activity for a specific cluster.
    """
    if 'mlm_estimate' not in adata.obsm:
        raise ValueError("Pathway activity not computed!")

    pathway_act = adata.obsm['mlm_estimate']
    pathway_pval = adata.obsm['mlm_pvals']

    mask = adata.obs[cluster_col] == cluster_id
    mean_act = pathway_act[mask].mean()
    mean_pval = pathway_pval[mask].mean()

    # Filter significant
    sig_mask = mean_pval < 0.05
    sig_pathways = mean_act[sig_mask].abs().nlargest(top_n)

    # Create DataFrame
    df = pd.DataFrame({
        'Pathway': sig_pathways.index,
        'Activity': mean_act[sig_pathways.index].values,
        'pval': mean_pval[sig_pathways.index].values
    })
    df = df.sort_values('Activity', ascending=True)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    colors = ['#d73027' if x > 0 else '#4575b4' for x in df['Activity']]
    bars = ax.barh(df['Pathway'], df['Activity'], color=colors, alpha=0.8)

    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('Pathway Activity Score')
    ax.set_title(f'Cluster {cluster_id} - Top Pathways')

    # Add significance stars
    for i, (idx, row) in enumerate(df.iterrows()):
        if row['pval'] < 0.001:
            star = '***'
        elif row['pval'] < 0.01:
            star = '**'
        elif row['pval'] < 0.05:
            star = '*'
        else:
            star = ''
        x_pos = row['Activity'] + 0.1 if row['Activity'] > 0 else row['Activity'] - 0.3
        ax.text(x_pos, i, star, va='center', fontsize=10)

    plt.tight_layout()

    if save_path:
        fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}.svg', bbox_inches='tight')

    plt.close(fig)
```

---

## 3. Trajectory Visualization

### Pseudotime UMAP

```python
def plot_pseudotime_umap(adata, figsize=(5, 5), save_path=None):
    """
    UMAP colored by pseudotime.
    """
    if 'pseudotime' not in adata.obs.columns:
        print("Pseudotime not computed, skipping plot")
        return

    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(
        adata,
        color='pseudotime',
        cmap='viridis',
        ax=ax,
        show=False,
        title='Pseudotime'
    )
    plt.tight_layout()

    if save_path:
        fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}.svg', bbox_inches='tight')

    plt.close(fig)
```

### Pseudotime Distribution per Cluster

```python
def plot_pseudotime_violin(adata, cluster_col, figsize=(8, 4), save_path=None):
    """
    Violin plot of pseudotime per cluster.
    """
    if 'pseudotime' not in adata.obs.columns:
        print("Pseudotime not computed, skipping plot")
        return

    fig, ax = plt.subplots(figsize=figsize)

    # Order clusters by mean pseudotime
    cluster_order = adata.obs.groupby(cluster_col)['pseudotime'].mean().sort_values().index

    sc.pl.violin(
        adata,
        keys='pseudotime',
        groupby=cluster_col,
        order=cluster_order,
        ax=ax,
        show=False,
        stripplot=False
    )

    ax.set_xlabel('Cluster')
    ax.set_ylabel('Pseudotime')
    ax.set_title('Pseudotime Distribution per Cluster')

    # Add horizontal lines for categories
    ax.axhline(0.33, color='gray', linestyle='--', alpha=0.5, label='Early/Mid')
    ax.axhline(0.67, color='gray', linestyle='--', alpha=0.5, label='Mid/Late')

    plt.tight_layout()

    if save_path:
        fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}.svg', bbox_inches='tight')

    plt.close(fig)
```

---

## 4. Combined Report Figure

### Tier 2 Combined (Markers + TF + Trajectory)

```python
def plot_tier2_combined(adata, annotation_col, marker_dict, annotation_order,
                        title_prefix, output_dir, top_tfs=10):
    """
    4-panel figure: UMAP, Dotplot, TF Heatmap, Pseudotime
    """
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, wspace=0.3, hspace=0.3)

    # 1. UMAP
    ax1 = fig.add_subplot(gs[0, 0])
    sc.pl.umap(adata, color=annotation_col, ax=ax1, show=False, title='Cell Types')

    # 2. Marker Dotplot
    ax2 = fig.add_subplot(gs[0, 1])
    markers = []
    for ct in annotation_order:
        markers.extend(marker_dict.get(ct, [])[:3])
    markers = list(dict.fromkeys(markers))[:20]  # Top 20 unique

    sc.pl.dotplot(
        adata, var_names=markers, groupby=annotation_col,
        categories_order=annotation_order, ax=ax2, show=False,
        standard_scale='var', title='Marker Expression'
    )
    plt.sca(ax2)
    plt.xticks(rotation=45, ha='right', fontsize=7)

    # 3. TF Activity Heatmap
    ax3 = fig.add_subplot(gs[1, 0])
    if 'ulm_estimate' in adata.obsm:
        tf_act = adata.obsm['ulm_estimate']
        mean_tf = pd.DataFrame(index=annotation_order)
        for ann in annotation_order:
            mask = adata.obs[annotation_col] == ann
            if mask.sum() > 0:
                mean_tf.loc[ann] = tf_act[mask].mean()
        top_tfs_list = mean_tf.var().nlargest(top_tfs).index.tolist()
        sns.heatmap(
            mean_tf[top_tfs_list].astype(float),
            cmap='RdBu_r', center=0, ax=ax3,
            xticklabels=True, yticklabels=True
        )
        ax3.set_title('TF Activity')
        plt.sca(ax3)
        plt.xticks(rotation=45, ha='right', fontsize=8)
    else:
        ax3.text(0.5, 0.5, 'TF Activity\nNot Computed', ha='center', va='center')
        ax3.axis('off')

    # 4. Pseudotime
    ax4 = fig.add_subplot(gs[1, 1])
    if 'pseudotime' in adata.obs.columns and not adata.obs['pseudotime'].isna().all():
        cluster_order = adata.obs.groupby(annotation_col)['pseudotime'].mean().reindex(annotation_order).sort_values().index
        sc.pl.violin(
            adata, keys='pseudotime', groupby=annotation_col,
            order=cluster_order, ax=ax4, show=False, stripplot=False
        )
        ax4.set_title('Pseudotime Distribution')
        ax4.axhline(0.33, color='gray', linestyle='--', alpha=0.5)
        ax4.axhline(0.67, color='gray', linestyle='--', alpha=0.5)
    else:
        ax4.text(0.5, 0.5, 'Pseudotime\nNot Computed', ha='center', va='center')
        ax4.axis('off')

    fig.suptitle(f'{title_prefix} - Tier 2 Analysis', fontsize=12, fontweight='bold')

    save_path = f'{output_dir}/{title_prefix.replace(" ", "_")}_tier2_combined'
    fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{save_path}.svg', bbox_inches='tight')
    print(f"Saved: {save_path}.png/.svg")
    plt.close(fig)
```

### Tier 3 Combined (Markers + Pathway)

```python
def plot_tier3_combined(adata, annotation_col, marker_dict, annotation_order,
                        title_prefix, output_dir):
    """
    3-panel figure: UMAP, Dotplot, Pathway Heatmap
    """
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 3, wspace=0.3)

    # 1. UMAP
    ax1 = fig.add_subplot(gs[0, 0])
    sc.pl.umap(adata, color=annotation_col, ax=ax1, show=False, title='Functional States')

    # 2. Marker Dotplot
    ax2 = fig.add_subplot(gs[0, 1])
    markers = []
    for ct in annotation_order:
        markers.extend(marker_dict.get(ct, [])[:3])
    markers = list(dict.fromkeys(markers))[:15]

    sc.pl.dotplot(
        adata, var_names=markers, groupby=annotation_col,
        categories_order=annotation_order, ax=ax2, show=False,
        standard_scale='var', title='Marker Expression'
    )
    plt.sca(ax2)
    plt.xticks(rotation=45, ha='right', fontsize=7)

    # 3. Pathway Activity
    ax3 = fig.add_subplot(gs[0, 2])
    if 'mlm_estimate' in adata.obsm:
        pathway_act = adata.obsm['mlm_estimate']
        mean_pw = pd.DataFrame(index=annotation_order)
        for ann in annotation_order:
            mask = adata.obs[annotation_col] == ann
            if mask.sum() > 0:
                mean_pw.loc[ann] = pathway_act[mask].mean()
        sns.heatmap(
            mean_pw.astype(float),
            cmap='RdBu_r', center=0, ax=ax3,
            xticklabels=True, yticklabels=True,
            annot=True, fmt='.1f', annot_kws={'size': 7}
        )
        ax3.set_title('Pathway Activity')
        plt.sca(ax3)
        plt.xticks(rotation=45, ha='right', fontsize=8)
    else:
        ax3.text(0.5, 0.5, 'Pathway Activity\nNot Computed', ha='center', va='center')
        ax3.axis('off')

    fig.suptitle(f'{title_prefix} - Tier 3 Analysis', fontsize=12, fontweight='bold')

    save_path = f'{output_dir}/{title_prefix.replace(" ", "_")}_tier3_combined'
    fig.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{save_path}.svg', bbox_inches='tight')
    print(f"Saved: {save_path}.png/.svg")
    plt.close(fig)
```

---

## 5. Output Checklist

```
Tier 2 Outputs:
- [ ] {prefix}_tier2_umap.png/.svg
- [ ] {prefix}_tier2_dotplot.png/.svg
- [ ] {prefix}_tier2_tf_heatmap.png/.svg
- [ ] {prefix}_tier2_pseudotime_umap.png/.svg (if computed)
- [ ] {prefix}_tier2_pseudotime_violin.png/.svg (if computed)
- [ ] {prefix}_tier2_combined.png/.svg

Tier 3 Outputs:
- [ ] {prefix}_tier3_umap.png/.svg
- [ ] {prefix}_tier3_dotplot.png/.svg
- [ ] {prefix}_tier3_pathway_heatmap.png/.svg
- [ ] {prefix}_tier3_combined.png/.svg
```
