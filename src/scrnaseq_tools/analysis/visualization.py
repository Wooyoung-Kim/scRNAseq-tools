"""
Visualization functions for QC and integration analysis.

Provides:
- QC metric distributions
- PCA variance plots
- UMAP visualizations
- Integration comparison plots
"""

import os
import logging
from typing import List, Optional, Dict, Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

from .config import PipelineConfig

logger = logging.getLogger(__name__)


def _setup_plot_style():
    """Setup matplotlib style for publication-quality plots."""
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['savefig.dpi'] = 150
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 10


def plot_qc_metrics(
    adata: sc.AnnData,
    output_path: str,
    title_prefix: str = "",
) -> str:
    """
    Generate comprehensive QC metric plots.
    
    Creates:
    - Histogram distributions for n_genes, n_counts, pct_counts_mt
    - Scatter plots showing relationships
    
    Args:
        adata: AnnData with QC metrics computed
        output_path: Path to save the plot
        title_prefix: Prefix for plot titles (e.g., "Before QC" or "After QC")
        
    Returns:
        Path to saved figure
    """
    _setup_plot_style()
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Row 1: Histograms
    ax = axes[0, 0]
    ax.hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax.set_xlabel('Genes per cell')
    ax.set_ylabel('Count')
    ax.set_title(f'{title_prefix}Genes per Cell')
    ax.axvline(adata.obs['n_genes_by_counts'].median(), color='red', linestyle='--',
               label=f"Median: {adata.obs['n_genes_by_counts'].median():.0f}")
    ax.legend()
    
    ax = axes[0, 1]
    ax.hist(adata.obs['total_counts'], bins=50, edgecolor='black', alpha=0.7, color='forestgreen')
    ax.set_xlabel('Total counts')
    ax.set_ylabel('Count')
    ax.set_title(f'{title_prefix}Total Counts per Cell')
    ax.axvline(adata.obs['total_counts'].median(), color='red', linestyle='--',
               label=f"Median: {adata.obs['total_counts'].median():.0f}")
    ax.legend()
    
    ax = axes[0, 2]
    ax.hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7, color='coral')
    ax.set_xlabel('MT %')
    ax.set_ylabel('Count')
    ax.set_title(f'{title_prefix}Mitochondrial %')
    ax.axvline(adata.obs['pct_counts_mt'].median(), color='red', linestyle='--',
               label=f"Median: {adata.obs['pct_counts_mt'].median():.1f}%")
    ax.legend()
    
    # Row 2: Scatter plots
    ax = axes[1, 0]
    ax.scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], 
               alpha=0.3, s=1, c='steelblue')
    ax.set_xlabel('Total counts')
    ax.set_ylabel('Genes per cell')
    ax.set_title('Counts vs Genes')
    
    ax = axes[1, 1]
    scatter = ax.scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'],
                        c=adata.obs['n_genes_by_counts'], alpha=0.3, s=1, cmap='viridis')
    ax.set_xlabel('Total counts')
    ax.set_ylabel('MT %')
    ax.set_title('Counts vs MT%')
    plt.colorbar(scatter, ax=ax, label='Genes')
    
    ax = axes[1, 2]
    ax.scatter(adata.obs['n_genes_by_counts'], adata.obs['pct_counts_mt'],
               alpha=0.3, s=1, c='coral')
    ax.set_xlabel('Genes per cell')
    ax.set_ylabel('MT %')
    ax.set_title('Genes vs MT%')
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"QC metrics plot saved to {output_path}")
    return output_path


def plot_pca_variance(
    adata: sc.AnnData,
    output_path: str,
) -> str:
    """
    Generate PCA variance explained plot (scree plot and cumulative).
    
    Args:
        adata: AnnData with PCA computed
        output_path: Path to save the plot
        
    Returns:
        Path to saved figure
    """
    _setup_plot_style()
    
    if 'pca' not in adata.uns or 'variance_ratio' not in adata.uns['pca']:
        logger.warning("PCA not computed. Skipping variance plot.")
        return None
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    variance_ratio = adata.uns['pca']['variance_ratio']
    n_pcs = len(variance_ratio)
    cum_var = np.cumsum(variance_ratio)
    
    # Get elbow point if available
    elbow_point = None
    n_pcs_selected = None
    if 'pc_selection' in adata.uns:
        elbow_point = adata.uns['pc_selection'].get('elbow_point')
        n_pcs_selected = adata.uns['pc_selection'].get('n_pcs_selected')
    
    # Scree plot
    ax = axes[0]
    ax.bar(range(1, n_pcs + 1), variance_ratio * 100, alpha=0.7, color='steelblue')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance Explained (%)')
    ax.set_title('Scree Plot')
    if elbow_point:
        ax.axvline(elbow_point, color='red', linestyle='--', label=f'Elbow: PC{elbow_point}')
        ax.legend()
    
    # Cumulative variance
    ax = axes[1]
    ax.plot(range(1, n_pcs + 1), cum_var * 100, 'b-', linewidth=2)
    ax.scatter(range(1, n_pcs + 1), cum_var * 100, color='steelblue', s=20)
    ax.set_xlabel('Number of Principal Components')
    ax.set_ylabel('Cumulative Variance Explained (%)')
    ax.set_title('Cumulative Variance Explained')
    ax.axhline(90, color='gray', linestyle=':', alpha=0.5, label='90% threshold')
    if n_pcs_selected and n_pcs_selected <= len(cum_var):
        ax.axvline(n_pcs_selected, color='red', linestyle='--',
                   label=f'Selected: {n_pcs_selected} PCs ({cum_var[n_pcs_selected-1]*100:.1f}%)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"PCA variance plot saved to {output_path}")
    return output_path


def plot_umap(
    adata: sc.AnnData,
    output_path: str,
    color_by: str = "leiden",
    title: str = "UMAP",
) -> str:
    """
    Generate UMAP visualization.
    
    Args:
        adata: AnnData with UMAP computed
        output_path: Path to save the plot
        color_by: Column to color points by
        title: Plot title
        
    Returns:
        Path to saved figure
    """
    _setup_plot_style()
    
    if 'X_umap' not in adata.obsm:
        logger.warning("UMAP not computed. Skipping UMAP plot.")
        return None
    
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=color_by, ax=ax, show=False, title=title)
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"UMAP plot saved to {output_path}")
    return output_path


def plot_integration_comparison(
    adata: sc.AnnData,
    output_path: str,
    batch_key: str,
    methods: List[str] = None,
    label_key: Optional[str] = None,
) -> str:
    """
    Generate side-by-side comparison of integration methods.
    
    Creates a grid of UMAP plots colored by batch (and optionally by label).
    
    Args:
        adata: AnnData with integration embeddings
        output_path: Path to save the plot
        batch_key: Column containing batch information
        methods: List of integration methods to compare
        label_key: Optional column for cell type labels
        
    Returns:
        Path to saved figure
    """
    _setup_plot_style()
    
    if methods is None:
        methods = []
        for key in ["X_scvi", "X_harmony", "X_scanorama", "X_pca"]:
            if key in adata.obsm:
                methods.append(key.replace("X_", ""))
    
    if not methods:
        logger.warning("No integration embeddings found. Skipping comparison plot.")
        return None
    
    n_methods = len(methods)
    n_rows = 2 if label_key else 1
    
    fig, axes = plt.subplots(n_rows, n_methods, figsize=(5 * n_methods, 5 * n_rows))

    # Handle different cases of axes shape from subplots
    # When n_rows=1 and n_methods=1: axes is a single Axes object
    # When n_rows=1 or n_methods=1 (but not both): axes is 1D array
    # When both > 1: axes is 2D array
    if n_rows == 1 and n_methods == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_methods == 1:
        axes = axes.reshape(-1, 1)
    
    for i, method in enumerate(methods):
        # Check for method-specific UMAP
        umap_key = f"X_umap_{method}" if f"X_umap_{method}" in adata.obsm else "X_umap"
        
        # Temporarily set the UMAP to use
        if f"X_umap_{method}" in adata.obsm:
            original_umap = adata.obsm.get("X_umap")
            adata.obsm["X_umap"] = adata.obsm[f"X_umap_{method}"]
        
        # Plot colored by batch
        ax = axes[0, i]
        sc.pl.umap(adata, color=batch_key, ax=ax, show=False, 
                   title=f"{method.upper()} - Batch")
        
        # Plot colored by label if provided
        if label_key and label_key in adata.obs.columns:
            ax = axes[1, i]
            sc.pl.umap(adata, color=label_key, ax=ax, show=False,
                       title=f"{method.upper()} - {label_key}")
        
        # Restore original UMAP
        if f"X_umap_{method}" in adata.obsm and original_umap is not None:
            adata.obsm["X_umap"] = original_umap
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Integration comparison plot saved to {output_path}")
    return output_path


def plot_benchmark_results(
    benchmark_scores: Dict[str, Dict[str, float]],
    output_path: str,
) -> str:
    """
    Generate bar chart of integration benchmark results.
    
    Args:
        benchmark_scores: Dictionary of {method: {metric: score}}
        output_path: Path to save the plot
        
    Returns:
        Path to saved figure
    """
    _setup_plot_style()
    
    if not benchmark_scores:
        logger.warning("No benchmark scores. Skipping benchmark plot.")
        return None
    
    df = pd.DataFrame(benchmark_scores).T
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(df.index))
    width = 0.8 / len(df.columns)
    
    for i, col in enumerate(df.columns):
        offset = (i - len(df.columns) / 2) * width + width / 2
        ax.bar(x + offset, df[col], width, label=col, alpha=0.8)
    
    ax.set_xlabel('Integration Method')
    ax.set_ylabel('Score')
    ax.set_title('Integration Benchmark Results')
    ax.set_xticks(x)
    ax.set_xticklabels([m.upper() for m in df.index])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.axhline(0, color='gray', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Benchmark results plot saved to {output_path}")
    return output_path


def generate_pipeline_plots(
    adata: sc.AnnData,
    output_dir: str,
    config: PipelineConfig,
) -> List[str]:
    """
    Generate all pipeline plots.
    
    Args:
        adata: AnnData with all computations done
        output_dir: Directory for output files
        config: Pipeline configuration
        
    Returns:
        List of paths to created plots
    """
    plots_created = []
    
    # QC plots
    qc_path = os.path.join(output_dir, "qc_metrics.png")
    if plot_qc_metrics(adata, qc_path, ""):
        plots_created.append(qc_path)
    
    # PCA variance plot
    pca_path = os.path.join(output_dir, "pca_variance.png")
    if plot_pca_variance(adata, pca_path):
        plots_created.append(pca_path)
    
    # UMAP by clusters
    cluster_key = 'leiden' if 'leiden' in adata.obs else 'louvain'
    if cluster_key in adata.obs:
        umap_cluster_path = os.path.join(output_dir, "umap_clusters.png")
        if plot_umap(adata, umap_cluster_path, color_by=cluster_key, title="UMAP - Clusters"):
            plots_created.append(umap_cluster_path)
    
    # UMAP by batch (if available)
    if config.batch_key and config.batch_key in adata.obs:
        umap_batch_path = os.path.join(output_dir, "umap_batch.png")
        if plot_umap(adata, umap_batch_path, color_by=config.batch_key, title="UMAP - Batch"):
            plots_created.append(umap_batch_path)
        
        # Integration comparison (if integration was run)
        methods_run = []
        for method in ["scvi", "harmony", "scanorama"]:
            if f"X_{method}" in adata.obsm:
                methods_run.append(method)
        
        if methods_run:
            comparison_path = os.path.join(output_dir, "integration_comparison.png")
            if plot_integration_comparison(adata, comparison_path, config.batch_key, methods_run):
                plots_created.append(comparison_path)
    
    logger.info(f"Generated {len(plots_created)} plots")
    return plots_created
