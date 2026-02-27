"""
Core pipeline for scRNA-seq QC to Clustering.

This module provides the main pipeline functions:
- load_data: Load h5ad/10x data
- run_qc: Adaptive QC with MAD-based thresholds
- normalize_and_scale: Size-factor normalization, HVG selection
- run_dimensionality_reduction: PCA + UMAP with elbow-based PC selection
- run_clustering: Leiden/Louvain clustering
- run_pipeline: Full pipeline orchestration
- run_agent_pipeline: Agent-friendly pipeline with auto-detection
"""

import os
import json
import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Union
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .config import PipelineConfig
from .stats import (
    compute_adaptive_qc_thresholds,
    apply_adaptive_qc,
    compute_qc_distributions,
    generate_qc_rationale,
    compute_pc_selection,
)

logger = logging.getLogger(__name__)


@dataclass
class PipelineResults:
    """Results from the pipeline run."""
    adata: Any = None
    qc_summary: pd.DataFrame = None
    clusters: pd.DataFrame = None
    hvg_list: List[str] = field(default_factory=list)
    run_metadata: Dict[str, Any] = field(default_factory=dict)
    output_dir: str = None
    files_created: list = field(default_factory=list)
    # Agent mode additions
    data_profile: Any = None  # DataProfile from profile module
    agent_report: Any = None  # AgentReport from report module
    doublet_result: Any = None  # DoubletResult from doublet module


def load_data(
    input_path: str,
    input_format: str = "auto"
) -> sc.AnnData:
    """
    Load single-cell data from various formats.

    Args:
        input_path: Path to input file or directory
        input_format: "h5ad", "10x", "10x_h5", or "auto" (auto-detect)

    Returns:
        AnnData object
    """
    input_path = Path(input_path)

    if input_format == "auto":
        if input_path.suffix == ".h5ad":
            input_format = "h5ad"
        elif input_path.suffix == ".h5":
            input_format = "10x_h5"
        elif input_path.is_dir():
            # Check for 10x files
            mtx_files = list(input_path.glob("*matrix*.mtx*"))
            if mtx_files:
                input_format = "10x"
            else:
                raise ValueError(f"Cannot auto-detect format for directory: {input_path}")
        else:
            raise ValueError(f"Cannot auto-detect format for: {input_path}")

    logger.info(f"Loading data from {input_path} (format: {input_format})")

    if input_format == "h5ad":
        adata = sc.read_h5ad(input_path)
    elif input_format == "10x_h5":
        # 10x CellRanger h5 format
        adata = sc.read_10x_h5(str(input_path))
    elif input_format == "10x":
        # Handle 10x CellRanger output (mtx format)
        if input_path.is_dir():
            # Check for filtered_feature_bc_matrix subdirectory
            filtered_dir = input_path / "filtered_feature_bc_matrix"
            if filtered_dir.exists():
                input_path = filtered_dir
            # Also check outs/filtered_feature_bc_matrix
            outs_filtered = input_path / "outs" / "filtered_feature_bc_matrix"
            if outs_filtered.exists():
                input_path = outs_filtered
        adata = sc.read_10x_mtx(str(input_path), var_names='gene_symbols', cache=True)
    else:
        raise ValueError(f"Unsupported format: {input_format}")

    # Ensure var_names are unique
    adata.var_names_make_unique()

    logger.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def run_qc(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> pd.DataFrame:
    """
    Run quality control filtering with statistically justified thresholds.

    Uses adaptive MAD-based thresholds when config.adaptive_qc=True:
    - Thresholds are computed as median ± nmads * MAD
    - MAD (Median Absolute Deviation) is robust to outliers
    - This approach adapts to dataset-specific distributions

    Args:
        adata: AnnData object (modified in place)
        config: Pipeline configuration
        run_metadata: Optional dict to store QC metadata

    Returns:
        DataFrame with QC summary including statistical rationale
    """
    logger.info("Running QC with statistically justified thresholds...")

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Compute distribution statistics for reporting
    distributions = compute_qc_distributions(adata)

    # Store pre-filter stats
    qc_data = {
        'metric': [],
        'value': [],
        'description': []
    }

    qc_data['metric'].append('n_cells_raw')
    qc_data['value'].append(n_cells_before)
    qc_data['description'].append('Total cells before QC')

    qc_data['metric'].append('n_genes_raw')
    qc_data['value'].append(n_genes_before)
    qc_data['description'].append('Total genes before QC')

    if config.adaptive_qc:
        # Compute adaptive thresholds using robust statistics
        thresholds = compute_adaptive_qc_thresholds(
            adata,
            nmads=config.qc_nmads,
            min_genes_floor=config.min_genes,
            max_mt_pct_ceiling=config.max_mt_pct,
            method=config.qc_method
        )

        # Apply adaptive filtering
        filter_stats = apply_adaptive_qc(adata, thresholds, config.min_cells)

        # Record thresholds
        qc_data['metric'].append('qc_method')
        qc_data['value'].append(config.qc_method)
        qc_data['description'].append('QC threshold selection method')

        qc_data['metric'].append('qc_nmads')
        qc_data['value'].append(config.qc_nmads)
        qc_data['description'].append('Number of MADs for outlier detection')

        qc_data['metric'].append('min_genes_threshold')
        qc_data['value'].append(thresholds.min_genes)
        qc_data['description'].append('Adaptive min genes threshold')

        qc_data['metric'].append('max_genes_threshold')
        qc_data['value'].append(thresholds.max_genes)
        qc_data['description'].append('Adaptive max genes threshold (doublet filter)')

        qc_data['metric'].append('max_mt_pct_threshold')
        qc_data['value'].append(thresholds.max_mt_pct)
        qc_data['description'].append('Adaptive max MT% threshold')

        # Store in run_metadata
        if run_metadata is not None:
            run_metadata['qc_thresholds'] = thresholds.to_dict()
            run_metadata['qc_distributions'] = distributions
            run_metadata['qc_rationale'] = generate_qc_rationale(thresholds, distributions)

    else:
        # Use fixed thresholds (legacy behavior)
        sc.pp.filter_cells(adata, min_genes=config.min_genes)
        sc.pp.filter_genes(adata, min_cells=config.min_cells)

        mt_mask = adata.obs['pct_counts_mt'] < config.max_mt_pct
        adata._inplace_subset_obs(mt_mask)

        qc_data['metric'].append('min_genes_threshold')
        qc_data['value'].append(config.min_genes)
        qc_data['description'].append('Fixed min genes threshold')

        qc_data['metric'].append('max_mt_pct_threshold')
        qc_data['value'].append(config.max_mt_pct)
        qc_data['description'].append('Fixed max MT% threshold')

    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars

    qc_data['metric'].append('n_cells_filtered')
    qc_data['value'].append(n_cells_after)
    qc_data['description'].append('Cells after QC')

    qc_data['metric'].append('n_genes_filtered')
    qc_data['value'].append(n_genes_after)
    qc_data['description'].append('Genes after QC')

    qc_data['metric'].append('retention_rate')
    qc_data['value'].append(n_cells_after / n_cells_before * 100)
    qc_data['description'].append('Cell retention rate (%)')

    qc_summary = pd.DataFrame(qc_data)

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['cells_before_qc'] = n_cells_before
        run_metadata['cells_after_qc'] = n_cells_after
        run_metadata['genes_before_qc'] = n_genes_before
        run_metadata['genes_after_qc'] = n_genes_after

    logger.info(f"QC complete: {n_cells_before} -> {n_cells_after} cells "
                f"({n_cells_after/n_cells_before*100:.1f}% retained), "
                f"{n_genes_before} -> {n_genes_after} genes")

    return qc_summary


def normalize_and_scale(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> List[str]:
    """
    Normalize, log-transform, and identify highly variable genes.

    Normalization approach:
    1. Size-factor normalization: counts / total_counts * target_sum
    2. Log1p transformation: log(x + 1)
    3. HVG selection using Seurat v3 method

    Args:
        adata: AnnData object (modified in place)
        config: Pipeline configuration
        run_metadata: Optional dict to store normalization metadata

    Returns:
        List of highly variable gene names
    """
    logger.info("Normalizing with size-factor correction and log1p transformation...")

    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Size-factor normalization
    sc.pp.normalize_total(adata, target_sum=config.target_sum)

    # Log1p transformation
    sc.pp.log1p(adata)

    # Store normalized data for later use
    adata.raw = adata.copy()

    # Find highly variable genes
    logger.info(f"Selecting HVGs using {config.hvg_flavor} method...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=config.n_top_genes,
        flavor=config.hvg_flavor,
        layer='counts'
    )

    # Get HVG list
    hvg_list = adata.var_names[adata.var['highly_variable']].tolist()
    n_hvg = len(hvg_list)

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['normalization_method'] = 'size_factor_log1p'
        run_metadata['target_sum'] = config.target_sum
        run_metadata['hvg_method'] = config.hvg_flavor
        run_metadata['n_hvgs'] = n_hvg

    logger.info(f"Identified {n_hvg} highly variable genes using {config.hvg_flavor} method")

    return hvg_list


def run_dimensionality_reduction(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> int:
    """
    Run PCA and UMAP with model-driven PC selection.

    Args:
        adata: AnnData object (modified in place)
        config: Pipeline configuration
        run_metadata: Optional dict to store dimensionality reduction metadata

    Returns:
        Number of PCs selected for downstream analysis
    """
    logger.info("Running dimensionality reduction with model-driven PC selection...")

    # Set random seed for reproducibility
    np.random.seed(config.seed)

    # PCA on HVGs
    logger.info(f"Computing PCA on highly variable genes (max {config.n_pcs} components)...")
    sc.tl.pca(adata, n_comps=config.n_pcs, use_highly_variable=True, random_state=config.seed)

    # Model-driven PC selection
    pc_result = compute_pc_selection(
        adata,
        method=config.pc_selection_method,
        variance_threshold=config.variance_threshold,
        min_pcs=10,
        max_pcs=config.n_pcs
    )

    n_pcs_use = pc_result.n_pcs_selected

    logger.info(f"PCA complete: {adata.obsm['X_pca'].shape}")
    logger.info(f"Selected {n_pcs_use} PCs using {config.pc_selection_method} method "
                f"(explains {pc_result.cumulative_variance[n_pcs_use-1]*100:.1f}% variance)")

    # Build kNN graph from PCA space
    logger.info(f"Building kNN graph (k={config.n_neighbors}) from {n_pcs_use} PCs...")
    sc.pp.neighbors(
        adata,
        n_neighbors=config.n_neighbors,
        n_pcs=n_pcs_use,
        random_state=config.seed
    )

    # UMAP
    logger.info("Computing UMAP embedding...")
    sc.tl.umap(adata, random_state=config.seed)

    # Store PC selection results in adata.uns
    adata.uns['pc_selection'] = {
        'method': config.pc_selection_method,
        'n_pcs_selected': n_pcs_use,
        'variance_explained': pc_result.variance_explained[:n_pcs_use],
        'cumulative_variance': pc_result.cumulative_variance[:n_pcs_use],
        'elbow_point': pc_result.elbow_point
    }

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['pc_selection'] = pc_result.to_dict()
        run_metadata['n_pcs_used'] = n_pcs_use
        run_metadata['n_neighbors'] = config.n_neighbors

    logger.info(f"UMAP complete: {adata.obsm['X_umap'].shape}")

    return n_pcs_use


def run_clustering(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> pd.DataFrame:
    """
    Run graph-based clustering on kNN graph.

    Uses community detection algorithms on the kNN graph built from PCA:
    - Leiden: Improved modularity optimization with guaranteed connectivity
    - Louvain: Classic modularity-based community detection

    Args:
        adata: AnnData with kNN graph computed
        config: Pipeline configuration
        run_metadata: Optional dict to store clustering metadata

    Returns:
        DataFrame with cluster assignments
    """
    logger.info(f"Running {config.clustering_method} clustering "
                f"(resolution={config.resolution})...")

    np.random.seed(config.seed)

    if config.clustering_method == "leiden":
        sc.tl.leiden(adata, resolution=config.resolution, random_state=config.seed)
        cluster_key = 'leiden'
    else:
        sc.tl.louvain(adata, resolution=config.resolution, random_state=config.seed)
        cluster_key = 'louvain'

    # Create clusters dataframe
    clusters_df = pd.DataFrame({
        'cell_id': adata.obs_names,
        'cluster': adata.obs[cluster_key].values
    })

    n_clusters = adata.obs[cluster_key].nunique()

    # Calculate cluster statistics
    cluster_sizes = adata.obs[cluster_key].value_counts()

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['clustering_method'] = config.clustering_method
        run_metadata['resolution'] = config.resolution
        run_metadata['n_clusters'] = n_clusters
        run_metadata['cluster_sizes'] = cluster_sizes.to_dict()

    logger.info(f"Found {n_clusters} clusters (sizes range: "
                f"{cluster_sizes.min()} - {cluster_sizes.max()})")

    return clusters_df


def run_multi_resolution_clustering(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> Dict[str, pd.DataFrame]:
    """
    Run clustering at multiple resolutions for hierarchical annotation.

    This function runs clustering at low, mid, and high resolutions,
    storing results for downstream annotation agent to use.

    Args:
        adata: AnnData with kNN graph computed
        config: Pipeline configuration with multi_resolution settings
        run_metadata: Optional dict to store clustering metadata

    Returns:
        Dict mapping resolution keys to cluster DataFrames
    """
    multi_res_config = config.multi_resolution

    if not multi_res_config.enabled:
        logger.info("Multi-resolution clustering disabled. Using single resolution.")
        return {}

    logger.info(f"Running multi-resolution clustering: {multi_res_config.resolutions}")

    np.random.seed(config.seed)

    results = {}
    cluster_info = {}

    for resolution in multi_res_config.resolutions:
        res_key = f"{multi_res_config.cluster_key_prefix}_{resolution}"

        logger.info(f"Clustering at resolution {resolution}...")

        if config.clustering_method == "leiden":
            sc.tl.leiden(adata, resolution=resolution, random_state=config.seed, key_added=res_key)
        else:
            sc.tl.louvain(adata, resolution=resolution, random_state=config.seed, key_added=res_key)

        # Convert to 1-based cluster numbering (for annotation agent)
        cluster_values = adata.obs[res_key].astype(int) + 1
        adata.obs[res_key] = cluster_values.astype(str).astype('category')

        n_clusters = adata.obs[res_key].nunique()
        cluster_sizes = adata.obs[res_key].value_counts()

        # Create clusters dataframe
        clusters_df = pd.DataFrame({
            'cell_id': adata.obs_names,
            'cluster': adata.obs[res_key].values
        })
        results[res_key] = clusters_df

        cluster_info[res_key] = {
            'resolution': resolution,
            'n_clusters': n_clusters,
            'cluster_sizes': cluster_sizes.to_dict(),
            'min_cluster_size': int(cluster_sizes.min()),
            'max_cluster_size': int(cluster_sizes.max()),
        }

        logger.info(f"  Resolution {resolution}: {n_clusters} clusters "
                   f"(sizes: {cluster_sizes.min()} - {cluster_sizes.max()})")

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['multi_resolution_clustering'] = {
            'enabled': True,
            'resolutions': multi_res_config.resolutions,
            'clustering_method': config.clustering_method,
            'cluster_info': cluster_info,
        }

    # Add convenience column: default to lowest resolution for Round 1
    lowest_res = min(multi_res_config.resolutions)
    lowest_res_key = f"{multi_res_config.cluster_key_prefix}_{lowest_res}"
    adata.obs['cluster'] = adata.obs[lowest_res_key].copy()

    logger.info(f"Default 'cluster' column set to {lowest_res_key} (resolution={lowest_res})")

    return results


def run_pipeline(
    input_path: str,
    output_dir: str,
    config: Optional[PipelineConfig] = None,
    input_format: str = "auto"
) -> PipelineResults:
    """
    Run the complete scRNA-seq QC to Clustering pipeline.

    This pipeline implements:
    1. Adaptive QC thresholds using robust statistics (MAD)
    2. Size-factor normalization with HVG selection
    3. PCA with elbow-based PC selection
    4. Leiden/Louvain clustering
    5. Optional batch integration

    Args:
        input_path: Path to input data (h5ad file or 10x directory)
        output_dir: Directory for output files
        config: Pipeline configuration (uses defaults if None)
        input_format: "h5ad", "10x", or "auto"

    Returns:
        PipelineResults with all outputs and metadata
    """
    if config is None:
        config = PipelineConfig()

    # Setup output directory
    os.makedirs(output_dir, exist_ok=True)

    results = PipelineResults(output_dir=output_dir)

    # Initialize run metadata
    run_metadata = {
        'pipeline_version': '1.0',
        'timestamp': datetime.now().isoformat(),
        'input_path': str(input_path),
        'input_format': input_format,
        'seed': config.seed,
    }

    logger.info("=" * 60)
    logger.info("scrna_qc: QC to Integration Pipeline")
    logger.info("=" * 60)
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 60)

    # 1. Load data
    adata = load_data(input_path, input_format)
    run_metadata['input_cells'] = adata.n_obs
    run_metadata['input_genes'] = adata.n_vars

    # 2. QC with adaptive thresholds
    qc_summary = run_qc(adata, config, run_metadata)
    qc_path = os.path.join(output_dir, "qc_summary.csv")
    qc_summary.to_csv(qc_path, index=False)
    results.qc_summary = qc_summary
    results.files_created.append(qc_path)

    # 3. Normalize and scale (with HVG selection)
    hvg_list = normalize_and_scale(adata, config, run_metadata)
    results.hvg_list = hvg_list

    # Save HVG list
    hvg_path = os.path.join(output_dir, "hvg_list.txt")
    with open(hvg_path, 'w') as f:
        f.write('\n'.join(hvg_list))
    results.files_created.append(hvg_path)

    # 4. Dimensionality reduction
    n_pcs_used = run_dimensionality_reduction(adata, config, run_metadata)

    # 5. Clustering (single resolution)
    clusters = run_clustering(adata, config, run_metadata)
    clusters_path = os.path.join(output_dir, "clusters.csv")
    clusters.to_csv(clusters_path, index=False)
    results.clusters = clusters
    results.files_created.append(clusters_path)

    # 5.1 Multi-resolution clustering (for annotation agent)
    if config.multi_resolution.enabled:
        multi_res_clusters = run_multi_resolution_clustering(adata, config, run_metadata)
        for res_key, res_df in multi_res_clusters.items():
            res_path = os.path.join(output_dir, f"{res_key}.csv")
            res_df.to_csv(res_path, index=False)
            results.files_created.append(res_path)

    # 6. Integration (if batch_key is specified)
    if config.batch_key is not None and config.batch_key in adata.obs.columns:
        logger.info(f"Batch key '{config.batch_key}' found. Running integration...")
        from .integration import run_integration
        
        integration_result = run_integration(
            adata,
            batch_key=config.batch_key,
            methods=config.integration.get_methods_list(),
            config=config.integration,
            output_dir=output_dir,
        )
        
        run_metadata['integration'] = {
            'methods': integration_result.methods_run,
            'best_method': integration_result.best_method,
            'benchmark_scores': integration_result.benchmark_scores,
        }
        results.files_created.extend(integration_result.files_created)

    # 7. Generate plots (if enabled)
    if config.generate_plots:
        from .visualization import generate_pipeline_plots
        plot_paths = generate_pipeline_plots(adata, output_dir, config)
        results.files_created.extend(plot_paths)

    # 8. Save processed data
    h5ad_path = os.path.join(output_dir, "processed.h5ad")
    adata.write_h5ad(h5ad_path)
    results.files_created.append(h5ad_path)

    # 9. Save run metadata
    run_metadata['summary'] = {
        'cells_analyzed': adata.n_obs,
        'genes_analyzed': adata.n_vars,
        'n_hvgs': len(hvg_list),
        'n_pcs_used': n_pcs_used,
        'n_clusters': run_metadata.get('n_clusters', 0),
    }

    # Clean metadata for JSON serialization
    run_metadata_clean = {}
    for key, value in run_metadata.items():
        if isinstance(value, (str, int, float, bool, type(None))):
            run_metadata_clean[key] = value
        elif isinstance(value, dict):
            run_metadata_clean[key] = {
                k: v for k, v in value.items()
                if isinstance(v, (str, int, float, bool, type(None), list, dict))
            }
        elif isinstance(value, list):
            run_metadata_clean[key] = value

    metadata_path = os.path.join(output_dir, "run_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(run_metadata_clean, f, indent=2, default=str)
    results.files_created.append(metadata_path)
    results.run_metadata = run_metadata

    results.adata = adata

    logger.info("=" * 60)
    logger.info("Pipeline complete!")
    logger.info(f"Files created: {len(results.files_created)}")
    logger.info("=" * 60)

    return results


def run_agent_pipeline(
    input_path: str,
    output_dir: str,
    config: Optional[PipelineConfig] = None,
    dry_run: bool = False,
) -> PipelineResults:
    """
    Run the Agent-friendly scRNA-seq pipeline with automatic detection.

    This pipeline is designed for AI agents and includes:
    1. Data profiling and format detection
    2. Smart loading (single or multi-sample)
    3. Automatic batch key detection
    4. Automatic MT gene pattern detection
    5. Doublet detection
    6. Adaptive QC
    7. Normalization and HVG selection
    8. Dimensionality reduction
    9. Clustering
    10. Optional integration
    11. Structured report generation

    Args:
        input_path: Path to input data (file or directory)
        output_dir: Directory for output files
        config: Pipeline configuration (uses defaults if None)
        dry_run: If True, only profile data without processing

    Returns:
        PipelineResults with all outputs, profile, and agent report
    """
    from .profile import profile_data, DataProfile
    from .loader import smart_load, MultiSampleConfig
    from .doublet import run_doublet_detection, filter_doublets, DoubletConfig as DoubletDetectionConfig
    from .report import generate_agent_report, save_agent_report

    if config is None:
        config = PipelineConfig()

    # Setup output directory
    os.makedirs(output_dir, exist_ok=True)

    results = PipelineResults(output_dir=output_dir)

    # Initialize run metadata
    run_metadata = {
        'pipeline_version': '1.0',
        'agent_mode': True,
        'timestamp': datetime.now().isoformat(),
        'input_path': str(input_path),
        'seed': config.seed,
    }

    logger.info("=" * 60)
    logger.info("scrna_qc: Agent-Friendly Pipeline")
    logger.info("=" * 60)
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 60)

    # ==================== Phase 1: Data Profiling ====================
    logger.info("Phase 1: Data Profiling...")

    # First, do a quick profile to understand the data
    data_profile = profile_data(input_path)
    results.data_profile = data_profile

    # Save profile
    profile_path = os.path.join(output_dir, "data_profile.json")
    with open(profile_path, 'w') as f:
        json.dump(data_profile.to_dict(), f, indent=2, default=str)
    results.files_created.append(profile_path)

    logger.info(f"Profile: {data_profile.n_cells:,} cells x {data_profile.n_genes:,} genes")
    logger.info(f"Format: {data_profile.file_detection.format if data_profile.file_detection else 'unknown'}")

    if dry_run:
        logger.info("Dry run mode - stopping after profiling")
        logger.info("\n" + data_profile.summary())
        results.run_metadata = run_metadata
        return results

    # ==================== Phase 2: Smart Loading ====================
    logger.info("Phase 2: Smart Loading...")

    load_config = MultiSampleConfig(
        sample_key="sample",
        min_cells_per_sample=10,
    )
    load_result = smart_load(input_path, config=load_config, profile_first=False)

    if load_result.load_errors:
        for error in load_result.load_errors:
            logger.error(f"Load error: {error}")
        raise RuntimeError(f"Failed to load data: {load_result.load_errors}")

    adata = load_result.adata
    run_metadata['load_result'] = load_result.to_dict()
    run_metadata['input_cells'] = adata.n_obs
    run_metadata['input_genes'] = adata.n_vars

    logger.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    if load_result.n_samples > 1:
        logger.info(f"Merged {load_result.n_samples} samples")

    # ==================== Phase 3: Auto-detection ====================
    logger.info("Phase 3: Auto-detection...")

    # Auto-detect batch key if enabled and not specified
    if config.agent.auto_detect_batch_key and config.batch_key is None:
        if data_profile.batch_key_detection and data_profile.batch_key_detection.detected_key:
            config.batch_key = data_profile.batch_key_detection.detected_key
            logger.info(f"Auto-detected batch key: {config.batch_key}")
            run_metadata['auto_detected_batch_key'] = config.batch_key
        elif load_result.n_samples > 1:
            config.batch_key = load_result.sample_key
            logger.info(f"Using sample key as batch key: {config.batch_key}")
            run_metadata['auto_detected_batch_key'] = config.batch_key

    # Auto-detect MT gene pattern
    if config.agent.auto_detect_mt_pattern:
        if data_profile.mt_gene_detection and data_profile.mt_gene_detection.detected:
            mt_pattern = data_profile.mt_gene_detection.pattern
            logger.info(f"Auto-detected MT pattern: {mt_pattern}")
            run_metadata['auto_detected_mt_pattern'] = mt_pattern

    # ==================== Phase 4: Doublet Detection ====================
    if config.doublet.enabled:
        logger.info("Phase 4: Doublet Detection...")

        doublet_config = DoubletDetectionConfig(
            method=config.doublet.method,
            expected_doublet_rate=config.doublet.expected_doublet_rate,
            n_prin_comps=config.doublet.n_prin_comps,
            random_state=config.seed,
        )

        try:
            doublet_result = run_doublet_detection(
                adata,
                batch_key=config.batch_key,
                config=doublet_config,
                inplace=True,
            )
            results.doublet_result = doublet_result
            run_metadata['doublet_detection'] = doublet_result.to_dict()

            logger.info(f"Doublets detected: {doublet_result.n_doublets:,} "
                       f"({doublet_result.doublet_rate:.1%})")

            if config.doublet.filter_doublets:
                adata = filter_doublets(adata, inplace=False)
                logger.info(f"After doublet removal: {adata.n_obs:,} cells")

        except ImportError:
            logger.warning("Scrublet not installed. Skipping doublet detection.")
            logger.warning("Install with: pip install scrublet")
        except Exception as e:
            logger.warning(f"Doublet detection failed: {e}")

    # ==================== Phase 5: QC ====================
    logger.info("Phase 5: Quality Control...")

    qc_summary = run_qc(adata, config, run_metadata)
    qc_path = os.path.join(output_dir, "qc_summary.csv")
    qc_summary.to_csv(qc_path, index=False)
    results.qc_summary = qc_summary
    results.files_created.append(qc_path)

    # ==================== Phase 6: Normalization ====================
    logger.info("Phase 6: Normalization...")

    hvg_list = normalize_and_scale(adata, config, run_metadata)
    results.hvg_list = hvg_list

    hvg_path = os.path.join(output_dir, "hvg_list.txt")
    with open(hvg_path, 'w') as f:
        f.write('\n'.join(hvg_list))
    results.files_created.append(hvg_path)

    # ==================== Phase 7: Dimensionality Reduction ====================
    logger.info("Phase 7: Dimensionality Reduction...")

    n_pcs_used = run_dimensionality_reduction(adata, config, run_metadata)

    # ==================== Phase 8: Clustering ====================
    logger.info("Phase 8: Clustering...")

    clusters = run_clustering(adata, config, run_metadata)
    clusters_path = os.path.join(output_dir, "clusters.csv")
    clusters.to_csv(clusters_path, index=False)
    results.clusters = clusters
    results.files_created.append(clusters_path)

    # ==================== Phase 9: Integration (if applicable) ====================
    if config.batch_key is not None and config.batch_key in adata.obs.columns:
        n_batches = adata.obs[config.batch_key].nunique()
        if n_batches > 1:
            logger.info(f"Phase 9: Integration ({n_batches} batches)...")
            from .integration import run_integration

            try:
                integration_result = run_integration(
                    adata,
                    batch_key=config.batch_key,
                    methods=config.integration.get_methods_list(),
                    config=config.integration,
                    output_dir=output_dir,
                )

                run_metadata['integration'] = {
                    'methods': integration_result.methods_run,
                    'best_method': integration_result.best_method,
                    'benchmark_scores': integration_result.benchmark_scores,
                }
                run_metadata['batch_key'] = config.batch_key
                results.files_created.extend(integration_result.files_created)

                # Re-run multi-resolution clustering on integrated embedding
                if config.multi_resolution.enabled:
                    logger.info("Phase 9.1: Multi-resolution clustering on integrated data...")
                    multi_res_clusters = run_multi_resolution_clustering(adata, config, run_metadata)
                    for res_key, res_df in multi_res_clusters.items():
                        res_path = os.path.join(output_dir, f"{res_key}.csv")
                        res_df.to_csv(res_path, index=False)
                        results.files_created.append(res_path)

            except Exception as e:
                logger.warning(f"Integration failed: {e}")

    # Multi-resolution clustering (if not done after integration)
    if config.multi_resolution.enabled and 'multi_resolution_clustering' not in run_metadata:
        logger.info("Phase 9.1: Multi-resolution clustering...")
        multi_res_clusters = run_multi_resolution_clustering(adata, config, run_metadata)
        for res_key, res_df in multi_res_clusters.items():
            res_path = os.path.join(output_dir, f"{res_key}.csv")
            res_df.to_csv(res_path, index=False)
            results.files_created.append(res_path)

    # ==================== Phase 10: Visualization ====================
    if config.generate_plots:
        logger.info("Phase 10: Generating Plots...")
        from .visualization import generate_pipeline_plots
        plot_paths = generate_pipeline_plots(adata, output_dir, config)
        results.files_created.extend(plot_paths)

    # ==================== Phase 11: Save Data ====================
    logger.info("Phase 11: Saving Processed Data...")

    h5ad_path = os.path.join(output_dir, "processed.h5ad")
    adata.write_h5ad(h5ad_path)
    results.files_created.append(h5ad_path)

    # ==================== Phase 12: Agent Report ====================
    if config.agent.generate_structured_report:
        logger.info("Phase 12: Generating Agent Report...")

        run_metadata['summary'] = {
            'cells_analyzed': adata.n_obs,
            'genes_analyzed': adata.n_vars,
            'n_hvgs': len(hvg_list),
            'n_pcs_used': n_pcs_used,
            'n_clusters': run_metadata.get('n_clusters', 0),
        }

        agent_report = generate_agent_report(
            adata,
            run_metadata,
            input_path=str(input_path),
            output_dir=output_dir,
        )
        agent_report.files_created = results.files_created.copy()
        results.agent_report = agent_report

        report_paths = save_agent_report(
            agent_report,
            output_dir,
            formats=config.agent.report_formats,
        )
        results.files_created.extend(report_paths)

    # Save run metadata
    run_metadata_clean = _clean_metadata_for_json(run_metadata)
    metadata_path = os.path.join(output_dir, "run_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(run_metadata_clean, f, indent=2, default=str)
    results.files_created.append(metadata_path)
    results.run_metadata = run_metadata

    results.adata = adata

    logger.info("=" * 60)
    logger.info("Agent Pipeline Complete!")
    logger.info(f"Files created: {len(results.files_created)}")
    if results.agent_report:
        logger.info(f"Quality score: {results.agent_report.data_quality_score:.2f}/1.0")
    logger.info("=" * 60)

    return results


def _clean_metadata_for_json(run_metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Clean metadata for JSON serialization."""
    run_metadata_clean = {}
    for key, value in run_metadata.items():
        if isinstance(value, (str, int, float, bool, type(None))):
            run_metadata_clean[key] = value
        elif isinstance(value, dict):
            run_metadata_clean[key] = {
                k: v for k, v in value.items()
                if isinstance(v, (str, int, float, bool, type(None), list, dict))
            }
        elif isinstance(value, list):
            run_metadata_clean[key] = value
        elif isinstance(value, np.ndarray):
            run_metadata_clean[key] = value.tolist()
    return run_metadata_clean
