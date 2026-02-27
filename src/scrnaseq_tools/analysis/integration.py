"""
Integration methods for scRNA-seq batch correction.

Supports multiple integration methods:
- scVI: Variational Autoencoder-based integration
- Harmony: Fast linear batch correction in PCA space
- Scanorama: Mutual Nearest Neighbors-based integration

Features:
- Run single or multiple methods
- Automatic benchmarking when multiple methods are run
- Auto-selection of best integration method
"""

import os
import logging
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scanpy as sc

from .config import IntegrationConfig

logger = logging.getLogger(__name__)


@dataclass
class IntegrationResult:
    """Results from integration."""
    methods_run: List[str] = field(default_factory=list)
    embeddings: Dict[str, str] = field(default_factory=dict)  # method -> obsm key
    benchmark_scores: Dict[str, Dict[str, float]] = field(default_factory=dict)
    best_method: Optional[str] = None
    files_created: List[str] = field(default_factory=list)


def prepare_for_integration(
    adata: sc.AnnData,
    batch_key: str,
    n_top_genes: int = 3000,
    layer: str = "counts",
) -> sc.AnnData:
    """
    Prepare data for integration.
    
    Ensures raw counts are available and selects HVGs per batch.
    
    Args:
        adata: AnnData object
        batch_key: Column containing batch information
        n_top_genes: Number of HVGs to select
        layer: Layer containing raw counts
        
    Returns:
        Prepared AnnData
    """
    logger.info(f"Preparing data for integration (batch_key={batch_key})...")
    
    # Ensure we have raw counts
    if layer in adata.layers:
        adata_prep = adata.copy()
        adata_prep.X = adata_prep.layers[layer].copy()
    elif adata.raw is not None:
        adata_prep = adata.raw.to_adata()
        adata_prep.obs = adata.obs.copy()
    else:
        adata_prep = adata.copy()
        logger.warning("No raw counts found. Using current .X matrix.")
    
    # Select HVGs considering batch
    logger.info(f"Selecting {n_top_genes} HVGs for integration...")
    sc.pp.highly_variable_genes(
        adata_prep,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
        batch_key=batch_key,
        subset=True,
    )
    
    logger.info(f"Prepared {adata_prep.n_obs} cells x {adata_prep.n_vars} genes")
    return adata_prep


def run_scvi(
    adata: sc.AnnData,
    batch_key: str,
    config: IntegrationConfig,
    layer: str = "counts",
) -> str:
    """
    Run scVI integration.
    
    scVI uses a variational autoencoder to learn batch-corrected latent representations.
    Best for complex batch effects and large datasets.
    
    Args:
        adata: AnnData object (modified in place)
        batch_key: Column containing batch information
        config: Integration configuration
        layer: Layer containing raw counts
        
    Returns:
        Key in obsm where embedding is stored ("X_scvi")
    """
    try:
        import scvi
    except ImportError:
        raise ImportError(
            "scvi-tools not installed. Install with: pip install scrna-qc[scvi]"
        )
    
    logger.info("Running scVI integration...")
    
    # Prepare data
    if layer in adata.layers:
        adata_scvi = adata.copy()
        adata_scvi.X = adata_scvi.layers[layer].copy()
    else:
        adata_scvi = adata.copy()
    
    # Setup scVI
    scvi.model.SCVI.setup_anndata(
        adata_scvi,
        layer=None,  # Use X
        batch_key=batch_key,
    )
    
    # Create and train model
    model = scvi.model.SCVI(
        adata_scvi,
        n_latent=config.scvi_n_latent,
        n_layers=config.scvi_n_layers,
        gene_likelihood=config.scvi_gene_likelihood,
    )
    
    model.train(
        max_epochs=config.scvi_max_epochs,
        early_stopping=config.scvi_early_stopping,
        early_stopping_patience=20,
    )
    
    # Get latent representation
    latent = model.get_latent_representation()
    adata.obsm["X_scvi"] = latent
    
    # Compute neighbors and UMAP from scVI embedding
    sc.pp.neighbors(adata, use_rep="X_scvi")
    sc.tl.umap(adata)
    adata.obsm["X_umap_scvi"] = adata.obsm["X_umap"].copy()
    
    logger.info(f"scVI complete: {latent.shape[1]} latent dimensions")
    
    return "X_scvi"


def run_harmony(
    adata: sc.AnnData,
    batch_key: str,
    config: IntegrationConfig,
) -> str:
    """
    Run Harmony integration.
    
    Harmony performs iterative batch correction in PCA space.
    Fast and effective for simpler batch effects.
    
    Args:
        adata: AnnData object (modified in place)
        batch_key: Column containing batch information
        config: Integration configuration
        
    Returns:
        Key in obsm where embedding is stored ("X_harmony")
    """
    try:
        import harmonypy as hm
    except ImportError:
        raise ImportError(
            "harmonypy not installed. Install with: pip install scrna-qc[harmony]"
        )
    
    logger.info("Running Harmony integration...")
    
    # Ensure PCA is computed
    if "X_pca" not in adata.obsm:
        logger.info("PCA not found. Computing PCA...")
        sc.tl.pca(adata, n_comps=config.harmony_n_pcs, use_highly_variable=True)
    
    # Run Harmony
    pca_key = "X_pca"
    n_pcs = min(config.harmony_n_pcs, adata.obsm[pca_key].shape[1])
    
    ho = hm.run_harmony(
        adata.obsm[pca_key][:, :n_pcs],
        adata.obs,
        batch_key,
        max_iter_harmony=config.harmony_max_iter,
    )
    
    # Store corrected embedding
    adata.obsm["X_harmony"] = ho.Z_corr.T
    
    # Compute neighbors and UMAP from Harmony embedding
    sc.pp.neighbors(adata, use_rep="X_harmony")
    sc.tl.umap(adata)
    adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
    
    logger.info(f"Harmony complete: {adata.obsm['X_harmony'].shape[1]} dimensions")
    
    return "X_harmony"


def run_scanorama(
    adata: sc.AnnData,
    batch_key: str,
    config: IntegrationConfig,
) -> str:
    """
    Run Scanorama integration.
    
    Scanorama uses mutual nearest neighbors for batch correction.
    Good for panoramic stitching of batches.
    
    Args:
        adata: AnnData object (modified in place)
        batch_key: Column containing batch information
        config: Integration configuration
        
    Returns:
        Key in obsm where embedding is stored ("X_scanorama")
    """
    try:
        import scanorama
    except ImportError:
        raise ImportError(
            "scanorama not installed. Install with: pip install scrna-qc[scanorama]"
        )
    
    logger.info("Running Scanorama integration...")
    
    # Split by batch
    batches = adata.obs[batch_key].unique().tolist()
    adatas = [adata[adata.obs[batch_key] == b].copy() for b in batches]
    
    # Run Scanorama
    scanorama.integrate_scanpy(adatas, dimred=50, knn=config.scanorama_knn)
    
    # Reconstruct merged embedding
    embeddings = [ad.obsm["X_scanorama"] for ad in adatas]
    
    # Create combined embedding with correct cell order
    combined = np.zeros((adata.n_obs, embeddings[0].shape[1]))
    for batch, ad, emb in zip(batches, adatas, embeddings):
        mask = adata.obs[batch_key] == batch
        combined[mask] = emb
    
    adata.obsm["X_scanorama"] = combined
    
    # Compute neighbors and UMAP
    sc.pp.neighbors(adata, use_rep="X_scanorama")
    sc.tl.umap(adata)
    adata.obsm["X_umap_scanorama"] = adata.obsm["X_umap"].copy()
    
    logger.info(f"Scanorama complete: {combined.shape[1]} dimensions")
    
    return "X_scanorama"


def benchmark_integration(
    adata: sc.AnnData,
    batch_key: str,
    embeddings: List[str],
    label_key: Optional[str] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Benchmark integration methods using scib-metrics.
    
    Computes:
    - Batch correction: silhouette_batch, iLISI
    - Bio conservation: silhouette_label, cLISI (if label_key provided)
    
    Args:
        adata: AnnData with integration embeddings
        batch_key: Column containing batch information
        embeddings: List of obsm keys to benchmark
        label_key: Optional column for biological conservation metrics
        
    Returns:
        Dictionary of {method: {metric: score}}
    """
    try:
        from scib_metrics.benchmark import Benchmarker
    except ImportError:
        logger.warning(
            "scib-metrics not installed. Skipping benchmarking. "
            "Install with: pip install scrna-qc[benchmark]"
        )
        # Return basic silhouette scores instead
        return _compute_basic_metrics(adata, batch_key, embeddings, label_key)
    
    logger.info(f"Benchmarking {len(embeddings)} integration methods...")
    
    # Run benchmarking
    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=embeddings,
        n_jobs=-1,
    )
    
    bm.benchmark()
    
    # Extract results
    results = bm.get_results(min_max_scale=False)
    scores = {}
    
    for embed_key in embeddings:
        method = embed_key.replace("X_", "")
        if embed_key in results.index:
            row = results.loc[embed_key]
            scores[method] = {col: float(row[col]) for col in results.columns}
    
    return scores


def _compute_basic_metrics(
    adata: sc.AnnData,
    batch_key: str,
    embeddings: List[str],
    label_key: Optional[str] = None,
) -> Dict[str, Dict[str, float]]:
    """Compute basic silhouette metrics without scib-metrics."""
    from sklearn.metrics import silhouette_score
    
    scores = {}
    batch_labels = adata.obs[batch_key].values
    
    for embed_key in embeddings:
        method = embed_key.replace("X_", "")
        embedding = adata.obsm[embed_key]
        
        # Batch silhouette (lower is better for batch mixing)
        try:
            batch_sil = silhouette_score(embedding, batch_labels, sample_size=5000)
        except:
            batch_sil = 0.0
        
        # Label silhouette (higher is better) if labels provided
        label_sil = 0.0
        if label_key and label_key in adata.obs.columns:
            try:
                label_sil = silhouette_score(
                    embedding, 
                    adata.obs[label_key].values, 
                    sample_size=5000
                )
            except:
                pass
        
        scores[method] = {
            "batch_silhouette": batch_sil,
            "label_silhouette": label_sil,
            # Negative batch silhouette is better (more mixing)
            "overall_score": label_sil - batch_sil,
        }
    
    return scores


def select_best_integration(
    scores: Dict[str, Dict[str, float]],
    priority: str = "balanced",
) -> str:
    """
    Select best integration method based on benchmark scores.
    
    Args:
        scores: Dictionary of {method: {metric: score}}
        priority: "batch" (batch correction), "bio" (bio conservation), 
                  or "balanced" (weighted average)
                  
    Returns:
        Name of best method
    """
    if not scores:
        return None
    
    method_scores = {}
    
    for method, metrics in scores.items():
        if priority == "batch":
            # Lower batch silhouette is better
            score = -metrics.get("batch_silhouette", 0)
        elif priority == "bio":
            score = metrics.get("label_silhouette", 0)
        else:  # balanced
            score = metrics.get("overall_score", 0)
        
        method_scores[method] = score
    
    best = max(method_scores, key=method_scores.get)
    logger.info(f"Best integration method: {best}")
    
    return best


def run_integration(
    adata: sc.AnnData,
    batch_key: str,
    methods: List[str] = None,
    config: IntegrationConfig = None,
    output_dir: str = "./",
) -> IntegrationResult:
    """
    Run one or multiple integration methods.
    
    When multiple methods are specified:
    1. Run each method and store embeddings (X_scvi, X_harmony, etc.)
    2. Benchmark all methods if more than one
    3. If auto_select_best=True, set "best" embedding as default
    4. Generate comparison plots
    
    Args:
        adata: AnnData object (modified in place)
        batch_key: Column containing batch information
        methods: List of methods to run (["scvi"], ["harmony"], ["scvi", "harmony"], ["all"])
        config: Integration configuration
        output_dir: Directory for output files
        
    Returns:
        IntegrationResult with embeddings, benchmark scores, and best method
    """
    if config is None:
        config = IntegrationConfig()
    
    if methods is None:
        methods = config.methods
    
    # Expand "all" to full method list
    if "all" in methods:
        methods = ["scvi", "harmony", "scanorama"]
    
    result = IntegrationResult()
    embeddings_run = []
    
    # Run each integration method
    for method in methods:
        try:
            if method == "scvi":
                embed_key = run_scvi(adata, batch_key, config)
            elif method == "harmony":
                embed_key = run_harmony(adata, batch_key, config)
            elif method == "scanorama":
                embed_key = run_scanorama(adata, batch_key, config)
            else:
                logger.warning(f"Unknown integration method: {method}")
                continue
            
            result.methods_run.append(method)
            result.embeddings[method] = embed_key
            embeddings_run.append(embed_key)
            
        except ImportError as e:
            logger.warning(f"Could not run {method}: {e}")
        except Exception as e:
            logger.error(f"Error running {method}: {e}")
    
    # Benchmark if multiple methods were run
    if len(embeddings_run) > 1:
        logger.info("Benchmarking integration methods...")
        result.benchmark_scores = benchmark_integration(
            adata,
            batch_key,
            embeddings_run,
            label_key=config.benchmark_label_key,
        )
        
        # Save benchmark results
        if result.benchmark_scores:
            benchmark_df = pd.DataFrame(result.benchmark_scores).T
            benchmark_path = os.path.join(output_dir, "integration_benchmark.csv")
            benchmark_df.to_csv(benchmark_path)
            result.files_created.append(benchmark_path)
            
            logger.info("\nBenchmark Results:")
            logger.info(benchmark_df.to_string())
    
    # Select best method
    if result.benchmark_scores and config.auto_select_best:
        result.best_method = select_best_integration(result.benchmark_scores)
        
        # Set best embedding as default
        if result.best_method in result.embeddings:
            best_key = result.embeddings[result.best_method]
            adata.obsm["X_integrated"] = adata.obsm[best_key].copy()
            
            # Update neighbors and UMAP to use best
            sc.pp.neighbors(adata, use_rep="X_integrated")
            sc.tl.umap(adata)
    
    elif len(result.methods_run) == 1:
        result.best_method = result.methods_run[0]
        best_key = result.embeddings[result.best_method]
        adata.obsm["X_integrated"] = adata.obsm[best_key].copy()
    
    logger.info(f"Integration complete. Methods run: {result.methods_run}")
    if result.best_method:
        logger.info(f"Best method: {result.best_method}")
    
    return result
