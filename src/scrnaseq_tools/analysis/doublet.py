"""
Doublet detection module for AI Agent-friendly scRNA-seq analysis.

This module provides automatic doublet detection:
- Scrublet-based doublet detection
- Batch-aware processing
- Automatic parameter selection
- Doublet filtering

Usage:
    from scrnaseq_tools.analysis.doublet import run_doublet_detection

    result = run_doublet_detection(adata, batch_key="sample")
    adata_filtered = result.adata_filtered
"""

import logging
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


# Empirical doublet rates based on 10x Genomics data
# Source: 10x Genomics technical documentation
DOUBLET_RATE_TABLE = {
    500: 0.004,
    1000: 0.008,
    2000: 0.016,
    3000: 0.023,
    4000: 0.031,
    5000: 0.039,
    6000: 0.046,
    7000: 0.054,
    8000: 0.061,
    9000: 0.069,
    10000: 0.076,
    15000: 0.115,
    20000: 0.153,
}


@dataclass
class DoubletConfig:
    """Configuration for doublet detection."""

    method: str = "scrublet"  # Currently only "scrublet" supported
    expected_doublet_rate: Optional[float] = None  # Auto-estimated if None
    n_prin_comps: int = 30
    min_counts: int = 2
    min_cells: int = 3
    min_gene_variability_pctl: float = 85.0
    n_neighbors: Optional[int] = None  # Auto-selected if None
    threshold: Optional[float] = None  # Auto-selected if None
    random_state: int = 42

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DoubletResult:
    """Result of doublet detection."""

    method: str = "scrublet"
    n_cells_total: int = 0
    n_doublets: int = 0
    n_singlets: int = 0
    doublet_rate: float = 0.0

    # Per-batch results (if batch_key provided)
    batch_results: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # Thresholds used
    threshold_used: float = 0.0
    expected_doublet_rate: float = 0.0

    # Scores (added to adata.obs)
    score_column: str = "doublet_score"
    prediction_column: str = "predicted_doublet"

    # Filtered data
    adata_filtered: Any = None  # AnnData with doublets removed

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "method": self.method,
            "n_cells_total": self.n_cells_total,
            "n_doublets": self.n_doublets,
            "n_singlets": self.n_singlets,
            "doublet_rate": self.doublet_rate,
            "threshold_used": self.threshold_used,
            "expected_doublet_rate": self.expected_doublet_rate,
            "score_column": self.score_column,
            "prediction_column": self.prediction_column,
        }
        if self.batch_results:
            result["batch_results"] = self.batch_results
        return result

    def summary(self) -> str:
        lines = [
            f"Doublet Detection Results ({self.method})",
            f"  Total cells: {self.n_cells_total:,}",
            f"  Doublets detected: {self.n_doublets:,} ({self.doublet_rate:.1%})",
            f"  Singlets: {self.n_singlets:,}",
            f"  Threshold: {self.threshold_used:.3f}",
        ]
        if self.batch_results:
            lines.append("  Per-batch results:")
            for batch, stats in self.batch_results.items():
                lines.append(f"    {batch}: {stats['n_doublets']}/{stats['n_cells']} "
                           f"({stats['doublet_rate']:.1%})")
        return "\n".join(lines)


def estimate_doublet_rate(n_cells: int) -> float:
    """
    Estimate expected doublet rate based on cell count.

    Uses empirical data from 10x Genomics. Linear interpolation
    between known points.

    Args:
        n_cells: Number of cells loaded

    Returns:
        Estimated doublet rate (0-1)
    """
    if n_cells <= 0:
        return 0.0

    # Sort the table by cell count
    sorted_counts = sorted(DOUBLET_RATE_TABLE.keys())

    # Handle edge cases
    if n_cells <= sorted_counts[0]:
        return DOUBLET_RATE_TABLE[sorted_counts[0]]
    if n_cells >= sorted_counts[-1]:
        # Extrapolate linearly for very large datasets
        # Approximately 0.008 per 1000 cells
        return min(0.3, n_cells * 0.008 / 1000)

    # Linear interpolation
    for i in range(len(sorted_counts) - 1):
        if sorted_counts[i] <= n_cells < sorted_counts[i + 1]:
            x1, x2 = sorted_counts[i], sorted_counts[i + 1]
            y1, y2 = DOUBLET_RATE_TABLE[x1], DOUBLET_RATE_TABLE[x2]
            rate = y1 + (y2 - y1) * (n_cells - x1) / (x2 - x1)
            return rate

    return 0.05  # Default fallback


def auto_select_doublet_params(
    n_cells: int,
    n_genes: int,
) -> Dict[str, Any]:
    """
    Automatically select doublet detection parameters.

    Args:
        n_cells: Number of cells
        n_genes: Number of genes

    Returns:
        Dictionary of recommended parameters
    """
    params = {
        "expected_doublet_rate": estimate_doublet_rate(n_cells),
    }

    # Adjust n_neighbors based on cell count
    if n_cells < 1000:
        params["n_neighbors"] = max(10, n_cells // 50)
    elif n_cells < 10000:
        params["n_neighbors"] = 30
    else:
        params["n_neighbors"] = min(50, n_cells // 200)

    # Adjust PCs based on data size
    if n_genes < 2000:
        params["n_prin_comps"] = min(20, n_genes // 100)
    else:
        params["n_prin_comps"] = 30

    logger.info(f"Auto-selected doublet params: {params}")
    return params


def run_scrublet(
    adata: sc.AnnData,
    config: Optional[DoubletConfig] = None,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Run Scrublet doublet detection.

    Args:
        adata: AnnData object (raw counts preferred)
        config: DoubletConfig with parameters

    Returns:
        Tuple of (doublet_scores, predicted_doublets, threshold)
    """
    try:
        import scrublet as scr
    except ImportError:
        raise ImportError(
            "scrublet is required for doublet detection. "
            "Install with: pip install scrublet"
        )

    if config is None:
        config = DoubletConfig()

    # Auto-select parameters if needed
    auto_params = auto_select_doublet_params(adata.n_obs, adata.n_vars)

    expected_rate = config.expected_doublet_rate or auto_params["expected_doublet_rate"]
    n_neighbors = config.n_neighbors or auto_params.get("n_neighbors", 30)
    n_prin_comps = config.n_prin_comps or auto_params.get("n_prin_comps", 30)

    logger.info(f"Running Scrublet (expected_rate={expected_rate:.3f}, "
                f"n_neighbors={n_neighbors}, n_pcs={n_prin_comps})")

    # Get count matrix
    if 'counts' in adata.layers:
        counts = adata.layers['counts']
    else:
        counts = adata.X

    # Ensure we have a matrix (not sparse for some scrublet versions)
    from scipy import sparse
    if sparse.issparse(counts):
        counts_dense = counts.toarray() if counts.shape[0] * counts.shape[1] < 1e8 else counts
    else:
        counts_dense = counts

    # Initialize Scrublet
    scrub = scr.Scrublet(
        counts_dense,
        expected_doublet_rate=expected_rate,
        random_state=config.random_state,
    )

    # Run doublet detection
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=config.min_counts,
        min_cells=config.min_cells,
        min_gene_variability_pctl=config.min_gene_variability_pctl,
        n_prin_comps=min(n_prin_comps, adata.n_vars - 1, adata.n_obs - 1),
    )

    threshold = scrub.threshold_

    return doublet_scores, predicted_doublets, threshold


def run_doublet_detection(
    adata: sc.AnnData,
    batch_key: Optional[str] = None,
    config: Optional[DoubletConfig] = None,
    inplace: bool = True,
) -> DoubletResult:
    """
    Run doublet detection, optionally per batch.

    This is the main entry point for doublet detection. If batch_key is
    provided, runs detection separately for each batch to avoid
    cross-batch artifacts.

    Args:
        adata: AnnData object
        batch_key: Column in obs containing batch information
        config: DoubletConfig with parameters
        inplace: Whether to modify adata in place

    Returns:
        DoubletResult with detection results
    """
    if config is None:
        config = DoubletConfig()

    result = DoubletResult(
        method=config.method,
        n_cells_total=adata.n_obs,
    )

    if not inplace:
        adata = adata.copy()

    # Initialize columns
    adata.obs[result.score_column] = 0.0
    adata.obs[result.prediction_column] = False

    if batch_key is not None and batch_key in adata.obs.columns:
        # Per-batch processing
        batches = adata.obs[batch_key].unique()
        logger.info(f"Running doublet detection per batch ({len(batches)} batches)")

        all_scores = np.zeros(adata.n_obs)
        all_predictions = np.zeros(adata.n_obs, dtype=bool)
        total_doublets = 0

        for batch in batches:
            batch_mask = adata.obs[batch_key] == batch
            batch_adata = adata[batch_mask].copy()

            if batch_adata.n_obs < 50:
                logger.warning(f"Batch '{batch}' has only {batch_adata.n_obs} cells. "
                             "Skipping doublet detection for this batch.")
                result.batch_results[str(batch)] = {
                    "n_cells": batch_adata.n_obs,
                    "n_doublets": 0,
                    "doublet_rate": 0.0,
                    "skipped": True,
                }
                continue

            try:
                scores, predictions, threshold = run_scrublet(batch_adata, config)

                # Store results
                batch_indices = np.where(batch_mask)[0]
                all_scores[batch_indices] = scores
                all_predictions[batch_indices] = predictions

                n_doublets = predictions.sum()
                total_doublets += n_doublets

                result.batch_results[str(batch)] = {
                    "n_cells": batch_adata.n_obs,
                    "n_doublets": int(n_doublets),
                    "doublet_rate": float(n_doublets / batch_adata.n_obs),
                    "threshold": float(threshold),
                }

                logger.info(f"  Batch '{batch}': {n_doublets}/{batch_adata.n_obs} "
                          f"doublets ({n_doublets/batch_adata.n_obs:.1%})")

            except Exception as e:
                logger.error(f"Doublet detection failed for batch '{batch}': {e}")
                result.batch_results[str(batch)] = {
                    "n_cells": batch_adata.n_obs,
                    "n_doublets": 0,
                    "doublet_rate": 0.0,
                    "error": str(e),
                }

        adata.obs[result.score_column] = all_scores
        adata.obs[result.prediction_column] = all_predictions
        result.n_doublets = total_doublets
        result.threshold_used = np.mean([
            b["threshold"] for b in result.batch_results.values()
            if "threshold" in b
        ]) if result.batch_results else 0.0

    else:
        # Single run for all data
        logger.info("Running doublet detection on all cells")

        try:
            scores, predictions, threshold = run_scrublet(adata, config)

            adata.obs[result.score_column] = scores
            adata.obs[result.prediction_column] = predictions
            result.n_doublets = int(predictions.sum())
            result.threshold_used = float(threshold)

        except Exception as e:
            logger.error(f"Doublet detection failed: {e}")
            raise

    # Calculate final stats
    result.n_singlets = result.n_cells_total - result.n_doublets
    result.doublet_rate = result.n_doublets / result.n_cells_total if result.n_cells_total > 0 else 0.0
    result.expected_doublet_rate = config.expected_doublet_rate or estimate_doublet_rate(result.n_cells_total)

    logger.info(f"Doublet detection complete: {result.n_doublets:,} doublets "
                f"({result.doublet_rate:.1%})")

    return result


def filter_doublets(
    adata: sc.AnnData,
    prediction_column: str = "predicted_doublet",
    inplace: bool = False,
) -> sc.AnnData:
    """
    Remove predicted doublets from AnnData.

    Args:
        adata: AnnData object with doublet predictions
        prediction_column: Column containing boolean doublet predictions
        inplace: Whether to filter in place

    Returns:
        Filtered AnnData
    """
    if prediction_column not in adata.obs.columns:
        raise ValueError(f"Column '{prediction_column}' not found in adata.obs. "
                        "Run doublet detection first.")

    singlet_mask = ~adata.obs[prediction_column].astype(bool)
    n_before = adata.n_obs
    n_doublets = (~singlet_mask).sum()

    if inplace:
        adata._inplace_subset_obs(singlet_mask)
        filtered_adata = adata
    else:
        filtered_adata = adata[singlet_mask].copy()

    logger.info(f"Filtered doublets: {n_before:,} -> {filtered_adata.n_obs:,} cells "
                f"(removed {n_doublets:,})")

    return filtered_adata
