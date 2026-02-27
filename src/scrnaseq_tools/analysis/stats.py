"""
Statistical methods for scRNA-seq QC and dimensionality reduction.

Provides:
- Adaptive QC threshold selection using robust statistics (MAD)
- Explained variance-based PC selection
"""

import logging
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

logger = logging.getLogger(__name__)


@dataclass
class QCThresholds:
    """Data class for QC thresholds with statistical rationale."""
    min_genes: int
    max_genes: int
    min_counts: int
    max_counts: int
    max_mt_pct: float

    # Statistical basis for thresholds
    n_genes_median: float = 0.0
    n_genes_mad: float = 0.0
    n_counts_median: float = 0.0
    n_counts_mad: float = 0.0
    mt_pct_median: float = 0.0
    mt_pct_mad: float = 0.0

    # Method used
    method: str = "mad"
    nmads: float = 3.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PCSelectionResult:
    """Results from PC selection analysis."""
    n_pcs_selected: int
    method: str  # "elbow" or "variance_ratio"
    variance_explained: List[float] = field(default_factory=list)
    cumulative_variance: List[float] = field(default_factory=list)
    elbow_point: Optional[int] = None
    variance_threshold: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def compute_mad(x: np.ndarray) -> float:
    """
    Compute Median Absolute Deviation (MAD).

    MAD = median(|X - median(X)|)

    This is a robust measure of spread that is resistant to outliers.
    For normally distributed data, MAD ≈ 0.6745 * σ
    """
    median = np.median(x)
    return np.median(np.abs(x - median))


def mad_to_std(mad: float) -> float:
    """
    Convert MAD to standard deviation equivalent.

    For normal distribution: σ ≈ 1.4826 * MAD
    """
    return 1.4826 * mad


def compute_adaptive_qc_thresholds(
    adata: sc.AnnData,
    nmads: float = 3.0,
    min_genes_floor: int = 200,
    max_mt_pct_ceiling: float = 20.0,
    method: str = "mad"
) -> QCThresholds:
    """
    Compute QC thresholds using robust statistics.

    Uses Median Absolute Deviation (MAD) for outlier detection:
    - Lower threshold: median - nmads * MAD
    - Upper threshold: median + nmads * MAD

    This is more robust than using mean ± SD because:
    1. Median is resistant to outliers
    2. MAD provides robust spread estimation
    3. Works well for skewed distributions common in scRNA-seq

    Args:
        adata: AnnData with QC metrics computed
        nmads: Number of MADs from median for threshold (default: 3)
        min_genes_floor: Minimum allowed value for min_genes
        max_mt_pct_ceiling: Maximum allowed value for max_mt_pct
        method: "mad" for MAD-based, "quantile" for quantile-based

    Returns:
        QCThresholds with computed values and rationale
    """
    # Ensure QC metrics are computed
    if 'n_genes_by_counts' not in adata.obs.columns:
        raise ValueError("QC metrics not computed. Run sc.pp.calculate_qc_metrics first.")

    n_genes = adata.obs['n_genes_by_counts'].values
    total_counts = adata.obs['total_counts'].values
    mt_pct = adata.obs['pct_counts_mt'].values

    if method == "mad":
        # Compute MAD-based thresholds
        n_genes_median = np.median(n_genes)
        n_genes_mad = compute_mad(n_genes)

        n_counts_median = np.median(total_counts)
        n_counts_mad = compute_mad(total_counts)

        mt_pct_median = np.median(mt_pct)
        mt_pct_mad = compute_mad(mt_pct)

        # Compute thresholds
        min_genes = max(
            int(n_genes_median - nmads * mad_to_std(n_genes_mad)),
            min_genes_floor
        )
        max_genes = int(n_genes_median + nmads * mad_to_std(n_genes_mad))

        min_counts = max(
            int(n_counts_median - nmads * mad_to_std(n_counts_mad)),
            500  # Reasonable floor
        )
        max_counts = int(n_counts_median + nmads * mad_to_std(n_counts_mad))

        # MT percentage: only upper threshold matters
        max_mt_pct = min(
            mt_pct_median + nmads * mad_to_std(mt_pct_mad),
            max_mt_pct_ceiling
        )

        return QCThresholds(
            min_genes=min_genes,
            max_genes=max_genes,
            min_counts=min_counts,
            max_counts=max_counts,
            max_mt_pct=max_mt_pct,
            n_genes_median=n_genes_median,
            n_genes_mad=n_genes_mad,
            n_counts_median=n_counts_median,
            n_counts_mad=n_counts_mad,
            mt_pct_median=mt_pct_median,
            mt_pct_mad=mt_pct_mad,
            method=method,
            nmads=nmads
        )

    elif method == "quantile":
        # Quantile-based thresholds (more conservative)
        min_genes = max(int(np.percentile(n_genes, 1)), min_genes_floor)
        max_genes = int(np.percentile(n_genes, 99))

        min_counts = max(int(np.percentile(total_counts, 1)), 500)
        max_counts = int(np.percentile(total_counts, 99))

        max_mt_pct = min(np.percentile(mt_pct, 95), max_mt_pct_ceiling)

        return QCThresholds(
            min_genes=min_genes,
            max_genes=max_genes,
            min_counts=min_counts,
            max_counts=max_counts,
            max_mt_pct=max_mt_pct,
            n_genes_median=np.median(n_genes),
            n_genes_mad=compute_mad(n_genes),
            n_counts_median=np.median(total_counts),
            n_counts_mad=compute_mad(total_counts),
            mt_pct_median=np.median(mt_pct),
            mt_pct_mad=compute_mad(mt_pct),
            method=method,
            nmads=nmads
        )

    else:
        raise ValueError(f"Unknown method: {method}")


def apply_adaptive_qc(
    adata: sc.AnnData,
    thresholds: QCThresholds,
    min_cells: int = 3
) -> Dict[str, int]:
    """
    Apply adaptive QC thresholds to filter cells and genes.

    Args:
        adata: AnnData object (modified in place)
        thresholds: QCThresholds computed by compute_adaptive_qc_thresholds
        min_cells: Minimum cells per gene

    Returns:
        Dictionary with filtering statistics
    """
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filter by genes per cell (lower bound)
    mask_min_genes = adata.obs['n_genes_by_counts'] >= thresholds.min_genes

    # Filter by genes per cell (upper bound - doublet detection)
    mask_max_genes = adata.obs['n_genes_by_counts'] <= thresholds.max_genes

    # Filter by total counts (lower bound)
    mask_min_counts = adata.obs['total_counts'] >= thresholds.min_counts

    # Filter by total counts (upper bound)
    mask_max_counts = adata.obs['total_counts'] <= thresholds.max_counts

    # Filter by MT percentage
    mask_mt = adata.obs['pct_counts_mt'] <= thresholds.max_mt_pct

    # Combined mask
    cell_mask = mask_min_genes & mask_max_genes & mask_min_counts & mask_max_counts & mask_mt

    # Apply cell filter
    adata._inplace_subset_obs(cell_mask)

    # Filter genes by minimum cells
    sc.pp.filter_genes(adata, min_cells=min_cells)

    stats = {
        'cells_before': n_cells_before,
        'cells_after': adata.n_obs,
        'genes_before': n_genes_before,
        'genes_after': adata.n_vars,
        'cells_removed': n_cells_before - adata.n_obs,
        'genes_removed': n_genes_before - adata.n_vars,
        'retention_rate': adata.n_obs / n_cells_before * 100,
    }

    logger.info(
        f"QC filtering: {n_cells_before} -> {adata.n_obs} cells "
        f"({stats['retention_rate']:.1f}% retained), "
        f"{n_genes_before} -> {adata.n_vars} genes"
    )

    return stats


def select_pcs_elbow(
    variance_ratio: np.ndarray,
    min_pcs: int = 10,
    max_pcs: int = 50
) -> int:
    """
    Select number of PCs using the elbow method.

    The elbow point is where the rate of decrease in variance explained
    slows down significantly. We use the point of maximum curvature.

    Args:
        variance_ratio: Explained variance ratio per PC
        min_pcs: Minimum PCs to consider
        max_pcs: Maximum PCs to consider

    Returns:
        Optimal number of PCs
    """
    n_pcs = min(len(variance_ratio), max_pcs)

    # Cumulative variance
    cum_var = np.cumsum(variance_ratio[:n_pcs])

    # Normalized coordinates for elbow detection
    x = np.arange(n_pcs)
    y = cum_var

    # Line from first to last point
    line_vec = np.array([x[-1] - x[0], y[-1] - y[0]])
    line_vec = line_vec / np.linalg.norm(line_vec)

    # Distance from each point to the line
    distances = []
    for i in range(n_pcs):
        point_vec = np.array([x[i] - x[0], y[i] - y[0]])
        # Project onto perpendicular
        proj = point_vec - np.dot(point_vec, line_vec) * line_vec
        distances.append(np.linalg.norm(proj))

    # Elbow is where distance is maximum
    elbow_idx = np.argmax(distances)

    # Ensure within bounds
    n_pcs_selected = max(min_pcs, min(elbow_idx + 1, max_pcs))

    return n_pcs_selected


def select_pcs_variance(
    variance_ratio: np.ndarray,
    variance_threshold: float = 0.90,
    min_pcs: int = 10,
    max_pcs: int = 50
) -> int:
    """
    Select number of PCs based on cumulative variance threshold.

    Args:
        variance_ratio: Explained variance ratio per PC
        variance_threshold: Target cumulative variance (default: 90%)
        min_pcs: Minimum PCs to select
        max_pcs: Maximum PCs to select

    Returns:
        Number of PCs needed to reach threshold
    """
    cum_var = np.cumsum(variance_ratio)

    # Find first PC where cumulative variance exceeds threshold
    n_pcs = np.searchsorted(cum_var, variance_threshold) + 1

    # Ensure within bounds
    n_pcs = max(min_pcs, min(n_pcs, max_pcs, len(variance_ratio)))

    return n_pcs


def compute_pc_selection(
    adata: sc.AnnData,
    method: str = "elbow",
    variance_threshold: float = 0.90,
    min_pcs: int = 10,
    max_pcs: int = 50
) -> PCSelectionResult:
    """
    Select optimal number of PCs using specified method.

    Args:
        adata: AnnData with PCA computed
        method: "elbow" or "variance_ratio"
        variance_threshold: For variance_ratio method
        min_pcs: Minimum PCs
        max_pcs: Maximum PCs

    Returns:
        PCSelectionResult with analysis
    """
    if 'X_pca' not in adata.obsm:
        raise ValueError("PCA not computed. Run sc.tl.pca first.")

    variance_ratio = adata.uns['pca']['variance_ratio']
    cum_var = np.cumsum(variance_ratio)

    if method == "elbow":
        n_pcs = select_pcs_elbow(variance_ratio, min_pcs, max_pcs)
        elbow_point = n_pcs
        var_thresh = None
    elif method == "variance_ratio":
        n_pcs = select_pcs_variance(variance_ratio, variance_threshold, min_pcs, max_pcs)
        elbow_point = select_pcs_elbow(variance_ratio, min_pcs, max_pcs)
        var_thresh = variance_threshold
    else:
        raise ValueError(f"Unknown method: {method}")

    result = PCSelectionResult(
        n_pcs_selected=n_pcs,
        method=method,
        variance_explained=variance_ratio.tolist(),
        cumulative_variance=cum_var.tolist(),
        elbow_point=elbow_point,
        variance_threshold=var_thresh
    )

    logger.info(
        f"PC selection ({method}): {n_pcs} PCs, "
        f"explaining {cum_var[n_pcs-1]*100:.1f}% variance"
    )

    return result


def compute_qc_distributions(adata: sc.AnnData) -> Dict[str, Dict[str, float]]:
    """
    Compute distribution statistics for QC metrics.

    Returns comprehensive statistics for reporting:
    - Mean, median, SD, MAD
    - Quantiles (5th, 25th, 50th, 75th, 95th)
    """
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    distributions = {}

    for metric in metrics:
        if metric not in adata.obs.columns:
            continue

        values = adata.obs[metric].values

        distributions[metric] = {
            'mean': float(np.mean(values)),
            'median': float(np.median(values)),
            'std': float(np.std(values)),
            'mad': float(compute_mad(values)),
            'min': float(np.min(values)),
            'max': float(np.max(values)),
            'q05': float(np.percentile(values, 5)),
            'q25': float(np.percentile(values, 25)),
            'q50': float(np.percentile(values, 50)),
            'q75': float(np.percentile(values, 75)),
            'q95': float(np.percentile(values, 95)),
        }

    return distributions


def generate_qc_rationale(
    thresholds: QCThresholds,
    distributions: Dict[str, Dict[str, float]]
) -> str:
    """
    Generate human-readable rationale for QC thresholds.
    """
    lines = [
        "### QC Threshold Rationale",
        "",
        f"**Method**: {thresholds.method.upper()}-based outlier detection "
        f"(±{thresholds.nmads} MADs from median)",
        "",
        "#### Genes per Cell",
        f"- Distribution: median={thresholds.n_genes_median:.0f}, MAD={thresholds.n_genes_mad:.0f}",
        f"- Lower threshold: {thresholds.min_genes} (removes low-quality/empty droplets)",
        f"- Upper threshold: {thresholds.max_genes} (removes potential doublets)",
        "",
        "#### UMI Counts per Cell",
        f"- Distribution: median={thresholds.n_counts_median:.0f}, MAD={thresholds.n_counts_mad:.0f}",
        f"- Lower threshold: {thresholds.min_counts} (removes low-quality cells)",
        f"- Upper threshold: {thresholds.max_counts} (removes potential doublets)",
        "",
        "#### Mitochondrial Percentage",
        f"- Distribution: median={thresholds.mt_pct_median:.1f}%, MAD={thresholds.mt_pct_mad:.1f}%",
        f"- Upper threshold: {thresholds.max_mt_pct:.1f}% (removes stressed/dying cells)",
    ]

    return "\n".join(lines)
