"""
Data profiling module for AI Agent-friendly scRNA-seq analysis.

This module provides automatic detection and profiling of scRNA-seq data:
- File format detection (h5ad, 10x_mtx, 10x_h5, multi-sample)
- Batch key auto-detection
- MT gene pattern detection
- Existing QC status detection
- Comprehensive data profiling for agent consumption

Usage:
    from scrnaseq_tools.analysis.profile import profile_data

    profile = profile_data("data.h5ad")
    print(profile.to_dict())  # Agent-friendly structured output
"""

import os
import re
import logging
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any, Tuple

import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)


@dataclass
class FileDetectionResult:
    """Result of file format detection."""

    format: str  # "h5ad", "10x_mtx", "10x_h5", "multi_sample"
    path: str
    is_directory: bool
    detected_files: List[str] = field(default_factory=list)
    n_samples: int = 1
    confidence: float = 1.0
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class BatchKeyDetection:
    """Result of batch key auto-detection."""

    detected_key: Optional[str] = None
    candidates: List[Dict[str, Any]] = field(default_factory=list)
    confidence: float = 0.0
    reasoning: str = ""
    n_batches: int = 0
    batch_sizes: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class MTGeneDetection:
    """Result of mitochondrial gene pattern detection."""

    pattern: str = ""  # e.g., "MT-", "mt-", "Mt-"
    n_mt_genes: int = 0
    mt_gene_examples: List[str] = field(default_factory=list)
    mt_pct_mean: float = 0.0
    mt_pct_median: float = 0.0
    detected: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ExistingQCStatus:
    """Detection of existing QC application status."""

    qc_metrics_present: bool = False
    appears_filtered: bool = False
    normalized: bool = False
    log_transformed: bool = False
    hvgs_selected: bool = False
    pca_computed: bool = False
    neighbors_computed: bool = False
    umap_computed: bool = False
    clustered: bool = False

    indicators: Dict[str, Any] = field(default_factory=dict)
    recommendation: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DataProfile:
    """Comprehensive data profile for AI Agent consumption."""

    # Basic info
    input_path: str = ""
    n_cells: int = 0
    n_genes: int = 0

    # Detection results
    file_detection: FileDetectionResult = None
    batch_key_detection: BatchKeyDetection = None
    mt_gene_detection: MTGeneDetection = None
    qc_status: ExistingQCStatus = None

    # Data characteristics
    sparsity: float = 0.0
    median_genes_per_cell: float = 0.0
    median_counts_per_cell: float = 0.0

    # Available annotations
    obs_columns: List[str] = field(default_factory=list)
    var_columns: List[str] = field(default_factory=list)
    layers: List[str] = field(default_factory=list)
    obsm_keys: List[str] = field(default_factory=list)
    uns_keys: List[str] = field(default_factory=list)

    # Recommendations
    recommended_batch_key: Optional[str] = None
    recommended_mt_pattern: Optional[str] = None
    recommended_actions: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {
            "input_path": self.input_path,
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "sparsity": self.sparsity,
            "median_genes_per_cell": self.median_genes_per_cell,
            "median_counts_per_cell": self.median_counts_per_cell,
            "obs_columns": self.obs_columns,
            "var_columns": self.var_columns,
            "layers": self.layers,
            "obsm_keys": self.obsm_keys,
            "uns_keys": self.uns_keys,
            "recommended_batch_key": self.recommended_batch_key,
            "recommended_mt_pattern": self.recommended_mt_pattern,
            "recommended_actions": self.recommended_actions,
            "warnings": self.warnings,
        }

        if self.file_detection:
            result["file_detection"] = self.file_detection.to_dict()
        if self.batch_key_detection:
            result["batch_key_detection"] = self.batch_key_detection.to_dict()
        if self.mt_gene_detection:
            result["mt_gene_detection"] = self.mt_gene_detection.to_dict()
        if self.qc_status:
            result["qc_status"] = self.qc_status.to_dict()

        return result

    def to_json(self) -> str:
        """Convert to JSON string."""
        import json
        return json.dumps(self.to_dict(), indent=2)

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "=" * 60,
            "Data Profile Summary",
            "=" * 60,
            f"Input: {self.input_path}",
            f"Cells: {self.n_cells:,}",
            f"Genes: {self.n_genes:,}",
            f"Sparsity: {self.sparsity:.1%}",
            "",
        ]

        if self.file_detection:
            lines.append(f"Format: {self.file_detection.format}")
            if self.file_detection.n_samples > 1:
                lines.append(f"Samples: {self.file_detection.n_samples}")

        if self.batch_key_detection and self.batch_key_detection.detected_key:
            lines.append(f"Detected batch key: {self.batch_key_detection.detected_key}")
            lines.append(f"  - {self.batch_key_detection.n_batches} batches")

        if self.mt_gene_detection and self.mt_gene_detection.detected:
            lines.append(f"MT genes: {self.mt_gene_detection.n_mt_genes} (pattern: {self.mt_gene_detection.pattern})")
            lines.append(f"  - Mean MT%: {self.mt_gene_detection.mt_pct_mean:.1f}%")

        if self.qc_status:
            status_items = []
            if self.qc_status.normalized:
                status_items.append("normalized")
            if self.qc_status.pca_computed:
                status_items.append("PCA")
            if self.qc_status.umap_computed:
                status_items.append("UMAP")
            if self.qc_status.clustered:
                status_items.append("clustered")
            if status_items:
                lines.append(f"Existing processing: {', '.join(status_items)}")

        if self.warnings:
            lines.append("")
            lines.append("Warnings:")
            for w in self.warnings:
                lines.append(f"  - {w}")

        if self.recommended_actions:
            lines.append("")
            lines.append("Recommended actions:")
            for i, action in enumerate(self.recommended_actions, 1):
                lines.append(f"  {i}. {action}")

        lines.append("=" * 60)
        return "\n".join(lines)


def detect_file_format(path: str) -> FileDetectionResult:
    """
    Detect the format of input data.

    Supports:
    - h5ad: AnnData HDF5 format
    - 10x_mtx: 10x Genomics MTX directory
    - 10x_h5: 10x Genomics H5 format
    - multi_sample: Directory with multiple samples

    Args:
        path: Path to file or directory

    Returns:
        FileDetectionResult with detected format
    """
    path = Path(path)
    result = FileDetectionResult(
        path=str(path),
        format="unknown",
        is_directory=path.is_dir(),
    )

    if not path.exists():
        result.notes.append(f"Path does not exist: {path}")
        result.confidence = 0.0
        return result

    # Single file detection
    if path.is_file():
        if path.suffix == ".h5ad":
            result.format = "h5ad"
            result.detected_files = [str(path)]
            result.confidence = 1.0
        elif path.suffix == ".h5":
            result.format = "10x_h5"
            result.detected_files = [str(path)]
            result.confidence = 0.9
            result.notes.append("Assuming 10x Genomics H5 format")
        else:
            result.notes.append(f"Unknown file extension: {path.suffix}")
            result.confidence = 0.0
        return result

    # Directory detection
    if path.is_dir():
        # Check for 10x MTX files
        mtx_files = list(path.glob("*matrix*.mtx*"))
        if mtx_files:
            result.format = "10x_mtx"
            result.detected_files = [str(f) for f in mtx_files]
            result.confidence = 0.95
            return result

        # Check for filtered_feature_bc_matrix subdirectory
        filtered_dir = path / "filtered_feature_bc_matrix"
        if filtered_dir.exists():
            result.format = "10x_mtx"
            result.detected_files = [str(filtered_dir)]
            result.notes.append("Found filtered_feature_bc_matrix subdirectory")
            result.confidence = 1.0
            return result

        # Check for outs/filtered_feature_bc_matrix (CellRanger output)
        outs_filtered = path / "outs" / "filtered_feature_bc_matrix"
        if outs_filtered.exists():
            result.format = "10x_mtx"
            result.detected_files = [str(outs_filtered)]
            result.notes.append("Found CellRanger outs directory")
            result.confidence = 1.0
            return result

        # Check for multiple h5ad files
        h5ad_files = list(path.glob("*.h5ad"))
        if len(h5ad_files) > 1:
            result.format = "multi_sample"
            result.detected_files = [str(f) for f in h5ad_files]
            result.n_samples = len(h5ad_files)
            result.confidence = 0.9
            return result
        elif len(h5ad_files) == 1:
            result.format = "h5ad"
            result.detected_files = [str(h5ad_files[0])]
            result.confidence = 1.0
            return result

        # Check for multiple sample directories
        subdirs = [d for d in path.iterdir() if d.is_dir()]
        sample_dirs = []
        for subdir in subdirs:
            # Check if subdir contains 10x data
            if list(subdir.glob("*matrix*.mtx*")) or \
               (subdir / "filtered_feature_bc_matrix").exists() or \
               list(subdir.glob("*.h5")):
                sample_dirs.append(subdir)

        if len(sample_dirs) > 1:
            result.format = "multi_sample"
            result.detected_files = [str(d) for d in sample_dirs]
            result.n_samples = len(sample_dirs)
            result.confidence = 0.85
            result.notes.append(f"Found {len(sample_dirs)} potential sample directories")
            return result

        # Check for h5 files
        h5_files = list(path.glob("*.h5"))
        if h5_files:
            if len(h5_files) > 1:
                result.format = "multi_sample"
                result.n_samples = len(h5_files)
            else:
                result.format = "10x_h5"
            result.detected_files = [str(f) for f in h5_files]
            result.confidence = 0.8
            return result

        result.notes.append("No recognizable data files found in directory")
        result.confidence = 0.0

    return result


def detect_batch_key(
    adata: sc.AnnData,
    candidate_patterns: Optional[List[str]] = None,
) -> BatchKeyDetection:
    """
    Auto-detect batch key from AnnData obs columns.

    Heuristics:
    1. Common batch key names (batch, sample, donor, patient, etc.)
    2. Categorical columns with reasonable number of categories
    3. Columns with "batch" or "sample" in the name

    Args:
        adata: AnnData object
        candidate_patterns: Optional list of patterns to match

    Returns:
        BatchKeyDetection with detected batch key
    """
    if candidate_patterns is None:
        candidate_patterns = [
            "batch", "sample", "donor", "patient", "subject",
            "library", "chip", "lane", "replicate", "condition",
            "sample_id", "batch_id", "donor_id", "patient_id",
        ]

    result = BatchKeyDetection()
    candidates = []

    for col in adata.obs.columns:
        col_lower = col.lower()

        # Skip QC metrics and numeric-only columns
        if col_lower in ['n_genes', 'n_counts', 'total_counts', 'pct_counts_mt',
                         'n_genes_by_counts', 'total_counts_mt']:
            continue

        # Check if column is categorical or can be treated as categorical
        try:
            if hasattr(adata.obs[col], 'cat'):
                n_categories = len(adata.obs[col].cat.categories)
            else:
                n_categories = adata.obs[col].nunique()
        except Exception:
            continue

        # Skip if too many or too few categories
        if n_categories < 2 or n_categories > adata.n_obs * 0.5:
            continue

        # Calculate score
        score = 0.0
        reasons = []

        # Check for pattern match
        for pattern in candidate_patterns:
            if pattern in col_lower:
                score += 0.5
                reasons.append(f"matches pattern '{pattern}'")
                break

        # Prefer reasonable number of categories (2-50)
        if 2 <= n_categories <= 50:
            score += 0.3
            reasons.append(f"{n_categories} categories")
        elif n_categories > 50:
            score += 0.1
            reasons.append(f"{n_categories} categories (many)")

        # Check for balanced distribution
        value_counts = adata.obs[col].value_counts()
        min_size = value_counts.min()
        max_size = value_counts.max()
        balance_ratio = min_size / max_size if max_size > 0 else 0

        if balance_ratio > 0.1:
            score += 0.2
            reasons.append("relatively balanced")

        if score > 0:
            candidates.append({
                "column": col,
                "score": score,
                "n_categories": n_categories,
                "reasons": reasons,
                "sizes": value_counts.to_dict(),
            })

    # Sort by score
    candidates.sort(key=lambda x: x["score"], reverse=True)
    result.candidates = candidates

    if candidates:
        best = candidates[0]
        result.detected_key = best["column"]
        result.confidence = min(best["score"], 1.0)
        result.n_batches = best["n_categories"]
        result.batch_sizes = best["sizes"]
        result.reasoning = f"Selected '{best['column']}' because: {', '.join(best['reasons'])}"
    else:
        result.reasoning = "No suitable batch key candidates found"

    return result


def detect_mt_genes(adata: sc.AnnData) -> MTGeneDetection:
    """
    Detect mitochondrial gene pattern.

    Checks for common patterns: MT-, mt-, Mt-

    Args:
        adata: AnnData object

    Returns:
        MTGeneDetection with detected pattern
    """
    result = MTGeneDetection()

    patterns = [
        ("MT-", r"^MT-"),
        ("mt-", r"^mt-"),
        ("Mt-", r"^Mt-"),
        ("MT_", r"^MT_"),
        ("mt_", r"^mt_"),
    ]

    gene_names = adata.var_names.tolist()

    for pattern_name, regex in patterns:
        matches = [g for g in gene_names if re.match(regex, g)]
        if matches:
            result.pattern = pattern_name
            result.n_mt_genes = len(matches)
            result.mt_gene_examples = matches[:5]
            result.detected = True
            break

    # Calculate MT percentage if detected
    if result.detected and result.n_mt_genes > 0:
        mt_mask = adata.var_names.str.match(f"^{re.escape(result.pattern)}")

        # Need to calculate QC metrics if not present
        if 'pct_counts_mt' not in adata.obs.columns:
            # Temporarily calculate
            adata_copy = adata.copy()
            adata_copy.var['mt'] = mt_mask
            sc.pp.calculate_qc_metrics(
                adata_copy,
                qc_vars=['mt'],
                percent_top=None,
                log1p=False,
                inplace=True
            )
            result.mt_pct_mean = float(adata_copy.obs['pct_counts_mt'].mean())
            result.mt_pct_median = float(adata_copy.obs['pct_counts_mt'].median())
        else:
            result.mt_pct_mean = float(adata.obs['pct_counts_mt'].mean())
            result.mt_pct_median = float(adata.obs['pct_counts_mt'].median())

    return result


def detect_qc_status(adata: sc.AnnData) -> ExistingQCStatus:
    """
    Detect if QC and processing have already been applied.

    Checks for:
    - QC metrics in obs
    - Signs of filtering
    - Normalization
    - HVG selection
    - PCA/UMAP
    - Clustering

    Args:
        adata: AnnData object

    Returns:
        ExistingQCStatus with detection results
    """
    result = ExistingQCStatus()
    indicators = {}

    # Check for QC metrics
    qc_columns = ['n_genes', 'n_counts', 'n_genes_by_counts', 'total_counts',
                  'pct_counts_mt', 'total_counts_mt']
    present_qc = [c for c in qc_columns if c in adata.obs.columns]
    result.qc_metrics_present = len(present_qc) > 0
    indicators['qc_columns_present'] = present_qc

    # Check if data appears filtered (low gene counts suggest filtering)
    if 'n_genes' in adata.obs.columns or 'n_genes_by_counts' in adata.obs.columns:
        gene_col = 'n_genes' if 'n_genes' in adata.obs.columns else 'n_genes_by_counts'
        min_genes = adata.obs[gene_col].min()
        if min_genes >= 100:  # Typical filtering threshold
            result.appears_filtered = True
            indicators['min_genes_per_cell'] = int(min_genes)

    # Check for normalization
    # If max value is around log scale (~10-15), likely log-transformed
    if hasattr(adata.X, 'max'):
        try:
            max_val = float(adata.X.max())
        except Exception:
            from scipy import sparse
            if sparse.issparse(adata.X):
                max_val = float(adata.X.max())
            else:
                max_val = float(np.max(adata.X))
    else:
        max_val = float(np.max(adata.X))

    indicators['max_expression_value'] = max_val

    if max_val < 20:  # Log-transformed data typically has max ~10-15
        result.log_transformed = True
        result.normalized = True
    elif 'counts' in adata.layers:
        # If counts layer exists, main matrix is likely normalized
        result.normalized = True

    # Check for HVG selection
    if 'highly_variable' in adata.var.columns:
        result.hvgs_selected = True
        indicators['n_hvgs'] = int(adata.var['highly_variable'].sum())

    # Check for PCA
    if 'X_pca' in adata.obsm:
        result.pca_computed = True
        indicators['n_pcs'] = adata.obsm['X_pca'].shape[1]

    # Check for neighbors
    if 'neighbors' in adata.uns or 'connectivities' in adata.obsp:
        result.neighbors_computed = True

    # Check for UMAP
    if 'X_umap' in adata.obsm:
        result.umap_computed = True

    # Check for clustering
    cluster_columns = [c for c in adata.obs.columns
                       if c in ['leiden', 'louvain', 'cluster', 'clusters']]
    if cluster_columns:
        result.clustered = True
        indicators['cluster_columns'] = cluster_columns

    result.indicators = indicators

    # Generate recommendation
    if result.clustered:
        result.recommendation = "Data appears fully processed. Consider re-analysis only if needed."
    elif result.pca_computed:
        result.recommendation = "Data has PCA. Can proceed with clustering and integration."
    elif result.normalized:
        result.recommendation = "Data is normalized. Can proceed with HVG selection and dimensionality reduction."
    elif result.qc_metrics_present and result.appears_filtered:
        result.recommendation = "Data has QC metrics and appears filtered. Can proceed with normalization."
    else:
        result.recommendation = "Data appears raw. Full QC pipeline recommended."

    return result


def profile_data(
    input_path: str,
    adata: Optional[sc.AnnData] = None,
) -> DataProfile:
    """
    Generate comprehensive data profile for AI Agent.

    This is the main entry point for data profiling. It combines all
    detection functions to create a complete profile.

    Args:
        input_path: Path to input data
        adata: Optional pre-loaded AnnData (will load if not provided)

    Returns:
        DataProfile with all detection results and recommendations
    """
    logger.info(f"Profiling data: {input_path}")

    profile = DataProfile(input_path=str(input_path))

    # Detect file format
    profile.file_detection = detect_file_format(input_path)
    logger.info(f"Detected format: {profile.file_detection.format}")

    # Load data if not provided
    if adata is None:
        from .pipeline import load_data
        try:
            adata = load_data(input_path)
        except Exception as e:
            profile.warnings.append(f"Failed to load data: {str(e)}")
            return profile

    # Basic stats
    profile.n_cells = adata.n_obs
    profile.n_genes = adata.n_vars

    # Calculate sparsity
    from scipy import sparse
    if sparse.issparse(adata.X):
        n_nonzero = adata.X.nnz
    else:
        n_nonzero = np.count_nonzero(adata.X)
    total_elements = adata.n_obs * adata.n_vars
    profile.sparsity = 1.0 - (n_nonzero / total_elements) if total_elements > 0 else 0.0

    # Available annotations
    profile.obs_columns = list(adata.obs.columns)
    profile.var_columns = list(adata.var.columns)
    profile.layers = list(adata.layers.keys())
    profile.obsm_keys = list(adata.obsm.keys())
    profile.uns_keys = list(adata.uns.keys())

    # Detect batch key
    profile.batch_key_detection = detect_batch_key(adata)
    if profile.batch_key_detection.detected_key:
        profile.recommended_batch_key = profile.batch_key_detection.detected_key
        logger.info(f"Detected batch key: {profile.recommended_batch_key}")

    # Detect MT genes
    profile.mt_gene_detection = detect_mt_genes(adata)
    if profile.mt_gene_detection.detected:
        profile.recommended_mt_pattern = profile.mt_gene_detection.pattern
        logger.info(f"Detected MT pattern: {profile.recommended_mt_pattern}")

    # Detect QC status
    profile.qc_status = detect_qc_status(adata)
    logger.info(f"QC status: {profile.qc_status.recommendation}")

    # Calculate median stats
    if 'n_genes' in adata.obs.columns:
        profile.median_genes_per_cell = float(adata.obs['n_genes'].median())
    elif 'n_genes_by_counts' in adata.obs.columns:
        profile.median_genes_per_cell = float(adata.obs['n_genes_by_counts'].median())

    if 'total_counts' in adata.obs.columns:
        profile.median_counts_per_cell = float(adata.obs['total_counts'].median())
    elif 'n_counts' in adata.obs.columns:
        profile.median_counts_per_cell = float(adata.obs['n_counts'].median())

    # Generate recommendations
    _generate_recommendations(profile, adata)

    return profile


def _generate_recommendations(profile: DataProfile, adata: sc.AnnData) -> None:
    """Generate recommended actions based on profile."""

    qc_status = profile.qc_status

    # Warnings
    if profile.n_cells < 100:
        profile.warnings.append(f"Very few cells ({profile.n_cells}). Results may be unreliable.")

    if profile.n_cells > 500000:
        profile.warnings.append(f"Large dataset ({profile.n_cells:,} cells). Consider subsampling or using GPU.")

    if profile.sparsity > 0.99:
        profile.warnings.append(f"Very sparse data ({profile.sparsity:.1%}). Check data quality.")

    if profile.mt_gene_detection and profile.mt_gene_detection.mt_pct_mean > 20:
        profile.warnings.append(f"High MT% mean ({profile.mt_gene_detection.mt_pct_mean:.1f}%). May indicate poor cell quality.")

    # Recommended actions
    if not qc_status.qc_metrics_present:
        profile.recommended_actions.append("Calculate QC metrics (n_genes, n_counts, pct_counts_mt)")

    if not qc_status.appears_filtered:
        profile.recommended_actions.append("Apply QC filtering")

    if not qc_status.normalized:
        profile.recommended_actions.append("Normalize and log-transform")

    if not qc_status.hvgs_selected:
        profile.recommended_actions.append("Select highly variable genes")

    if not qc_status.pca_computed:
        profile.recommended_actions.append("Run PCA")

    if not qc_status.neighbors_computed:
        profile.recommended_actions.append("Compute nearest neighbors graph")

    if not qc_status.umap_computed:
        profile.recommended_actions.append("Compute UMAP embedding")

    if not qc_status.clustered:
        profile.recommended_actions.append("Run clustering")

    # Batch integration recommendation
    if profile.batch_key_detection and profile.batch_key_detection.n_batches > 1:
        if profile.batch_key_detection.n_batches >= 2:
            profile.recommended_actions.append(
                f"Consider batch integration (detected {profile.batch_key_detection.n_batches} batches)"
            )

    # Doublet detection recommendation
    if profile.n_cells > 500 and not qc_status.appears_filtered:
        profile.recommended_actions.append("Run doublet detection")
