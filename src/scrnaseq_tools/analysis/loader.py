"""
Smart data loading module for AI Agent-friendly scRNA-seq analysis.

This module provides intelligent data loading capabilities:
- Automatic sample discovery in directories
- Multiple format support (h5ad, 10x MTX, 10x H5)
- Multi-sample loading and merging
- Sample metadata extraction

Usage:
    from scrnaseq_tools.analysis.loader import smart_load

    result = smart_load("/path/to/data")
    adata = result.adata
"""

import os
import re
import logging
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any, Union

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


@dataclass
class SampleInfo:
    """Information about a single sample."""

    sample_id: str
    path: str
    format: str  # "h5ad", "10x_mtx", "10x_h5"
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class LoadResult:
    """Result of data loading operation."""

    adata: Any  # sc.AnnData
    samples: List[SampleInfo] = field(default_factory=list)
    n_samples: int = 1
    sample_key: str = "sample"
    total_cells: int = 0
    total_genes: int = 0
    load_warnings: List[str] = field(default_factory=list)
    load_errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_samples": self.n_samples,
            "sample_key": self.sample_key,
            "total_cells": self.total_cells,
            "total_genes": self.total_genes,
            "samples": [s.to_dict() for s in self.samples],
            "load_warnings": self.load_warnings,
            "load_errors": self.load_errors,
        }

    def summary(self) -> str:
        lines = [
            f"Loaded {self.n_samples} sample(s)",
            f"Total cells: {self.total_cells:,}",
            f"Total genes: {self.total_genes:,}",
        ]
        if self.n_samples > 1:
            lines.append(f"Sample key: '{self.sample_key}'")
            for sample in self.samples:
                lines.append(f"  - {sample.sample_id}: {sample.n_cells:,} cells")
        if self.load_warnings:
            lines.append("Warnings:")
            for w in self.load_warnings:
                lines.append(f"  - {w}")
        return "\n".join(lines)


@dataclass
class MultiSampleConfig:
    """Configuration for multi-sample loading."""

    sample_key: str = "sample"  # Column name for sample IDs in merged data
    min_cells_per_sample: int = 10  # Minimum cells to include a sample
    extract_metadata_from_path: bool = True  # Extract sample info from path
    metadata_pattern: Optional[str] = None  # Regex pattern for metadata extraction
    join_vars: str = "inner"  # How to join var (genes): "inner" or "outer"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def discover_samples(
    directory: str,
    recursive: bool = True,
    max_depth: int = 3,
) -> List[SampleInfo]:
    """
    Discover samples in a directory.

    Searches for:
    - .h5ad files
    - 10x MTX directories (containing matrix.mtx)
    - 10x H5 files (.h5)

    Args:
        directory: Directory to search
        recursive: Whether to search recursively
        max_depth: Maximum recursion depth

    Returns:
        List of SampleInfo objects
    """
    directory = Path(directory)
    samples = []

    if not directory.exists():
        logger.warning(f"Directory does not exist: {directory}")
        return samples

    def _search_dir(path: Path, depth: int = 0):
        if depth > max_depth:
            return

        # Check for h5ad files
        for h5ad_file in path.glob("*.h5ad"):
            sample_id = h5ad_file.stem
            samples.append(SampleInfo(
                sample_id=sample_id,
                path=str(h5ad_file),
                format="h5ad",
            ))

        # Check for h5 files (10x format)
        for h5_file in path.glob("*.h5"):
            if h5_file.stem.endswith("_filtered_feature_bc_matrix") or \
               h5_file.stem.endswith("_raw_feature_bc_matrix"):
                # CellRanger output
                sample_id = h5_file.stem.replace("_filtered_feature_bc_matrix", "") \
                                        .replace("_raw_feature_bc_matrix", "")
            else:
                sample_id = h5_file.stem
            samples.append(SampleInfo(
                sample_id=sample_id,
                path=str(h5_file),
                format="10x_h5",
            ))

        # Check for 10x MTX directories
        mtx_files = list(path.glob("*matrix*.mtx*"))
        if mtx_files:
            # This directory contains MTX files
            sample_id = path.name
            samples.append(SampleInfo(
                sample_id=sample_id,
                path=str(path),
                format="10x_mtx",
            ))
        else:
            # Check for filtered_feature_bc_matrix subdirectory
            for subdir_name in ["filtered_feature_bc_matrix", "raw_feature_bc_matrix"]:
                subdir = path / subdir_name
                if subdir.exists() and list(subdir.glob("*matrix*.mtx*")):
                    sample_id = path.name
                    samples.append(SampleInfo(
                        sample_id=sample_id,
                        path=str(subdir),
                        format="10x_mtx",
                    ))
                    break

        # Recurse into subdirectories
        if recursive and depth < max_depth:
            for subdir in path.iterdir():
                if subdir.is_dir() and not subdir.name.startswith('.'):
                    # Skip common non-data directories
                    if subdir.name in ["outs", "filtered_feature_bc_matrix",
                                       "raw_feature_bc_matrix"]:
                        continue
                    _search_dir(subdir, depth + 1)

    _search_dir(directory)

    # Remove duplicates (same path)
    seen_paths = set()
    unique_samples = []
    for sample in samples:
        if sample.path not in seen_paths:
            seen_paths.add(sample.path)
            unique_samples.append(sample)

    logger.info(f"Discovered {len(unique_samples)} samples in {directory}")
    return unique_samples


def load_single_sample(
    sample_info: SampleInfo,
) -> sc.AnnData:
    """
    Load a single sample.

    Args:
        sample_info: SampleInfo with path and format

    Returns:
        AnnData object
    """
    path = Path(sample_info.path)
    logger.info(f"Loading sample '{sample_info.sample_id}' from {path}")

    if sample_info.format == "h5ad":
        adata = sc.read_h5ad(path)
    elif sample_info.format == "10x_h5":
        adata = sc.read_10x_h5(str(path))
    elif sample_info.format == "10x_mtx":
        adata = sc.read_10x_mtx(str(path), var_names='gene_symbols', cache=True)
    else:
        raise ValueError(f"Unsupported format: {sample_info.format}")

    # Ensure unique var_names
    adata.var_names_make_unique()

    # Update sample info
    sample_info.n_cells = adata.n_obs
    sample_info.n_genes = adata.n_vars

    return adata


def load_multiple_samples(
    samples: List[SampleInfo],
    config: Optional[MultiSampleConfig] = None,
) -> LoadResult:
    """
    Load and merge multiple samples.

    Args:
        samples: List of SampleInfo objects
        config: MultiSampleConfig for merging options

    Returns:
        LoadResult with merged AnnData
    """
    if config is None:
        config = MultiSampleConfig()

    result = LoadResult(
        adata=None,
        sample_key=config.sample_key,
    )

    if not samples:
        result.load_errors.append("No samples provided")
        return result

    adatas = []
    valid_samples = []

    for sample_info in samples:
        try:
            adata = load_single_sample(sample_info)

            # Check minimum cells
            if adata.n_obs < config.min_cells_per_sample:
                result.load_warnings.append(
                    f"Sample '{sample_info.sample_id}' has only {adata.n_obs} cells "
                    f"(min: {config.min_cells_per_sample}). Skipping."
                )
                continue

            # Add sample ID to obs
            adata.obs[config.sample_key] = sample_info.sample_id

            adatas.append(adata)
            valid_samples.append(sample_info)

        except Exception as e:
            result.load_errors.append(f"Failed to load '{sample_info.sample_id}': {str(e)}")
            logger.error(f"Failed to load sample {sample_info.sample_id}: {e}")

    if not adatas:
        result.load_errors.append("No samples were successfully loaded")
        return result

    # Merge samples
    logger.info(f"Merging {len(adatas)} samples...")

    if len(adatas) == 1:
        merged = adatas[0]
    else:
        # Use scanpy's concatenate
        merged = sc.concat(
            adatas,
            join=config.join_vars,
            label=config.sample_key,
            keys=[s.sample_id for s in valid_samples],
            index_unique="-",
        )

    # Ensure unique obs_names
    merged.obs_names_make_unique()

    result.adata = merged
    result.samples = valid_samples
    result.n_samples = len(valid_samples)
    result.total_cells = merged.n_obs
    result.total_genes = merged.n_vars

    logger.info(f"Merged data: {result.total_cells:,} cells x {result.total_genes:,} genes")

    return result


def smart_load(
    input_path: str,
    config: Optional[MultiSampleConfig] = None,
    profile_first: bool = True,
) -> LoadResult:
    """
    Smart loading: automatically detect format and load data.

    This is the main entry point for data loading. It automatically:
    1. Detects the input format (single file, multi-sample directory)
    2. Discovers samples if needed
    3. Loads and merges data

    Args:
        input_path: Path to file or directory
        config: Optional MultiSampleConfig for multi-sample loading
        profile_first: Whether to profile the data first

    Returns:
        LoadResult with loaded data
    """
    if config is None:
        config = MultiSampleConfig()

    path = Path(input_path)
    result = LoadResult(adata=None, sample_key=config.sample_key)

    if not path.exists():
        result.load_errors.append(f"Path does not exist: {input_path}")
        return result

    logger.info(f"Smart loading: {input_path}")

    # Profile first to understand the data
    if profile_first:
        from .profile import detect_file_format
        file_detection = detect_file_format(input_path)
        logger.info(f"Detected format: {file_detection.format}")
    else:
        file_detection = None

    # Single file
    if path.is_file():
        if path.suffix == ".h5ad":
            sample_info = SampleInfo(
                sample_id=path.stem,
                path=str(path),
                format="h5ad",
            )
        elif path.suffix == ".h5":
            sample_info = SampleInfo(
                sample_id=path.stem,
                path=str(path),
                format="10x_h5",
            )
        else:
            result.load_errors.append(f"Unsupported file format: {path.suffix}")
            return result

        try:
            adata = load_single_sample(sample_info)
            result.adata = adata
            result.samples = [sample_info]
            result.n_samples = 1
            result.total_cells = adata.n_obs
            result.total_genes = adata.n_vars
        except Exception as e:
            result.load_errors.append(f"Failed to load file: {str(e)}")

        return result

    # Directory
    if path.is_dir():
        # Check if it's a single 10x directory
        mtx_files = list(path.glob("*matrix*.mtx*"))
        filtered_dir = path / "filtered_feature_bc_matrix"
        outs_filtered = path / "outs" / "filtered_feature_bc_matrix"

        if mtx_files:
            # Single 10x MTX directory
            sample_info = SampleInfo(
                sample_id=path.name,
                path=str(path),
                format="10x_mtx",
            )
            try:
                adata = load_single_sample(sample_info)
                result.adata = adata
                result.samples = [sample_info]
                result.n_samples = 1
                result.total_cells = adata.n_obs
                result.total_genes = adata.n_vars
            except Exception as e:
                result.load_errors.append(f"Failed to load 10x data: {str(e)}")
            return result

        elif filtered_dir.exists():
            # CellRanger output with filtered_feature_bc_matrix
            sample_info = SampleInfo(
                sample_id=path.name,
                path=str(filtered_dir),
                format="10x_mtx",
            )
            try:
                adata = load_single_sample(sample_info)
                result.adata = adata
                result.samples = [sample_info]
                result.n_samples = 1
                result.total_cells = adata.n_obs
                result.total_genes = adata.n_vars
            except Exception as e:
                result.load_errors.append(f"Failed to load 10x data: {str(e)}")
            return result

        elif outs_filtered.exists():
            # CellRanger outs directory
            sample_info = SampleInfo(
                sample_id=path.name,
                path=str(outs_filtered),
                format="10x_mtx",
            )
            try:
                adata = load_single_sample(sample_info)
                result.adata = adata
                result.samples = [sample_info]
                result.n_samples = 1
                result.total_cells = adata.n_obs
                result.total_genes = adata.n_vars
            except Exception as e:
                result.load_errors.append(f"Failed to load 10x data: {str(e)}")
            return result

        # Check for single h5ad file in directory
        h5ad_files = list(path.glob("*.h5ad"))
        if len(h5ad_files) == 1:
            sample_info = SampleInfo(
                sample_id=h5ad_files[0].stem,
                path=str(h5ad_files[0]),
                format="h5ad",
            )
            try:
                adata = load_single_sample(sample_info)
                result.adata = adata
                result.samples = [sample_info]
                result.n_samples = 1
                result.total_cells = adata.n_obs
                result.total_genes = adata.n_vars
            except Exception as e:
                result.load_errors.append(f"Failed to load h5ad: {str(e)}")
            return result

        # Multi-sample directory - discover and load
        samples = discover_samples(str(path))

        if not samples:
            result.load_errors.append(f"No samples found in directory: {input_path}")
            return result

        if len(samples) == 1:
            # Single sample found
            try:
                adata = load_single_sample(samples[0])
                result.adata = adata
                result.samples = samples
                result.n_samples = 1
                result.total_cells = adata.n_obs
                result.total_genes = adata.n_vars
            except Exception as e:
                result.load_errors.append(f"Failed to load sample: {str(e)}")
            return result

        # Multiple samples - load and merge
        return load_multiple_samples(samples, config)

    result.load_errors.append(f"Unable to process input: {input_path}")
    return result


def extract_sample_metadata(
    sample_id: str,
    pattern: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Extract metadata from sample ID using regex pattern.

    Example patterns:
    - "(?P<patient>P\\d+)_(?P<timepoint>T\\d+)" matches "P001_T1"
    - "(?P<condition>\\w+)_(?P<replicate>\\d+)" matches "control_1"

    Args:
        sample_id: Sample ID string
        pattern: Regex pattern with named groups

    Returns:
        Dictionary of extracted metadata
    """
    if pattern is None:
        return {}

    try:
        match = re.match(pattern, sample_id)
        if match:
            return match.groupdict()
    except re.error as e:
        logger.warning(f"Invalid regex pattern: {e}")

    return {}
