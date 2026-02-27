"""
Agent-friendly report generation module for scRNA-seq analysis.

This module generates structured reports optimized for AI agent consumption:
- JSON-serializable report sections
- Markdown formatting for human readability
- LLM-friendly summarization
- Actionable next steps

Usage:
    from scrnaseq_tools.analysis.report import generate_agent_report, save_agent_report

    report = generate_agent_report(adata, run_metadata)
    save_agent_report(report, "output/")
"""

import json
import logging
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


@dataclass
class QCReportSection:
    """QC section of the agent report."""

    cells_before: int = 0
    cells_after: int = 0
    genes_before: int = 0
    genes_after: int = 0
    retention_rate: float = 0.0

    # Thresholds used
    min_genes_threshold: float = 0.0
    max_genes_threshold: Optional[float] = None
    max_mt_pct_threshold: float = 0.0
    min_cells_per_gene: int = 0

    # QC method info
    qc_method: str = "mad"
    nmads: float = 3.0

    # Distribution stats
    median_genes_per_cell: float = 0.0
    median_counts_per_cell: float = 0.0
    median_mt_pct: float = 0.0

    # Doublet detection (if performed)
    doublets_detected: int = 0
    doublet_rate: float = 0.0
    doublet_method: Optional[str] = None

    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class NormalizationReportSection:
    """Normalization section of the agent report."""

    method: str = "size_factor_log1p"
    target_sum: float = 10000.0

    # HVG selection
    hvg_method: str = "seurat_v3"
    n_hvgs: int = 0
    hvg_top_genes: List[str] = field(default_factory=list)

    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DimensionalityReportSection:
    """Dimensionality reduction section of the agent report."""

    # PCA
    n_pcs_computed: int = 0
    n_pcs_selected: int = 0
    pc_selection_method: str = "elbow"
    variance_explained_by_selected: float = 0.0
    elbow_point: Optional[int] = None

    # Neighbors
    n_neighbors: int = 15

    # UMAP
    umap_computed: bool = False

    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ClusteringReportSection:
    """Clustering section of the agent report."""

    method: str = "leiden"
    resolution: float = 0.5
    n_clusters: int = 0

    # Cluster statistics
    cluster_sizes: Dict[str, int] = field(default_factory=dict)
    min_cluster_size: int = 0
    max_cluster_size: int = 0
    median_cluster_size: float = 0.0

    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class IntegrationReportSection:
    """Integration section of the agent report."""

    batch_key: str = ""
    n_batches: int = 0
    batch_sizes: Dict[str, int] = field(default_factory=dict)

    # Methods run
    methods_run: List[str] = field(default_factory=list)
    best_method: Optional[str] = None

    # Benchmark scores
    benchmark_scores: Dict[str, Dict[str, float]] = field(default_factory=dict)

    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class AgentReport:
    """Complete agent-friendly report for scRNA-seq analysis."""

    # Metadata
    timestamp: str = ""
    pipeline_version: str = "1.0"
    input_path: str = ""
    output_dir: str = ""

    # Data overview
    final_cells: int = 0
    final_genes: int = 0

    # Report sections
    qc: Optional[QCReportSection] = None
    normalization: Optional[NormalizationReportSection] = None
    dimensionality: Optional[DimensionalityReportSection] = None
    clustering: Optional[ClusteringReportSection] = None
    integration: Optional[IntegrationReportSection] = None

    # Agent-specific fields
    data_quality_score: float = 0.0  # 0-1 score
    quality_flags: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    next_steps: List[Dict[str, Any]] = field(default_factory=list)

    # Files created
    files_created: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {
            "metadata": {
                "timestamp": self.timestamp,
                "pipeline_version": self.pipeline_version,
                "input_path": self.input_path,
                "output_dir": self.output_dir,
            },
            "summary": {
                "final_cells": self.final_cells,
                "final_genes": self.final_genes,
                "data_quality_score": self.data_quality_score,
                "quality_flags": self.quality_flags,
            },
            "warnings": self.warnings,
            "next_steps": self.next_steps,
            "files_created": self.files_created,
        }

        if self.qc:
            result["qc"] = self.qc.to_dict()
        if self.normalization:
            result["normalization"] = self.normalization.to_dict()
        if self.dimensionality:
            result["dimensionality"] = self.dimensionality.to_dict()
        if self.clustering:
            result["clustering"] = self.clustering.to_dict()
        if self.integration:
            result["integration"] = self.integration.to_dict()

        return result

    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent, default=str)

    def to_markdown(self) -> str:
        """Convert to Markdown format."""
        lines = [
            "# scRNA-seq Analysis Report",
            "",
            f"**Generated:** {self.timestamp}",
            f"**Input:** `{self.input_path}`",
            f"**Output:** `{self.output_dir}`",
            "",
            "## Summary",
            "",
            f"- **Final cells:** {self.final_cells:,}",
            f"- **Final genes:** {self.final_genes:,}",
            f"- **Data quality score:** {self.data_quality_score:.2f}/1.0",
            "",
        ]

        if self.quality_flags:
            lines.append("### Quality Flags")
            for flag in self.quality_flags:
                lines.append(f"- {flag}")
            lines.append("")

        if self.warnings:
            lines.append("### Warnings")
            for warning in self.warnings:
                lines.append(f"- ⚠️ {warning}")
            lines.append("")

        # QC Section
        if self.qc:
            lines.extend([
                "## Quality Control",
                "",
                f"- **Method:** {self.qc.qc_method} (n_mads={self.qc.nmads})",
                f"- **Cells:** {self.qc.cells_before:,} → {self.qc.cells_after:,} "
                f"({self.qc.retention_rate:.1%} retained)",
                f"- **Genes:** {self.qc.genes_before:,} → {self.qc.genes_after:,}",
                "",
                "### Thresholds",
                f"- Min genes: {self.qc.min_genes_threshold:.0f}",
            ])
            if self.qc.max_genes_threshold:
                lines.append(f"- Max genes: {self.qc.max_genes_threshold:.0f}")
            lines.extend([
                f"- Max MT%: {self.qc.max_mt_pct_threshold:.1f}%",
                "",
            ])

            if self.qc.doublet_method:
                lines.extend([
                    "### Doublet Detection",
                    f"- **Method:** {self.qc.doublet_method}",
                    f"- **Doublets detected:** {self.qc.doublets_detected:,} "
                    f"({self.qc.doublet_rate:.1%})",
                    "",
                ])

        # Normalization Section
        if self.normalization:
            lines.extend([
                "## Normalization",
                "",
                f"- **Method:** {self.normalization.method}",
                f"- **Target sum:** {self.normalization.target_sum:,.0f}",
                f"- **HVG method:** {self.normalization.hvg_method}",
                f"- **HVGs selected:** {self.normalization.n_hvgs:,}",
                "",
            ])

        # Dimensionality Section
        if self.dimensionality:
            lines.extend([
                "## Dimensionality Reduction",
                "",
                f"- **PCs computed:** {self.dimensionality.n_pcs_computed}",
                f"- **PCs selected:** {self.dimensionality.n_pcs_selected} "
                f"({self.dimensionality.pc_selection_method})",
                f"- **Variance explained:** {self.dimensionality.variance_explained_by_selected:.1%}",
                f"- **Neighbors:** {self.dimensionality.n_neighbors}",
                f"- **UMAP:** {'Yes' if self.dimensionality.umap_computed else 'No'}",
                "",
            ])

        # Clustering Section
        if self.clustering:
            lines.extend([
                "## Clustering",
                "",
                f"- **Method:** {self.clustering.method}",
                f"- **Resolution:** {self.clustering.resolution}",
                f"- **Clusters found:** {self.clustering.n_clusters}",
                f"- **Cluster sizes:** {self.clustering.min_cluster_size} - "
                f"{self.clustering.max_cluster_size} (median: {self.clustering.median_cluster_size:.0f})",
                "",
            ])

        # Integration Section
        if self.integration and self.integration.methods_run:
            lines.extend([
                "## Integration",
                "",
                f"- **Batch key:** `{self.integration.batch_key}`",
                f"- **Batches:** {self.integration.n_batches}",
                f"- **Methods run:** {', '.join(self.integration.methods_run)}",
            ])
            if self.integration.best_method:
                lines.append(f"- **Best method:** {self.integration.best_method}")
            lines.append("")

        # Next Steps
        if self.next_steps:
            lines.extend([
                "## Recommended Next Steps",
                "",
            ])
            for i, step in enumerate(self.next_steps, 1):
                lines.append(f"{i}. **{step['action']}**")
                if 'description' in step:
                    lines.append(f"   {step['description']}")
                if 'code_example' in step:
                    lines.append(f"   ```python")
                    lines.append(f"   {step['code_example']}")
                    lines.append(f"   ```")
            lines.append("")

        # Files Created
        if self.files_created:
            lines.extend([
                "## Files Created",
                "",
            ])
            for f in self.files_created:
                lines.append(f"- `{f}`")
            lines.append("")

        return "\n".join(lines)


def generate_agent_report(
    adata: sc.AnnData,
    run_metadata: Dict[str, Any],
    input_path: str = "",
    output_dir: str = "",
) -> AgentReport:
    """
    Generate comprehensive agent report from pipeline results.

    Args:
        adata: Processed AnnData object
        run_metadata: Metadata dict from pipeline run
        input_path: Original input path
        output_dir: Output directory

    Returns:
        AgentReport with all sections populated
    """
    report = AgentReport(
        timestamp=datetime.now().isoformat(),
        pipeline_version=run_metadata.get("pipeline_version", "1.0"),
        input_path=input_path or run_metadata.get("input_path", ""),
        output_dir=output_dir,
        final_cells=adata.n_obs,
        final_genes=adata.n_vars,
    )

    # QC Section
    qc_section = QCReportSection()
    qc_section.cells_before = run_metadata.get("cells_before_qc", run_metadata.get("input_cells", 0))
    qc_section.cells_after = run_metadata.get("cells_after_qc", adata.n_obs)
    qc_section.genes_before = run_metadata.get("genes_before_qc", run_metadata.get("input_genes", 0))
    qc_section.genes_after = run_metadata.get("genes_after_qc", adata.n_vars)

    if qc_section.cells_before > 0:
        qc_section.retention_rate = qc_section.cells_after / qc_section.cells_before

    if "qc_thresholds" in run_metadata:
        thresholds = run_metadata["qc_thresholds"]
        qc_section.min_genes_threshold = thresholds.get("min_genes", 0)
        qc_section.max_genes_threshold = thresholds.get("max_genes")
        qc_section.max_mt_pct_threshold = thresholds.get("max_mt_pct", 0)

    # Doublet info
    if "doublet_score" in adata.obs.columns:
        qc_section.doublet_method = "scrublet"
        if "predicted_doublet" in adata.obs.columns:
            qc_section.doublets_detected = int(adata.obs["predicted_doublet"].sum())
            qc_section.doublet_rate = qc_section.doublets_detected / adata.n_obs

    report.qc = qc_section

    # Normalization Section
    norm_section = NormalizationReportSection()
    norm_section.method = run_metadata.get("normalization_method", "size_factor_log1p")
    norm_section.target_sum = run_metadata.get("target_sum", 10000)
    norm_section.hvg_method = run_metadata.get("hvg_method", "seurat_v3")
    norm_section.n_hvgs = run_metadata.get("n_hvgs", 0)

    if "highly_variable" in adata.var.columns:
        hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()
        norm_section.n_hvgs = len(hvg_genes)
        norm_section.hvg_top_genes = hvg_genes[:10]

    report.normalization = norm_section

    # Dimensionality Section
    dim_section = DimensionalityReportSection()
    if "X_pca" in adata.obsm:
        dim_section.n_pcs_computed = adata.obsm["X_pca"].shape[1]

    if "pc_selection" in run_metadata:
        pc_info = run_metadata["pc_selection"]
        dim_section.n_pcs_selected = pc_info.get("n_pcs_selected", dim_section.n_pcs_computed)
        dim_section.pc_selection_method = pc_info.get("method", "elbow")
        dim_section.elbow_point = pc_info.get("elbow_point")

        cum_var = pc_info.get("cumulative_variance", [])
        if cum_var and len(cum_var) >= dim_section.n_pcs_selected:
            dim_section.variance_explained_by_selected = cum_var[dim_section.n_pcs_selected - 1]

    dim_section.n_neighbors = run_metadata.get("n_neighbors", 15)
    dim_section.umap_computed = "X_umap" in adata.obsm

    report.dimensionality = dim_section

    # Clustering Section
    cluster_section = ClusteringReportSection()
    cluster_section.method = run_metadata.get("clustering_method", "leiden")
    cluster_section.resolution = run_metadata.get("resolution", 0.5)
    cluster_section.n_clusters = run_metadata.get("n_clusters", 0)

    cluster_key = "leiden" if "leiden" in adata.obs.columns else "louvain"
    if cluster_key in adata.obs.columns:
        cluster_sizes = adata.obs[cluster_key].value_counts()
        cluster_section.cluster_sizes = cluster_sizes.to_dict()
        cluster_section.n_clusters = len(cluster_sizes)
        cluster_section.min_cluster_size = int(cluster_sizes.min())
        cluster_section.max_cluster_size = int(cluster_sizes.max())
        cluster_section.median_cluster_size = float(cluster_sizes.median())

    report.clustering = cluster_section

    # Integration Section
    if "integration" in run_metadata:
        int_info = run_metadata["integration"]
        int_section = IntegrationReportSection()
        int_section.methods_run = int_info.get("methods", [])
        int_section.best_method = int_info.get("best_method")
        int_section.benchmark_scores = int_info.get("benchmark_scores", {})

        if "batch_key" in run_metadata:
            int_section.batch_key = run_metadata["batch_key"]
            if int_section.batch_key in adata.obs.columns:
                batch_counts = adata.obs[int_section.batch_key].value_counts()
                int_section.n_batches = len(batch_counts)
                int_section.batch_sizes = batch_counts.to_dict()

        report.integration = int_section

    # Calculate quality score
    report.data_quality_score = _calculate_quality_score(report)

    # Generate quality flags
    report.quality_flags = _generate_quality_flags(report)

    # Generate warnings
    report.warnings = _generate_warnings(report, adata)

    # Generate next steps
    report.next_steps = suggest_next_steps(report, adata)

    return report


def _calculate_quality_score(report: AgentReport) -> float:
    """Calculate overall data quality score (0-1)."""
    scores = []

    # Retention rate (higher is better, but too high might indicate insufficient QC)
    if report.qc and report.qc.retention_rate > 0:
        retention = report.qc.retention_rate
        if 0.5 <= retention <= 0.95:
            scores.append(1.0)
        elif retention > 0.95:
            scores.append(0.8)  # Might be under-filtered
        elif retention > 0.3:
            scores.append(0.6)
        else:
            scores.append(0.3)

    # Cell count
    if report.final_cells >= 1000:
        scores.append(1.0)
    elif report.final_cells >= 500:
        scores.append(0.7)
    elif report.final_cells >= 100:
        scores.append(0.4)
    else:
        scores.append(0.2)

    # Clustering quality (balanced cluster sizes)
    if report.clustering and report.clustering.n_clusters > 0:
        if report.clustering.min_cluster_size >= 10:
            scores.append(1.0)
        elif report.clustering.min_cluster_size >= 5:
            scores.append(0.7)
        else:
            scores.append(0.4)

    # Variance explained by PCs
    if report.dimensionality and report.dimensionality.variance_explained_by_selected > 0:
        var_exp = report.dimensionality.variance_explained_by_selected
        if var_exp >= 0.8:
            scores.append(1.0)
        elif var_exp >= 0.6:
            scores.append(0.8)
        else:
            scores.append(0.5)

    return np.mean(scores) if scores else 0.5


def _generate_quality_flags(report: AgentReport) -> List[str]:
    """Generate quality flags for quick assessment."""
    flags = []

    if report.final_cells < 500:
        flags.append("LOW_CELL_COUNT")
    elif report.final_cells > 100000:
        flags.append("LARGE_DATASET")

    if report.qc:
        if report.qc.retention_rate < 0.5:
            flags.append("HIGH_CELL_LOSS")
        if report.qc.doublet_rate > 0.15:
            flags.append("HIGH_DOUBLET_RATE")

    if report.clustering:
        if report.clustering.min_cluster_size < 10:
            flags.append("SMALL_CLUSTERS")
        if report.clustering.n_clusters > 50:
            flags.append("MANY_CLUSTERS")

    if not flags:
        flags.append("GOOD_QUALITY")

    return flags


def _generate_warnings(report: AgentReport, adata: sc.AnnData) -> List[str]:
    """Generate warning messages."""
    warnings = []

    if report.final_cells < 100:
        warnings.append("Very few cells. Results may be unreliable.")

    if report.qc and report.qc.retention_rate < 0.3:
        warnings.append(f"Only {report.qc.retention_rate:.1%} of cells passed QC. "
                       "Consider relaxing thresholds.")

    if report.qc and report.qc.doublet_rate > 0.2:
        warnings.append(f"High doublet rate ({report.qc.doublet_rate:.1%}). "
                       "Data quality may be affected.")

    if report.clustering and report.clustering.min_cluster_size < 5:
        warnings.append(f"Some clusters have very few cells ({report.clustering.min_cluster_size}). "
                       "Consider lowering resolution.")

    return warnings


def suggest_next_steps(
    report: AgentReport,
    adata: sc.AnnData,
) -> List[Dict[str, Any]]:
    """
    Suggest next analysis steps based on current state.

    Args:
        report: AgentReport with current analysis state
        adata: AnnData object

    Returns:
        List of suggested next steps with descriptions and code examples
    """
    steps = []

    # Cell type annotation
    cluster_key = "leiden" if "leiden" in adata.obs.columns else "louvain"
    if cluster_key in adata.obs.columns and "cell_type" not in adata.obs.columns:
        steps.append({
            "action": "Cell Type Annotation",
            "description": "Annotate cell types using marker genes or reference data",
            "priority": "high",
            "code_example": "sc.tl.rank_genes_groups(adata, cluster_key, method='wilcoxon')",
        })

    # Marker genes
    marker_key = f"rank_genes_groups"
    if cluster_key in adata.obs.columns and marker_key not in adata.uns:
        steps.append({
            "action": "Find Marker Genes",
            "description": "Identify differentially expressed genes for each cluster",
            "priority": "high",
            "code_example": f"sc.tl.rank_genes_groups(adata, '{cluster_key}', method='wilcoxon')",
        })

    # Integration suggestion
    if report.integration is None or not report.integration.methods_run:
        # Check if there are batches
        potential_batch_keys = [c for c in adata.obs.columns
                               if c in ["batch", "sample", "donor", "patient"]]
        if potential_batch_keys:
            batch_key = potential_batch_keys[0]
            n_batches = adata.obs[batch_key].nunique()
            if n_batches > 1:
                steps.append({
                    "action": "Batch Integration",
                    "description": f"Integrate {n_batches} batches using scVI or Harmony",
                    "priority": "medium",
                    "code_example": f"from scrnaseq_tools.analysis import run_integration\n"
                                   f"run_integration(adata, batch_key='{batch_key}', methods=['scvi'])",
                })

    # Trajectory analysis
    if "X_umap" in adata.obsm and "dpt_pseudotime" not in adata.obs.columns:
        steps.append({
            "action": "Trajectory Analysis",
            "description": "Infer developmental trajectories if applicable",
            "priority": "low",
            "code_example": "sc.tl.paga(adata)\nsc.pl.paga(adata)",
        })

    # Export for downstream analysis
    steps.append({
        "action": "Export Results",
        "description": "Export processed data for downstream analysis",
        "priority": "low",
        "code_example": "adata.write_h5ad('processed.h5ad')",
    })

    return steps


def save_agent_report(
    report: AgentReport,
    output_dir: str,
    formats: Optional[List[str]] = None,
) -> List[str]:
    """
    Save agent report in multiple formats.

    Args:
        report: AgentReport to save
        output_dir: Output directory
        formats: List of formats ("json", "markdown", "txt"). Default: all

    Returns:
        List of created file paths
    """
    if formats is None:
        formats = ["json", "markdown"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files_created = []

    if "json" in formats:
        json_path = output_dir / "agent_report.json"
        with open(json_path, 'w') as f:
            f.write(report.to_json())
        files_created.append(str(json_path))
        logger.info(f"Saved JSON report: {json_path}")

    if "markdown" in formats:
        md_path = output_dir / "agent_report.md"
        with open(md_path, 'w') as f:
            f.write(report.to_markdown())
        files_created.append(str(md_path))
        logger.info(f"Saved Markdown report: {md_path}")

    return files_created


def format_for_llm(
    report: AgentReport,
    max_tokens: int = 2000,
) -> str:
    """
    Format report for LLM context window.

    Creates a condensed version of the report optimized for
    inclusion in LLM prompts.

    Args:
        report: AgentReport to format
        max_tokens: Approximate maximum tokens (rough estimate)

    Returns:
        Formatted string for LLM context
    """
    lines = [
        "=== scRNA-seq Analysis Summary ===",
        f"Cells: {report.final_cells:,} | Genes: {report.final_genes:,}",
        f"Quality: {report.data_quality_score:.2f}/1.0",
    ]

    if report.quality_flags:
        lines.append(f"Flags: {', '.join(report.quality_flags)}")

    if report.qc:
        lines.append(f"QC: {report.qc.cells_before:,} → {report.qc.cells_after:,} "
                    f"({report.qc.retention_rate:.1%})")
        if report.qc.doublet_method:
            lines.append(f"Doublets: {report.qc.doublets_detected:,} ({report.qc.doublet_rate:.1%})")

    if report.clustering:
        lines.append(f"Clusters: {report.clustering.n_clusters} "
                    f"(sizes: {report.clustering.min_cluster_size}-{report.clustering.max_cluster_size})")

    if report.integration and report.integration.methods_run:
        lines.append(f"Integration: {', '.join(report.integration.methods_run)}")
        if report.integration.best_method:
            lines.append(f"Best method: {report.integration.best_method}")

    if report.warnings:
        lines.append("Warnings: " + "; ".join(report.warnings[:3]))

    if report.next_steps:
        lines.append("Next steps: " + ", ".join(s["action"] for s in report.next_steps[:3]))

    return "\n".join(lines)
