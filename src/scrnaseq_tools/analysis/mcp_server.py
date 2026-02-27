"""
MCP (Model Context Protocol) Server for scrna_qc.

This server exposes scrna_qc functionality as tools that can be called
by Claude Code or other MCP-compatible clients.

Tools provided:
- profile_data: Profile scRNA-seq data (format detection, batch key, MT genes)
- smart_load: Load data with automatic format detection
- run_qc: Run quality control filtering
- run_doublet_detection: Detect and filter doublets
- normalize_data: Normalize and select HVGs
- run_dimensionality_reduction: PCA + UMAP
- run_clustering: Leiden/Louvain clustering
- run_integration: Batch integration
- run_full_pipeline: Run complete agent pipeline
- get_analysis_status: Check current analysis state
- suggest_next_steps: Get AI-recommended next steps

Usage:
    # Run as standalone server
    python -m scrna_qc.mcp_server

    # Or via entry point
    scrna-qc-mcp
"""

import asyncio
import json
import logging
import os
import sys
import traceback
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

# MCP imports
try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import (
        Tool,
        TextContent,
        ImageContent,
        EmbeddedResource,
    )
except ImportError:
    print("MCP package not installed. Install with: pip install mcp", file=sys.stderr)
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger("scrna_qc.mcp")

# Create MCP server
server = Server("scrna-qc")

# Global state for the current analysis session
class AnalysisSession:
    """Maintains state across tool calls."""

    def __init__(self):
        self.adata = None
        self.data_profile = None
        self.run_metadata = {}
        self.output_dir = None
        self.input_path = None
        self.config = None

    def reset(self):
        self.adata = None
        self.data_profile = None
        self.run_metadata = {}
        self.output_dir = None
        self.input_path = None
        self.config = None

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "has_data": self.adata is not None,
            "input_path": self.input_path,
            "output_dir": self.output_dir,
        }
        if self.adata is not None:
            result["n_cells"] = self.adata.n_obs
            result["n_genes"] = self.adata.n_vars
            result["obs_columns"] = list(self.adata.obs.columns)
            result["obsm_keys"] = list(self.adata.obsm.keys())
        if self.data_profile is not None:
            result["profile"] = self.data_profile.to_dict()
        return result


session = AnalysisSession()


def json_response(data: Any) -> list:
    """Convert data to JSON text content."""
    if isinstance(data, str):
        return [TextContent(type="text", text=data)]
    return [TextContent(type="text", text=json.dumps(data, indent=2, default=str))]


def error_response(error: str) -> list:
    """Create error response."""
    return [TextContent(type="text", text=json.dumps({"error": error}, indent=2))]


# ==================== Tool Definitions ====================

@server.list_tools()
async def list_tools() -> list[Tool]:
    """List available tools."""
    return [
        Tool(
            name="scrna_profile_data",
            description="""Profile scRNA-seq data to understand its structure and characteristics.

Returns:
- File format detection (h5ad, 10x_mtx, 10x_h5, multi_sample)
- Number of cells and genes
- Detected batch key candidates
- MT gene pattern detection
- Existing QC/processing status
- Recommended actions

Use this FIRST before running any analysis to understand the data.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "input_path": {
                        "type": "string",
                        "description": "Path to input file (h5ad) or directory (10x/multi-sample)"
                    }
                },
                "required": ["input_path"]
            }
        ),
        Tool(
            name="scrna_load_data",
            description="""Load scRNA-seq data with automatic format detection.

Supports:
- Single h5ad file
- 10x Genomics MTX directory
- 10x Genomics H5 file
- Multi-sample directory (auto-merge)

The loaded data is stored in the session for subsequent operations.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "input_path": {
                        "type": "string",
                        "description": "Path to input file or directory"
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Output directory for results"
                    }
                },
                "required": ["input_path", "output_dir"]
            }
        ),
        Tool(
            name="scrna_run_qc",
            description="""Run quality control filtering on loaded data.

Uses adaptive MAD-based thresholds to filter:
- Low-quality cells (low gene count)
- Potential doublets (high gene count)
- Dying cells (high MT%)
- Lowly-expressed genes

Returns QC summary with thresholds used and cell retention rate.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "adaptive_qc": {
                        "type": "boolean",
                        "description": "Use adaptive MAD-based thresholds (default: true)",
                        "default": True
                    },
                    "nmads": {
                        "type": "number",
                        "description": "Number of MADs for outlier detection (default: 3.0)",
                        "default": 3.0
                    },
                    "min_genes": {
                        "type": "integer",
                        "description": "Minimum genes per cell floor (default: 200)",
                        "default": 200
                    },
                    "max_mt_pct": {
                        "type": "number",
                        "description": "Maximum MT% ceiling (default: 20.0)",
                        "default": 20.0
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="scrna_run_doublet_detection",
            description="""Detect doublets using Scrublet.

Doublets are artifacts where two cells are captured together.
Detection is performed per-batch if batch_key is specified.

Returns doublet scores and predictions. Optionally filters doublets.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "batch_key": {
                        "type": "string",
                        "description": "Column for batch-aware detection (optional)"
                    },
                    "filter_doublets": {
                        "type": "boolean",
                        "description": "Remove detected doublets (default: true)",
                        "default": True
                    },
                    "expected_doublet_rate": {
                        "type": "number",
                        "description": "Expected doublet rate (auto-estimated if not specified)"
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="scrna_normalize",
            description="""Normalize data and select highly variable genes.

Steps:
1. Size-factor normalization (counts per 10k)
2. Log1p transformation
3. HVG selection using Seurat v3 method

Returns list of selected HVGs.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "target_sum": {
                        "type": "number",
                        "description": "Normalization target sum (default: 10000)",
                        "default": 10000
                    },
                    "n_top_genes": {
                        "type": "integer",
                        "description": "Number of HVGs to select (default: 2000)",
                        "default": 2000
                    },
                    "hvg_flavor": {
                        "type": "string",
                        "description": "HVG method: seurat_v3, seurat, cell_ranger (default: seurat_v3)",
                        "default": "seurat_v3"
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="scrna_reduce_dimensions",
            description="""Run dimensionality reduction (PCA + UMAP).

Steps:
1. PCA on highly variable genes
2. Elbow-based PC selection
3. kNN graph construction
4. UMAP embedding

Returns number of PCs selected and variance explained.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "n_pcs": {
                        "type": "integer",
                        "description": "Maximum PCs to compute (default: 50)",
                        "default": 50
                    },
                    "n_neighbors": {
                        "type": "integer",
                        "description": "kNN neighbors (default: 15)",
                        "default": 15
                    },
                    "pc_selection_method": {
                        "type": "string",
                        "description": "PC selection: elbow, variance_ratio (default: elbow)",
                        "default": "elbow"
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="scrna_cluster",
            description="""Run graph-based clustering.

Uses Leiden or Louvain algorithm on the kNN graph.
Resolution controls the granularity of clusters.

Returns cluster assignments and statistics.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "resolution": {
                        "type": "number",
                        "description": "Clustering resolution (default: 0.5, higher = more clusters)",
                        "default": 0.5
                    },
                    "method": {
                        "type": "string",
                        "description": "Clustering method: leiden, louvain (default: leiden)",
                        "default": "leiden"
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="scrna_integrate",
            description="""Run batch integration.

Integrates data across batches using scVI, Harmony, or Scanorama.
Benchmarks methods if multiple are specified.

Returns integration results and best method recommendation.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "batch_key": {
                        "type": "string",
                        "description": "Column containing batch information (required)"
                    },
                    "methods": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Integration methods: scvi, harmony, scanorama (default: [scvi])",
                        "default": ["scvi"]
                    }
                },
                "required": ["batch_key"]
            }
        ),
        Tool(
            name="scrna_run_full_pipeline",
            description="""Run the complete agent-friendly pipeline.

Automatically:
1. Profiles data
2. Detects batch key and MT pattern
3. Loads and merges samples
4. Runs doublet detection
5. Applies QC filtering
6. Normalizes and selects HVGs
7. Reduces dimensions
8. Clusters cells
9. Integrates batches (if applicable)
10. Generates structured report

Returns comprehensive analysis report.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "input_path": {
                        "type": "string",
                        "description": "Path to input file or directory"
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Output directory for results"
                    },
                    "batch_key": {
                        "type": "string",
                        "description": "Batch key (auto-detected if not specified)"
                    },
                    "resolution": {
                        "type": "number",
                        "description": "Clustering resolution (default: 0.5)",
                        "default": 0.5
                    },
                    "dry_run": {
                        "type": "boolean",
                        "description": "Only profile data, don't process (default: false)",
                        "default": False
                    }
                },
                "required": ["input_path", "output_dir"]
            }
        ),
        Tool(
            name="scrna_get_status",
            description="""Get current analysis session status.

Returns information about:
- Loaded data (cells, genes)
- Available annotations
- Processing state
- Files created""",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="scrna_suggest_next_steps",
            description="""Get AI-recommended next analysis steps.

Based on current data state, suggests:
- Cell type annotation
- Marker gene analysis
- Trajectory inference
- Integration recommendations""",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="scrna_save_data",
            description="""Save current data to h5ad file.

Saves the processed AnnData object to the specified path.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "output_path": {
                        "type": "string",
                        "description": "Output file path (default: {output_dir}/processed.h5ad)"
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="scrna_find_markers",
            description="""Find marker genes for each cluster.

Uses Wilcoxon rank-sum test to identify differentially expressed genes.

Returns top marker genes per cluster.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "groupby": {
                        "type": "string",
                        "description": "Column to group by (default: leiden or louvain)"
                    },
                    "n_genes": {
                        "type": "integer",
                        "description": "Number of top genes per cluster (default: 10)",
                        "default": 10
                    },
                    "method": {
                        "type": "string",
                        "description": "Test method: wilcoxon, t-test (default: wilcoxon)",
                        "default": "wilcoxon"
                    }
                },
                "required": []
            }
        ),
    ]


# ==================== Tool Implementations ====================

@server.call_tool()
async def call_tool(name: str, arguments: dict) -> Sequence[TextContent | ImageContent | EmbeddedResource]:
    """Handle tool calls."""
    try:
        if name == "scrna_profile_data":
            return await tool_profile_data(arguments)
        elif name == "scrna_load_data":
            return await tool_load_data(arguments)
        elif name == "scrna_run_qc":
            return await tool_run_qc(arguments)
        elif name == "scrna_run_doublet_detection":
            return await tool_run_doublet_detection(arguments)
        elif name == "scrna_normalize":
            return await tool_normalize(arguments)
        elif name == "scrna_reduce_dimensions":
            return await tool_reduce_dimensions(arguments)
        elif name == "scrna_cluster":
            return await tool_cluster(arguments)
        elif name == "scrna_integrate":
            return await tool_integrate(arguments)
        elif name == "scrna_run_full_pipeline":
            return await tool_run_full_pipeline(arguments)
        elif name == "scrna_get_status":
            return await tool_get_status(arguments)
        elif name == "scrna_suggest_next_steps":
            return await tool_suggest_next_steps(arguments)
        elif name == "scrna_save_data":
            return await tool_save_data(arguments)
        elif name == "scrna_find_markers":
            return await tool_find_markers(arguments)
        else:
            return error_response(f"Unknown tool: {name}")
    except Exception as e:
        logger.error(f"Tool {name} failed: {e}")
        logger.error(traceback.format_exc())
        return error_response(f"Tool failed: {str(e)}")


async def tool_profile_data(args: dict) -> list:
    """Profile data."""
    from .profile import profile_data

    input_path = args["input_path"]
    logger.info(f"Profiling: {input_path}")

    profile = profile_data(input_path)
    session.data_profile = profile
    session.input_path = input_path

    return json_response(profile.to_dict())


async def tool_load_data(args: dict) -> list:
    """Load data."""
    from .loader import smart_load

    input_path = args["input_path"]
    output_dir = args["output_dir"]

    logger.info(f"Loading: {input_path}")
    os.makedirs(output_dir, exist_ok=True)

    result = smart_load(input_path)

    if result.load_errors:
        return error_response(f"Load failed: {result.load_errors}")

    session.adata = result.adata
    session.input_path = input_path
    session.output_dir = output_dir
    session.run_metadata = {
        "input_path": input_path,
        "output_dir": output_dir,
        "input_cells": result.adata.n_obs,
        "input_genes": result.adata.n_vars,
    }

    return json_response({
        "status": "success",
        "n_cells": result.adata.n_obs,
        "n_genes": result.adata.n_vars,
        "n_samples": result.n_samples,
        "sample_key": result.sample_key,
        "obs_columns": list(result.adata.obs.columns),
    })


async def tool_run_qc(args: dict) -> list:
    """Run QC."""
    from .pipeline import run_qc
    from .config import PipelineConfig

    if session.adata is None:
        return error_response("No data loaded. Use scrna_load_data first.")

    config = PipelineConfig(
        adaptive_qc=args.get("adaptive_qc", True),
        qc_nmads=args.get("nmads", 3.0),
        min_genes=args.get("min_genes", 200),
        max_mt_pct=args.get("max_mt_pct", 20.0),
    )

    cells_before = session.adata.n_obs
    qc_summary = run_qc(session.adata, config, session.run_metadata)

    # Save QC summary
    if session.output_dir:
        qc_path = os.path.join(session.output_dir, "qc_summary.csv")
        qc_summary.to_csv(qc_path, index=False)

    return json_response({
        "status": "success",
        "cells_before": cells_before,
        "cells_after": session.adata.n_obs,
        "retention_rate": session.adata.n_obs / cells_before,
        "thresholds": session.run_metadata.get("qc_thresholds", {}),
    })


async def tool_run_doublet_detection(args: dict) -> list:
    """Run doublet detection."""
    from .doublet import run_doublet_detection, filter_doublets, DoubletConfig

    if session.adata is None:
        return error_response("No data loaded. Use scrna_load_data first.")

    batch_key = args.get("batch_key")
    filter_flag = args.get("filter_doublets", True)

    config = DoubletConfig(
        expected_doublet_rate=args.get("expected_doublet_rate"),
    )

    try:
        result = run_doublet_detection(
            session.adata,
            batch_key=batch_key,
            config=config,
            inplace=True,
        )

        if filter_flag:
            session.adata = filter_doublets(session.adata, inplace=False)

        session.run_metadata["doublet_detection"] = result.to_dict()

        return json_response({
            "status": "success",
            "n_doublets": result.n_doublets,
            "doublet_rate": result.doublet_rate,
            "cells_after_filtering": session.adata.n_obs if filter_flag else "not filtered",
            "batch_results": result.batch_results,
        })
    except ImportError:
        return error_response("Scrublet not installed. Install with: pip install scrublet")


async def tool_normalize(args: dict) -> list:
    """Normalize data."""
    from .pipeline import normalize_and_scale
    from .config import PipelineConfig

    if session.adata is None:
        return error_response("No data loaded. Use scrna_load_data first.")

    config = PipelineConfig(
        target_sum=args.get("target_sum", 10000),
        n_top_genes=args.get("n_top_genes", 2000),
        hvg_flavor=args.get("hvg_flavor", "seurat_v3"),
    )

    hvg_list = normalize_and_scale(session.adata, config, session.run_metadata)

    # Save HVG list
    if session.output_dir:
        hvg_path = os.path.join(session.output_dir, "hvg_list.txt")
        with open(hvg_path, 'w') as f:
            f.write('\n'.join(hvg_list))

    return json_response({
        "status": "success",
        "n_hvgs": len(hvg_list),
        "top_hvgs": hvg_list[:20],
        "normalization_method": "size_factor_log1p",
        "target_sum": config.target_sum,
    })


async def tool_reduce_dimensions(args: dict) -> list:
    """Reduce dimensions."""
    from .pipeline import run_dimensionality_reduction
    from .config import PipelineConfig

    if session.adata is None:
        return error_response("No data loaded. Use scrna_load_data first.")

    config = PipelineConfig(
        n_pcs=args.get("n_pcs", 50),
        n_neighbors=args.get("n_neighbors", 15),
        pc_selection_method=args.get("pc_selection_method", "elbow"),
    )

    n_pcs_used = run_dimensionality_reduction(session.adata, config, session.run_metadata)

    return json_response({
        "status": "success",
        "n_pcs_computed": config.n_pcs,
        "n_pcs_selected": n_pcs_used,
        "pc_selection_method": config.pc_selection_method,
        "n_neighbors": config.n_neighbors,
        "umap_computed": "X_umap" in session.adata.obsm,
    })


async def tool_cluster(args: dict) -> list:
    """Run clustering."""
    from .pipeline import run_clustering
    from .config import PipelineConfig

    if session.adata is None:
        return error_response("No data loaded. Use scrna_load_data first.")

    if "neighbors" not in session.adata.uns:
        return error_response("Neighbors not computed. Use scrna_reduce_dimensions first.")

    config = PipelineConfig(
        resolution=args.get("resolution", 0.5),
        clustering_method=args.get("method", "leiden"),
    )

    clusters = run_clustering(session.adata, config, session.run_metadata)

    # Save clusters
    if session.output_dir:
        clusters_path = os.path.join(session.output_dir, "clusters.csv")
        clusters.to_csv(clusters_path, index=False)

    cluster_key = config.clustering_method
    cluster_sizes = session.adata.obs[cluster_key].value_counts().to_dict()

    return json_response({
        "status": "success",
        "n_clusters": len(cluster_sizes),
        "cluster_sizes": cluster_sizes,
        "resolution": config.resolution,
        "method": config.clustering_method,
    })


async def tool_integrate(args: dict) -> list:
    """Run integration."""
    from .integration import run_integration
    from .config import IntegrationConfig

    if session.adata is None:
        return error_response("No data loaded. Use scrna_load_data first.")

    batch_key = args["batch_key"]
    methods = args.get("methods", ["scvi"])

    if batch_key not in session.adata.obs.columns:
        return error_response(f"Batch key '{batch_key}' not found in data.")

    config = IntegrationConfig(methods=methods)

    try:
        result = run_integration(
            session.adata,
            batch_key=batch_key,
            methods=methods,
            config=config,
            output_dir=session.output_dir,
        )

        session.run_metadata["integration"] = {
            "methods": result.methods_run,
            "best_method": result.best_method,
            "benchmark_scores": result.benchmark_scores,
        }

        return json_response({
            "status": "success",
            "methods_run": result.methods_run,
            "best_method": result.best_method,
            "benchmark_scores": result.benchmark_scores,
        })
    except Exception as e:
        return error_response(f"Integration failed: {str(e)}")


async def tool_run_full_pipeline(args: dict) -> list:
    """Run full pipeline."""
    from .pipeline import run_agent_pipeline
    from .config import PipelineConfig

    input_path = args["input_path"]
    output_dir = args["output_dir"]
    dry_run = args.get("dry_run", False)

    config = PipelineConfig(
        batch_key=args.get("batch_key"),
        resolution=args.get("resolution", 0.5),
    )

    results = run_agent_pipeline(
        input_path=input_path,
        output_dir=output_dir,
        config=config,
        dry_run=dry_run,
    )

    session.adata = results.adata
    session.data_profile = results.data_profile
    session.run_metadata = results.run_metadata
    session.output_dir = output_dir
    session.input_path = input_path

    if results.agent_report:
        return json_response(results.agent_report.to_dict())
    elif results.data_profile:
        return json_response({
            "status": "dry_run_complete",
            "profile": results.data_profile.to_dict(),
        })
    else:
        return json_response({
            "status": "complete",
            "n_cells": results.adata.n_obs if results.adata else 0,
            "files_created": results.files_created,
        })


async def tool_get_status(args: dict) -> list:
    """Get session status."""
    return json_response(session.to_dict())


async def tool_suggest_next_steps(args: dict) -> list:
    """Suggest next steps."""
    from .report import suggest_next_steps, AgentReport

    if session.adata is None:
        return json_response({
            "next_steps": [
                {"action": "Load Data", "description": "Use scrna_load_data to load your data"},
                {"action": "Profile Data", "description": "Use scrna_profile_data to understand your data first"},
            ]
        })

    # Create a minimal report for suggestions
    report = AgentReport(
        final_cells=session.adata.n_obs,
        final_genes=session.adata.n_vars,
    )

    steps = suggest_next_steps(report, session.adata)

    return json_response({
        "current_state": {
            "n_cells": session.adata.n_obs,
            "has_pca": "X_pca" in session.adata.obsm,
            "has_umap": "X_umap" in session.adata.obsm,
            "has_clusters": "leiden" in session.adata.obs.columns or "louvain" in session.adata.obs.columns,
        },
        "next_steps": steps,
    })


async def tool_save_data(args: dict) -> list:
    """Save data."""
    if session.adata is None:
        return error_response("No data loaded.")

    output_path = args.get("output_path")
    if output_path is None:
        if session.output_dir:
            output_path = os.path.join(session.output_dir, "processed.h5ad")
        else:
            return error_response("No output path specified and no output_dir set.")

    session.adata.write_h5ad(output_path)

    return json_response({
        "status": "success",
        "output_path": output_path,
        "n_cells": session.adata.n_obs,
        "n_genes": session.adata.n_vars,
    })


async def tool_find_markers(args: dict) -> list:
    """Find marker genes."""
    import scanpy as sc

    if session.adata is None:
        return error_response("No data loaded.")

    # Determine groupby column
    groupby = args.get("groupby")
    if groupby is None:
        if "leiden" in session.adata.obs.columns:
            groupby = "leiden"
        elif "louvain" in session.adata.obs.columns:
            groupby = "louvain"
        else:
            return error_response("No clustering found. Run scrna_cluster first or specify groupby.")

    n_genes = args.get("n_genes", 10)
    method = args.get("method", "wilcoxon")

    # Run marker gene analysis
    sc.tl.rank_genes_groups(session.adata, groupby, method=method)

    # Extract top markers per cluster
    markers = {}
    for cluster in session.adata.obs[groupby].unique():
        cluster_markers = sc.get.rank_genes_groups_df(
            session.adata, group=str(cluster)
        ).head(n_genes)
        markers[str(cluster)] = cluster_markers[["names", "scores", "pvals_adj", "logfoldchanges"]].to_dict(orient="records")

    # Save markers
    if session.output_dir:
        markers_path = os.path.join(session.output_dir, "marker_genes.json")
        with open(markers_path, 'w') as f:
            json.dump(markers, f, indent=2)

    return json_response({
        "status": "success",
        "groupby": groupby,
        "n_clusters": len(markers),
        "markers": markers,
    })


# ==================== Main ====================

async def main():
    """Run the MCP server."""
    logger.info("Starting scrna-qc MCP server...")

    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options(),
        )


def run_server():
    """Entry point for the MCP server."""
    asyncio.run(main())


if __name__ == "__main__":
    run_server()
