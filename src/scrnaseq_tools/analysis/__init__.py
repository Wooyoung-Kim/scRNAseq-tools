"""
scRNAseq-tools analysis subpackage.

Provides scRNA-seq QC-to-Integration pipeline functionality.
Requires optional dependencies: pip install 'scRNAseq-tools[analysis]'

Public API:
    - run_pipeline: Standard QC→Clustering pipeline
    - run_agent_pipeline: Agent-friendly pipeline with auto-detection
    - load_data, run_qc, normalize_and_scale, etc.: Individual steps
    - PipelineConfig, IntegrationConfig: Configuration classes
    - profile_data: Data profiling for AI agents
    - smart_load: Smart data loading
"""

__version__ = "1.1.0"

from .config import (
    PipelineConfig,
    IntegrationConfig,
    DoubletConfig,
    AgentConfig,
    MultiResolutionConfig,
)

from .pipeline import (
    run_pipeline,
    run_agent_pipeline,
    load_data,
    run_qc,
    normalize_and_scale,
    run_dimensionality_reduction,
    run_clustering,
    run_multi_resolution_clustering,
    PipelineResults,
)

__all__ = [
    # Pipeline
    "run_pipeline",
    "run_agent_pipeline",
    "load_data",
    "run_qc",
    "normalize_and_scale",
    "run_dimensionality_reduction",
    "run_clustering",
    "run_multi_resolution_clustering",
    "PipelineResults",
    # Config
    "PipelineConfig",
    "IntegrationConfig",
    "DoubletConfig",
    "AgentConfig",
    "MultiResolutionConfig",
    # Version
    "__version__",
]


def run_mcp_server():
    """Run the MCP server (requires mcp package)."""
    try:
        from .mcp_server import run_server
        run_server()
    except ImportError:
        raise ImportError(
            "MCP package not installed. Install with: pip install 'scRNAseq-tools[mcp]'"
        )
