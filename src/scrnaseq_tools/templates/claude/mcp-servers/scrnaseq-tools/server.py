#!/usr/bin/env python3
"""
scRNAseq-tools Unified MCP Server

Basic tools (always available):
- list_dir: list directory contents
- scan_data: find common scRNA-seq data files
- summarize: summarize h5ad or table files
- env: show installed analysis packages

Analysis tools (available with pip install 'scRNAseq-tools[analysis]'):
- scrna_profile_data, scrna_load_data, scrna_run_qc, scrna_run_doublet_detection
- scrna_normalize, scrna_reduce_dimensions, scrna_cluster, scrna_integrate
- scrna_run_full_pipeline, scrna_get_status, scrna_suggest_next_steps
- scrna_save_data, scrna_find_markers
"""

import asyncio
import json
from pathlib import Path
from typing import Iterable, Optional

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool

DATA_EXTENSIONS = (
    ".h5ad", ".loom", ".mtx", ".mtx.gz",
    ".h5", ".hdf5", ".csv", ".tsv", ".txt",
)

server = Server("scrnaseq-tools")

# ==================== Check for analysis subpackage ====================

_HAS_ANALYSIS = False
try:
    from scrnaseq_tools.analysis.mcp_server import (
        list_tools as analysis_list_tools,
        call_tool as analysis_call_tool,
    )
    _HAS_ANALYSIS = True
except ImportError:
    pass


# ==================== Basic Tool Definitions ====================

BASIC_TOOLS = [
    Tool(
        name="list_dir",
        description="List directory contents.",
        inputSchema={
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Directory path (default: .)"}
            },
        },
    ),
    Tool(
        name="scan_data",
        description="Find data files in a directory.",
        inputSchema={
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to scan (default: .)"},
                "recursive": {"type": "boolean", "description": "Recursive search", "default": False},
            },
        },
    ),
    Tool(
        name="summarize",
        description="Summarize a file (h5ad/csv/tsv/txt).",
        inputSchema={
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "File path"},
            },
            "required": ["path"],
        },
    ),
    Tool(
        name="env",
        description="Show installed analysis packages.",
        inputSchema={"type": "object", "properties": {}},
    ),
]


# ==================== Tool Listing ====================

@server.list_tools()
async def list_tools() -> list[Tool]:
    tools = list(BASIC_TOOLS)
    if _HAS_ANALYSIS:
        # Get analysis tools from merged server
        analysis_tools = await analysis_list_tools()
        tools.extend(analysis_tools)
    return tools


# ==================== Basic Tool Implementations ====================

BASIC_TOOL_NAMES = {"list_dir", "scan_data", "summarize", "env"}


@server.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    if name in BASIC_TOOL_NAMES:
        # Handle basic tools
        if name == "list_dir":
            result = list_dir_impl(arguments.get("path"))
        elif name == "scan_data":
            result = scan_data_impl(arguments.get("path"), arguments.get("recursive", False))
        elif name == "summarize":
            result = summarize_impl(Path(arguments["path"]).expanduser())
        elif name == "env":
            result = env_info_impl()
        else:
            result = {"error": f"Unknown tool: {name}"}
        return [TextContent(type="text", text=json.dumps(result, indent=2))]

    # Delegate to analysis tools
    if _HAS_ANALYSIS:
        return await analysis_call_tool(name, arguments)

    return [TextContent(
        type="text",
        text=json.dumps({"error": f"Unknown tool: {name}. Analysis tools require: pip install 'scRNAseq-tools[analysis]'"}, indent=2)
    )]


def list_dir_impl(path: Optional[str]) -> dict:
    target = Path(path or ".").expanduser()
    if not target.exists():
        return {"error": f"not found: {target}"}
    if target.is_file():
        return {"path": str(target), "items": [target.name]}
    items = []
    for entry in sorted(target.iterdir(), key=lambda p: (p.is_file(), p.name.lower())):
        items.append({"name": entry.name, "type": "dir" if entry.is_dir() else "file"})
    return {"path": str(target), "items": items}


def scan_data_impl(path: Optional[str], recursive: bool) -> dict:
    target = Path(path or ".").expanduser()
    if not target.exists():
        return {"error": f"not found: {target}"}
    if target.is_file():
        return {"path": str(target), "matches": [str(target)]}
    iterator: Iterable[Path]
    iterator = target.rglob("*") if recursive else target.iterdir()
    matches = [
        str(p)
        for p in iterator
        if p.is_file() and p.name.lower().endswith(DATA_EXTENSIONS)
    ]
    return {"path": str(target), "recursive": recursive, "matches": sorted(matches)}


def summarize_impl(path: Path) -> dict:
    if not path.exists():
        return {"error": f"not found: {path}"}

    # Try rich profiling from analysis if available
    if _HAS_ANALYSIS and path.suffix.lower() == ".h5ad":
        try:
            from scrnaseq_tools.analysis.profile import profile_data
            profile = profile_data(str(path))
            return profile.to_dict()
        except Exception:
            pass

    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        return summarize_h5ad_impl(path)
    if suffix in (".csv", ".tsv", ".txt"):
        return summarize_table_impl(path)
    return {"path": str(path), "size_bytes": path.stat().st_size}


def summarize_h5ad_impl(path: Path) -> dict:
    try:
        import anndata as ad
    except ImportError:
        return {"error": "anndata is not installed; cannot read .h5ad"}
    adata = ad.read_h5ad(path, backed="r")
    try:
        return {
            "path": str(path),
            "shape": [adata.n_obs, adata.n_vars],
            "obs_columns": list(adata.obs.keys())[:20],
            "var_columns": list(adata.var.keys())[:20],
            "layers": list(adata.layers.keys())[:20],
            "obsm": list(adata.obsm.keys())[:20],
            "uns": list(adata.uns.keys())[:20],
        }
    finally:
        backing = getattr(adata, "file", None)
        close = getattr(backing, "close", None)
        if callable(close):
            close()


def summarize_table_impl(path: Path) -> dict:
    try:
        import pandas as pd
    except ImportError:
        return {"error": "pandas is not installed; cannot read table"}
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(path, sep=sep, nrows=5)
    return {
        "path": str(path),
        "columns": list(df.columns),
        "head": df.head().to_dict(orient="records"),
    }


def env_info_impl() -> dict:
    def get_version(module_name: str) -> str:
        try:
            module = __import__(module_name)
        except ImportError:
            return "not installed"
        return getattr(module, "__version__", "unknown")

    info = {
        "anndata": get_version("anndata"),
        "scanpy": get_version("scanpy"),
        "pandas": get_version("pandas"),
        "analysis_available": _HAS_ANALYSIS,
    }
    if _HAS_ANALYSIS:
        from scrnaseq_tools.analysis import __version__ as analysis_version
        info["analysis_version"] = analysis_version
        for pkg in ("scvi", "harmonypy", "scanorama", "scrublet"):
            info[pkg] = get_version(pkg)
    return info


async def main() -> None:
    async with stdio_server() as (read, write):
        await server.run(read, write)


if __name__ == "__main__":
    asyncio.run(main())
