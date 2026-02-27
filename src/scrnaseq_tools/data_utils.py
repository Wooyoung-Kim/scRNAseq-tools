import importlib
from pathlib import Path
from typing import Any, Dict, Iterable, List


DATA_EXTENSIONS = (
    ".h5ad",
    ".loom",
    ".mtx",
    ".mtx.gz",
    ".h5",
    ".hdf5",
    ".csv",
    ".tsv",
    ".txt",
)

TABLE_EXTENSIONS = (".csv", ".tsv", ".txt")


def scan_data(path: Path, recursive: bool) -> List[Path]:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.is_file():
        return [path]
    iterator: Iterable[Path] = path.rglob("*") if recursive else path.iterdir()
    matches = [p for p in iterator if p.is_file() and p.name.lower().endswith(DATA_EXTENSIONS)]
    return sorted(matches)


def summarize_file(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".h5ad":
        return summarize_h5ad(path)
    if path.suffix.lower() in TABLE_EXTENSIONS:
        return summarize_table(path)
    return {
        "path": str(path),
        "type": "file",
        "size_bytes": path.stat().st_size,
    }


def summarize_h5ad(path: Path) -> Dict[str, Any]:
    # Try rich profiling from analysis subpackage first
    try:
        from scrnaseq_tools.analysis.profile import profile_data
        profile = profile_data(str(path))
        return profile.to_dict()
    except ImportError:
        pass

    # Fallback: lightweight summary (anndata only)
    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError("anndata is not installed; cannot read .h5ad") from exc

    adata = ad.read_h5ad(path, backed="r")
    try:
        return {
            "path": str(path),
            "type": "h5ad",
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
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


def summarize_table(path: Path) -> Dict[str, Any]:
    try:
        import pandas as pd
    except ImportError as exc:
        raise RuntimeError("pandas is not installed; cannot read table") from exc

    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(path, sep=sep, nrows=5)
    preview_rows = []
    for _, row in df.iterrows():
        preview_rows.append(
            {
                column: (None if pd.isna(value) else str(value))
                for column, value in row.items()
            }
        )
    return {
        "path": str(path),
        "type": "table",
        "columns": [str(col) for col in df.columns],
        "preview_rows": preview_rows,
        "preview_text": df.head().to_string(index=False),
    }


def env_info() -> Dict[str, str]:
    info = {
        "anndata": _module_version("anndata"),
        "scanpy": _module_version("scanpy"),
        "pandas": _module_version("pandas"),
    }
    # Show additional packages when analysis is available
    try:
        from scrnaseq_tools.analysis import __version__ as analysis_version
        info["scrnaseq_tools.analysis"] = analysis_version
        for pkg in ("scvi", "harmonypy", "scanorama", "scrublet"):
            info[pkg] = _module_version(pkg)
    except ImportError:
        pass
    return info


def render_env_text(data: Dict[str, str]) -> str:
    return "\n".join(
        [
            f"anndata: {data['anndata']}",
            f"scanpy: {data['scanpy']}",
            f"pandas: {data['pandas']}",
        ]
    )


def render_summary_text(data: Dict[str, Any]) -> str:
    kind = data.get("type")
    if kind == "h5ad":
        lines = [
            f"path: {data['path']}",
            f"shape: {data['n_obs']} cells x {data['n_vars']} genes",
            f"obs columns: {', '.join(data['obs_columns'])}",
            f"var columns: {', '.join(data['var_columns'])}",
        ]
        if data["layers"]:
            lines.append(f"layers: {', '.join(data['layers'])}")
        if data["obsm"]:
            lines.append(f"obsm: {', '.join(data['obsm'])}")
        if data["uns"]:
            lines.append(f"uns: {', '.join(data['uns'])}")
        return "\n".join(lines)

    if kind == "table":
        return "\n".join(
            [
                f"path: {data['path']}",
                f"columns: {', '.join(data['columns'])}",
                data["preview_text"],
            ]
        )

    return "\n".join(
        [
            f"file: {data['path']}",
            f"size: {data['size_bytes']} bytes",
        ]
    )


def _module_version(module_name: str) -> str:
    try:
        module = importlib.import_module(module_name)
    except ImportError:
        return "not installed"
    return getattr(module, "__version__", "unknown")
