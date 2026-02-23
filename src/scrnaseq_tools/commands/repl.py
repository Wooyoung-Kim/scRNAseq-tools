import os
import shlex
from pathlib import Path
from typing import Iterable, List, Optional


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


def register(subparsers) -> None:
    parser = subparsers.add_parser("repl", help="Local analysis REPL (no external LLM)")
    parser.set_defaults(func=run)


def list_dir(path: Path) -> None:
    if not path.exists():
        print(f"not found: {path}")
        return
    if path.is_file():
        print(path.name)
        return
    for entry in sorted(path.iterdir(), key=lambda p: (p.is_file(), p.name.lower())):
        suffix = "/" if entry.is_dir() else ""
        print(f"{entry.name}{suffix}")


def scan_data(path: Path, recursive: bool) -> List[Path]:
    if not path.exists():
        print(f"not found: {path}")
        return []
    if path.is_file():
        return [path]
    iterator: Iterable[Path]
    iterator = path.rglob("*") if recursive else path.iterdir()
    matches = [p for p in iterator if p.is_file() and p.name.lower().endswith(DATA_EXTENSIONS)]
    return sorted(matches)


def summarize_file(path: Path) -> None:
    if not path.exists():
        print(f"not found: {path}")
        return
    if path.suffix.lower() == ".h5ad":
        summarize_h5ad(path)
        return
    if path.suffix.lower() in (".csv", ".tsv", ".txt"):
        summarize_table(path)
        return
    print(f"file: {path}")
    print(f"size: {path.stat().st_size} bytes")


def summarize_h5ad(path: Path) -> None:
    try:
        import anndata as ad
    except ImportError:
        print("anndata is not installed; cannot read .h5ad")
        return
    adata = ad.read_h5ad(path, backed="r")
    print(f"path: {path}")
    print(f"shape: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"obs columns: {', '.join(list(adata.obs.keys())[:20])}")
    print(f"var columns: {', '.join(list(adata.var.keys())[:20])}")
    if adata.layers:
        print(f"layers: {', '.join(list(adata.layers.keys())[:20])}")
    if adata.obsm:
        print(f"obsm: {', '.join(list(adata.obsm.keys())[:20])}")
    if adata.uns:
        print(f"uns: {', '.join(list(adata.uns.keys())[:20])}")


def summarize_table(path: Path) -> None:
    try:
        import pandas as pd
    except ImportError:
        print("pandas is not installed; cannot read table")
        return
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(path, sep=sep, nrows=5)
    print(f"path: {path}")
    print(f"columns: {', '.join(df.columns)}")
    print(df.head().to_string(index=False))


def show_env() -> None:
    try:
        import anndata
        ad_ver = anndata.__version__
    except ImportError:
        ad_ver = "not installed"
    try:
        import scanpy
        sc_ver = scanpy.__version__
    except ImportError:
        sc_ver = "not installed"
    try:
        import pandas
        pd_ver = pandas.__version__
    except ImportError:
        pd_ver = "not installed"
    print(f"anndata: {ad_ver}")
    print(f"scanpy: {sc_ver}")
    print(f"pandas: {pd_ver}")


def print_help() -> None:
    print("commands:")
    print("  help                         show this help")
    print("  exit | quit                  exit REPL")
    print("  pwd                          print working directory")
    print("  cd <path>                    change directory")
    print("  ls [path]                    list directory contents")
    print("  scan [path] [-r]             list data files")
    print("  summary <path>               summarize a file (h5ad/csv/tsv)")
    print("  env                          show installed analysis packages")


def run(_args) -> int:
    print("scRNAseq-tools REPL (local analysis only). Type 'help' for commands.")
    while True:
        try:
            line = input("scrna> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        if not line:
            continue
        try:
            tokens = shlex.split(line)
        except ValueError as exc:
            print(f"parse error: {exc}")
            continue
        cmd = tokens[0].lower()
        args = tokens[1:]

        if cmd in ("exit", "quit"):
            return 0
        if cmd in ("help", "?"):
            print_help()
            continue
        if cmd == "pwd":
            print(os.getcwd())
            continue
        if cmd == "cd":
            if not args:
                print("usage: cd <path>")
                continue
            new_dir = Path(args[0]).expanduser()
            if not new_dir.is_dir():
                print(f"not a directory: {new_dir}")
                continue
            os.chdir(new_dir)
            continue
        if cmd == "ls":
            target = Path(args[0]).expanduser() if args else Path(os.getcwd())
            list_dir(target)
            continue
        if cmd == "scan":
            recursive = False
            target: Optional[Path] = None
            for item in args:
                if item in ("-r", "--recursive"):
                    recursive = True
                elif target is None:
                    target = Path(item).expanduser()
            if target is None:
                target = Path(os.getcwd())
            matches = scan_data(target, recursive)
            if not matches:
                print("no matching files found")
                continue
            for match in matches:
                print(match)
            continue
        if cmd == "summary":
            if not args:
                print("usage: summary <path>")
                continue
            summarize_file(Path(args[0]).expanduser())
            continue
        if cmd == "env":
            show_env()
            continue

        print(f"unknown command: {cmd}")
        print("type 'help' to see available commands")
