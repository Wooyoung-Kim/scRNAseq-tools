import os
import shlex
from pathlib import Path
from typing import Optional

from ..data_utils import env_info, render_env_text, render_summary_text, scan_data, summarize_file


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
            try:
                matches = scan_data(target, recursive)
            except FileNotFoundError:
                print(f"not found: {target}")
                continue
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
            try:
                data = summarize_file(Path(args[0]).expanduser())
            except FileNotFoundError:
                print(f"not found: {Path(args[0]).expanduser()}")
                continue
            except RuntimeError as exc:
                print(str(exc))
                continue
            print(render_summary_text(data))
            continue
        if cmd == "env":
            print(render_env_text(env_info()))
            continue

        print(f"unknown command: {cmd}")
        print("type 'help' to see available commands")
