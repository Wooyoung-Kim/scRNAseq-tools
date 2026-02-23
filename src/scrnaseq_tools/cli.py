import argparse
import sys

from .commands import register_subcommands
from .manifest import write_run_manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scRNAseq-tools",
        description="CLI toolkit for scRNA-seq workflows.",
    )
    subparsers = parser.add_subparsers(dest="command")
    register_subcommands(subparsers)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    exit_code = 0
    error = None
    try:
        exit_code = args.func(args)
        if exit_code is None:
            exit_code = 0
        return exit_code
    except Exception as exc:
        error = f"{exc.__class__.__name__}: {exc}"
        exit_code = 1
        raise
    finally:
        write_run_manifest(
            getattr(args, "command", "unknown"),
            list(argv) if argv is not None else sys.argv[1:],
            exit_code,
            error,
        )
