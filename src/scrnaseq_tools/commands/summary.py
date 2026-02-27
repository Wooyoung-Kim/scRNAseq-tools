import json
from pathlib import Path

from ..data_utils import render_summary_text, summarize_file


def register(subparsers) -> None:
    parser = subparsers.add_parser("summary", help="Summarize a data file")
    parser.add_argument("path", help="Path to file (.h5ad/.csv/.tsv/.txt or other)")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output summary as JSON",
    )
    parser.set_defaults(func=run)


def run(args) -> int:
    target = Path(args.path).expanduser()
    try:
        data = summarize_file(target)
    except FileNotFoundError:
        print(f"not found: {target}")
        return 1
    except RuntimeError as exc:
        print(str(exc))
        return 1

    if args.json:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        print(render_summary_text(data))
    return 0
