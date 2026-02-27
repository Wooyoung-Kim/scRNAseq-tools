import json
from pathlib import Path

from ..data_utils import scan_data


def register(subparsers) -> None:
    parser = subparsers.add_parser("scan", help="Scan directory for data files")
    parser.add_argument(
        "path",
        nargs="?",
        default=".",
        help="Directory or file to scan (default: current directory)",
    )
    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help="Scan subdirectories recursively",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON",
    )
    parser.set_defaults(func=run)


def run(args) -> int:
    target = Path(args.path).expanduser()
    try:
        matches = scan_data(target, args.recursive)
    except FileNotFoundError:
        print(f"not found: {target}")
        return 1

    if args.json:
        payload = {
            "path": str(target),
            "recursive": bool(args.recursive),
            "count": len(matches),
            "files": [str(path) for path in matches],
        }
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0

    if not matches:
        print("no matching files found")
        return 0

    for path in matches:
        print(path)
    return 0
