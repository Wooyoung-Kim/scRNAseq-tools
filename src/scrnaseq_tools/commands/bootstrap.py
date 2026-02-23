import os
import subprocess
import sys
from pathlib import Path

from .init import install_claude, resolve_source_dir


def register(subparsers) -> None:
    parser = subparsers.add_parser(
        "bootstrap",
        help="Initialize .claude and install MCP server dependencies",
    )
    parser.add_argument(
        "--source",
        default=None,
        help="Path to source .claude directory",
    )
    parser.add_argument(
        "--target",
        default=os.getcwd(),
        help="Target project directory (defaults to current directory)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing .claude directory",
    )
    parser.add_argument(
        "--no-mcp-install",
        action="store_true",
        help="Skip MCP server dependency installation",
    )
    parser.set_defaults(func=run)


def install_with_pip(server_dir: Path) -> bool:
    cmd = [sys.executable, "-m", "pip", "install", "-e", str(server_dir)]
    return subprocess.call(cmd) == 0


def install_requirements(req_file: Path) -> bool:
    cmd = [sys.executable, "-m", "pip", "install", "-r", str(req_file)]
    return subprocess.call(cmd) == 0


def run(args) -> int:
    target_dir = Path(args.target).expanduser()
    source_dir = resolve_source_dir(args.source)

    dest_dir = target_dir / ".claude"
    if not dest_dir.exists() or args.force:
        try:
            dest_dir = install_claude(source_dir, target_dir, args.force)
        except FileExistsError as exc:
            print(f"Target already has .claude: {exc}")
            print("Use --force to overwrite.")
            return 1
        except FileNotFoundError as exc:
            print(f"Source .claude not found: {exc}")
            return 1

    if args.no_mcp_install:
        print("MCP install: skipped")
        return 0

    mcp_root = dest_dir / "mcp-servers"
    if not mcp_root.is_dir():
        print("MCP install: no mcp-servers directory found")
        return 0

    failures = 0
    for server_dir in sorted(mcp_root.iterdir()):
        if not server_dir.is_dir():
            continue
        pyproject = server_dir / "pyproject.toml"
        requirements = server_dir / "requirements.txt"
        print(f"MCP install: {server_dir.name}")
        if pyproject.is_file():
            ok = install_with_pip(server_dir)
        elif requirements.is_file():
            ok = install_requirements(requirements)
        else:
            print("  skip: no pyproject.toml or requirements.txt")
            continue
        if not ok:
            print("  failed: check pip output and network access")
            failures += 1

    if failures:
        print(f"MCP install: {failures} failure(s)")
        return 1

    print("MCP install: complete")
    return 0

