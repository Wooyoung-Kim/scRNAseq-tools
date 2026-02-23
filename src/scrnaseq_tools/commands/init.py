import os
import shutil
from importlib import resources
from pathlib import Path
from typing import Optional

DEFAULT_CLAUDE_SOURCE = "/home/kwy7605/data_61/Vaccine_V2/.claude"


def register(subparsers) -> None:
    parser = subparsers.add_parser("init", help="Install .claude into a project")
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
    parser.set_defaults(func=run)


def copy_resource_tree(source, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for entry in source.iterdir():
        target = destination / entry.name
        if entry.is_dir():
            copy_resource_tree(entry, target)
        else:
            with entry.open("rb") as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)


def resolve_source_dir(source_arg: Optional[str]) -> Optional[Path]:
    source_env = os.environ.get("SCRNASEQ_TOOLS_CLAUDE_SOURCE")
    default_source_dir = Path(DEFAULT_CLAUDE_SOURCE)

    if source_arg:
        return Path(source_arg).expanduser()
    if source_env:
        return Path(source_env).expanduser()
    if default_source_dir.is_dir():
        return default_source_dir
    return None


def install_claude(source_dir: Optional[Path], target_dir: Path, force: bool) -> Path:
    target_dir.mkdir(parents=True, exist_ok=True)
    dest_dir = target_dir / ".claude"

    if dest_dir.exists():
        if not force:
            raise FileExistsError(dest_dir)
        shutil.rmtree(dest_dir)

    if source_dir is not None:
        if not source_dir.is_dir():
            raise FileNotFoundError(source_dir)
        shutil.copytree(source_dir, dest_dir)
    else:
        source = resources.files("scrnaseq_tools").joinpath("templates/claude")
        copy_resource_tree(source, dest_dir)
    return dest_dir


PROJECT_ROOT_FILES = ["annotation_gate.py", "phase0_config.py"]


def install_project_root_files(dest_claude_dir: Path, target_dir: Path) -> list:
    """Move annotation_gate.py and phase0_config.py from .claude/ to project root."""
    installed = []
    for fname in PROJECT_ROOT_FILES:
        src = dest_claude_dir / fname
        dst = target_dir / fname
        if src.is_file():
            shutil.copy2(src, dst)
            src.unlink()  # remove from .claude/
            installed.append(str(dst))
        else:
            # Try embedded template
            try:
                template = resources.files("scrnaseq_tools").joinpath(f"templates/claude/{fname}")
                if template.is_file():
                    with template.open("rb") as s, dst.open("wb") as d:
                        shutil.copyfileobj(s, d)
                    installed.append(str(dst))
            except Exception:
                pass
    return installed


def install_claude_md(target_dir: Path) -> Optional[Path]:
    """Generate CLAUDE.md at project root if not present."""
    claude_md = target_dir / "CLAUDE.md"
    if claude_md.exists():
        return None  # Don't overwrite user's CLAUDE.md

    content = """\
# Project Rules

## ABSOLUTE RULE: Phase 0 Configuration (CODE-LEVEL ENFORCEMENT)

Two Python files enforce Phase 0 at the code level:

| File | Purpose |
|------|---------|
| `phase0_config.py` | Interactive script — generates `annotation_config.json` |
| `annotation_gate.py` | Validation module — `sys.exit(1)` if config missing |

### MANDATORY: Every annotation script MUST start with

```python
from annotation_gate import load_config
config = load_config()   # sys.exit(1) if config missing
```

### Rules

- **NEVER** hardcode tools, methods, databases, or marker lists
- **NEVER** skip Phase 0
- **NEVER** write a tier script without `from annotation_gate import load_config`
- **ALWAYS** run `python phase0_config.py` before annotation
"""
    claude_md.write_text(content)
    return claude_md


def run(args) -> int:
    source_dir = resolve_source_dir(args.source)
    target_dir = Path(args.target).expanduser()

    try:
        dest_dir = install_claude(source_dir, target_dir, args.force)
    except FileExistsError as exc:
        print(f"Target already has .claude: {exc}")
        print("Use --force to overwrite.")
        return 1
    except FileNotFoundError as exc:
        print(f"Source .claude not found: {exc}")
        return 1

    print(f"Installed .claude to: {dest_dir}")

    # Install project-root files (annotation_gate.py, phase0_config.py)
    root_files = install_project_root_files(dest_dir, target_dir)
    for f in root_files:
        print(f"Installed: {f}")

    # Generate CLAUDE.md if not present
    claude_md = install_claude_md(target_dir)
    if claude_md:
        print(f"Generated: {claude_md}")

    return 0
