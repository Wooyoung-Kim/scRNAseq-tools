import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

from . import __version__


def write_run_manifest(command: str, argv: List[str], exit_code: int, error: Optional[str]) -> None:
    if os.environ.get("SCRNASEQ_MANIFEST", "1").lower() in ("0", "false", "no"):
        return

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "argv": argv,
        "exit_code": exit_code,
        "error": error,
        "cwd": os.getcwd(),
        "user": os.environ.get("USER") or os.environ.get("USERNAME") or "unknown",
        "host": platform.node(),
        "python": platform.python_version(),
        "scrnaseq_tools_version": __version__,
    }

    try:
        path = Path(os.environ.get("SCRNASEQ_MANIFEST_PATH", "run_manifest.json"))
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    except OSError:
        # Do not fail CLI execution due to manifest write errors.
        return
