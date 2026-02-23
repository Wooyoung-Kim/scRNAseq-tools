"""
Annotation Gate — HARD validation for Phase 0 configuration.

EVERY annotation script MUST call `load_config()` before doing any analysis.
If annotation_config.json does not exist, the script exits immediately.

Usage in any tier script:
    from annotation_gate import load_config
    config = load_config()   # SystemExit if config missing
"""

import json
import os
import sys

CONFIG_FILENAME = "annotation_config.json"


def _find_config(search_dirs=None):
    """Search for annotation_config.json in common locations."""
    if search_dirs is None:
        search_dirs = [
            ".",
            "annotation_output",
            os.environ.get("ANNOTATION_CONFIG_DIR", ""),
        ]

    for d in search_dirs:
        if not d:
            continue
        path = os.path.join(d, CONFIG_FILENAME)
        if os.path.isfile(path):
            return os.path.abspath(path)
    return None


def load_config(config_path=None):
    """
    Load and validate Phase 0 configuration.

    Raises SystemExit if config is missing or invalid.
    This is intentionally a HARD exit — no try/except can recover from missing config.

    Returns
    -------
    dict
        The validated configuration dictionary.
    """
    if config_path is None:
        config_path = _find_config()

    if config_path is None:
        msg = (
            "\n" + "=" * 70 + "\n"
            "HARD FAIL: annotation_config.json NOT FOUND\n" +
            "=" * 70 + "\n\n"
            "Phase 0 configuration has not been completed.\n\n"
            "You MUST run Phase 0 first:\n"
            "  python phase0_config.py\n\n"
            "This will interactively ask you to select:\n"
            "  1. TF Activity     -- tool / database / method\n"
            "  2. Pathway Activity -- tool / database / method\n"
            "  3. Trajectory       -- tool / root cell strategy\n"
            "  4. Re-clustering    -- strategy (A / B / C)\n"
            "  5. Resolution       -- fixed / multi-resolution scan\n\n"
            "After completion, annotation_config.json will be created\n"
            "and annotation can proceed.\n" +
            "=" * 70
        )
        print(msg, file=sys.stderr)
        sys.exit(1)

    with open(config_path) as f:
        config = json.load(f)

    # Validate required top-level keys
    required = ["tf_activity", "pathway_activity", "trajectory", "reclustering"]
    missing = [k for k in required if k not in config]
    if missing:
        print(
            f"\nHARD FAIL: annotation_config.json is INCOMPLETE.\n"
            f"Missing keys: {missing}\n"
            f"Re-run: python phase0_config.py\n",
            file=sys.stderr,
        )
        sys.exit(1)

    # Validate sub-keys
    tf = config["tf_activity"]
    for k in ["tool", "database", "method"]:
        if k not in tf or tf[k] is None:
            print(f"\nHARD FAIL: tf_activity.{k} is not set. Re-run Phase 0.", file=sys.stderr)
            sys.exit(1)

    pw = config["pathway_activity"]
    for k in ["tool", "database", "method"]:
        if k not in pw or pw[k] is None:
            print(f"\nHARD FAIL: pathway_activity.{k} is not set. Re-run Phase 0.", file=sys.stderr)
            sys.exit(1)

    traj = config["trajectory"]
    if "tool" not in traj or traj["tool"] is None:
        print("\nHARD FAIL: trajectory.tool is not set. Re-run Phase 0.", file=sys.stderr)
        sys.exit(1)

    rc = config["reclustering"]
    if "strategy" not in rc or rc["strategy"] is None:
        print("\nHARD FAIL: reclustering.strategy is not set. Re-run Phase 0.", file=sys.stderr)
        sys.exit(1)

    print(f"Phase 0 config loaded from: {config_path}")
    print(f"  TF Activity:  {tf['tool']} / {tf['database']} / {tf['method']}")
    print(f"  Pathway:      {pw['tool']} / {pw['database']} / {pw['method']}")
    print(f"  Trajectory:   {traj['tool']} / root={traj.get('root_cell', 'automatic')}")
    print(f"  Re-clustering: Strategy {rc['strategy']}")
    print(f"  Resolution:   {rc.get('resolution', {}).get('strategy', 'multi')}")
    return config


def store_config_in_adata(adata, config):
    """Store Phase 0 config in adata.uns for reproducibility."""
    adata.uns["annotation_config"] = config
    print("Config stored in adata.uns['annotation_config']")


if __name__ == "__main__":
    # Quick test: just try to load
    cfg = load_config()
    print("\nConfig is valid:")
    print(json.dumps(cfg, indent=2))
