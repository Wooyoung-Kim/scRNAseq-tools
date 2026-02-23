#!/usr/bin/env python3
"""
Phase 0: Annotation Configuration with Preset Profiles.

Usage:
    python phase0_config.py --list                    # Show all presets
    python phase0_config.py --preset standard          # Apply preset directly
    python phase0_config.py --preset custom             # Interactive mode
    python phase0_config.py --preset standard -o out/   # Custom output dir
"""

import json
import os
import argparse
import sys

# ═══════════════════════════════════════════════════════════════════
#  PRESET DEFINITIONS — all valid config combinations
# ═══════════════════════════════════════════════════════════════════

PRESETS = {
    "standard": {
        "label": "Standard (Recommended)",
        "description": (
            "decoupler CollecTRI/ULM + PROGENy/MLM + Palantir + "
            "Harmony reuse (A) + multi-resolution silhouette"
        ),
        "config": {
            "tf_activity": {"tool": "decoupler", "database": "collectri", "method": "ulm"},
            "pathway_activity": {"tool": "decoupler", "database": "progeny_500", "method": "mlm"},
            "trajectory": {"tool": "palantir", "use": True, "root_cell": "automatic"},
            "reclustering": {
                "strategy": "A",
                "neighbors": {"n_neighbors": 15, "n_pcs": 50},
                "resolution": {
                    "strategy": "multi",
                    "multi": {"range": [0.3, 0.5, 0.8, 1.0, 1.2, 1.5], "metric": "silhouette"},
                },
            },
        },
    },
    "thorough": {
        "label": "Thorough",
        "description": (
            "decoupler CollecTRI/ULM + PROGENy/MLM + Palantir + "
            "Full re-Harmony (C) + multi-resolution silhouette"
        ),
        "config": {
            "tf_activity": {"tool": "decoupler", "database": "collectri", "method": "ulm"},
            "pathway_activity": {"tool": "decoupler", "database": "progeny_500", "method": "mlm"},
            "trajectory": {"tool": "palantir", "use": True, "root_cell": "automatic"},
            "reclustering": {
                "strategy": "C",
                "hvg": {"flavor": "seurat_v3", "n_top_genes": 3000, "batch_key": "batch"},
                "integration": {"tool": "harmony", "batch_key": "batch"},
                "neighbors": {"n_neighbors": 15, "n_pcs": 50},
                "resolution": {
                    "strategy": "multi",
                    "multi": {"range": [0.3, 0.5, 0.8, 1.0, 1.2, 1.5], "metric": "silhouette"},
                },
            },
        },
    },
    "quick": {
        "label": "Quick",
        "description": (
            "decoupler CollecTRI/ULM + PROGENy/MLM + DPT (fast) + "
            "Harmony reuse (A) + fixed resolution"
        ),
        "config": {
            "tf_activity": {"tool": "decoupler", "database": "collectri", "method": "ulm"},
            "pathway_activity": {"tool": "decoupler", "database": "progeny_500", "method": "mlm"},
            "trajectory": {"tool": "dpt", "use": True, "root_cell": "automatic"},
            "reclustering": {
                "strategy": "A",
                "neighbors": {"n_neighbors": 15, "n_pcs": 50},
                "resolution": {
                    "strategy": "fixed",
                    "fixed": {"tier2": 0.8, "tier3": 1.0},
                },
            },
        },
    },
    "no_trajectory": {
        "label": "No Trajectory",
        "description": (
            "decoupler CollecTRI/ULM + PROGENy/MLM + Skip trajectory + "
            "Harmony reuse (A) + multi-resolution"
        ),
        "config": {
            "tf_activity": {"tool": "decoupler", "database": "collectri", "method": "ulm"},
            "pathway_activity": {"tool": "decoupler", "database": "progeny_500", "method": "mlm"},
            "trajectory": {"tool": "skip", "use": False, "root_cell": "none"},
            "reclustering": {
                "strategy": "A",
                "neighbors": {"n_neighbors": 15, "n_pcs": 50},
                "resolution": {
                    "strategy": "multi",
                    "multi": {"range": [0.3, 0.5, 0.8, 1.0, 1.2, 1.5], "metric": "silhouette"},
                },
            },
        },
    },
    "dorothea_hallmark": {
        "label": "DoRothEA + MSigDB Hallmark",
        "description": (
            "decoupler DoRothEA(A,B)/ULM + MSigDB Hallmark/MLM + Palantir + "
            "Harmony reuse (A) + multi-resolution"
        ),
        "config": {
            "tf_activity": {"tool": "decoupler", "database": "dorothea_ab", "method": "ulm"},
            "pathway_activity": {"tool": "decoupler", "database": "msigdb_hallmark", "method": "mlm"},
            "trajectory": {"tool": "palantir", "use": True, "root_cell": "automatic"},
            "reclustering": {
                "strategy": "A",
                "neighbors": {"n_neighbors": 15, "n_pcs": 50},
                "resolution": {
                    "strategy": "multi",
                    "multi": {"range": [0.3, 0.5, 0.8, 1.0, 1.2, 1.5], "metric": "silhouette"},
                },
            },
        },
    },
    "gsea_based": {
        "label": "GSEA-based Pathway",
        "description": (
            "decoupler CollecTRI/ULM + GSEApy MSigDB Hallmark/prerank + Palantir + "
            "Harmony reuse (A) + multi-resolution"
        ),
        "config": {
            "tf_activity": {"tool": "decoupler", "database": "collectri", "method": "ulm"},
            "pathway_activity": {"tool": "gseapy", "database": "msigdb_hallmark", "method": "gsea_prerank"},
            "trajectory": {"tool": "palantir", "use": True, "root_cell": "automatic"},
            "reclustering": {
                "strategy": "A",
                "neighbors": {"n_neighbors": 15, "n_pcs": 50},
                "resolution": {
                    "strategy": "multi",
                    "multi": {"range": [0.3, 0.5, 0.8, 1.0, 1.2, 1.5], "metric": "silhouette"},
                },
            },
        },
    },
    "new_hvg": {
        "label": "New HVG per subset (no re-integration)",
        "description": (
            "decoupler CollecTRI/ULM + PROGENy/MLM + Palantir + "
            "New HVG+PCA (B) + multi-resolution"
        ),
        "config": {
            "tf_activity": {"tool": "decoupler", "database": "collectri", "method": "ulm"},
            "pathway_activity": {"tool": "decoupler", "database": "progeny_500", "method": "mlm"},
            "trajectory": {"tool": "palantir", "use": True, "root_cell": "automatic"},
            "reclustering": {
                "strategy": "B",
                "hvg": {"flavor": "seurat_v3", "n_top_genes": 3000, "batch_key": "batch"},
                "neighbors": {"n_neighbors": 15, "n_pcs": 50},
                "resolution": {
                    "strategy": "multi",
                    "multi": {"range": [0.3, 0.5, 0.8, 1.0, 1.2, 1.5], "metric": "silhouette"},
                },
            },
        },
    },
}

# ═══════════════════════════════════════════════════════════════════
#  Interactive custom mode (for terminal use)
# ═══════════════════════════════════════════════════════════════════

def ask_choice(prompt, options, default=None):
    """Ask user to select from numbered options."""
    print(f"\n  {prompt}")
    for i, (label, desc) in enumerate(options, 1):
        marker = " *" if default and label == default else ""
        print(f"    [{i}] {label}{marker}")
        if desc:
            print(f"        {desc}")
    while True:
        choice = input(f"  Select [1-{len(options)}] (default={default}): ").strip()
        if not choice and default:
            for i, (label, _) in enumerate(options):
                if label == default:
                    return label
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(options):
                return options[idx][0]
        except ValueError:
            pass


def ask_value(prompt, default=None):
    val = input(f"  {prompt} (default={default}): ").strip()
    return val if val else default


def build_custom_config():
    """Interactive custom configuration builder."""
    config = {}

    print("\n>>> TF Activity (Tier 2)")
    tf_tool = ask_choice("Tool:", [
        ("decoupler", "Linear model-based"), ("pySCENIC", "Network inference"),
    ], default="decoupler")
    tf_db = ask_choice("Database:", [
        ("collectri", "Curated TF-targets"), ("dorothea_ab", "DoRothEA A+B"),
        ("dorothea_abc", "DoRothEA A+B+C"),
    ], default="collectri") if tf_tool == "decoupler" else "pyscenic"
    tf_method = ask_choice("Method:", [
        ("ulm", "Univariate Linear Model"), ("mlm", "Multivariate Linear Model"),
        ("viper", "VIPER"), ("aucell", "AUCell"),
    ], default="ulm") if tf_tool == "decoupler" else "pyscenic"
    config["tf_activity"] = {"tool": tf_tool, "database": tf_db, "method": tf_method}

    print("\n>>> Pathway Activity (Tier 3)")
    pw_tool = ask_choice("Tool:", [
        ("decoupler", "Linear model-based"), ("gseapy", "GSEA-based"),
    ], default="decoupler")
    if pw_tool == "decoupler":
        pw_db = ask_choice("Database:", [
            ("progeny_500", "PROGENy top 500"), ("progeny_100", "PROGENy top 100"),
            ("msigdb_hallmark", "MSigDB Hallmark"),
        ], default="progeny_500")
        pw_method = ask_choice("Method:", [
            ("mlm", "MLM"), ("ulm", "ULM"), ("aucell", "AUCell"),
        ], default="mlm")
    else:
        pw_db = ask_choice("Database:", [
            ("msigdb_hallmark", "Hallmark"), ("kegg", "KEGG"),
            ("reactome", "Reactome"), ("go_bp", "GO BP"),
        ], default="msigdb_hallmark")
        pw_method = ask_choice("Method:", [
            ("gsea_prerank", "Prerank"), ("ora", "ORA"), ("ssgsea", "ssGSEA"),
        ], default="gsea_prerank")
    config["pathway_activity"] = {"tool": pw_tool, "database": pw_db, "method": pw_method}

    print("\n>>> Trajectory (Tier 2)")
    traj = ask_choice("Tool:", [
        ("palantir", "Diffusion pseudotime + fate"), ("dpt", "Simple DPT"),
        ("skip", "Skip trajectory"),
    ], default="palantir")
    root = "none" if traj == "skip" else ask_choice("Root cell:", [
        ("automatic", "Auto from DE"), ("manual", "Manual"),
    ], default="automatic")
    config["trajectory"] = {"tool": traj, "use": traj != "skip", "root_cell": root}

    print("\n>>> Re-clustering")
    strategy = ask_choice("Strategy:", [
        ("A", "Reuse Harmony"), ("B", "New HVG+PCA"), ("C", "New HVG+PCA+Harmony"),
    ], default="A")
    rc = {"strategy": strategy, "neighbors": {"n_neighbors": 15, "n_pcs": 50}}
    if strategy in ("B", "C"):
        rc["hvg"] = {"flavor": "seurat_v3", "n_top_genes": 3000, "batch_key": "batch"}
    if strategy == "C":
        rc["integration"] = {"tool": "harmony", "batch_key": "batch"}

    print("\n>>> Resolution")
    res = ask_choice("Strategy:", [
        ("multi", "Multi-resolution scan"), ("fixed", "Fixed"),
    ], default="multi")
    if res == "multi":
        rc["resolution"] = {
            "strategy": "multi",
            "multi": {"range": [0.3, 0.5, 0.8, 1.0, 1.2, 1.5], "metric": "silhouette"},
        }
    else:
        rc["resolution"] = {
            "strategy": "fixed",
            "fixed": {"tier2": 0.8, "tier3": 1.0},
        }
    config["reclustering"] = rc
    return config


# ═══════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════

def list_presets():
    """Print all available presets."""
    print("\n" + "=" * 70)
    print("  Available Annotation Config Presets")
    print("=" * 70)
    for i, (name, preset) in enumerate(PRESETS.items(), 1):
        print(f"\n  [{i}] {name}")
        print(f"      {preset['label']}")
        print(f"      {preset['description']}")
        c = preset["config"]
        print(f"      TF:          {c['tf_activity']['tool']}/{c['tf_activity']['database']}/{c['tf_activity']['method']}")
        print(f"      Pathway:     {c['pathway_activity']['tool']}/{c['pathway_activity']['database']}/{c['pathway_activity']['method']}")
        print(f"      Trajectory:  {c['trajectory']['tool']}")
        print(f"      Recluster:   Strategy {c['reclustering']['strategy']}")
        print(f"      Resolution:  {c['reclustering']['resolution']['strategy']}")
    print("\n" + "=" * 70)
    print(f"  Total: {len(PRESETS)} presets")
    print(f"  Usage: python phase0_config.py --preset <name>")
    print(f"         python phase0_config.py --preset custom  (interactive)")
    print("=" * 70)


def save_config(config, output_path):
    """Save config JSON and print summary."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(config, f, indent=2)

    print(f"\n{'='*70}")
    print(f"  Config saved: {os.path.abspath(output_path)}")
    print(f"{'='*70}")
    print(f"  TF Activity:   {config['tf_activity']['tool']} / "
          f"{config['tf_activity']['database']} / {config['tf_activity']['method']}")
    print(f"  Pathway:       {config['pathway_activity']['tool']} / "
          f"{config['pathway_activity']['database']} / {config['pathway_activity']['method']}")
    print(f"  Trajectory:    {config['trajectory']['tool']} / "
          f"root={config['trajectory']['root_cell']}")
    print(f"  Re-clustering: Strategy {config['reclustering']['strategy']}")
    print(f"  Resolution:    {config['reclustering']['resolution']['strategy']}")
    print(f"{'='*70}")
    print("  Annotation scripts can now proceed.")


def main():
    parser = argparse.ArgumentParser(
        description="Phase 0: Annotation Configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python phase0_config.py --list                 # Show all presets\n"
            "  python phase0_config.py --preset standard       # Use standard preset\n"
            "  python phase0_config.py --preset custom          # Interactive mode\n"
            "  python phase0_config.py --select                 # Pick preset by number\n"
        ),
    )
    parser.add_argument("--list", action="store_true", help="List all preset profiles")
    parser.add_argument("--preset", type=str, help="Apply a preset (name or 'custom')")
    parser.add_argument("--select", action="store_true", help="Interactively select a preset by number")
    parser.add_argument("-o", "--output", default="annotation_config.json", help="Output JSON path")
    args = parser.parse_args()

    if args.list:
        list_presets()
        return

    if args.select:
        list_presets()
        names = list(PRESETS.keys()) + ["custom"]
        print(f"\n  [  ] custom — Build your own interactively")
        while True:
            choice = input(f"\n  Select preset [1-{len(names)}] or name: ").strip()
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(PRESETS):
                    preset_name = list(PRESETS.keys())[idx]
                    break
                elif idx == len(PRESETS):
                    preset_name = "custom"
                    break
            except ValueError:
                if choice in PRESETS or choice == "custom":
                    preset_name = choice
                    break
            print("  Invalid choice.")

        if preset_name == "custom":
            config = build_custom_config()
        else:
            config = PRESETS[preset_name]["config"]
            print(f"\n  Selected: {preset_name} — {PRESETS[preset_name]['label']}")
        save_config(config, args.output)
        return

    if args.preset:
        if args.preset == "custom":
            config = build_custom_config()
        elif args.preset in PRESETS:
            config = PRESETS[args.preset]["config"]
            print(f"\n  Applying preset: {args.preset} — {PRESETS[args.preset]['label']}")
        else:
            print(f"  Unknown preset: {args.preset}", file=sys.stderr)
            print(f"  Available: {', '.join(PRESETS.keys())}, custom", file=sys.stderr)
            sys.exit(1)
        save_config(config, args.output)
        return

    # No arguments — show help
    parser.print_help()
    print("\n  Tip: Use --list to see presets, --preset <name> to apply one.")


if __name__ == "__main__":
    main()
