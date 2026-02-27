"""CLI command: profile — Profile scRNA-seq data for agent consumption."""

import json


def register(subparsers) -> None:
    parser = subparsers.add_parser("profile", help="Profile scRNA-seq data (agent-friendly)")
    parser.add_argument("input", help="Input data path (h5ad, 10x directory)")
    parser.add_argument("--json", action="store_true", dest="output_json", help="Output as JSON")
    parser.add_argument("-o", "--output", default=None, help="Save profile to file")
    parser.set_defaults(func=run)


def run(args) -> int:
    try:
        from scrnaseq_tools.analysis.profile import profile_data
    except ImportError:
        print("Analysis package not available.")
        print("Install with: pip install 'scRNAseq-tools[analysis]'")
        return 1

    print(f"Profiling: {args.input}")
    profile = profile_data(args.input)

    if args.output_json:
        output = profile.to_json()
        print(output)
    else:
        print(profile.summary())

    if args.output:
        with open(args.output, "w") as f:
            f.write(profile.to_json())
        print(f"Profile saved: {args.output}")

    return 0
