"""CLI command: integrate — Run batch integration."""

import os


def register(subparsers) -> None:
    parser = subparsers.add_parser("integrate", help="Run batch integration methods")
    parser.add_argument("input", help="Input h5ad file (processed)")
    parser.add_argument("-o", "--output", default="integration_output", help="Output directory")
    parser.add_argument("--batch-key", required=True, help="Batch key column name")
    parser.add_argument(
        "--methods",
        nargs="+",
        default=["harmony"],
        choices=["scvi", "harmony", "scanorama", "all"],
        help="Integration methods to run",
    )
    parser.add_argument("--benchmark", action="store_true", help="Run integration benchmarking")
    parser.set_defaults(func=run)


def run(args) -> int:
    try:
        from scrnaseq_tools.analysis.integration import run_integration
        from scrnaseq_tools.analysis.config import IntegrationConfig
    except ImportError:
        print("Analysis package not available.")
        print("Install with: pip install 'scRNAseq-tools[analysis]'")
        return 1

    try:
        import scanpy as sc
    except ImportError:
        print("scanpy not installed. Install with: pip install 'scRNAseq-tools[analysis]'")
        return 1

    methods = args.methods
    if "all" in methods:
        methods = ["scvi", "harmony", "scanorama"]

    config = IntegrationConfig(
        methods=methods,
        auto_select_best=args.benchmark,
    )

    os.makedirs(args.output, exist_ok=True)

    print(f"Loading: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    print(f"Running integration: {', '.join(methods)}")
    result = run_integration(adata, args.batch_key, config)

    out_path = os.path.join(args.output, "integrated.h5ad")
    result.write(out_path)
    print(f"  Saved: {out_path}")

    if args.benchmark:
        try:
            from scrnaseq_tools.analysis.integration import benchmark_integration
            embeddings = [f"X_{m}" for m in methods if f"X_{m}" in result.obsm]
            if embeddings:
                scores = benchmark_integration(result, args.batch_key, embeddings)
                print(f"  Benchmark scores: {scores}")
        except ImportError:
            print("  Benchmarking requires: pip install 'scRNAseq-tools[all]'")

    return 0
