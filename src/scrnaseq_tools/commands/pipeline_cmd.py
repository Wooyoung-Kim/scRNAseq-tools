"""CLI command: pipeline — Run the full scRNA-seq pipeline."""

import os


def register(subparsers) -> None:
    parser = subparsers.add_parser("pipeline", help="Run full QC-to-clustering pipeline")
    parser.add_argument("input", help="Input data path (h5ad, 10x directory)")
    parser.add_argument("-o", "--output", default="pipeline_output", help="Output directory")
    parser.add_argument("--batch-key", default=None, help="Batch key column name")
    parser.add_argument("--config", default=None, help="YAML config file")
    parser.add_argument("--resolution", type=float, default=1.0, help="Clustering resolution")
    parser.add_argument("--n-top-genes", type=int, default=2000, help="Number of HVGs")
    parser.add_argument("--agent", action="store_true", help="Run agent-friendly pipeline")
    parser.set_defaults(func=run)


def run(args) -> int:
    try:
        from scrnaseq_tools.analysis.pipeline import run_pipeline, run_agent_pipeline
        from scrnaseq_tools.analysis.config import PipelineConfig
    except ImportError:
        print("Analysis package not available.")
        print("Install with: pip install 'scRNAseq-tools[analysis]'")
        return 1

    if args.config:
        config = PipelineConfig.from_yaml(args.config)
    else:
        config = PipelineConfig(
            resolution=args.resolution,
            n_top_genes=args.n_top_genes,
        )

    if args.batch_key:
        config.batch_key = args.batch_key

    os.makedirs(args.output, exist_ok=True)

    if args.agent:
        print(f"Running agent pipeline: {args.input}")
        results = run_agent_pipeline(args.input, args.output, config=config)
    else:
        print(f"Running pipeline: {args.input}")
        results = run_pipeline(args.input, args.output, config=config)

    if results and results.adata is not None:
        print(f"Pipeline complete: {results.adata.n_obs} cells")
        print(f"Output: {args.output}")
        if results.files_created:
            for f in results.files_created:
                print(f"  Created: {f}")
    else:
        print("Pipeline failed. Check logs for details.")
        return 1

    return 0
