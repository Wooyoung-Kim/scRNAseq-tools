"""CLI command: qc — Run scRNA-seq quality control."""

import os


def register(subparsers) -> None:
    parser = subparsers.add_parser("qc", help="Run scRNA-seq quality control")
    parser.add_argument("input", help="Input data path (h5ad, 10x directory)")
    parser.add_argument("-o", "--output", default="qc_output", help="Output directory")
    parser.add_argument("--adaptive", action="store_true", default=True, help="Use adaptive QC (MAD-based)")
    parser.add_argument("--nmads", type=float, default=3.0, help="Number of MADs for thresholds")
    parser.add_argument("--max-mt-pct", type=float, default=None, help="Max mitochondrial percentage")
    parser.add_argument("--min-genes", type=int, default=200, help="Min genes per cell")
    parser.add_argument("--min-cells", type=int, default=3, help="Min cells per gene")
    parser.set_defaults(func=run)


def run(args) -> int:
    try:
        from scrnaseq_tools.analysis.pipeline import load_data, run_qc
        from scrnaseq_tools.analysis.config import PipelineConfig
    except ImportError:
        print("Analysis package not available.")
        print("Install with: pip install 'scRNAseq-tools[analysis]'")
        return 1

    config = PipelineConfig(
        adaptive_qc=args.adaptive,
        qc_nmads=args.nmads,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
    )
    if args.max_mt_pct is not None:
        config.max_pct_mt = args.max_mt_pct

    os.makedirs(args.output, exist_ok=True)

    print(f"Loading data: {args.input}")
    adata = load_data(args.input)
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    print("Running QC...")
    adata, qc_summary = run_qc(adata, config)
    print(f"  After QC: {adata.n_obs} cells × {adata.n_vars} genes")

    out_path = os.path.join(args.output, "qc_filtered.h5ad")
    adata.write(out_path)
    print(f"  Saved: {out_path}")

    if qc_summary is not None:
        csv_path = os.path.join(args.output, "qc_summary.csv")
        qc_summary.to_csv(csv_path)
        print(f"  Summary: {csv_path}")

    return 0
