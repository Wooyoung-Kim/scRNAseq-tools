# Reproducibility Guide

This project aims to provide deterministic, reproducible analyses. Use the steps below to
create a stable environment and record every run.

## 1) Environment setup (conda)

```bash
conda env create -f environment.yml
conda activate scrnaseq-tools
python -m pip install -e .
```

## 2) Environment setup (Docker)

```bash
docker build -t scrnaseq-tools:latest .
docker run --rm -it -v "$PWD:/workspace" scrnaseq-tools:latest scRNAseq-tools --help
```

## 3) Run manifest

Every CLI command writes `run_manifest.json` in the current directory by default.
You can disable it or change the output path:

```bash
export SCRNASEQ_MANIFEST=0
export SCRNASEQ_MANIFEST_PATH=/path/to/run_manifest.json
```

The manifest records:
- timestamp (UTC)
- command and argv
- exit code and error (if any)
- working directory, user, host
- Python and scRNAseq-tools versions

## 4) Randomness control

Set a fixed seed in your analysis scripts and keep it in your run notes. For Scanpy:

```python
import numpy as np
import scanpy as sc

np.random.seed(0)
sc.settings.set_figure_params(dpi=100)
```

## 5) Suggested citation

See `CITATION.cff` for citation metadata.
