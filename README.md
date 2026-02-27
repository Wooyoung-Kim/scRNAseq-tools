# scRNAseq-tools

CLI-first toolkit for single-cell RNA-seq workflows. The default mode does not use any external LLM API.

## Quick start

```bash
python -m pip install -e .                # Basic install (init/chat/scan)
python -m pip install -e ".[analysis]"    # + QC/pipeline/integration
python -m pip install -e ".[all]"         # + scVI/Harmony/Scanorama/doublet/MCP
scRNAseq-tools --help
```

## Analysis pipeline (v0.2.0+)

Integrated scRNA-seq QC-to-Integration pipeline. Requires `[analysis]` extras.

### Run QC only
```bash
scRNAseq-tools qc data.h5ad -o output/ --nmads 3.0
```

### Run full pipeline
```bash
scRNAseq-tools pipeline data.h5ad -o output/ --batch-key sample --resolution 0.5
scRNAseq-tools pipeline data.h5ad -o output/ --agent   # Agent-friendly mode
```

### Batch integration
```bash
scRNAseq-tools integrate processed.h5ad -o output/ --batch-key sample --methods harmony scvi
```

### Data profiling (agent-friendly)
```bash
scRNAseq-tools profile data.h5ad --json
```

### Optional dependency groups
| Group | Contents |
|-------|----------|
| `analysis` | scanpy, anndata, numpy, scipy, pandas, matplotlib, seaborn, leidenalg |
| `scvi` | analysis + scvi-tools |
| `harmony` | analysis + harmonypy |
| `scanorama` | analysis + scanorama |
| `doublet` | analysis + scrublet |
| `benchmark` | analysis + scib-metrics |
| `mcp` | mcp |
| `all` | Everything above |

## Step-by-step guides

### Initialize a project (.claude only)
1) `python -m pip install -e .`
2) `scRNAseq-tools init --target /path/to/project`
3) `cd /path/to/project` and launch `claude` or `codex`

### Bootstrap a project (includes MCP deps)
1) `scRNAseq-tools bootstrap --target /path/to/project --force`
2) `cd /path/to/project` and launch `claude` or `codex`
3) Use `mcp__scrnaseq-tools__scan_data` or `mcp__scrnaseq-tools__summarize`

### Chat with a local LLM (no API key)
1) Start a local server (Ollama, vLLM, or LM Studio)
2) Run one of:
   `scRNAseq-tools chat --provider ollama`
   `scRNAseq-tools chat --provider vllm --base-url http://localhost:8000 --model local-model`
   `scRNAseq-tools chat --provider lmstudio --base-url http://localhost:1234 --model local-model`
3) Type `/exit` to leave

### Chat with API keys (usage-based billing)
1) `export OPENAI_API_KEY=...` or `export ANTHROPIC_API_KEY=...`
2) Run one of:
   `scRNAseq-tools chat --provider openai --model gpt-4.1-mini`
   `scRNAseq-tools chat --provider anthropic --model claude-3-5-sonnet-20240620`
3) Type `/exit` to leave

### Chat via CLI login (Claude Code / Codex)
1) Ensure `claude` or `codex` CLI is installed and logged in
2) Run:
   `scRNAseq-tools chat --provider claude-cli`
   `scRNAseq-tools chat --provider codex-cli`
3) The CLI uses the current directory to load `.claude`

### Local REPL (no LLM)
1) `scRNAseq-tools repl`
2) Use `scan`, `summary`, and `env`

### Non-interactive data commands
1) `scRNAseq-tools scan .`
2) `scRNAseq-tools scan /path/to/data -r --json`
3) `scRNAseq-tools summary /path/to/file.h5ad`

### Reproducible setup (conda/Docker)
1) `conda env create -f environment.yml`
2) `conda activate scrnaseq-tools`
3) `python -m pip install -e .`

Docker:
```bash
docker build -t scrnaseq-tools:latest .
docker run --rm -it -v "$PWD:/workspace" scrnaseq-tools:latest scRNAseq-tools --help
```

## LLM usage (optional)

By default, no API calls are made. To enable an external model, set a provider and key:

```bash
export LLM_PROVIDER=openai
export OPENAI_API_KEY=your_key_here
export LLM_MODEL=gpt-4.1-mini
```

Supported providers (built-in): `none`, `openai`, `anthropic`, `ollama`.

### LLM configuration combinations

Provider resolution order:
1) `LLM_PROVIDER`
2) `OPENAI_API_KEY` (implies `openai`)
3) `ANTHROPIC_API_KEY` (implies `anthropic`)
4) default `none`

Model resolution order:
1) `LLM_MODEL`
2) `MODEL`
3) `OPENAI_MODEL`
4) `ANTHROPIC_MODEL`

Common combinations:
- `LLM_PROVIDER=none` (or unset) and no API keys: LLM disabled.
- `OPENAI_API_KEY=...` with no `LLM_PROVIDER`: provider auto-detects as `openai`.
- `ANTHROPIC_API_KEY=...` with no `LLM_PROVIDER`: provider auto-detects as `anthropic`.
- `LLM_PROVIDER=ollama`: uses local Ollama, no API key required.
- `LLM_PROVIDER=openai` but no `OPENAI_API_KEY`: provider set, key missing.
- `LLM_PROVIDER=anthropic` but no `ANTHROPIC_API_KEY`: provider set, key missing.

## Project bootstrap

Install the default `.claude` bundle into a project directory:

```bash
scRNAseq-tools init --target /path/to/project
```

If the default source path is not found, an embedded `.claude` template is used.

Bootstrap a project and install MCP server dependencies:

```bash
scRNAseq-tools bootstrap --target /path/to/project
```

`bootstrap` uses `pip` and may require network access depending on server dependencies.

Override the source directory:

```bash
scRNAseq-tools init --source /path/to/.claude --target /path/to/project
```

You can also set `SCRNASEQ_TOOLS_CLAUDE_SOURCE` to change the default source path.

## MCP integration

The embedded `.claude` template includes an MCP server named `scrnaseq-tools` that exposes
local analysis helpers (`list_dir`, `scan_data`, `summarize`, `env`). Run:

```bash
scRNAseq-tools bootstrap --target /path/to/project
```

Then open Claude Code/Codex in that project to access `mcp__scrnaseq-tools__*` tools.

### CLI command combinations

General commands:
```bash
scRNAseq-tools --help
scRNAseq-tools info
scRNAseq-tools doctor
scRNAseq-tools chat
scRNAseq-tools repl
scRNAseq-tools scan
scRNAseq-tools summary
scRNAseq-tools version
```

`init` options and behavior:
- `--target DIR`: install `.claude` into `DIR/.claude` (default: current directory).
- `--source PATH`: explicit source `.claude` directory.
- `--force`: overwrite existing `DIR/.claude`.
- Source resolution order when `--source` is not provided:
  1) `SCRNASEQ_TOOLS_CLAUDE_SOURCE`
  2) `/home/kwy7605/data_61/Vaccine_V2/.claude` (if it exists)
  3) embedded template inside the package
- If `DIR/.claude` already exists and `--force` is not set, `init` exits with an error.

Examples:
```bash
scRNAseq-tools init --target /path/to/project
scRNAseq-tools init --target /path/to/project --force
scRNAseq-tools init --source /path/to/.claude --target /path/to/project
SCRNASEQ_TOOLS_CLAUDE_SOURCE=/path/to/.claude scRNAseq-tools init --target /path/to/project
```

`bootstrap` options and behavior:
- `--target DIR`: same target semantics as `init`.
- `--source PATH`: same source semantics as `init`.
- `--force`: overwrites existing `.claude` before installing MCP dependencies.
- `--no-mcp-install`: skip MCP dependency installation.
- If `DIR/.claude` exists and `--force` is not set, `bootstrap` uses the existing directory.

Examples:
```bash
scRNAseq-tools bootstrap --target /path/to/project
scRNAseq-tools bootstrap --target /path/to/project --force
scRNAseq-tools bootstrap --target /path/to/project --no-mcp-install
scRNAseq-tools bootstrap --source /path/to/.claude --target /path/to/project
```

## Chat mode

`chat` prompts for a provider first, then starts a REPL session. Press Enter to use the default (`claude-cli`):

```bash
scRNAseq-tools chat
```

Skip the prompt with `--provider`:

```bash
scRNAseq-tools chat --provider ollama
scRNAseq-tools chat --provider vllm --base-url http://localhost:8000 --model local-model
scRNAseq-tools chat --provider lmstudio --base-url http://localhost:1234 --model local-model
scRNAseq-tools chat --provider openai --model gpt-4.1-mini
scRNAseq-tools chat --provider anthropic --model claude-3-5-sonnet-20240620
```

Local providers:
- `ollama`: uses `OLLAMA_HOST` (default `http://localhost:11434`) and `OLLAMA_MODEL`.
- `vllm`: OpenAI-compatible; uses `VLLM_BASE_URL` (default `http://localhost:8000`) and `VLLM_MODEL`.
- `lmstudio`: OpenAI-compatible; uses `LMSTUDIO_BASE_URL` (default `http://localhost:1234`) and `LMSTUDIO_MODEL`.

Cloud providers:
- `openai`: requires `OPENAI_API_KEY`, optional `OPENAI_MODEL` or `LLM_MODEL`.
- `anthropic`: requires `ANTHROPIC_API_KEY`, optional `ANTHROPIC_MODEL` or `LLM_MODEL`.

CLI login reuse:
- `codex-cli`: runs `codex` (override with `SCRNASEQ_CODEX_CMD` or `--cli-cmd`)
- `claude-cli`: runs `claude` (override with `SCRNASEQ_CLAUDE_CMD` or `--cli-cmd`)
- Use `--cli-args` to pass extra arguments to the CLI command.
- If Codex is detected behind `claude`, an `AGENTS.md` file is auto-generated when `.claude` exists.
- Disable auto-detect with `SCRNASEQ_DISABLE_CLI_DETECT=1`.
- Disable AGENTS bootstrap with `SCRNASEQ_DISABLE_AGENTS_BOOTSTRAP=1`.
- When `claude-cli` runs and `.claude` exists, a `CLAUDE.md` file is auto-generated if missing.
- Disable CLAUDE.md bootstrap with `SCRNASEQ_DISABLE_CLAUDE_BOOTSTRAP=1`.

Chat commands:
- `/exit` or `/quit`: exit
- `/reset`: clear history
- `/help`: show commands

## Local REPL

`repl` runs a built-in local analysis console without external LLMs:

```bash
scRNAseq-tools repl
```

Commands:
- `help`: show commands
- `pwd`: print working directory
- `cd <path>`: change directory
- `ls [path]`: list directory contents
- `scan [path] [-r]`: list data files (`.h5ad`, `.mtx`, `.csv`, ...)
- `summary <path>`: quick summary (uses `anndata` for `.h5ad`, `pandas` for tables)
- `env`: show installed analysis packages

## Data Commands

Use non-interactive commands for scripts and pipelines:

```bash
scRNAseq-tools scan [path] [-r] [--json]
scRNAseq-tools summary <path> [--json]
```

`scan`:
- Lists matching data files (`.h5ad`, `.loom`, `.mtx`, `.csv`, `.tsv`, ...).
- Returns JSON with `--json` (includes file count and file list).

`summary`:
- Summarizes `.h5ad` (shape, obs/var/layers keys), table files (`.csv/.tsv/.txt`), and generic files.
- Returns JSON with `--json`.

## Reproducibility

Every CLI command writes `run_manifest.json` in the current directory by default.
You can disable or redirect it:

```bash
export SCRNASEQ_MANIFEST=0
export SCRNASEQ_MANIFEST_PATH=/path/to/run_manifest.json
```

See `reproducibility.md` for the full guide.

## Add new commands

Create a new module in `src/scrnaseq_tools/commands/` and register it in
`src/scrnaseq_tools/commands/__init__.py`.
