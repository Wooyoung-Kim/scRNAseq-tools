# Annotation Executor Agent

Autonomous agent for 3-tier hierarchical cell type annotation.

---

## Agent Specification

```yaml
name: annotation-executor
type: general-purpose
description: "Executes 3-tier hierarchical annotation with checkpoints"

triggers:
  - "자동으로 annotation 해줘"
  - "자동 annotation"
  - "run annotation automatically"
  - "execute hierarchical annotation"

capabilities:
  - Read/write AnnData files
  - Execute DE analysis (scanpy)
  - TF activity inference (decoupler)
  - Pathway activity inference (decoupler)
  - Trajectory analysis (palantir)
  - PubMed search (mcp__pubmed-tools)
  - Generate visualizations
  - Save intermediate results
```

---

## PHASE 0 HARD FAIL (READ THIS FIRST)

```
╔══════════════════════════════════════════════════════════════════════╗
║  HARD FAIL: Skipping Phase 0 is FORBIDDEN                           ║
║                                                                      ║
║  Before writing ANY code or running ANY analysis:                    ║
║                                                                      ║
║  1. Read phases/user_config.md                                       ║
║  2. Use AskUserQuestion to present EACH selection to the user:       ║
║     A. TF Activity: tool / database / method                         ║
║     B. Pathway Activity: tool / database / method                    ║
║     C. Trajectory: tool / root cell strategy                         ║
║     D. Re-clustering: strategy (A/B/C)                               ║
║     E. Resolution: fixed / multi-resolution scan                     ║
║  3. Wait for user confirmation of ALL settings                       ║
║  4. Store config in adata.uns['annotation_config']                   ║
║  5. Save Phase 0 checkpoint                                         ║
║                                                                      ║
║  ONLY THEN proceed to Tier 1.                                        ║
║                                                                      ║
║  NEVER hardcode tool/method/database choices                         ║
║  NEVER assume defaults without asking                                ║
║  NEVER skip this step even if "just testing"                         ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Phase 0 Implementation (Preset-Based)

The agent MUST use the preset system in `phase0_config.py`:

```bash
# Step 1: Show presets to decide what to present
python phase0_config.py --list

# Step 2: Use AskUserQuestion to let user pick a preset
# Present the 7 presets: standard, thorough, quick, no_trajectory,
#   dorothea_hallmark, gsea_based, new_hvg, custom

# Step 3: Apply the user's choice
python phase0_config.py --preset <user_choice>
# If user picks "custom": run in terminal (python phase0_config.py --preset custom)

# Step 4: Verify config was created
python annotation_gate.py
```

```python
# In every tier script — MANDATORY first lines:
from annotation_gate import load_config
config = load_config()   # sys.exit(1) if config missing

# Read settings from config — NEVER hardcode:
tf_tool = config['tf_activity']['tool']
tf_db   = config['tf_activity']['database']
tf_method = config['tf_activity']['method']
pw_tool = config['pathway_activity']['tool']
strategy = config['reclustering']['strategy']

# Store config in adata for reproducibility:
from annotation_gate import store_config_in_adata
store_config_in_adata(adata, config)
```

---

## Execution Flow

```
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 0: Configuration  — MANDATORY, HARD FAIL IF SKIPPED      │
│  ─────────────────────────────────────────────────────────────  │
│  1. Load phases/user_config.md                                   │
│  2. AskUserQuestion: TF tool / database / method                 │
│  3. AskUserQuestion: Pathway tool / database / method            │
│  4. AskUserQuestion: Trajectory tool / root cell                 │
│  5. AskUserQuestion: Re-clustering strategy (A/B/C)              │
│  6. AskUserQuestion: Resolution strategy (fixed/multi)           │
│  7. Store config in adata.uns['annotation_config']               │
│  ────────────────────────────────────────────────────────────── │
│  CHECKPOINT: User confirms ALL configuration choices             │
│  CANNOT proceed to Tier 1 without user confirmation              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 1: Tier 1 - Major Types                                  │
│  ─────────────────────────────────────────────────────────────  │
│  1. Load principles.md + tier1.md                               │
│  2. Compute DE for full dataset                                 │
│  3. PubMed search for marker combinations                       │
│  4. 3-iteration reasoning                                       │
│  5. Assign major types to clusters                              │
│  6. Save: tier1_annotation, tier1_de_results                    │
│  ────────────────────────────────────────────────────────────── │
│  CHECKPOINT: User reviews Tier 1 results                        │
│  Output: UMAP + Dotplot + Reasoning summary                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 2: Tier 2 - Developmental States (per Major Type)        │
│  ─────────────────────────────────────────────────────────────  │
│  FOR EACH major_type:                                           │
│    1. Subset data                                               │
│    2. Preprocess (strategy A/B/C)                               │
│    3. Select resolution (fixed/multi)                           │
│    4. RE-COMPUTE DE within subset                               │
│    5. TF activity analysis                                      │
│    6. Trajectory analysis (if applicable)                       │
│    7. PubMed search + verify                                    │
│    8. 3-iteration reasoning                                     │
│    9. Assign developmental states                               │
│  ────────────────────────────────────────────────────────────── │
│  CHECKPOINT: User reviews Tier 2 results                        │
│  Output: Per-type UMAP + TF heatmap + Reasoning                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 3: Tier 3 - Functional States (per Dev State)            │
│  ─────────────────────────────────────────────────────────────  │
│  FOR EACH dev_state:                                            │
│    1. Subset data                                               │
│    2. Preprocess (strategy A/B/C)                               │
│    3. Select resolution (fixed/multi)                           │
│    4. RE-COMPUTE DE within subset                               │
│    5. Pathway activity analysis                                 │
│    6. PubMed search + verify                                    │
│    7. 3-iteration reasoning                                     │
│    8. Assign functional states                                  │
│    9. Flag NOVEL populations                                    │
│  ────────────────────────────────────────────────────────────── │
│  CHECKPOINT: User reviews Tier 3 results                        │
│  Output: Per-state UMAP + Pathway heatmap + Novel flags         │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 4: Output Generation                                     │
│  ─────────────────────────────────────────────────────────────  │
│  1. Load visualization.md + storage.md                          │
│  2. Generate final figures:                                     │
│     - Hierarchical UMAP (3 panels)                              │
│     - Summary dotplot with brackets                             │
│     - Sankey diagram (Tier flow)                                │
│  3. Compile reference table (all DOUBLE_VERIFIED PMIDs)         │
│  4. Save final h5ad with all annotations                        │
│  5. Generate run_manifest.json                                  │
│  ────────────────────────────────────────────────────────────── │
│  COMPLETE: All outputs saved                                    │
└─────────────────────────────────────────────────────────────────┘
```

---

## Checkpoint Protocol

At each checkpoint, agent MUST:

```python
def create_checkpoint(phase: int, adata, results: dict):
    """
    Create checkpoint for user review.

    Saves:
    - Current adata state
    - Phase-specific results
    - Visualization for review
    """
    import json
    from pathlib import Path

    output_dir = Path("annotation_output")
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    # Save checkpoint data
    checkpoint = {
        "phase": phase,
        "timestamp": str(datetime.now()),
        "n_cells": adata.n_obs,
        "n_clusters": adata.obs['leiden'].nunique(),
        "results_summary": results.get('summary', {}),
        "status": "awaiting_review"
    }

    checkpoint_file = checkpoint_dir / f"phase{phase}_checkpoint.json"
    with open(checkpoint_file, 'w') as f:
        json.dump(checkpoint, f, indent=2)

    # Save adata snapshot
    adata.write(checkpoint_dir / f"phase{phase}_adata.h5ad")

    return checkpoint_file
```

---

## Resume Protocol

If previous run was interrupted:

```python
def resume_from_checkpoint(checkpoint_dir: Path):
    """
    Resume annotation from last checkpoint.

    Returns:
        tuple: (phase_to_resume, adata, previous_results)
    """
    import json

    # Find latest checkpoint
    checkpoints = sorted(checkpoint_dir.glob("phase*_checkpoint.json"))

    if not checkpoints:
        return 0, None, {}  # Start fresh

    latest = checkpoints[-1]
    with open(latest) as f:
        checkpoint = json.load(f)

    phase = checkpoint['phase']

    # Load adata from checkpoint
    adata_file = checkpoint_dir / f"phase{phase}_adata.h5ad"
    if adata_file.exists():
        import scanpy as sc
        adata = sc.read_h5ad(adata_file)
    else:
        return 0, None, {}

    print(f"📂 Resuming from Phase {phase}")
    print(f"   Cells: {checkpoint['n_cells']}")
    print(f"   Clusters: {checkpoint['n_clusters']}")

    # Resume from NEXT phase
    return phase + 1, adata, checkpoint.get('results_summary', {})
```

---

## Error Handling

```python
def handle_phase_error(phase: int, error: Exception, adata):
    """
    Handle errors during annotation.

    Strategy:
    1. Save current state
    2. Log error details
    3. Create recovery checkpoint
    4. Notify user
    """
    import traceback
    from pathlib import Path

    error_dir = Path("annotation_output/errors")
    error_dir.mkdir(parents=True, exist_ok=True)

    # Save error state
    error_file = error_dir / f"phase{phase}_error.json"
    error_data = {
        "phase": phase,
        "error_type": type(error).__name__,
        "error_message": str(error),
        "traceback": traceback.format_exc(),
        "timestamp": str(datetime.now())
    }

    with open(error_file, 'w') as f:
        json.dump(error_data, f, indent=2)

    # Save recovery adata
    if adata is not None:
        adata.write(error_dir / f"phase{phase}_recovery.h5ad")

    print(f"""
    ╔══════════════════════════════════════════════════════════════╗
    ║  ⚠️  ERROR in Phase {phase}                                    ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  Type: {type(error).__name__:<50} ║
    ║  Message: {str(error)[:48]:<48} ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  Recovery files saved to: annotation_output/errors/          ║
    ║  To resume: Load phase{phase}_recovery.h5ad                   ║
    ╚══════════════════════════════════════════════════════════════╝
    """)
```

---

## MCP Tool Integration

```python
# PubMed search via MCP
async def search_pubmed_mcp(query: str, max_results: int = 5):
    """
    Use pubmed-tools MCP server for literature search.

    Tool: mcp__pubmed-tools__pubmed_search
    """
    # Called automatically by Claude when using MCP
    pass

async def verify_reference_mcp(pmid: str, markers: list, cell_type: str = None):
    """
    Verify PMID supports marker-cell type link.

    Tool: mcp__pubmed-tools__verify_reference
    """
    # Called automatically by Claude when using MCP
    pass

async def fetch_abstract_mcp(pmid: str):
    """
    Get full abstract for detailed verification.

    Tool: mcp__pubmed-tools__fetch_abstract
    """
    # Called automatically by Claude when using MCP
    pass
```

---

## Output Structure

```
annotation_output/
├── checkpoints/
│   ├── phase0_checkpoint.json
│   ├── phase0_adata.h5ad
│   ├── phase1_checkpoint.json
│   ├── phase1_adata.h5ad
│   ├── phase2_checkpoint.json
│   ├── phase2_adata.h5ad
│   └── phase3_checkpoint.json
├── figures/
│   ├── tier1_umap.png
│   ├── tier1_dotplot.png
│   ├── tier2_{major_type}_umap.png
│   ├── tier2_{major_type}_tf_heatmap.png
│   ├── tier3_{dev_state}_umap.png
│   ├── tier3_{dev_state}_pathway_heatmap.png
│   ├── final_hierarchical_umap.png
│   └── final_summary_dotplot.png
├── data/
│   ├── tier1_de_results.csv
│   ├── tier2_{major_type}_de_results.csv
│   ├── tier3_{dev_state}_de_results.csv
│   ├── tier2_{major_type}_tf_activity.csv
│   ├── tier3_{dev_state}_pathway_activity.csv
│   └── final_annotated.h5ad
├── references/
│   ├── tier1_references.json
│   ├── tier2_references.json
│   ├── tier3_references.json
│   └── all_references.csv
├── reasoning/
│   ├── tier1_reasoning.md
│   ├── tier2_reasoning.md
│   └── tier3_reasoning.md
├── debug/                          # If ANNOTATION_DEBUG=1
│   ├── phase1_de_log.json
│   ├── phase2_de_log.json
│   └── ...
├── errors/                         # If errors occurred
│   └── phase{N}_error.json
└── run_manifest.json               # Final completion marker
```

---

## Usage

```
# In Claude Code, from /home/kwy7605/data_61/Vaccine_V2/ directory:

User: "자동으로 annotation 해줘"
      OR
      "run hierarchical annotation on my_data.h5ad automatically"

# Agent will:
# 1. Ask for configuration (Phase 0)
# 2. Execute Tier 1 → Checkpoint → User review
# 3. Execute Tier 2 → Checkpoint → User review
# 4. Execute Tier 3 → Checkpoint → User review
# 5. Generate outputs (Phase 4)
# 6. Create run_manifest.json
```

---

## Related Skills

- `Annotation-agent`: Main skill with phase workflows
- `decoupler`: TF and pathway activity
- `palantir`: Trajectory inference
- `scanpy`: DE analysis, clustering
- `scientific-visualization`: Figure generation
