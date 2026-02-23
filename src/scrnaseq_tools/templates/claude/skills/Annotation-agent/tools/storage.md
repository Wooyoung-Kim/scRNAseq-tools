# Data Storage Code

Data persistence strategy for hierarchical annotation.

---

## MANDATORY: run_manifest.json (Completion Contract)

```
╔══════════════════════════════════════════════════════════════════════╗
║  Task is INCOMPLETE until run_manifest.json exists and validates     ║
║  Console output "completed" means NOTHING without valid manifest     ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Manifest Schema

```json
{
  "version": "v3",
  "timestamp": "2024-01-15T10:30:00Z",
  "status": "completed",

  "input": {
    "file_path": "/path/to/input.h5ad",
    "file_hash": "sha256:abc123...",
    "n_cells": 50000,
    "n_genes": 20000
  },

  "config": {
    "tf_activity": {"tool": "decoupler", "database": "dorothea", "method": "ulm"},
    "pathway_activity": {"tool": "decoupler", "database": "progeny", "method": "mlm"},
    "trajectory": {"tool": "palantir", "use": true},
    "marker_thresholds": {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05}
  },

  "outputs": {
    "columns": {
      "tier1": "tier1_annotation",
      "tier2": "tier2_annotation",
      "tier3": "tier3_annotation",
      "final": "final_annotation",
      "confidence": "annotation_confidence"
    },
    "subsets": {
      "tier1": "subsets/tier1_full.h5ad",
      "tier2": ["subsets/tier2/T_cells.h5ad", "..."],
      "tier3": ["subsets/tier3/T_cells_Effector.h5ad", "..."]
    },
    "figures": [
      "figures/tier1_umap.png",
      "figures/tier1_dotplot.png",
      "..."
    ],
    "final_h5ad": "final_annotated.h5ad",
    "evidence_json": "references/annotation_evidence.json"
  },

  "summary": {
    "n_tier1_types": 6,
    "n_tier2_states": 18,
    "n_tier3_states": 45,
    "n_novel_populations": 2,
    "evidence_entries": 127
  },

  "de_tracking": {
    "tier1": {"subset_id": "full_dataset", "de_table": "tier1_de.csv"},
    "tier2": [
      {"subset_id": "T_cells", "de_table": "tier2_T_cells_de.csv"},
      "..."
    ],
    "tier3": [
      {"subset_id": "T_cells_Effector", "de_table": "tier3_T_cells_Effector_de.csv"},
      "..."
    ]
  }
}
```

### Manifest Generation Code

```python
import json
import hashlib
from datetime import datetime

def generate_manifest(adata_original, config, output_dir, all_outputs, all_evidence):
    """Generate run_manifest.json - REQUIRED for completion."""

    # Calculate input hash
    import tempfile
    with tempfile.NamedTemporaryFile(suffix='.h5ad') as tmp:
        adata_original.write(tmp.name)
        with open(tmp.name, 'rb') as f:
            file_hash = hashlib.sha256(f.read()).hexdigest()[:16]

    manifest = {
        "version": "v3",
        "timestamp": datetime.now().isoformat(),
        "status": "completed",

        "input": {
            "file_path": str(config.get('input_path', 'unknown')),
            "file_hash": f"sha256:{file_hash}",
            "n_cells": adata_original.n_obs,
            "n_genes": adata_original.n_vars
        },

        "config": {
            "tf_activity": config.get('tf_activity', {}),
            "pathway_activity": config.get('pathway_activity', {}),
            "trajectory": config.get('trajectory', {}),
            "marker_thresholds": {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05}
        },

        "outputs": all_outputs,

        "summary": {
            "n_tier1_types": len(adata_original.obs['tier1_annotation'].unique()),
            "n_tier2_states": len(adata_original.obs['tier2_annotation'].unique()),
            "n_tier3_states": len(adata_original.obs['tier3_annotation'].unique()),
            "n_novel_populations": sum(1 for x in adata_original.obs['tier3_annotation'].unique()
                                       if 'Novel' in str(x)),
            "evidence_entries": len(all_evidence)
        },

        "de_tracking": config.get('de_tracking', {})
    }

    manifest_path = f'{output_dir}/run_manifest.json'
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"✅ Manifest saved: {manifest_path}")
    return manifest


def validate_manifest(output_dir):
    """Validate manifest and all referenced files exist."""
    import os

    manifest_path = f'{output_dir}/run_manifest.json'
    if not os.path.exists(manifest_path):
        raise RuntimeError("❌ INCOMPLETE: run_manifest.json not found")

    with open(manifest_path) as f:
        manifest = json.load(f)

    errors = []

    # Check final h5ad
    final_path = f"{output_dir}/{manifest['outputs']['final_h5ad']}"
    if not os.path.exists(final_path):
        errors.append(f"Missing: {final_path}")

    # Check evidence json
    evidence_path = f"{output_dir}/{manifest['outputs']['evidence_json']}"
    if not os.path.exists(evidence_path):
        errors.append(f"Missing: {evidence_path}")

    # Check figures
    for fig in manifest['outputs'].get('figures', []):
        fig_path = f"{output_dir}/{fig}"
        if not os.path.exists(fig_path):
            errors.append(f"Missing figure: {fig_path}")

    if errors:
        raise RuntimeError(f"❌ INCOMPLETE: {len(errors)} missing files\n" + "\n".join(errors))

    print(f"✅ VALIDATED: All {len(manifest['outputs']['figures'])} figures + data files present")
    return True
```

---

## Output Structure

```
annotation_output/
├── subsets/
│   ├── tier1_full.h5ad              # Full data + Tier1
│   ├── tier2/
│   │   ├── T_cells.h5ad             # Per major type + Tier2
│   │   ├── B_cells.h5ad
│   │   └── ...
│   └── tier3/
│       ├── T_cells_Naive.h5ad       # Per dev state + Tier3
│       ├── T_cells_Effector.h5ad
│       └── ...
├── figures/
│   ├── tier1_*.png/.svg
│   ├── tier2_*.png/.svg
│   └── tier3_*.png/.svg
├── reports/
│   ├── annotation_report.md         # ← PRIMARY: consolidated markdown (format below)
│   └── annotation_evidence.json     #   machine-readable evidence
├── references/
│   └── annotation_evidence.json
└── final_annotated.h5ad             # Original + all annotations
```

---

## Report Format (`annotation_report.md`)

PRIMARY deliverable. Generated by `save_markdown_report()` in `tools/visualization.md`.

```markdown
# Title
## Marker Genes, Functions, and Literature Evidence

**Species**: Ferret (*Mustela putorius furo*)
**Date**: 2026-02-09
**PMID Verification**: 125/125 VERIFIED via NCBI Entrez

---

# I. B Lineage (153,307 cells, 12 original types)

**Source**: `SP-sbcl_ct.Bl.h5ad`
**Clean output**: 134,948 cells (removed 18,359 = 12%)
**Reclassifications**: 3 (NB3->Pre_B, MB1->ABC, MB2->Follicular_B)
**Removals**: 2 clusters (PB2=T/NK contamination, PC2=Low quality)

---

## 1. ImmB -> Immature_B (Pro/Pre-B) | 3,046 cells | CONFIRMED

### Marker Genes and Functions

| Gene | pct | Enrichment | Biological Function |
|------|-----|-----------|-------------------|
| **DNTT** (TdT) | 95.4% | 211.8x | Terminal deoxynucleotidyl transferase. ... |
| **RAG1** | 90.8% | 11.6x | Recombination-activating gene 1. ... |

### Key TFs: PAX5 (5.43), MYB (4.43), EBF1 (2.48), FOXO1 (2.24)

### Literature Evidence

| PMID | First Author | Year | Journal | Title |
|------|-------------|------|---------|-------|
| [35354960](https://pubmed.ncbi.nlm.nih.gov/35354960/) | Klein F | 2022 | *Nature Immunology* | Dntt expression reveals... |

---
```

Structure rules:
- **Header**: title, species, date, `N/N VERIFIED via NCBI Entrez`
- **Lineage section**: `# I. Lineage (N cells, N original types)` + source/clean/reclassifications/removals
- **Per cell type**: `## N. original -> Verified | N cells | STATUS`
- **Marker table**: `| **GENE** (alias) | pct | Enrichment | Biological Function |`
- **Key TFs**: `### Key TFs: TF1 (score), TF2 (score), ...`
- **Literature**: `| [PMID](hyperlink) | First Author | Year | *Journal* | Title |`
- **Machine-readable companion**: `annotation_evidence.json` (same data, JSON format)

---

## 1. Save Tier Subsets

```python
import os
import scanpy as sc

def save_tier_subsets(adata, output_dir='./annotation_output/subsets/'):
    """Save subset AnnData at each tier."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f'{output_dir}/tier2', exist_ok=True)
    os.makedirs(f'{output_dir}/tier3', exist_ok=True)

    # Tier 1 (full data)
    adata.write(f'{output_dir}/tier1_full.h5ad')
    print(f"Saved: tier1_full.h5ad ({adata.n_obs} cells)")

    # Tier 2 (per major type)
    for major_type in adata.obs['tier1_annotation'].unique():
        subset = adata[adata.obs['tier1_annotation'] == major_type].copy()
        safe_name = major_type.replace(' ', '_').replace('/', '_')
        subset.write(f'{output_dir}/tier2/{safe_name}.h5ad')
        print(f"Saved: tier2/{safe_name}.h5ad ({subset.n_obs} cells)")

    # Tier 3 (per dev state)
    if 'tier2_annotation' in adata.obs.columns:
        for major_type in adata.obs['tier1_annotation'].unique():
            major_subset = adata[adata.obs['tier1_annotation'] == major_type]
            for dev_state in major_subset.obs['tier2_annotation'].unique():
                subset = major_subset[major_subset.obs['tier2_annotation'] == dev_state].copy()
                safe_name = f"{major_type}_{dev_state}".replace(' ', '_').replace('/', '_')
                subset.write(f'{output_dir}/tier3/{safe_name}.h5ad')
                print(f"Saved: tier3/{safe_name}.h5ad ({subset.n_obs} cells)")
```

---

## 2. Write Final Annotations to Original

```python
def write_final_annotations(adata_original, annotation_results):
    """Write all tier annotations back to ORIGINAL adata.

    Parameters
    ----------
    adata_original : AnnData
        Original full dataset
    annotation_results : dict
        {cell_barcode: {'tier1': ..., 'tier2': ..., 'tier3': ..., 'final': ..., 'confidence': ...}}

    Adds columns
    ------------
    tier1_annotation, tier2_annotation, tier3_annotation,
    final_annotation, annotation_confidence
    """
    # Initialize
    for col in ['tier1_annotation', 'tier2_annotation', 'tier3_annotation',
                'final_annotation', 'annotation_confidence']:
        adata_original.obs[col] = 'Unknown'

    # Map annotations
    for bc, info in annotation_results.items():
        if bc in adata_original.obs_names:
            adata_original.obs.loc[bc, 'tier1_annotation'] = info.get('tier1', 'Unknown')
            adata_original.obs.loc[bc, 'tier2_annotation'] = info.get('tier2', 'Unknown')
            adata_original.obs.loc[bc, 'tier3_annotation'] = info.get('tier3', 'Unknown')
            adata_original.obs.loc[bc, 'final_annotation'] = info.get('final', 'Unknown')
            adata_original.obs.loc[bc, 'annotation_confidence'] = info.get('confidence', 'Low')

    # Convert to categorical
    for col in ['tier1_annotation', 'tier2_annotation', 'tier3_annotation',
                'final_annotation', 'annotation_confidence']:
        adata_original.obs[col] = adata_original.obs[col].astype('category')

    return adata_original
```

---

## 3. Propagate from Subset Files

```python
import glob

def propagate_from_subsets(adata_original, subset_dir='./annotation_output/subsets/'):
    """Propagate annotations from saved subsets back to original."""
    all_annotations = {}

    # From Tier 1
    tier1 = sc.read_h5ad(f'{subset_dir}/tier1_full.h5ad')
    for bc in tier1.obs_names:
        all_annotations[bc] = {'tier1': tier1.obs.loc[bc, 'tier1_annotation']}

    # From Tier 2
    for f in glob.glob(f'{subset_dir}/tier2/*.h5ad'):
        subset = sc.read_h5ad(f)
        for bc in subset.obs_names:
            if bc in all_annotations:
                all_annotations[bc]['tier2'] = subset.obs.loc[bc, 'tier2_annotation']

    # From Tier 3
    for f in glob.glob(f'{subset_dir}/tier3/*.h5ad'):
        subset = sc.read_h5ad(f)
        for bc in subset.obs_names:
            if bc in all_annotations:
                t1 = all_annotations[bc].get('tier1', 'Unknown')
                t2 = all_annotations[bc].get('tier2', 'Unknown')
                t3 = subset.obs.loc[bc, 'tier3_annotation']
                all_annotations[bc]['tier3'] = t3
                all_annotations[bc]['final'] = f"{t1}_{t2}_{t3}"
                all_annotations[bc]['confidence'] = subset.obs.loc[bc].get(
                    'annotation_confidence', 'Medium')

    return write_final_annotations(adata_original, all_annotations)
```

---

## 4. Complete Storage Pipeline

```python
import json

def complete_storage_pipeline(adata_original, tier1_result, tier2_results,
                              tier3_results, all_evidence,
                              output_dir='./annotation_output/'):
    """Complete pipeline with storage at each step."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f'{output_dir}/subsets/tier2', exist_ok=True)
    os.makedirs(f'{output_dir}/subsets/tier3', exist_ok=True)
    os.makedirs(f'{output_dir}/figures', exist_ok=True)
    os.makedirs(f'{output_dir}/references', exist_ok=True)

    adata_work = adata_original.copy()

    # Tier 1
    adata_work.obs['tier1_annotation'] = tier1_result
    adata_work.write(f'{output_dir}/subsets/tier1_full.h5ad')

    # Tier 2 (per major type)
    tier2_all = {}
    for major_type, annotations in tier2_results.items():
        subset = adata_work[adata_work.obs['tier1_annotation'] == major_type].copy()
        subset.obs['tier2_annotation'] = annotations
        safe_name = major_type.replace(' ', '_')
        subset.write(f'{output_dir}/subsets/tier2/{safe_name}.h5ad')
        for bc in subset.obs_names:
            tier2_all[bc] = subset.obs.loc[bc, 'tier2_annotation']

    adata_work.obs['tier2_annotation'] = adata_work.obs_names.map(tier2_all).fillna('Unknown')

    # Tier 3 (per dev state)
    tier3_all = {}
    for key, annotations in tier3_results.items():
        # key = "{major_type}_{dev_state}"
        parts = key.split('_', 1)
        major_type, dev_state = parts[0], parts[1] if len(parts) > 1 else 'Unknown'

        major_subset = adata_work[adata_work.obs['tier1_annotation'] == major_type]
        subset = major_subset[major_subset.obs['tier2_annotation'] == dev_state].copy()
        subset.obs['tier3_annotation'] = annotations
        safe_name = key.replace(' ', '_')
        subset.write(f'{output_dir}/subsets/tier3/{safe_name}.h5ad')

        for bc in subset.obs_names:
            t3 = subset.obs.loc[bc, 'tier3_annotation']
            tier3_all[bc] = t3

    adata_work.obs['tier3_annotation'] = adata_work.obs_names.map(tier3_all).fillna('Unknown')

    # Final annotation
    adata_work.obs['final_annotation'] = (
        adata_work.obs['tier1_annotation'].astype(str) + '_' +
        adata_work.obs['tier2_annotation'].astype(str) + '_' +
        adata_work.obs['tier3_annotation'].astype(str)
    )

    # Copy to original
    for col in ['tier1_annotation', 'tier2_annotation', 'tier3_annotation',
                'final_annotation', 'annotation_confidence']:
        if col in adata_work.obs.columns:
            adata_original.obs[col] = adata_work.obs[col]

    # Save evidence
    with open(f'{output_dir}/references/annotation_evidence.json', 'w') as f:
        json.dump(all_evidence, f, indent=2, ensure_ascii=False)

    # Save final
    adata_original.write(f'{output_dir}/final_annotated.h5ad')
    print(f"Final annotated data saved: {output_dir}/final_annotated.h5ad")

    return adata_original, all_evidence
```

---

## Completion Checklist

```
╔══════════════════════════════════════════════════════════════════════╗
║  COMPLETION IS DETERMINED BY run_manifest.json ONLY                  ║
║  Do NOT report "completed" until validate_manifest() passes          ║
╚══════════════════════════════════════════════════════════════════════╝

REQUIRED files:
- [ ] run_manifest.json (MANDATORY - defines completion)
- [ ] tier1_full.h5ad
- [ ] tier2/*.h5ad (per major type)
- [ ] tier3/*.h5ad (per dev state)
- [ ] final_annotated.h5ad
- [ ] annotation_evidence.json
- [ ] annotation_report.md (consolidated markdown report)
- [ ] All figures listed in manifest

REQUIRED columns in final_annotated.h5ad:
- [ ] tier1_annotation
- [ ] tier2_annotation
- [ ] tier3_annotation
- [ ] final_annotation
- [ ] annotation_confidence

VALIDATION:
>>> validate_manifest('./annotation_output/')
✅ VALIDATED: All files present
>>> # Only NOW can you report completion
```
