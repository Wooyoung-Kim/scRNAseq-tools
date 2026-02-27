"""
Storage Template
================
Code reference for tools/storage.md rules.
Read storage.md for output structure/schema FIRST.
"""
import json
import hashlib
import os
import glob
import scanpy as sc
from datetime import datetime

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
            "marker_thresholds": {
                "tier1_cluster_de": {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                                     "note": "Tier1 초기 클러스터 reasoning용 (Steps 3-6.5)"},
                "tier1_dotplot":    {"pct_min": 0.40, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                                     "hk_removal": True, "dotplot_highest": True,
                                     "pmid_required": True, "final_n": "3-5",
                                     "note": "Tier1 annotation DE 기반 pool → PMID 확인된 gene만 DotPlot"},
                "tier2_cluster_de": {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                                     "note": "Tier2 클러스터 reasoning용 (Steps 7-10)"},
                "tier2_dotplot":    {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                                     "hk_removal": True, "dotplot_highest": True,
                                     "pmid_required": True, "final_n": "3-5",
                                     "note": "Tier2 annotation DE 기반 pool → PMID 확인된 gene만 DotPlot (pct 0.25 유지)"},
                "tier3_group_de":   {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                                     "note": "Tier3 tier3_group reasoning용 (Steps 7-10)"},
                "tier3_dotplot":    {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                                     "hk_removal": True, "dotplot_highest": True,
                                     "pmid_required": True, "final_n": "3-5",
                                     "note": "Tier3 annotation DE 기반 pool → PMID 확인된 gene만 DotPlot (pct 0.25 유지)"}
            }
        },

        "outputs": all_outputs,

        "summary": {
            "n_tier1_types": len(adata_original.obs['tier1_annotation'].unique()),
            "n_tier2_states": len(adata_original.obs['tier2_annotation'].unique())
                              if 'tier2_annotation' in adata_original.obs.columns else 0,
            "n_tier3_states": len(adata_original.obs['tier3_annotation'].unique())
                              if 'tier3_annotation' in adata_original.obs.columns else 0,
            "n_novel_populations": sum(1 for x in adata_original.obs.get('tier3_annotation',
                                       []).unique() if 'Novel' in str(x)),
            "evidence_entries": len(all_evidence),
            "skipped": config.get('skipped', {})   # {major_type: 'n_clusters=N < 4'}
        },

        "de_tracking": config.get('de_tracking', {})
    }

    manifest_path = f'{output_dir}/run_manifest.json'
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"✅ Manifest saved: {manifest_path}")
    return manifest


def validate_evidence_schema(all_evidence):
    """annotation_evidence 리스트의 per-marker 스키마를 검증한다.

    각 entry의 markers[]에 gene / pct_in / log2fc / pmid / title / reasoning
    필드가 모두 있는지 확인. reasoning 누락 시 경고, pmid 누락 시 경고.

    Args:
        all_evidence: list of annotation_evidence dicts

    Returns:
        bool: True if valid, raises ValueError if critical fields missing
    """
    required_fields = {'gene', 'pct_in', 'log2fc', 'pmid', 'title', 'reasoning'}
    errors, warnings = [], []

    for entry in all_evidence:
        annotation = entry.get('annotation', '?')
        markers = entry.get('markers', [])

        if not markers:
            warnings.append(f"[{annotation}] markers[] is empty")
            continue

        for m in markers:
            gene = m.get('gene', '?')
            missing = required_fields - set(m.keys())
            if missing:
                errors.append(f"[{annotation}] {gene}: missing fields {missing}")
                continue
            if not m.get('pmid'):
                warnings.append(f"[{annotation}] {gene}: pmid is empty")
            if not m.get('reasoning'):
                warnings.append(f"[{annotation}] {gene}: reasoning is empty")

    for w in warnings:
        print(f"   ⚠️ {w}")
    if errors:
        raise ValueError(f"❌ Evidence schema errors ({len(errors)}):\n" +
                         "\n".join(f"  - {e}" for e in errors))

    print(f"   ✅ Evidence schema valid: {len(all_evidence)} entries, "
          f"{sum(len(e.get('markers',[])) for e in all_evidence)} markers")
    return True


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

    # Validate + save evidence
    validate_evidence_schema(all_evidence)
    with open(f'{output_dir}/references/annotation_evidence.json', 'w') as f:
        json.dump(all_evidence, f, indent=2, ensure_ascii=False)

    # Save final
    adata_original.write(f'{output_dir}/final_annotated.h5ad')
    print(f"Final annotated data saved: {output_dir}/final_annotated.h5ad")

    return adata_original, all_evidence
