# Data Storage Rules

> 📄 Code reference: `storage_template.py` (same directory)

---

## MANDATORY: run_manifest.json

Every completed annotation run MUST generate `run_manifest.json`.

```json
{
  "version": "v3",
  "timestamp": "ISO8601",
  "status": "completed",
  "input": { "file_path": "...", "file_hash": "sha256:...", "n_cells": N, "n_genes": N },
  "config": {
    "marker_thresholds": {
      "tier1_cluster_de": {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                           "note": "Tier1 초기 클러스터 reasoning용 (Steps 3-6.5)"},
      "tier1_dotplot":    {"pct_min": 0.40, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                           "hk_removal": true, "dotplot_highest": true,
                           "pmid_required": true, "final_n": "3-5",
                           "note": "Tier1 annotation DE 기반 pool → PMID 확인된 gene만 DotPlot"},
      "tier2_cluster_de": {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                           "note": "Tier2 클러스터 reasoning용 (Steps 7-10)"},
      "tier2_dotplot":    {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                           "hk_removal": true, "dotplot_highest": true,
                           "pmid_required": true, "final_n": "3-5",
                           "note": "Tier2 annotation DE 기반 pool → PMID 확인된 gene만 DotPlot (pct 0.25 유지)"},
      "tier3_group_de":   {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                           "note": "Tier3 tier3_group reasoning용 (Steps 7-10)"},
      "tier3_dotplot":    {"pct_min": 0.25, "lfc_min": 1.0, "padj_max": 0.05, "top_n": 50,
                           "hk_removal": true, "dotplot_highest": true,
                           "pmid_required": true, "final_n": "3-5",
                           "note": "Tier3 annotation DE 기반 pool → PMID 확인된 gene만 DotPlot (pct 0.25 유지)"}
    }
  },
  "outputs": { "report_md": "...", "evidence_json": "..." },
  "summary": { "n_tier1_types": N, "n_tier2_states": N, "n_tier3_states": N,
               "skipped": {"B_cells": "n_clusters=2 < 4"} },
  "de_tracking": { "tier1": {...}, "tier2": [...], "tier3": [...] }
}
```

### Marker Threshold 요약

| 단계 | 용도 | pct | HK 제거 | Dotplot-Highest | PMID 필수 |
|------|------|-----|---------|----------------|----------|
| Tier1 클러스터 DE | Reasoning | ≥ 0.25 | — | — | — |
| Tier1 DotPlot | 시각화 | ≥ 0.40 | ✅ | ✅ | ✅ |
| Tier2 클러스터 DE | Reasoning | ≥ 0.25 | — | — | — |
| Tier2 DotPlot | 시각화 | ≥ 0.25 | ✅ | ✅ | ✅ |
| Tier3 group DE | Reasoning | ≥ 0.25 | — | — | — |
| Tier3 DotPlot | 시각화 | ≥ 0.25 | ✅ | ✅ | ✅ |

---

## Output Structure

```
annotation_output/
├── subsets/
│   ├── tier1_full.h5ad
│   ├── tier2/
│   │   ├── T_cells.h5ad
│   │   └── B_cells.h5ad
│   └── tier3/
│       ├── T_cells_Naive.h5ad
│       └── T_cells_Effector.h5ad
├── figures/
│   ├── tier1/
│   │   ├── tier1_umap.png/.svg
│   │   └── tier1_annotation_marker_dotplot.png/.svg
│   ├── tier2/
│   │   ├── T_cells_umap.png/.svg
│   │   ├── T_cells_annotation_marker_dotplot.png/.svg
│   │   ├── T_cells_tf_heatmap.png/.svg
│   │   ├── B_cells_umap.png/.svg
│   │   ├── B_cells_annotation_marker_dotplot.png/.svg
│   │   └── B_cells_tf_heatmap.png/.svg
│   └── tier3/
│       ├── T_cells_Naive_umap.png/.svg
│       ├── T_cells_Naive_annotation_marker_dotplot.png/.svg
│       ├── T_cells_Naive_pathway_heatmap.png/.svg
│       ├── T_cells_Effector_umap.png/.svg
│       ├── T_cells_Effector_annotation_marker_dotplot.png/.svg
│       └── T_cells_Effector_pathway_heatmap.png/.svg
├── references/
│   ├── annotation_evidence.json
│   └── tier1_de.csv
└── run_manifest.json  ← MANDATORY
```

---

## Key Functions (in storage_template.py)

| Function | Purpose |
|----------|---------|
| `generate_manifest()` | Create run_manifest.json (marker_thresholds 포함) |
| `validate_manifest()` | Verify all referenced files exist |
| `validate_evidence_schema()` | per-marker pmid/title/reasoning 필드 검증 |
| `save_tier_subsets()` | Save per-tier .h5ad subsets |
| `write_final_annotations()` | Write annotations back to original adata |
| `propagate_from_subsets()` | Merge subset annotations into full dataset |
| `complete_storage_pipeline()` | Full pipeline: annotate → evidence validate → save |

---

## annotation_evidence.json Schema

`markers[]`는 per-marker 구조 — gene별 pmid + title + reasoning 필수.

```json
{
  "annotation": "Naive_T",
  "n_cells": 1234,
  "confidence_level": "HIGH",
  "markers": [
    {
      "gene":      "Ccr7",
      "pct_in":    0.78,
      "log2fc":    2.1,
      "pmid":      "12345678",
      "title":     "CCR7 controls naive T cell homing...",
      "reasoning": "Naive T cell homing receptor, Dotplot-highest 통과."
    }
  ],
  "tf_activity": [...],
  "references":  [...]
}
```

```
❌ FORBIDDEN: markers[]에 pmid/title/reasoning 없이 gene + pct만 저장
✅ REQUIRED:  선정된 3-5개 마커 각각에 pmid + title + reasoning 포함
```

---

## Report Format (PRIMARY deliverable)

```markdown
# Report Title
## Marker Genes, Functions, and Literature Evidence

**Species**: ... | **Date**: YYYY-MM-DD | **PMID Verification**: N/N markers VERIFIED (per-marker)

## 1. original -> Verified_Name | N cells | STATUS

### Marker Evidence
| Gene | pct | LFC | PMID (hyperlinked) | Title | Reasoning |

### Key TFs: TF1 (score), ...
```

---

## Checklist

```
- [ ] run_manifest.json generated
- [ ] All tier subsets saved (.h5ad)
- [ ] annotation_evidence.json saved
- [ ] report.md saved (PRIMARY)
- [ ] All figures saved (UMAP + Annotation Marker DotPlot)
- [ ] Final annotations written to original adata
- [ ] validate_manifest() PASSED
```
