---
name: Annotation-agent
description: "3-tier hierarchical annotation with MANDATORY functional analysis (TF/Pathway/Trajectory). All analyses integrated into reasoning."
---

# Annotation Agent v2 (Hierarchical Cell Type Annotation)

**Version 2.0**: Functional analysis (TF, Pathway, Trajectory) is MANDATORY and integrated into reasoning.

---

## TIER GOALS (핵심 목표)

```
╔══════════════════════════════════════════════════════════════════════╗
║  Tier 1: WHAT lineage?   → "이 세포가 어떤 계통인가?"                 ║
║          예: T cells, B lineage, Myeloid, NK cells                   ║
║          도구: DE markers only                                       ║
╠══════════════════════════════════════════════════════════════════════╣
║  Tier 2: WHEN in development? → "분화 과정에서 어느 단계인가?"        ║
║          예: Naive, Effector, Memory, Pro_B, Plasma                  ║
║          도구: DE + TF Activity + Trajectory                         ║
╠══════════════════════════════════════════════════════════════════════╣
║  Tier 3: WHAT is it doing? → "현재 무엇을 하고 있는가?"               ║
║          예: Cytotoxic, Exhausted, TRM, Light_Zone, M1_like          ║
║          도구: DE + Pathway Activity                                 ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Quick Reference

| Tier | 질문 | 분류 예시 | 분석 도구 |
|------|------|-----------|-----------|
| **1** | What lineage? | T cells, B cells | DE |
| **2** | When in development? | Naive, Effector | DE + TF + Trajectory |
| **3** | What is it doing? | Cytotoxic, Exhausted | DE + Pathway |

**자세한 목표 정의**: `core/tier_goals.md` 참조

---

## PREREQUISITE: QC-pipeline

```
╔══════════════════════════════════════════════════════════════════════╗
║  ⚠️ QC-pipeline MUST be completed before running Annotation-agent    ║
╚══════════════════════════════════════════════════════════════════════╝
```

**Required QC Output:**

```python
# Before starting annotation, verify:
assert 'X_pca_harmony' in adata.obsm or 'X_pca' in adata.obsm, "PCA not computed!"
assert 'X_umap' in adata.obsm, "UMAP not computed!"
assert 'doublet_score' in adata.obs.columns, "Doublet detection not run!"
assert adata.obs['pct_counts_mt'].max() < 25, "High mito cells present!"
assert 'highly_variable' in adata.var.columns, "HVGs not selected!"
assert adata.raw is not None, "Raw counts not preserved!"
```

**Pipeline Flow:**

```
Raw Data → [QC-pipeline] → qc_complete.h5ad → [Annotation-agent] → annotated.h5ad
                ↓
     ├── Cell filtering
     ├── Gene filtering
     ├── Doublet detection
     ├── Normalization + HVG
     └── Batch integration (Harmony)
```

**Reference:** `../QC-pipeline/SKILL.md`

---

## ABSOLUTE RULES (NEVER VIOLATE)

```
╔══════════════════════════════════════════════════════════════════════╗
║  0. PHASE 0 CONFIGURATION IS MANDATORY (HARD FAIL)                   ║
║     - MUST ask user for ALL configuration choices BEFORE any tier     ║
║     - MUST use AskUserQuestion tool to present options interactively  ║
║     - MUST store config in adata.uns['annotation_config']            ║
║     - If annotation_config NOT in adata.uns → HARD FAIL              ║
║     - NEVER hardcode tools/methods/markers without user selection     ║
║                                                                      ║
║     Required selections (phases/user_config.md):                     ║
║       A. TF Activity: tool / database / method                       ║
║       B. Pathway Activity: tool / database / method                  ║
║       C. Trajectory: tool / root cell strategy                       ║
║       D. Re-clustering: strategy (A/B/C) per tier                    ║
║       E. Resolution: fixed / multi-resolution scan                   ║
║                                                                      ║
║     User MUST confirm all settings before proceeding to Tier 1.      ║
╠══════════════════════════════════════════════════════════════════════╣
║  1. FUNCTIONAL ANALYSIS IS MANDATORY (NEW in v2)                     ║
║     - Tier 2: TF Activity + Trajectory (MUST compute)                ║
║     - Tier 3: Pathway Activity (MUST compute)                        ║
║     - Results MUST be in obsm before annotation proceeds             ║
╠══════════════════════════════════════════════════════════════════════╣
║  2. DE RE-COMPUTATION ENFORCEMENT (HARD FAIL)                        ║
║     - Each tier MUST have unique subset_id in DE results             ║
║     - If marker source matches previous tier's DE table → HARD FAIL  ║
╠══════════════════════════════════════════════════════════════════════╣
║  3. PMID VERIFICATION MANDATORY                                      ║
║     - ALL literature citations MUST include PMID                     ║
║     - References WITHOUT PMID are INVALID                            ║
╠══════════════════════════════════════════════════════════════════════╣
║  4. INTEGRATED REASONING                                             ║
║     - Markers + TF/Pathway + Trajectory = Combined evidence          ║
║     - Cannot annotate without functional analysis results            ║
╠══════════════════════════════════════════════════════════════════════╣
║  5. DYNAMIC KNOWLEDGE BASE (NEW - MANDATORY)                         ║
║     - NO hardcoded cell type lists or marker sets                    ║
║     - Cell type discovery via literature search (dynamic_knowledge)  ║
║     - tier_goals.md is REFERENCE ONLY (NOT constraints)              ║
║     - Novel population detection is EXPECTED and ENCOURAGED          ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Critical Rules (ALWAYS IN CONTEXT)

```
0. PHASE 0:         User MUST select all config options BEFORE Tier 1 (HARD FAIL)
1. MARKERS:         pct >= 25%, LFC >= 1, padj < 0.05 (Top 50)
2. DE:              RE-COMPUTED at EACH tier subset
3. TF:              MANDATORY at Tier 2 (tool/db/method from Phase 0 config)
4. PATHWAY:         MANDATORY at Tier 3 (tool/db/method from Phase 0 config)
5. TRAJECTORY:      MANDATORY at Tier 2 if >= 2000 cells (tool from Phase 0 config)
6. PMID:            ALL literature citations MUST include PMID
7. COMBINATION:     2-5 markers per cell type (NEVER single marker)
8. NO HARDCODE:     NEVER hardcode tool/method/markers — use Phase 0 config values
9. DYNAMIC SEARCH:  ALWAYS use dynamic_knowledge.py for cell type discovery
10. NO CONSTRAINTS: tier_goals.md is REFERENCE ONLY, NOT a fixed classification system
```

### Phase 0 Config Validation (MUST pass before ANY tier)

```python
def validate_phase0_config(adata):
    """
    HARD FAIL if Phase 0 configuration is missing.
    MUST call this before starting ANY tier.
    """
    if 'annotation_config' not in adata.uns:
        raise AssertionError(
            "\n" + "=" * 60 + "\n"
            "HARD FAIL: Phase 0 Configuration NOT found!\n"
            "=" * 60 + "\n"
            "adata.uns['annotation_config'] is MISSING.\n\n"
            "You MUST complete Phase 0 FIRST:\n"
            "  1. Ask user: TF Activity (tool/database/method)\n"
            "  2. Ask user: Pathway Activity (tool/database/method)\n"
            "  3. Ask user: Trajectory (tool/root cell)\n"
            "  4. Ask user: Re-clustering strategy (A/B/C)\n"
            "  5. Ask user: Resolution strategy (fixed/multi)\n"
            "  6. Store in adata.uns['annotation_config']\n"
            "  7. User confirms -> THEN proceed to Tier 1\n"
            "\nCANNOT PROCEED WITHOUT USER CONFIGURATION!\n"
            "=" * 60
        )

    config = adata.uns['annotation_config']
    required_keys = ['tf_activity', 'pathway_activity', 'trajectory', 'reclustering']
    missing = [k for k in required_keys if k not in config]
    if missing:
        raise AssertionError(
            f"HARD FAIL: annotation_config missing keys: {missing}\n"
            "Re-run Phase 0 to collect all user selections."
        )

    print("Phase 0 configuration validated")
    return config
```

---

## Functional Analysis Requirements

### Tier 2: TF Activity (MANDATORY)

```python
# MUST execute before annotation
import decoupler as dc

net = dc.op.collectri(organism='human')

# decoupler v2 API — use dc.mt.ulm(), NOT dc.run_ulm()
dc.mt.ulm(subset, net, verbose=True, raw=False)

# decoupler v2 stores as 'score_ulm'/'padj_ulm'
# Create legacy aliases for backward compatibility
if 'score_ulm' in subset.obsm:
    subset.obsm['ulm_estimate'] = subset.obsm['score_ulm']
    subset.obsm['ulm_pvals'] = subset.obsm['padj_ulm']

# VERIFY: Must exist before proceeding (check both v2 and legacy keys)
assert 'score_ulm' in subset.obsm or 'ulm_estimate' in subset.obsm, "❌ TF activity not computed!"
```

### Tier 2: Trajectory (MANDATORY if >= 2000 cells)

```python
import palantir

if subset.n_obs >= 2000:
    # Run diffusion maps + multiscale space
    palantir.utils.run_diffusion_maps(subset, n_components=10, pca_key='X_pca')
    palantir.utils.determine_multiscale_space(subset)

    # Run Palantir (results stored directly in AnnData)
    start_cell = find_earliest_cell(subset)  # User-defined or auto
    palantir.core.run_palantir(subset, start_cell, num_waypoints=500)

    subset.obs['pseudotime'] = subset.obs['palantir_pseudotime']

    # VERIFY
    assert 'palantir_pseudotime' in subset.obs.columns, "❌ Pseudotime not computed!"
```

### Tier 3: Pathway Activity (MANDATORY)

```python
net = dc.op.progeny(organism='human', top=500)

# decoupler v2 API — use dc.mt.mlm(), NOT dc.run_mlm()
dc.mt.mlm(subset, net, verbose=True, raw=False)

# decoupler v2 stores as 'score_mlm'/'padj_mlm'
if 'score_mlm' in subset.obsm:
    subset.obsm['mlm_estimate'] = subset.obsm['score_mlm']
    subset.obsm['mlm_pvals'] = subset.obsm['padj_mlm']

# VERIFY (check both v2 and legacy keys)
assert 'score_mlm' in subset.obsm or 'mlm_estimate' in subset.obsm, "❌ Pathway activity not computed!"
```

---

## Phase Flow

| Phase | REQUIRED Analysis | Marker Criteria | Output |
|-------|-------------------|-----------------|--------|
| **0** | User config | — | Config stored |
| **1** | DE only | pct ≥ 25%, LFC ≥ 1, padj < 0.05, Top 50, **2-5개 마커 조합 (단일 마커 금지)** | Tier 1 labels |
| **2** | DE + TF + Trajectory | pct ≥ 25%, LFC ≥ 1, padj < 0.05, Top 50, **2-5개 마커 조합 (단일 마커 금지)** | ① Subset UMAP ② Annotation Marker DotPlot (상단 브라켓) ③ Reasoning + 논문 PMID/제목 |
| **3** | DE + Pathway | pct ≥ 25%, LFC ≥ 1, padj < 0.05, Top 50, **2-5개 마커 조합 (단일 마커 금지)** | ① Subset UMAP ② Annotation Marker DotPlot (상단 브라켓) ③ 기능 상태 Reasoning + PMID/제목 |
| **4** | Visualization + Report | — | **report.md** (PRIMARY) |

### Tier 3 입력 구성: Tier 2 Annotation + Cluster 결합

Tier 3의 입력은 Tier 2 annotation 결과와 cluster 번호를 **결합**하여 생성합니다:

```python
# Tier 2 결과: cluster 1, 2 → Memory_B / cluster 3 → Naive_B 등
# Tier 3 입력: Memory_B_1, Memory_B_2, Naive_B_3 형태로 결합

subset.obs['tier3_group'] = (
    subset.obs['tier2_annotation'].astype(str) + '_' +
    subset.obs['tier2_cluster'].astype(str)
)

# DE는 이 결합 레이블로 실행
sc.tl.rank_genes_groups(subset, groupby='tier3_group', method='wilcoxon')
de = sc.get.rank_genes_groups_df(subset, group=None)
de = de[(de['pcts'] >= 0.25) & (de['logfoldchanges'] >= 1) & (de['pvals_adj'] < 0.05)]
de = de.groupby('group').head(50)

# 목적: Memory_B_1 vs Memory_B_2의 기능적 차이를 파악
# → Memory_B_1 = Resting Memory? Memory_B_2 = Activated Memory?
```

> **핵심**: 같은 cell type 내에서도 cluster별 기능적 차이를 DE + Pathway로 밝히는 것이 Tier 3의 목적

---

## Report Format (PRIMARY Output)

Generated by `save_markdown_report()` in `tools/visualization.md`. This is the PRIMARY deliverable of Phase 4.

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
- **Machine-readable companion**: `annotation_evidence.json`

### Integrated Evidence Collection (Tier 2)

During annotation (before report), collect evidence per cluster:

```markdown
**Evidence for Cluster X**:

## DE Markers (Top 50)
| Gene | LFC | pct_in | padj |
|------|-----|--------|------|
| ... |

## TF Activity (Top 5)
| TF | Activity Score | p-value | Known Function |
|----|----------------|---------|----------------|
| ... |

## Pseudotime Position
- Mean: 0.35
- Category: EARLY (0-0.3: Early, 0.3-0.7: Mid, 0.7-1.0: Late)

## Combined Evidence
- Markers suggest: [candidate1, candidate2]
- TF activity supports: [candidate1] (TCF7 active → Naive/Memory)
- Trajectory position: EARLY → Supports progenitor/naive state
```

### Integrated Evidence Collection (Tier 3)

```markdown
**Evidence for Cluster X**:

## DE Markers (Top 50)
| Gene | LFC | pct_in | padj |
|------|-----|--------|------|
| ... |

## Pathway Activity (Top 5)
| Pathway | Activity Score | p-value | Interpretation |
|---------|----------------|---------|----------------|
| TNFa | 2.3 | 0.001 | Pro-inflammatory |
| Hypoxia | 1.8 | 0.01 | Tissue adaptation |
| ... |

## Combined Evidence
- Markers suggest: [func_state1, func_state2]
- Pathway activity supports: [func_state1] (TNFa + NFkB high → Activated/Effector)
```

---

## Confidence Scoring (Updated)

| Criteria | High (3 pts) | Medium (2 pts) | Low (1 pt) | None (0 pts) |
|----------|--------------|----------------|------------|--------------|
| **Markers** | >= 3 meeting criteria | 2 markers | 1 marker | 0 |
| **References** | >= 2 DOUBLE_VERIFIED | 1 DOUBLE or 2 VERIFIED | 1 VERIFIED | None |
| **TF/Pathway** | Strong match (>= 2 support) | Partial (1 supports) | Weak/conflicting | Not computed |
| **Trajectory** | Matches category | Borderline | Conflicts | Not computed |

**Total >= 10**: HIGH | **7-9**: MEDIUM | **4-6**: LOW | **< 4**: INSUFFICIENT

---

## Output Verification Checklist

Before completing each tier:

```
TIER 2 CHECKLIST:
- [ ] DE computed (rank_genes_groups)
- [ ] TF activity computed (score_ulm / ulm_estimate in obsm)
- [ ] Trajectory computed if >= 2000 cells (pseudotime in obs)
- [ ] Reasoning includes TF + trajectory evidence
- [ ] PMID verified for all references

TIER 3 CHECKLIST:
- [ ] DE RE-computed within subset
- [ ] Pathway activity computed (score_mlm / mlm_estimate in obsm)
- [ ] Reasoning includes pathway evidence
- [ ] PMID verified for all references

PHASE 4 CHECKLIST (Visualization + Report):
- [ ] UMAP figures saved (adaptive dot size, legend_loc='on data')
- [ ] Annotation Marker Dotplot figures saved (DotPlot OOP API, dot_max=1.0, largest_dot=60)
- [ ] {prefix}_report.md saved (PRIMARY deliverable — structured markdown)
- [ ] annotation_evidence.json saved (machine-readable)
- [ ] run_manifest.json saved and validated
```

---

## Quick Start

```python
# Tier 2 workflow
subset = adata[adata.obs['tier1_annotation'] == major_type].copy()

# 1. Re-cluster
sc.tl.leiden(subset, resolution=0.8, key_added='tier2_cluster')

# 2. DE (MANDATORY) — pct >= 25%, LFC >= 1, padj < 0.05, Top 50
sc.tl.rank_genes_groups(subset, groupby='tier2_cluster', method='wilcoxon')
de = sc.get.rank_genes_groups_df(subset, group=None)
de = de[(de['pcts'] >= 0.25) & (de['logfoldchanges'] >= 1) & (de['pvals_adj'] < 0.05)]
de = de.groupby('group').head(50)

# 3. TF Activity (MANDATORY) — decoupler v2 API
import decoupler as dc
net = dc.op.collectri(organism='human')
dc.mt.ulm(subset, net, verbose=True, raw=False)
# Create legacy aliases
if 'score_ulm' in subset.obsm:
    subset.obsm['ulm_estimate'] = subset.obsm['score_ulm']
    subset.obsm['ulm_pvals'] = subset.obsm['padj_ulm']

# 4. Trajectory (MANDATORY if >= 2000 cells)
if subset.n_obs >= 2000:
    import palantir
    palantir.utils.run_diffusion_maps(subset, pca_key='X_pca')
    palantir.utils.determine_multiscale_space(subset)
    start_cell = find_earliest_cell(subset, naive_markers)  # Data-driven, not hardcoded
    palantir.core.run_palantir(subset, start_cell)
    subset.obs['pseudotime'] = subset.obs['palantir_pseudotime']

# 5. VERIFY before annotation (check both v2 and legacy keys)
assert 'score_ulm' in subset.obsm or 'ulm_estimate' in subset.obsm
if subset.n_obs >= 2000:
    assert 'pseudotime' in subset.obs.columns

# 6. Now proceed to annotation with integrated reasoning
```

---

## Related Files

| File | Purpose | When to Load |
|------|---------|--------------|
| **`../QC-pipeline/SKILL.md`** | **QC prerequisite** | **Before annotation** |
| `../QC-pipeline/connection/to_annotation.md` | QC → Annotation handoff | Before Tier 1 |
| **`core/tier_goals.md`** | **Tier 목표 정의 (REFERENCE ONLY, NOT constraints)** | **최초 읽기 (예시 참고용)** |
| **`tools/dynamic_knowledge.py`** | **데이터 드라이븐 세포 타입 발견 (NO hardcoding)** | **Tier 1, 2, 3 (annotation)** |
| **`tools/dynamic_knowledge.md`** | **Dynamic knowledge 사용 가이드** | **Tier 1, 2, 3 (annotation)** |
| `core/principles.md` | Critical rules | Every tier |
| `core/functional_analysis.md` | TF/Pathway/Trajectory code | Tier 2, 3 |
| **`phases/user_config.md`** | **Phase 0 사용자 설정 (TF/Pathway/Trajectory/Re-clustering)** | **Before Tier 1 (Phase 0)** |
| `phases/tier1.md` | Major cell type workflow | Phase 1 |
| `phases/tier2.md` | Developmental state workflow | Phase 2 |
| `phases/tier3.md` | Functional state workflow | Phase 3 |
| **`reasoning/agent_format.md`** | **3-iteration 추론 템플릿 + evidence.json 스키마 매핑** | **All tiers (reasoning)** |
| `reasoning/integrated_format.md` | Combined evidence reasoning | All tiers |
| `tools/functional_plots.md` | TF/Pathway/Trajectory 시각화 | Phase 4 |
| **`tools/visualization.md`** | **DotPlot OOP API (dot_max=1.0, largest_dot=60), adaptive UMAP sizing, save_markdown_report()** | **Phase 4** |
| **`tools/literature_verification.md`** | **PubMed API 설정, PMID 검증 워크플로우** | **All tiers (literature)** |
| **`tools/storage.md`** | **run_manifest.json 스키마, 저장 파이프라인, 완료 검증** | **Every tier (save) + Final** |
| `decision/matrices.md` | Confidence & Novel criteria | When deciding |
| `examples/workflow_example.md` | 실제 작동 예제 | Reference |
| `examples/compact_examples.md` | Compact 예제 (Novel/Low Quality 포함) | Reference |
