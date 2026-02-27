# Phase 2: Tier 2 - Developmental States (v2)

> ⚠️ TF Activity + Trajectory are **MANDATORY**. Cannot annotate without functional analysis.
> 📄 Code reference: `tier2_template.py` (same directory)

---

## Goal

| 이 Tier에서 결정 | 이 Tier에서 결정 안함 |
|------------------|----------------------|
| 분화 단계 (Naive/Effector/Memory) | 기능 상태 (Tier 3) |
| Developmental trajectory position | Lineage identity (Tier 1) |

---

## Input / Output

- **Input**: Tier 1 subset (한 lineage만, ex: T cells)
- **Output**:
  - ① Subset UMAP
  - ② Annotation Marker DotPlot (상단 브라켓)
  - ③ Reasoning + 논문 PMID/제목

---

## Marker Criteria

| 단계 | 목적 | pct 기준 | 기타 |
|------|------|----------|------|
| Step 7 (클러스터 reasoning) | LLM reasoning용 후보 마커 | pct ≥ 25%, LFC ≥ 1, padj < 0.05 | Top 50 pool |
| Step 11b-11d (DotPlot 후보 pool 구성) | 코드: annotation DE에서 strict filter → pool | pct ≥ 25%, LFC ≥ 1, padj < 0.05 (annotation DE 기준) | HK 제거 → Dotplot-Highest → per-type Top 50 |
| Step 11.5 (최종 마커 선정) | Agent: PubMed 검색 → PMID 확인된 gene만 DotPlot | pool 내 gene 각각 PubMed 검색 | **PMID 없는 gene은 DotPlot 제외, 3-5개 최종** |

> ※ Tier 2는 lineage subset 내 발달 상태 구분 → Tier 1의 0.40과 달리 **pct 0.25 유지**
> (subset이 작고 발달 마커의 발현율이 상대적으로 낮음)

---

## Skip Condition

```
Re-cluster 결과 클러스터 수 < 4  →  Tier 2 / Tier 3 전체 생략
  - tier2_annotation = tier1_annotation 그대로 유지
  - tier3_annotation 생략
  - skip 사유를 run_manifest.json에 기록
```

> 클러스터가 4개 미만이면 발달 상태 구분이 불가능하거나 무의미함.
> `recluster()` 직후 `check_cluster_count()` 호출.

---

## Workflow (13 Steps)

```
Step 1:   Subset Data (from tier1_annotation)
Step 2:   Re-cluster (resolution 0.1-0.5, use harmony embedding)
          ⚠️ 클러스터 수 < 4 → SKIP (Tier 2/3 생략, manifest에 기록)
Step 3:   RE-COMPUTE DE (⚠️ MUST re-compute within subset)
Step 4:   TF Activity (MANDATORY — decoupler CollecTRI/ULM)
Step 5:   Trajectory (MANDATORY if >= 2000 cells — Palantir)
Step 5.5: Statistical Outlier Detection (TF z-scores)
Step 5.6: Evidence Conflict Detection
Step 6:   VERIFY Functional Analysis
Step 7:   Filter Valid Markers (Top 50 after pct/LFC/padj filter)
Step 8:   Integrated Pre-Analysis (markers + TF + trajectory + outliers)
Step 9:   MCP Verification (PMID Required, per-marker)
Step 10:  3-Iteration Reasoning (→ reasoning/integrated_format.md)
Step 11:  Assign Annotations
Step 11a: Annotation 기반 DE 재계산 (groupby=tier2_annotation)
Step 11b: Strict Filter — annotation DE에서 pct≥25%, LFC≥1, padj<0.05 → per-type Top 50 pool
Step 11c: Housekeeping 유전자 제거 (Rps/Rpl/mt-/Gm*/기타)
Step 11d: Dotplot-Highest 검증 (argmax pct == target annotation)
Step 11.5: [Agent] pool에서 PubMed 검색 → PMID 확인된 gene만 DotPlot 후보
           PMID 없는 gene은 DotPlot 제외 → 3-5개 최종 선정 (PMID + reasoning 필수)
Step 12:  Save Results + Visualizations (UMAP + Bracketed DotPlot + Report)
```

> ⚠️ **Step 11a-11.5은 필수**: `compute_de_by_annotation()` → `build_marker_pool()` → Agent PMID 검색 → `build_marker_dict_from_selections()` → `save_results(..., marker_dict)`.

> 📄 각 Step의 코드: `tier2_template.py` 참조

---

## MANDATORY Analysis

| Analysis | Tool | Key | Required |
|----------|------|-----|----------|
| TF Activity | decoupler CollecTRI/ULM | `obsm['ulm_estimate']` | ✅ Always |
| TF p-values | decoupler | `obsm['ulm_pvals']` | ✅ Always |
| Trajectory | Palantir pseudotime | `obs['pseudotime']` | ✅ if >= 2000 cells |

---

## Evidence Integration

Tier 2 reasoning **MUST** combine 3 evidence types:

1. **DE Markers** — Top 50 valid markers per cluster
2. **TF Activity** — Top 10 TFs + z-score outliers
3. **Trajectory** — Pseudotime category (Early/Mid/Late)

### Outlier Handling

| Outlier Count | Interpretation | Action |
|---------------|---------------|--------|
| 0 | Standard developmental state | Normal workflow |
| 1-2 | Minor variation | Check if markers also unusual |
| >= 3 | **HIGH PRIORITY** | Literature search for specialized subset |

### Confidence Scoring (12 pts max)

| Criteria | 3 pts (High) | 2 pts (Medium) | 1 pt (Low) | 0 pts |
|----------|-------------|----------------|------------|-------|
| Markers | >= 3 meeting criteria | 2 markers | 1 marker | 0 |
| References | >= 2 DOUBLE_VERIFIED | 1 DOUBLE or 2 VERIFIED | 1 VERIFIED | None |
| TF consistency | >= 2 TFs support | 1 TF supports | Weak/conflicting | N/A |
| Trajectory | Matches category | Borderline | Conflicts | N/A |

---

## Per-Marker Evidence Schema (Step 9 Required)

`annotation_evidence['markers']`에 각 marker별 **pmid + title + reasoning** 포함 필수.

```json
"markers": [
  {"gene": "Ccr7", "pct_in": 0.78, "log2fc": 2.1,
   "pmid": "12345678", "title": "CCR7 controls naive T cell homing...",
   "reasoning": "Naive T cell homing receptor, Dotplot-highest 통과."},
  {"gene": "Tcf7", "pct_in": 0.71, "log2fc": 1.9,
   "pmid": "23456789", "title": "TCF7 maintains naive T cell identity...",
   "reasoning": "Stemness TF, Naive/Memory 상태 구분의 핵심 마커."}
]
```

```
❌ FORBIDDEN: reasoning 없이 gene + pct만 저장
✅ REQUIRED: 선정된 3-5개 마커 각각에 pmid + title + reasoning (1-2문장)
```

→ `build_per_marker_evidence(markers, cell_type, reasonings={gene: reason})` 호출

---

## Pre-Analysis Template (Step 8)

```markdown
=== {Major Type} Cluster {X} Analysis ===

## DE Markers (Top 50 Valid)
| Rank | Gene | pct_in | log2FC | padj |
[show all 50]

## TF Activity (Top 10 + ALL outliers)
| TF | Activity | p-value | Z-score | Outlier? |

## Trajectory Position
- Mean pseudotime: {value}
- Category: {Early/Mid/Late}

## REASONING INSTRUCTION
- IF NO outliers: Standard developmental state reasoning
- IF outliers: Priority literature search → specialized subset candidates
```

---

## Checklist (Per Major Type)

```
MANDATORY:
- [ ] Subset extracted
- [ ] Re-clustered → check_cluster_count() 호출
- [ ] 클러스터 수 >= 4 확인 (미만이면 SKIP → manifest 기록)
- [ ] DE RE-COMPUTED within subset
- [ ] TF activity computed (ulm_estimate in obsm)
- [ ] Trajectory computed if >= 2000 cells
- [ ] Outliers detected (tf_outliers in uns)
- [ ] Functional analysis VERIFIED

ANNOTATION:
- [ ] Valid markers filtered (Top 50)
- [ ] Integrated Pre-Analysis (markers + TF + trajectory)
- [ ] MCP verification with PMID (per-marker)
- [ ] 3-iteration reasoning
- [ ] References double-verified
- [ ] Confidence scored (max 12 pts)

MARKER SELECTION & DOTPLOT:
- [ ] Step 11a: compute_de_by_annotation(subset) → annotation 단위 DE
- [ ] Step 11b-11d: build_marker_pool(ann_de_df, subset) → pct≥0.25, HK 제거, Dotplot-highest
- [ ] Step 11.5: [Agent] pool에서 PubMed 검색 → PMID 확인된 gene만 DotPlot 후보
               PMID 없는 gene 제외 → 3-5개 선정 (reasoning 필수)
- [ ] build_marker_dict_from_selections(agent_selections) → {annotation: [markers]} dict
- [ ] DotPlot with brackets (var_names=dict, use_raw=False, standard_scale='var')

OUTPUT:
- [ ] Annotations assigned
- [ ] Subset saved (.h5ad)
- [ ] UMAP saved
- [ ] Annotation Marker DotPlot saved (WITH brackets)
- [ ] TF heatmap saved
- [ ] annotation_evidence.json saved (markers[].pmid / .title / .reasoning 포함)
- [ ] report.md saved — Marker Evidence 표: Gene | pct | LFC | PMID | Title | Reasoning
```

---

## Troubleshooting

- **TF fails**: Check gene symbol overlap with CollecTRI network
- **Trajectory fails**: Need >= 2000 cells, check PCA exists
- **Memory issues**: Subsample to 20K for trajectory, impute back
