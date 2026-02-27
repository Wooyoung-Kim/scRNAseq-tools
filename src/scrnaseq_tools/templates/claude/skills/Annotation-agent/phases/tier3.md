# Phase 3: Tier 3 - Functional States (v2)

> ⚠️ Pathway Activity is **MANDATORY**. Cannot annotate without pathway analysis.
> 📄 Code reference: `tier3_template.py` (same directory)

---

## Goal

| 이 Tier에서 결정 | 이 Tier에서 결정 안함 |
|------------------|----------------------|
| 기능 상태 (Cytotoxic/Exhausted/TRM) | Lineage (Tier 1) |
| 같은 cell type 내 기능적 차이 | 발달 단계 (Tier 2) |

---

## Tier 3 Input: Annotation + Cluster 결합

Tier 3의 입력은 Tier 2 annotation + cluster 번호를 **결합**하여 생성:

```
Tier 2 결과: cluster 1,2 → Memory_B / cluster 3 → Naive_B
         ↓ 결합
Tier 3 입력: Memory_B_1, Memory_B_2, Naive_B_3
         ↓ DE + Pathway
목적: Memory_B_1 vs Memory_B_2의 기능적 차이 파악
```

> 📄 `tier3_template.py`의 `create_tier3_groups()` 함수 참조

---

## Input / Output

- **Input**: Tier 2 subset + annotation + cluster 결합 레이블
- **Output**:
  - ① Subset UMAP
  - ② Annotation Marker DotPlot (상단 브라켓)
  - ③ 기능 상태 Reasoning + PMID/제목

---

## Marker Criteria

| 단계 | 목적 | pct 기준 | 기타 |
|------|------|----------|------|
| Step 7 (tier3_group reasoning) | LLM reasoning용 후보 마커 | pct ≥ 25%, LFC ≥ 1, padj < 0.05 | Top 50 pool |
| Step 12b-12d (DotPlot 후보 pool 구성) | 코드: annotation DE에서 strict filter → pool | pct ≥ 25%, LFC ≥ 1, padj < 0.05 (annotation DE 기준) | HK 제거 → Dotplot-Highest → per-type Top 50 |
| Step 12.5 (최종 마커 선정) | Agent: PubMed 검색 → PMID 확인된 gene만 DotPlot | pool 내 gene 각각 PubMed 검색 | **PMID 없는 gene은 DotPlot 제외, 3-5개 최종** |

> ※ Tier 3은 developmental state 내 기능적 차이 구분 → **pct 0.25 유지**
> (functional marker는 세포 내 발현율이 낮을 수 있으며 subset도 작음)

---

## Workflow (13 Steps)

```
Step 1:   Subset Data (from Tier 2 annotation + cluster)
Step 1.5: Create tier3_groups (annotation + cluster 결합)
Step 2:   RE-COMPUTE DE (on tier3_group = annotation + cluster)
Step 3:   Pathway Activity (MANDATORY — decoupler PROGENy/MLM)
Step 4.5: Statistical Outlier Detection (pathway z-scores)
Step 5.6: Evidence Conflict Detection
Step 6:   VERIFY Functional Analysis
Step 7:   Filter Valid Markers (Top 50)
Step 8:   Integrated Pre-Analysis (markers + pathway + outliers)
Step 9:   MCP Verification (PMID Required, per-marker)
Step 10:  3-Iteration Reasoning (→ reasoning/integrated_format.md)
Step 11:  Detect Novel Populations
Step 12:  Assign Annotations
Step 12a: Annotation 기반 DE 재계산 (groupby=tier3_annotation)
Step 12b: Strict Filter — annotation DE에서 pct≥25%, LFC≥1, padj<0.05 → per-type Top 50 pool
Step 12c: Housekeeping 유전자 제거 (Rps/Rpl/mt-/Gm*/기타)
Step 12d: Dotplot-Highest 검증 (argmax pct == target annotation)
Step 12.5: [Agent] pool에서 PubMed 검색 → PMID 확인된 gene만 DotPlot 후보
           PMID 없는 gene은 DotPlot 제외 → 3-5개 최종 선정 (PMID + reasoning 필수)
Step 13:  Save Results + Visualizations (UMAP + Bracketed DotPlot + Report)
```

> ⚠️ **Re-clustering 필요 없음**: Tier 3은 tier3_group (annotation+cluster 결합)을 입력으로 사용. `create_tier3_groups()` 함수 참조.
> ⚠️ **Step 12a-12.5은 필수**: `compute_de_by_annotation()` → `build_marker_pool()` → Agent PMID 검색 → `build_marker_dict_from_selections()` → `save_results(..., marker_dict)`.

---

## MANDATORY Analysis

| Analysis | Tool | Key | Required |
|----------|------|-----|----------|
| Pathway Activity | decoupler PROGENy/MLM | `obsm['mlm_estimate']` | ✅ Always |
| Pathway p-values | decoupler | `obsm['mlm_pvals']` | ✅ Always |

---

## Pathway Interpretation

| Pathway | High → | Low → | Outlier → |
|---------|--------|-------|-----------|
| **TNFa** | Pro-inflammatory | Quiescent | Chronic inflammation |
| **NFkB** | Activation/survival | Suppressed | Dysfunction |
| **JAK-STAT** | Cytokine response | Unresponsive | Chronic stimulation |
| **TGFb** | Immunosuppression | Pro-inflammatory | Regulatory dysfunction |
| **WNT** | Stemness | Terminal diff. | Precursor/plastic |
| **Hypoxia** | Tissue adaptation | Normoxic | TRM / tumor |
| **Trail** | Cytotoxic | Non-cytotoxic | Atypical killing |

### Outlier Pathway Combinations (Specialized States)

| Combination | Interpretation | Keywords |
|-------------|---------------|----------|
| WNT high in Effector | Tpex / stem-like | "TCF7 effector", "Tpex" |
| TNFa low in Cytotoxic | Exhausted | "exhausted cytotoxic" |
| TGFb + GZMB high | Suppressed TIL | "TGFb cytotoxic" |
| Hypoxia outlier | TRM / tumor-infiltrating | "TRM", "tumor microenvironment" |

---

## Per-Marker Evidence Schema (Step 9 Required)

`annotation_evidence['markers']`에 각 marker별 **pmid + title + reasoning** 포함 필수.

```json
"markers": [
  {"gene": "Gzmb", "pct_in": 0.82, "log2fc": 3.1,
   "pmid": "12345678", "title": "Granzyme B mediates cytotoxic killing...",
   "reasoning": "Cytotoxic effector 핵심 마커, Trail pathway 활성과 일치."},
  {"gene": "Prf1", "pct_in": 0.74, "log2fc": 2.8,
   "pmid": "23456789", "title": "Perforin pore formation in target cells...",
   "reasoning": "Gzmb와 함께 세포독성 기능 직접 반영, Dotplot-highest 통과."}
]
```

```
❌ FORBIDDEN: reasoning 없이 gene + pct만 저장
✅ REQUIRED: 선정된 3-5개 마커 각각에 pmid + title + reasoning (1-2문장)
```

→ `build_per_marker_evidence(markers, cell_type, reasonings={gene: reason})` 호출

---

## Novel Population Detection

Flag as **NOVEL** if:
- No literature match after 3+ queries (including outlier-specific)
- Opposing markers present (e.g., GZMB + TCF7)
- Statistical outliers >= 3 with no literature explanation
- HIGH severity conflicts unresolved

Naming: `{ParentContext}_{Feature}` (e.g., `T_Effector_WNT_high`)

---

## Confidence Scoring (12 pts max)

| Criteria | 3 pts | 2 pts | 1 pt | 0 pts |
|----------|-------|-------|------|-------|
| Markers | >= 3 | 2 | 1 | 0 |
| References | >= 2 DOUBLE_VERIFIED | 1 DOUBLE | 1 VERIFIED | None |
| Pathway consistency | >= 2 pathways support | 1 supports | Weak | N/A |
| Literature/Conflict | Resolved | Partial | Unresolved | N/A |

---

## Checklist (Per Developmental State)

```
MANDATORY:
- [ ] Subset from Tier 2 (annotation + cluster combined)
- [ ] tier3_groups created (NO re-clustering)
- [ ] DE RE-COMPUTED on tier3_group
- [ ] Pathway activity computed (mlm_estimate)
- [ ] Outliers detected (pathway_outliers in uns)
- [ ] Functional analysis VERIFIED

ANNOTATION:
- [ ] Valid markers filtered (Top 50)
- [ ] Integrated Pre-Analysis (markers + pathways + outliers)
- [ ] MCP verification with PMID (per-marker)
- [ ] 3-iteration reasoning
- [ ] Novel populations checked
- [ ] Confidence scored (max 12 pts)

MARKER SELECTION & DOTPLOT:
- [ ] Step 12a: compute_de_by_annotation(subset) → annotation 단위 DE
- [ ] Step 12b-12d: build_marker_pool(ann_de_df, subset) → pct≥0.25, HK 제거, Dotplot-highest
- [ ] Step 12.5: [Agent] pool에서 PubMed 검색 → PMID 확인된 gene만 DotPlot 후보
               PMID 없는 gene 제외 → 3-5개 선정 (reasoning 필수)
- [ ] build_marker_dict_from_selections(agent_selections) → {annotation: [markers]} dict
- [ ] DotPlot with brackets (var_names=dict, use_raw=False, standard_scale='var')

OUTPUT:
- [ ] Tier 3 annotations assigned
- [ ] Final hierarchical labels: {Major}_{Dev}_{Func}
- [ ] Subset saved (.h5ad)
- [ ] UMAP saved
- [ ] Annotation Marker DotPlot saved (WITH brackets)
- [ ] Pathway heatmap saved (MANDATORY)
- [ ] annotation_evidence.json saved (markers[].pmid / .title / .reasoning 포함)
- [ ] report.md saved — Marker Evidence 표: Gene | pct | LFC | PMID | Title | Reasoning
```

---

## Output Format

```
{Major Type}_{Dev State}_{Functional State}

Examples:
- T cells_Effector_Cytotoxic
- T cells_Effector_Exhausted
- B cells_GC_Light_Zone
- T cells_Effector_Novel_GZMB+TCF7
```
