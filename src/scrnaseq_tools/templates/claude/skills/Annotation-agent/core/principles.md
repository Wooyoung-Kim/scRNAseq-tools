# Core Principles (v3)

> Critical rules for ALL Tiers. Code in `tier1/2/3_template.py`.

---

## 1. Marker Criteria (모든 Tier 공통)

| Criterion | Threshold |
|-----------|-----------|
| `pct_in` | >= 25% |
| `log2FC` | >= 1 |
| `padj` | < 0.05 |
| **Top N** | **50** (after filtering) |
| **조합** | **2-5개 마커 (단일 마커 금지)** |

---

## 2. DE Re-computation (HARD FAIL)

각 Tier에서 **반드시** subset 내에서 DE 재계산. 이전 Tier DE 사용 금지.

| Tier | Subset | subset_id |
|------|--------|-----------|
| 1 | Full dataset | `full_dataset` |
| 2 | Per major type | `{major_type}` |
| 3 | Per dev state | `{major_type}_{dev_state}` |

---

## 3. Functional Analysis (MANDATORY)

| Tier | 필수 분석 | Key |
|------|----------|-----|
| Tier 2 | TF Activity (decoupler ULM) | `obsm['ulm_estimate']` |
| Tier 2 | Trajectory (>= 2000 cells) | `obs['pseudotime']` |
| Tier 3 | Pathway Activity (decoupler MLM) | `obsm['mlm_estimate']` |

**기능 분석 없이 annotation 진행 → HARD FAIL**

---

## 4. PMID Requirement (Per-Marker 필수)

- 모든 문헌 참조에 PMID 필수
- PMID 없는 참조 = INVALID → REJECTED
- 검증: MCP `pubmed_search` → `verify_reference` → DOUBLE_VERIFIED
- **Per-Marker 검증 필수**: `pmid_verified`는 각 marker gene별 dict

```
❌ FORBIDDEN: "pmid_verified": [{"pmid": "..."}]        ← flat list
✅ REQUIRED:  "pmid_verified": {"Gene1": [...], "Gene2": [...]}  ← per-marker dict
```

각 선택된 marker gene에 대해 개별적으로 `search_marker_reference()` 호출 → 결과를 marker-keyed dict로 저장

---

## 5. NO Hardcoding (절대 규칙)

| ❌ FORBIDDEN | ✅ REQUIRED |
|-------------|------------|
| `t_markers = ['CD3D', 'CD3E']` | DE에서 추출 |
| `TF_ASSOCIATIONS = {'TBX21': 'Th1'}` | z-score + 문헌 검색 |
| `PATHWAY_ASSOCIATIONS = {'WNT': 'Stemness'}` | 통계적 outlier detection |
| 단일 마커 annotation | 2-5 마커 조합 |

---

## 6. Statistical Outlier Detection

- TF/Pathway activity z-score > 2.5 → outlier 플래그
- Outlier → 전문 문헌 검색 우선

| Outlier Count | Action |
|---------------|--------|
| 0 | Standard workflow |
| 1-2 | Note + investigate |
| >= 3 | **HIGH PRIORITY** — specialized subset likely |

---

## 7. Confidence Scoring (12 pts max)

| Criteria | 3 pts | 2 pts | 1 pt | 0 pts |
|----------|-------|-------|------|-------|
| Markers | >= 3 | 2 | 1 | 0 |
| References | >= 2 DOUBLE_VERIFIED | 1 DOUBLE | 1 VERIFIED | None |
| TF/Pathway consistency | >= 2 support | 1 supports | Weak | N/A |
| Trajectory/Literature | Aligned | Partial | Conflicting | N/A |

### Outlier Confidence Rules

```
RULE 1: Outlier + standard annotation → max 6-7 pts (PENALTY)
RULE 2: Outlier + literature explains → 12 pts (FULL)
RULE 3: >= 3 outliers unexplained → max 6 pts (SEVERE)
```

| Score | Level |
|-------|-------|
| >= 10 | HIGH ✅ |
| 7-9 | MEDIUM ⚠️ |
| 4-6 | LOW ⚠️ |
| < 4 | INSUFFICIENT ❌ |

---

## 8. NEVER DO

1. ❌ Skip functional analysis → annotate
2. ❌ Use upstream DE markers (이전 Tier 마커)
3. ❌ Single marker annotation
4. ❌ Ignore TF/Pathway conflicts
5. ❌ Skip PMID verification
6. ❌ Hardcode canonical markers
7. ❌ Ignore anti-enrichment (enrichment < 1.0)
8. ❌ Skip contamination check
9. ❌ Hardcode TF/Pathway associations
10. ❌ Ignore statistical outliers
11. ❌ Standard annotation with unexplained outliers (→ penalty)

---

## 9. Validation Checkpoints

### Before Tier 2 Annotation
- [ ] DE computed within subset
- [ ] TF activity computed (`ulm_estimate`)
- [ ] Trajectory if >= 2000 cells
- [ ] Outliers detected (`tf_outliers` in uns)

### Before Tier 3 Annotation
- [ ] DE RE-computed within subset
- [ ] Pathway activity (`mlm_estimate`)
- [ ] Outliers detected (`pathway_outliers` in uns)

### Before Final Assignment
- [ ] All evidence considered (markers + TF/pathway + outliers)
- [ ] IF outliers: specialized literature search done
- [ ] IF unexplained: confidence penalty applied
- [ ] >= 2 references DOUBLE_VERIFIED
- [ ] Novel populations flagged if applicable

---

## 10. Output Requirements

### Tier 2 (.h5ad must contain)
- `obs['tier2_annotation']`
- `obsm['ulm_estimate']`, `obsm['ulm_pvals']`
- `obs['pseudotime']` (if >= 2000 cells)
- `uns['tf_outliers']`, `uns['evidence_conflicts']`

### Tier 3 (.h5ad must contain)
- `obs['tier3_annotation']`, `obs['final_annotation']`
- `obsm['mlm_estimate']`, `obsm['mlm_pvals']`
- `uns['pathway_outliers']`, `uns['evidence_conflicts']`
