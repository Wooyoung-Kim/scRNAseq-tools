# Phase 1: Tier 1 - Major Cell Types (v3 - DATA-DRIVEN)

> ⚠️ NO hardcoded lineage marker lists. Uses dynamic knowledge base.
> 📄 Code reference: `tier1_template.py` (same directory)

---

## Goal

| 이 Tier에서 결정 | 이 Tier에서 결정 안함 |
|------------------|----------------------|
| 세포 계통 (T/B/Myeloid/NK) | 발달 단계 (Tier 2) |
| Lineage identity | 기능 상태 (Tier 3) |

---

## Input / Output

- **Input**: `phase0_adata.h5ad` — QC 완료, `adata.X` = log-normalized, `adata.raw.X` = raw counts (`X_pca` or `X_pca_harmony` 필수)
- **Output**: `adata.obs['tier1_annotation']` — Major cell type labels

---

## Marker Criteria

| 단계 | 목적 | pct 기준 | 기타 |
|------|------|----------|------|
| Step 3 (초기 클러스터 추론) | LLM reasoning용 후보 마커 | pct ≥ 25%, LFC ≥ 1, padj < 0.05 | Top 50 pool |
| Step 7b-7d (DotPlot 후보 pool 구성) | 코드: annotation DE에서 strict filter → pool | pct ≥ 40%, LFC ≥ 1, padj < 0.05 (annotation DE 기준) | HK 제거 → Dotplot-Highest → per-type Top 50 |
| Step 7e (최종 마커 선정) | Agent: PubMed 검색 → PMID 확인된 gene만 DotPlot | pool 내 gene 각각 PubMed 검색 | **PMID 없는 gene은 DotPlot 제외, 3-5개 최종** |

---

## Workflow (9 Steps)

```
Step 1:   Clustering (Leiden, resolution 0.5-1.0)
Step 2:   DE 계산 (rank_genes_groups, wilcoxon)
Step 2.5: Statistical Outlier Detection (z-scores, data-driven)
Step 2.6: Cross-Lineage Contamination Detection (data-driven)
Step 3:   Filter Valid Markers (Top 50)
Step 4:   Identify Cell Types (Dynamic Knowledge — NO hardcoding)
Step 5:   LLM Pre-Analysis (outlier-aware)
Step 6:   MCP Verification (PMID Required)
Step 6.5: 3-Iteration Reasoning
Step 7:   Assign Annotations
Step 7a:  Annotation 기반 DE 재계산 (groupby=tier1_annotation, 16 types)
Step 7b:  Strict Filter — annotation DE에서 pct≥40%, LFC≥1, padj<0.05 → per-type Top 50 pool
Step 7c:  Housekeeping 유전자 제거 (Rps/Rpl/mt-/Gm*/기타)
Step 7d:  Dotplot-Highest 검증 (argmax pct == target cell type)
Step 7e:  [Agent] pool에서 PubMed 검색 → PMID 확인된 gene만 DotPlot 후보
          PMID 없는 gene은 DotPlot 제외 → 3-5개 최종 선정 (PMID + reasoning 필수)
Step 7.5: Create annotation_evidence
Step 8:   Cross-Contamination Verification (post-annotation)
Step 9:   Save Results + Visualizations
```

---

## Key Rules

### Data-Driven Cell Type Discovery (NO Hardcoding)

```
❌ FORBIDDEN:
   t_markers = ['CD3D', 'CD3E', 'TRAC']  # 하드코딩
   b_markers = ['CD79A', 'CD79B']         # 하드코딩

✅ REQUIRED:
   1. Extract top DE markers from current data
   2. Search literature via dynamic_knowledge.py
   3. Get cell type candidates from PubMed
   4. Rank by confidence → Assign
   5. No match → Flag as NOVEL
```

> 📄 See `tools/dynamic_knowledge.md` and `tools/dynamic_knowledge.py` for implementation

### Statistical Outlier Detection

- Detect marker expression outliers: z-score > 2.5 and expression > 0.1
- Detect cross-lineage contamination: top markers in other clusters > 10%
- Results stored in `adata.uns['marker_outliers']` and `adata.uns['contamination_issues']`

### Pre-Analysis Template (Step 5)

```markdown
=== Cluster {X} Analysis (Tier 1) ===

## Top 50 DE Markers
| Rank | Gene | pct_in | log2FC | padj | Z-score | Outlier? |

## STATISTICAL OUTLIERS (if detected)
## CONTAMINATION ALERTS (if detected)

## REASONING INSTRUCTION
- Standard: Identify lineage → search literature → assign
- If outliers: Priority investigation of marker combinations
- If contamination: Check for doublets
```

---

## Per-Marker Evidence Schema (Step 6 / Step 7e Required)

`annotation_evidence['markers']`는 **각 marker gene별** PMID + Reasoning이 포함된 리스트여야 합니다.

```json
"markers": [
  {
    "gene":      "Cd79a",
    "pct_in":    0.92,
    "log2fc":    3.4,
    "pmid":      "12345678",
    "title":     "CD79A expression defines B lineage identity...",
    "reasoning": "Cd79a는 B cell receptor complex의 신호 전달 서브유닛으로, B lineage 전반에 걸쳐 발현되며 T/Myeloid 계통과 명확히 구분됨."
  },
  {
    "gene":      "Ebf1",
    "pct_in":    0.81,
    "log2fc":    2.9,
    "pmid":      "23456789",
    "title":     "EBF1 controls B lymphocyte fate...",
    "reasoning": "Ebf1은 B lineage commitment의 핵심 TF로, dotplot-highest 조건도 통과한 특이적 마커."
  }
]
```

```
❌ FORBIDDEN:
   markers에 pmid/title/reasoning 없이 gene + pct만 저장
   reasoning을 cell type 레벨로 통합 작성 (마커별 개별 서술 필수)

✅ REQUIRED:
   선정된 3-5개 마커 각각에 pmid + title + reasoning 작성
   reasoning = "왜 이 유전자가 이 세포 타입의 대표 마커인가" (1-2문장)
```

검증 절차 (Step 7e):
1. `build_marker_pool()` 결과에서 3-5개 선정
2. 각 gene에 대해 `pubmed_search(gene, cell_type)` → 최상위 PMID + title 저장
3. 각 gene에 대해 reasoning 1-2문장 작성
4. `build_per_marker_evidence(markers, cell_type, reasonings={gene: reason})` 호출

---

## Post-Annotation Verification

- Cross-type contamination must be < 5%
- Each annotation needs 2-4 lineage markers (from DE, NOT hardcoded)
- All annotations verified with PMID (per-marker)

---

## Checklist

```
ANALYSIS:
- [ ] Leiden clustering (resolution 0.5-1.0)
- [ ] DE computed
- [ ] Statistical outliers detected
- [ ] Cross-lineage contamination detected
- [ ] Valid markers filtered (Top 50)

ANNOTATION:
- [ ] Each cluster: 2-5 lineage markers (from DE)
- [ ] Dynamic knowledge search performed (NO hardcoding)
- [ ] MCP PMID verification (per-marker)
- [ ] 3-iteration reasoning
- [ ] Cross-contamination < 5%
- [ ] Annotation 기반 DE 재계산 (groupby=tier1_annotation, use_raw=True, pts=True)

MARKER SELECTION & DOTPLOT:
- [ ] Step 7b: build_marker_pool(ann_de_df, adata) → pct≥0.40, LFC≥1, padj<0.05, Top 50
- [ ] Step 7c: Housekeeping 유전자 제거 (자동 적용됨)
- [ ] Step 7d: Dotplot-Highest 검증 (argmax 세포 타입 일치, 자동 적용됨)
- [ ] Step 7e: [Agent] pool에서 각 gene PubMed 검색 → PMID 없는 gene 제외
              PMID 확인된 gene 중 3-5개 선정 (reasoning 필수)
- [ ] build_marker_dict_from_selections(agent_selections) → {annotation: [markers]} dict
- [ ] DotPlot with brackets (var_names=dict, use_raw=False, standard_scale='var')

OUTPUT:
- [ ] tier1_annotated.h5ad saved
- [ ] Per-type subsets saved (tier1/{type}.h5ad)
- [ ] UMAP saved
- [ ] Annotation Marker DotPlot saved (WITH brackets)
- [ ] annotation_evidence.json saved (markers[].pmid / .title / .reasoning 포함)
- [ ] report.md saved — Marker Evidence 표: Gene | pct | LFC | PMID | Title | Reasoning
```

---

## DotPlot 설정

```python
sc.pl.DotPlot(
    var_names      = marker_dict,   # {cell_type: [genes]} → 브라켓 자동 생성
    groupby        = 'tier1_annotation',
    use_raw        = False,         # adata.X (log-normalized) 사용
    standard_scale = 'var',         # 유전자별 0~1 정규화 → 컬러 균일
).style(largest_dot=60, dot_max=1.0)

# figsize = (total_genes × 0.28 + 1, n_types × 0.38 + 1.5)
# dpi     = 100
```

---

## 주의사항

```
❌ 하지 말 것:
- 발달 상태 포함 (예: "Naive T cells" → Tier 2에서)
- 기능 상태 포함 (예: "Cytotoxic T cells" → Tier 3에서)
- 단일 마커로 annotation (예: "CD3D+ cells")
- Hardcoded lineage marker lists 사용

✅ 해야 할 것:
- 계통만 결정 (예: "T cells", "B lineage")
- 2-5개 마커 조합 사용 (from current dataset DE)
- Outlier 기반 우선순위 설정
- 문헌 검색 (dynamic_knowledge.py)
```
