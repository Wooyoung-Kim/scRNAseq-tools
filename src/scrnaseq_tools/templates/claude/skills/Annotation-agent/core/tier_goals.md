# Tier Goals (REFERENCE GUIDE - NOT CONSTRAINTS)

⚠️ **CRITICAL: This is a REFERENCE GUIDE, NOT hardcoded constraints.**

```
╔══════════════════════════════════════════════════════════════════════╗
║  WARNING: DO NOT use this as a fixed classification system.         ║
║                                                                      ║
║  The categories and markers listed below are EXAMPLES ONLY.         ║
║  They represent COMMON cell types in HUMAN IMMUNE datasets.         ║
║                                                                      ║
║  ✓ ALWAYS prioritize:                                               ║
║    1. Data-driven DE results from current dataset                   ║
║    2. Statistical outlier detection (z-score > 2.5)                 ║
║    3. Independent PubMed literature search                          ║
║    4. Novel population detection                                    ║
║    5. Dynamic knowledge base (tools/dynamic_knowledge.py)           ║
║                                                                      ║
║  ✓ Use this guide ONLY to:                                          ║
║    - Understand the tier structure (Tier 1/2/3 goals)              ║
║    - Learn the TYPES of questions each tier asks                   ║
║    - See examples of marker combinations (NOT as ground truth)     ║
║                                                                      ║
║  ✗ DO NOT use this to:                                              ║
║    - Constrain cell type discovery to listed categories            ║
║    - Validate markers without PubMed verification                  ║
║    - Limit annotations to example markers                          ║
║    - Skip dynamic literature search                                ║
╚══════════════════════════════════════════════════════════════════════╝
```

**For data-driven annotation, see: `tools/dynamic_knowledge.md`**

각 Tier의 목표, Input, Output, 성공 기준을 명확히 정의합니다.

---

## Tier Overview

```
╔══════════════════════════════════════════════════════════════════════╗
║  Tier 1: WHAT cell type?     → Major lineage identity               ║
║  Tier 2: WHEN in development? → Developmental/differentiation stage ║
║  Tier 3: WHAT is it doing?   → Functional state/activity            ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Tier 1: Major Cell Types

### 목표 (Goal)

```
"이 세포가 어떤 계통(lineage)에 속하는가?"

면역세포, 상피세포, 기질세포 등 세포의 근본적인 정체성을 결정합니다.
이 단계에서는 발달 상태나 기능 상태를 구분하지 않습니다.
```

### 분류 기준 (EXAMPLE ONLY - NOT CONSTRAINTS)

> ⚠️ **WARNING: These are REFERENCE EXAMPLES from human immune datasets.**
>
> **DO NOT use as hardcoded constraints.**
> Your dataset may have DIFFERENT cell types (ILC, MAIT, γδ T, etc.).
>
> **ALWAYS**:
> 1. Use DE markers from YOUR data
> 2. Search PubMed for marker combinations
> 3. Use dynamic knowledge base (tools/dynamic_knowledge.py)
> 4. Flag novel populations if no match

| 카테고리 (Example) | 예시 | 핵심 마커 (Example, NOT definitive) |
|-------------------|------|-------------------------------------|
| **T cells** | 모든 T 세포 | CD3D, CD3E, CD3G, TRAC |
| **B cells / B lineage** | 모든 B 계통 | CD79A, CD79B, MS4A1, CD19 |
| **NK cells** | NK 세포 | NCAM1, NKG7, GNLY (CD3-) |
| **Myeloid** | 단핵구, 대식세포, DC | CD14, CD68, ITGAM, CD1C |
| **Epithelial** | 상피세포 | EPCAM, KRT |
| **Stromal** | 섬유아세포, 내피세포 | COL1A1, PECAM1 |

**Missing from this list (examples of what you might find in YOUR data):**
- ILC (Innate Lymphoid Cells)
- MAIT cells
- γδ T cells
- Mesothelial cells
- Pericytes
- Tissue-resident memory cells
- Disease-specific populations
- **→ Use dynamic_knowledge.py to discover these!**

### Input

```
- Full dataset (QC 완료)
- Leiden clustering (resolution 0.5-1.0)
- DE markers per cluster
```

### Output

```python
adata.obs['tier1_annotation']  # 예: "T cells", "B lineage", "Myeloid"
```

### 성공 기준

```
✓ 모든 클러스터가 major type으로 annotation됨
✓ 각 major type에 2-4개의 lineage-defining marker 확인
✓ Cross-contamination < 5% (예: T cell cluster에 CD79A+ 세포)
✓ 모든 annotation에 PMID 참조
```

### 실패 케이스

```
✗ 단일 마커로 annotation (예: "CD3D+ = T cells")
✗ 기능 상태를 포함 (예: "Exhausted T cells" - 이건 Tier 3)
✗ 발달 상태를 포함 (예: "Naive T cells" - 이건 Tier 2)
```

---

## Tier 2: Developmental States

### 목표 (Goal)

```
"이 세포가 분화 과정에서 어느 단계에 있는가?"

같은 계통 내에서 세포의 발달/분화 상태를 결정합니다.
세포가 현재 수행하는 기능이 아니라, "어디서 왔고 어디로 가는지"를 정의합니다.
```

### 분류 기준 (EXAMPLE ONLY - NOT CONSTRAINTS)

> ⚠️ **WARNING: These are REFERENCE EXAMPLES from human T/B cell studies.**
>
> **DO NOT use as fixed developmental states.**
> Your data may show:
> - Different intermediate states
> - Tissue-specific differentiation
> - Disease-specific trajectories
> - Alternative developmental paths
>
> **ALWAYS**:
> 1. Use DE + TF activity + Trajectory from YOUR data
> 2. Search PubMed for TF + marker combinations
> 3. Use dynamic knowledge for developmental states
> 4. Respect pseudotime order (data-driven)

**T cells (Example - NOT exhaustive):**
| 발달 상태 (Example) | 의미 | 핵심 마커 (Example) | Pseudotime |
|--------------------|------|---------------------|------------|
| **Naive** | 항원 미경험 | CCR7, TCF7, SELL, LEF1 | Early |
| **Activated** | 최근 활성화 | CD69, IL2RA, CD38 | Early-Mid |
| **Effector** | 분화된 효과 세포 | GZMB, PRF1, NKG7 | Mid-Late |
| **Memory** | 기억 세포 | IL7R, CCR7-, CD44 | Late |

**B cells (Example - NOT exhaustive):**
| 발달 상태 (Example) | 의미 | 핵심 마커 (Example, canonical) | Pseudotime |
|--------------------|------|--------------------------------|------------|
| **Pro_B** | 초기 발달 | DNTT, RAG1, ERG | Early |
| **Pre_B** | Pre-BCR 단계 | SOX4, BACH2, MYB | Early |
| **Transitional_B** | 골수 이탈 | IGHM, SOX5, EBF1 | Early-Mid |
| **Naive_B** | 성숙, 항원 미경험 | BANK1, BTLA, CD74 | Mid |
| **GC_B** | 배중심 반응 | BCL6, AICDA | Mid |
| **Memory_B** | 기억 B | IGHG1, LTB, BLK | Late |
| **Plasma_cells** | 항체 분비 | JCHAIN, XBP1, PRDM1 | Late |

> **⚠️ Dataset-specific pitfalls**: Canonical markers may fail in specific datasets.
> Always verify against DE results before using. Known issues:
> - **VPREB1**: May be inverted (higher in non-Pre_B). Use SOX4, BACH2 instead.
> - **BCL6/AICDA**: Often absent at low resolution (pct < 7%). GC_B may not form a distinct cluster.
> - **CD27**: Can be anti-enriched in Memory_B. Use IGHG1, LTB, BLK instead.
>
> **→ These pitfalls illustrate WHY hardcoding fails. Use dynamic_knowledge.py!**

### Input

```
- Tier 1 subset (예: T cells만)
- Re-computed DE within subset
- TF activity (MANDATORY)
- Trajectory/Pseudotime (MANDATORY if >= 2000 cells)
```

### Output

```python
adata.obs['tier2_annotation']  # 예: "Naive", "Effector", "Memory"
adata.obsm['ulm_estimate']     # TF activity scores
adata.obs['pseudotime']        # Differentiation trajectory
```

### 성공 기준

```
✓ Tier 1 마커가 아닌 NEW DE 마커 사용
✓ TF activity가 developmental state와 일치
  - Naive: TCF7+, LEF1+
  - Effector: TBX21+, EOMES+
  - Memory: TCF7+/-, BCL6-
✓ Trajectory position이 developmental order와 일치
  - Naive → Early pseudotime
  - Effector/Memory → Late pseudotime
✓ 모든 annotation에 PMID 참조
```

### 실패 케이스

```
✗ Tier 1 마커 사용 (예: "CD3D high" - 모든 T cell에서 높음)
✗ TF activity 무시하고 마커만으로 annotation
✗ Trajectory와 충돌하는 annotation
✗ 기능 상태 포함 (예: "Exhausted" - 이건 Tier 3)
```

---

## Tier 3: Functional States

### 목표 (Goal)

```
"이 세포가 현재 무엇을 하고 있는가?"

같은 developmental state 내에서 세포의 기능적 상태를 결정합니다.
활성화, 억제, 피로, 조절 등 현재 수행 중인 생물학적 기능을 정의합니다.
```

### 분류 기준 (EXAMPLE ONLY - NOT CONSTRAINTS)

> ⚠️ **WARNING: These are REFERENCE EXAMPLES from cancer immunology studies.**
>
> **DO NOT use as fixed functional states.**
> Your data may show:
> - Tissue-specific functional states
> - Disease-specific activity profiles
> - Novel functional subsets
> - Context-dependent activation
>
> **ALWAYS**:
> 1. Use DE + Pathway activity from YOUR data
> 2. Search PubMed for pathway + marker combinations
> 3. Use dynamic knowledge for functional states
> 4. Novel functions are COMMON (discover them!)

**T cells → Effector (Example - tissue/context dependent):**
| 기능 상태 (Example) | 의미 | 핵심 마커 (Example) | Key Pathway |
|-------------------|------|---------------------|-------------|
| **Cytotoxic** | 세포살해 활성 | GZMB, PRF1, NKG7 | TNFa, Trail |
| **Exhausted** | 기능 저하 | PDCD1, LAG3, TOX | (low activity) |
| **TRM** | 조직 상주 | ITGAE, CXCR6, ZNF683 | - |
| **Tpex** | 전구 피로 | TCF7, PDCD1, CXCR5 | WNT |

**B cells → GC_B (Example - secondary lymphoid organs):**
| 기능 상태 (Example) | 의미 | 핵심 마커 (Example) | Key Pathway |
|-------------------|------|---------------------|-------------|
| **Light_Zone** | 선택 단계 | CD83, BCL6 | - |
| **Dark_Zone** | 증식/변이 | AICDA, MKI67, CXCR4 | JAK-STAT |

**What's missing (examples of functional states you might find):**
- IFN-responsive states
- Hypoxia-adapted cells
- Metabolically distinct subsets
- Cycling vs. quiescent
- Stress response states
- **→ Use pathway activity + dynamic_knowledge.py to discover these!**

### Input

```
- Tier 2 subset (예: T cells → Effector만)
- Re-computed DE within subset
- Pathway activity (MANDATORY)
```

### Output

```python
adata.obs['tier3_annotation']   # 예: "Cytotoxic", "Exhausted"
adata.obs['final_annotation']   # 예: "T cells_Effector_Cytotoxic"
adata.obsm['mlm_estimate']      # Pathway activity scores
```

### 성공 기준

```
✓ Tier 2 마커가 아닌 NEW DE 마커 사용
✓ Pathway activity가 functional state와 일치
  - Cytotoxic: TNFa+, NFkB+, Trail+
  - Exhausted: Low pathway activity
  - Regulatory: TGFb+
✓ 마커 조합이 문헌에서 검증됨
✓ Novel population 적절히 flagging
```

### 실패 케이스

```
✗ Tier 2 마커 사용 (예: "GZMB high" - 모든 Effector에서 높음)
✗ Pathway activity 무시
✗ 단일 마커로 기능 상태 결정
```

---

## 계층 구조 예시

```
Full Dataset
├── T cells (Tier 1)
│   ├── Naive (Tier 2)
│   │   ├── Resting (Tier 3)
│   │   └── Activated_Naive (Tier 3)
│   ├── Effector (Tier 2)
│   │   ├── Cytotoxic (Tier 3)
│   │   ├── Exhausted (Tier 3)
│   │   └── TRM (Tier 3)
│   └── Memory (Tier 2)
│       ├── Central_Memory (Tier 3)
│       └── Effector_Memory (Tier 3)
│
├── B lineage (Tier 1)
│   ├── Pro_B (Tier 2)
│   │   ├── Cycling (Tier 3)
│   │   └── Resting (Tier 3)
│   ├── GC_B (Tier 2)
│   │   ├── Light_Zone (Tier 3)
│   │   └── Dark_Zone (Tier 3)
│   └── Plasma_cells (Tier 2)
│       ├── Short_lived (Tier 3)
│       └── Long_lived (Tier 3)
│
└── Myeloid (Tier 1)
    ├── Monocyte (Tier 2)
    ├── Macrophage (Tier 2)
    │   ├── M1_like (Tier 3)
    │   └── M2_like (Tier 3)
    └── DC (Tier 2)
```

---

## Quick Reference: 어떤 Tier인지 판단하기

| 질문 | Tier |
|------|------|
| "T cell인가 B cell인가?" | **Tier 1** |
| "Naive인가 Effector인가?" | **Tier 2** |
| "Cytotoxic인가 Exhausted인가?" | **Tier 3** |

| 특징 | Tier 1 | Tier 2 | Tier 3 |
|------|--------|--------|--------|
| **질문** | What lineage? | When in development? | What is it doing? |
| **한국어** | 어떤 계통? | 발달 단계? | 기능 상태? |
| **분석 도구** | DE markers | DE + TF + Trajectory | DE + Pathway |
| **마커 유형** | Lineage markers | Differentiation markers | Functional markers |
| **Pseudotime** | N/A | Early/Mid/Late | N/A |
| **Pathway** | N/A | N/A | TNFa, NFkB, etc. |

---

## 목표 달성 체크리스트

### Tier 1 완료 체크

```
[ ] 모든 클러스터가 major type으로 분류됨
[ ] Lineage-defining marker (2-4개) 확인됨
[ ] Cross-contamination 체크 완료
[ ] PMID 검증 완료
[ ] Subset 저장됨 (tier1/{major_type}.h5ad)
```

### Tier 2 완료 체크

```
[ ] DE RE-COMPUTED within major type
[ ] TF activity 계산됨 (ulm_estimate)
[ ] Trajectory 계산됨 (pseudotime) - if >= 2000 cells
[ ] TF + Trajectory가 developmental state와 일치
[ ] 모든 클러스터가 developmental state로 분류됨
[ ] PMID 검증 완료
[ ] Subset 저장됨 (tier2/{major_type}.h5ad)
```

### Tier 3 완료 체크

```
[ ] DE RE-COMPUTED within developmental state
[ ] Pathway activity 계산됨 (mlm_estimate)
[ ] Pathway가 functional state와 일치
[ ] 모든 클러스터가 functional state로 분류됨
[ ] Novel population 검토됨
[ ] final_annotation 생성됨
[ ] PMID 검증 완료
[ ] Subset 저장됨 (tier3/{major_type}_{dev_state}.h5ad)
```
