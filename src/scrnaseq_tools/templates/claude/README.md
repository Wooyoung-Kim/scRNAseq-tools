# .claude 디렉토리 구조

## 개요

이 디렉토리는 Claude Code의 프로젝트별 설정을 포함합니다.
Vaccine_V2 프로젝트의 single-cell RNA-seq 분석 및 hierarchical cell type annotation을 지원합니다.

## 구조

```
.claude/
├── agents/                    # Autonomous 실행용 agent 명세
│   └── annotation-executor.md
├── mcp-servers/               # 외부 도구 연동 MCP 서버
│   └── pubmed-tools/
│       ├── pyproject.toml
│       └── server.py
├── skills/                    # 재사용 가능한 도구/워크플로우
│   ├── Annotation-agent/
│   ├── scanpy/
│   ├── decoupler/
│   ├── palantir/
│   ├── scientific-visualization/
│   ├── anndata/
│   ├── biopython/
│   └── pptx/
├── settings.json              # 프로젝트 전역 설정
└── settings.local.json        # 로컬 개발 설정
```

---

## ⚠️ Skill Loading Rules (Context Optimization)

**DO NOT read all skill files at once.** Read ONLY the files needed for your CURRENT phase.  
Each `.md` contains rules/criteria. Code is in the matching `_template.py`.

| Phase | Read These Files | Code Reference |
|-------|-----------------|----------------|
| Phase 0 | `SKILL.md` + `phases/user_config.md` | — |
| Tier 1 | `phases/tier1.md` | `phases/tier1_template.py` |
| Tier 2 | `phases/tier2.md` | `phases/tier2_template.py` |
| Tier 3 | `phases/tier3.md` | `phases/tier3_template.py` |
| Phase 4 | `tools/visualization.md` | `tools/visualization_template.py` |
| As needed | `reasoning/integrated_format.md` | — |
| As needed | `tools/dynamic_knowledge.md` | `tools/dynamic_knowledge.py` |

---

## Approved Visualization Formats

### UMAP

Compact UMAP figures with cell-count-adaptive dot sizing and on-data legend.

```python
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=150, frameon=True, fontsize=5, figsize=(3, 3))

# Adaptive dot size based on cell count
n_cells = adata.n_obs
if n_cells < 50_000:
    size, legend_fs, legend_fo = 1.5, 5, 1.5
elif n_cells < 100_000:
    size, legend_fs, legend_fo = 1.0, 4, 1
elif n_cells < 200_000:
    size, legend_fs, legend_fo = 0.5, 4, 1
else:
    size, legend_fs, legend_fo = 0.3, 3.5, 1

# Single UMAP (3.35 x 3.0)
fig, ax = plt.subplots(figsize=(3.35, 3.0))
sc.pl.umap(adata, color='annotation', ax=ax, show=False,
           legend_loc='on data', legend_fontsize=legend_fs,
           legend_fontoutline=legend_fo, frameon=True, size=size)
fig.savefig('umap.png', dpi=150, bbox_inches='tight')

# 2-panel comparison (6.7 x 3.0)
fig, axes = plt.subplots(1, 2, figsize=(6.7, 3.0))
sc.pl.umap(adata, color='original', ax=axes[0], show=False,
           legend_loc='on data', legend_fontsize=legend_fs,
           legend_fontoutline=legend_fo, frameon=True, size=size,
           title='Original')
sc.pl.umap(adata, color='verified', ax=axes[1], show=False,
           legend_loc='on data', legend_fontsize=legend_fs,
           legend_fontoutline=legend_fo, frameon=True, size=size,
           title='Verified')
fig.suptitle('Verification', fontsize=7, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('umap_comparison.png', dpi=150, bbox_inches='tight')

# Multi-panel gene markers (4-col: 6.7 x 1.7*nrows, 3-col: 5.3 x 1.9*nrows)
markers = ['CD3E', 'CD14', 'MS4A1', 'JCHAIN']
ncols = 4
nrows = int(np.ceil(len(markers) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(6.7, 1.7 * nrows))
for i, gene in enumerate(markers):
    sc.pl.umap(adata, color=gene, ax=axes.flatten()[i], show=False,
               frameon=True, size=max(0.2, size * 0.5),
               title=gene, color_map='viridis', vmin=0)
    axes.flatten()[i].set_title(gene, fontsize=5, fontweight='bold')
    axes.flatten()[i].set_xlabel(''); axes.flatten()[i].set_ylabel('')
fig.savefig('umap_markers.png', dpi=150, bbox_inches='tight')
```

Key points:
- `legend_loc='on data'` with small fontsize + outline for readability
- Dot `size` scales inversely with cell count: 1.5 (<50K) → 0.3 (>200K)
- `frameon=True` for visible axes frame
- `dpi=150` for compact file size
- 2-panel: `(6.7, 3.0)`, single: `(3.35, 3.0)`, gene grid 4-col: `(6.7, 1.7*nrows)`
- `suptitle` fontsize=7 for compact titles

---

### Dotplot with Brackets

For cell type annotation dotplots, use grouped markers with brackets.
**Markers MUST be data-driven** — derived from DE results and evidence CSV of the current dataset.
NEVER hardcode canonical markers without verifying enrichment in the actual data.

```python
import scanpy as sc
import pandas as pd

# 1. Order annotations from actual adata.obs (no hardcoded list)
annotation_order = sorted(
    adata.obs['tier2_annotation'].unique(),
    key=lambda x: adata[adata.obs['tier2_annotation'] == x].obs['pseudotime'].median()
        if 'pseudotime' in adata.obs.columns else 0
)
adata.obs['tier2_annotation'] = pd.Categorical(
    adata.obs['tier2_annotation'],
    categories=annotation_order,
    ordered=True
)

# 2. Build marker_groups from evidence (e.g., tier2_annotation_evidence.csv)
#    Each marker must have pct >= 25% AND enrichment >= 1.5x vs others in the dataset.
#    If a cell type has NO specific markers (e.g., ribosomal-only), omit the bracket.
#    Example (must be rebuilt per dataset):
marker_groups = build_marker_groups_from_evidence(adata, evidence_csv_path)

# 3. Filter existing markers
filtered_marker_groups = {
    group: [m for m in markers if m in adata.var_names]
    for group, markers in marker_groups.items()
}
filtered_marker_groups = {k: v for k, v in filtered_marker_groups.items() if v}

# 4. Create dotplot with brackets using DotPlot OOP API
dp = sc.pl.DotPlot(
    adata,
    var_names=filtered_marker_groups,   # dict → automatic bracket grouping
    groupby='tier2_annotation',
    standard_scale='var',
    cmap='RdYlBu_r',
)
dp.style(
    dot_max=1.0,            # fraction legend 0-100%
    smallest_dot=0,         # 0% fraction → invisible
    largest_dot=60,         # physical max dot size (points^2)
    dot_edge_color='none',  # no edge border
    dot_edge_lw=0,
)
dp.var_group_rotation = 90  # bracket labels at 90°
dp.savefig('dotplot.png', dpi=300, bbox_inches='tight')
```

Key points:
- Use `sc.pl.DotPlot` OOP API (not `sc.pl.dotplot` functional) for precise dot size control
- `dot_max=1.0`: fraction legend always shows full 0-100% range
- `largest_dot=60`: controls physical max dot size to prevent overlap
- `dot_edge_color='none'`: no border around dots for cleaner look
- `var_group_rotation=90`: bracket labels at 90° (vertical)
- Use `pd.Categorical` with `ordered=True` for proper row ordering
- Pass `var_names` as dict for automatic bracket grouping
- `standard_scale='var'` for consistent scaling
- `cmap='RdYlBu_r'` for red-blue colormap
- **NEVER include markers that are anti-enriched in the target cell type** (e.g., CD27 in Memory_B if pct < 5%)
- If a cell type lacks specific DE markers, **omit the bracket** rather than using canonical markers that don't match the data

---

### Annotation Report (`save_markdown_report()`)

PRIMARY deliverable of Phase 4. Structured markdown with per-cell-type evidence.

```markdown
# Consolidated sbcl_ct Verification Report
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
| **DNTT** (TdT) | 95.4% | 211.8x | Terminal deoxynucleotidyl transferase. Adds non-templated nucleotides to V(D)J junctions. |
| **RAG1** | 90.8% | 11.6x | Recombination-activating gene 1. Initiates V(D)J recombination. |

### Key TFs: PAX5 (5.43), MYB (4.43), EBF1 (2.48), FOXO1 (2.24)

### Literature Evidence

| PMID | First Author | Year | Journal | Title |
|------|-------------|------|---------|-------|
| [35354960](https://pubmed.ncbi.nlm.nih.gov/35354960/) | Klein F | 2022 | *Nature Immunology* | Dntt expression reveals developmental hierarchy... |

---
```

Key points:
- Header: title, species, date, PMID verification count (`N/N VERIFIED via NCBI Entrez`)
- Lineage section: `# I. Lineage (N cells, N original types)` + source/clean/reclassifications/removals
- Per cell type: `## N. original -> Verified | N cells | STATUS`
- Marker table: `| **GENE** (alias) | pct | Enrichment | Biological Function |`
- Key TFs: `### Key TFs: TF1 (score), TF2 (score), ...`
- Literature table: `| [PMID](hyperlink) | First Author | Year | *Journal* | Title |`
- Use `save_markdown_report()` from `tools/visualization.md` with `report_kwargs` for metadata

---

## decoupler v2 Compatibility

decoupler v2.x changed its API. All code MUST use the v2 API:

| Old (v1) | New (v2) |
|----------|----------|
| `dc.run_ulm(mat=adata, net=net, source=..., target=..., weight=..., use_raw=False)` | `dc.mt.ulm(adata, net, verbose=True, raw=False)` |
| `dc.run_mlm(mat=adata, ...)` | `dc.mt.mlm(adata, net, verbose=True, raw=False)` |
| `adata.obsm['ulm_estimate']` | `adata.obsm['score_ulm']` (create alias for legacy code) |
| `adata.obsm['ulm_pvals']` | `adata.obsm['padj_ulm']` |
| `adata.obsm['mlm_estimate']` | `adata.obsm['score_mlm']` |

**Always** check for both key names (`score_ulm` or `ulm_estimate`) when reading results.

---

## Annotation Workflow Insights

### Claude Pre-Analysis Workflow
1. Show top 50 DE markers per cluster to Claude
2. Claude proposes cell type candidates with key marker combinations
3. MCP verifies only proposed combinations (reduces API calls)
4. Use simplified PubMed queries (2-4 key terms) for better results

### Contamination Handling
- **T cell contamination**: CD3E, CD3G, TRAC markers → Remove before Tier 3
- **Myeloid contamination**: S100A8, S100A9, LYZ markers → Flag and remove
- **Erythroid contamination**: HBA1, HBB, HBA2 markers → Flag and remove
- ~1-2% contamination is acceptable to flag and remove

### Data-Driven Marker Validation
Always validate canonical markers against actual DE results:
- Some markers may show unexpected patterns (e.g., VPREB1)
- Use multiple markers per cell type for robust identification
- Calculate enrichment ratio (cell type mean / max other) for specificity

---

## agents/

Autonomous 실행을 위한 agent 명세입니다.

| Agent | 설명 |
|-------|------|
| `annotation-executor.md` | 3-tier hierarchical annotation 자동 실행 |

## mcp-servers/

외부 도구 연동을 위한 MCP (Model Context Protocol) 서버입니다.

### pubmed-tools/

PubMed 검색 및 문헌 검증을 위한 MCP 서버입니다.

**제공 도구:**
- `pubmed_search`: PubMed 검색, PMID/제목/저자 반환
- `verify_reference`: PMID가 특정 마커-세포 유형 연관을 지지하는지 검증
- `fetch_abstract`: PMID의 전체 초록 가져오기

**API Key**: NCBI API key 등록됨 (rate limit 10 req/sec)

## skills/

재사용 가능한 도구 및 워크플로우입니다.

| Skill | 용도 |
|-------|------|
| **Annotation-agent** | 3-tier 계층적 세포 유형 annotation (Major → Developmental → Functional) |
| **scanpy** | Single-cell RNA-seq 분석 (QC, normalization, clustering, marker genes) |
| **decoupler** | TF/Pathway activity 추론 (DoRothEA, PROGENy, CollecTRI) |
| **palantir** | Trajectory 분석 및 pseudotime 추론 |
| **scientific-visualization** | Publication-quality 시각화 (matplotlib, seaborn) |
| **anndata** | AnnData 데이터 구조 처리 |
| **biopython** | 분자생물학 도구 (Bio.Entrez, BLAST 등) |
| **pptx** | PowerPoint 생성/편집 |

## 설정 파일

### settings.json

프로젝트 전역 설정입니다.

```json
{
  "mcpServers": {
    "pubmed-tools": { ... }
  },
  "permissions": {
    "allow": ["Bash(uv run*)", "Bash(python*)", "mcp__pubmed-tools__*"]
  }
}
```

### settings.local.json

로컬 개발 설정입니다. 필요에 따라 권한을 추가합니다.

## 데이터 흐름

```
Raw h5ad
    │
    ▼
[scanpy: QC/Preprocessing]
    │
    ▼
[Annotation-agent: Tier 1-3]
    │
    ├── Tier 1: Major Types (T cells, B cells, Myeloid, ...)
    │       │
    │       ▼
    ├── Tier 2: Developmental States (resolution 0.4-0.6)
    │       ├── [decoupler: TF Activity]
    │       └── [palantir: Trajectory]
    │       │
    │       ▼
    └── Tier 3: Functional States (resolution 0.1)
            └── [decoupler: Pathway Activity]
    │
    ▼
[scientific-visualization: Figures]
    │
    ▼
[pptx: Presentation]
```

## 관련 링크

- [Claude Code 공식 문서](https://docs.anthropic.com/claude-code)
- [MCP 프로토콜](https://modelcontextprotocol.io)
