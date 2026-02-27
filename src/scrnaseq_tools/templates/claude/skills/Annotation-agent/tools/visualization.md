# Visualization Rules

> 📄 Code reference: `visualization_template.py` (same directory)

---

## Output Files (Per Tier)

```
Independent:
  {prefix}_umap.png/.svg
  {prefix}_annotation_marker_dotplot.png/.svg

Combined:
  {prefix}_umap_dotplot.png/.svg (2-panel)
  {prefix}_full_report.png/.svg  (4-panel)

Reports (PRIMARY deliverable):
  {prefix}_report.md              ← Consolidated markdown report
  {prefix}_evidence.json          ← Machine-readable evidence
  {prefix}_report.pptx            (optional)
```

---

## Style Rules

- **Font**: Arial (sans-serif fallback: DejaVu Sans, Helvetica)
- **DPI**: 150 (save at 300)
- **Font sizes**: title=7, label=6, tick=5, legend=5
- **SVG**: fonttype='none'

---

## DotPlot Rules (Critical)

- **dot_max**: Always `1.0` (full 0-100% fraction range)
- **Dot size**: Controlled via `largest_dot` (40-100 based on gene count)
- **Brackets**: Use `sc.pl.DotPlot OOP API` with `var_names=dict` for automatic brackets
- **var_group_rotation**: 90° for vertical bracket labels
- **Layout decision**:
  - `<= 15 genes`: Standard dotplot
  - `16-20 genes`: Split into panels
  - `> 20 genes`: Swap axes (genes on y-axis)
- **NO manual bracket drawing** — scanpy DotPlot handles it

---

## UMAP Dot Size (Adaptive)

| Cell Count | Dot Size | Legend Font |
|------------|----------|-------------|
| < 50K | 1.5 | 5 |
| 50K-100K | 1.0 | 4 |
| 100K-200K | 0.5 | 4 |
| > 200K | 0.3 | 3.5 |

---

## Report (.md) Format

```markdown
# Title
## Marker Genes, Functions, and Literature Evidence

**Species**: ...
**Date**: YYYY-MM-DD
**PMID Verification**: N/N VERIFIED

---

## 1. original -> Verified_Name | N cells | STATUS

### Marker Genes and Functions
| Gene | pct | Enrichment | Biological Function |

### Key TFs: TF1 (score), TF2 (score), ...

### Literature Evidence
| PMID (hyperlinked) | First Author | Year | Journal | Title |
```

---

## Key Functions (in visualization_template.py)

| Function | Purpose |
|----------|---------|
| `setup_publication_style()` | Set Arial, DPI, sizes |
| `save_all_visualizations()` | Main entry: UMAP + dotplot + report |
| `save_dotplot()` | Router: standard/swapped/split |
| `save_dotplot_standard()` | DotPlot OOP API with brackets |
| `save_markdown_report()` | Generate .md report |
| `save_combined()` | 2-panel UMAP + dotplot |
| `save_full_report()` | 4-panel report figure |

---

## Checklist

```
Per tier:
- [ ] {prefix}_report.md (PRIMARY)
- [ ] {prefix}_evidence.json
- [ ] {prefix}_umap.png/.svg
- [ ] {prefix}_annotation_marker_dotplot.png/.svg
- [ ] {prefix}_umap_dotplot.png/.svg
- [ ] {prefix}_full_report.png/.svg
- [ ] Font=Arial, DPI=150
```
