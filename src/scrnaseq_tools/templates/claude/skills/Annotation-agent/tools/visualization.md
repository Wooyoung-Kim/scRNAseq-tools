# Visualization Code

Publication-ready visualization for hierarchical annotation.

---

## Output Files (Per Tier)

```
Independent files:
  {prefix}_umap.png/.svg
  {prefix}_dotplot.png/.svg

Combined files:
  {prefix}_umap_dotplot.png/.svg (2-panel)
  {prefix}_full_report.png/.svg  (4-panel)

Reports (PRIMARY deliverable):
  {prefix}_report.md              ← Consolidated markdown report
  {prefix}_evidence.json          ← Machine-readable evidence
  {prefix}_report.pptx            (optional)

Report format:
  Per cell type section with:
  - ## N. original -> Verified_Name | N cells | STATUS
  - Marker table: Gene | pct | Enrichment | Biological Function
  - Key TFs: TF1 (score), TF2 (score), ...
  - Literature table: PMID (hyperlinked) | First Author | Year | Journal | Title

Requirements: Font=Arial, DPI=150
```

---

## 1. Style Configuration

```python
import matplotlib.pyplot as plt
import scanpy as sc

def setup_publication_style():
    """Publication-ready style. Font: Arial."""
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
    plt.rcParams['font.size'] = 6
    plt.rcParams['axes.titlesize'] = 7
    plt.rcParams['axes.labelsize'] = 6
    plt.rcParams['xtick.labelsize'] = 5
    plt.rcParams['ytick.labelsize'] = 5
    plt.rcParams['legend.fontsize'] = 5
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['savefig.dpi'] = 150
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['svg.fonttype'] = 'none'
    sc.settings.set_figure_params(dpi=150, frameon=True, fontsize=5, figsize=(3, 3))

FIGURE_SIZES = {
    'single_column': (3.35, 3.0),        # single UMAP
    'single_column_tall': (3.35, 4.0),
    'double_column': (6.7, 3.0),          # 2-panel UMAP comparison
    'umap_dotplot': (6.7, 3.0),
    'marker_grid_3col': (5.3, None),      # per-row height = 1.9
    'marker_grid_4col': (6.7, None),      # per-row height = 1.7
    'full_report': (10.0, 8.0),
}


def calculate_umap_size(n_cells):
    """
    Calculate UMAP dot size based on cell count to prevent overlap.

    Returns dict with size (scatter point size) and legend fontsize.
    Larger datasets need smaller dots.
    """
    if n_cells < 50_000:
        return {'size': 1.5, 'legend_fontsize': 5, 'legend_fontoutline': 1.5}
    elif n_cells < 100_000:
        return {'size': 1.0, 'legend_fontsize': 4, 'legend_fontoutline': 1}
    elif n_cells < 200_000:
        return {'size': 0.5, 'legend_fontsize': 4, 'legend_fontoutline': 1}
    else:
        return {'size': 0.3, 'legend_fontsize': 3.5, 'legend_fontoutline': 1}
```

---

## 2. Main Save Function

```python
import scanpy as sc
import os
import json

def save_all_visualizations(adata, annotation_col, marker_dict,
                            evidence_list, annotation_order,
                            title_prefix, output_dir='./figures/',
                            report_kwargs=None):
    """
    Save all visualizations for a tier.

    Parameters
    ----------
    evidence_list : list[dict]
        List of evidence entries (from annotation_evidence.json schema).
        Each entry: annotation, markers[], tf_activity[], references[],
        confidence_level, n_cells, original_name, status, notes, etc.
    report_kwargs : dict, optional
        Extra kwargs for save_markdown_report():
        species, source_file, lineage_name, n_cells_total,
        n_original_types, reclassifications, removals, n_cells_clean.
    """
    os.makedirs(output_dir, exist_ok=True)
    setup_publication_style()
    prefix = title_prefix.replace(' ', '_').replace(':', '').replace('->', '_')

    # 1. UMAP (independent)
    umap_config = calculate_umap_size(adata.n_obs)
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single_column'])
    sc.pl.umap(adata, color=annotation_col, title=title_prefix,
               legend_loc='on data',
               legend_fontsize=umap_config['legend_fontsize'],
               legend_fontoutline=umap_config['legend_fontoutline'],
               size=umap_config['size'],
               frameon=True, ax=ax, show=False)
    save_figure(fig, f'{output_dir}/{prefix}_umap')
    plt.close(fig)

    # 2. Dotplot (independent)
    save_dotplot(adata, marker_dict, annotation_col, annotation_order,
                 f'{output_dir}/{prefix}_dotplot')

    # 3. Evidence report (markdown + json)
    rpt_kw = report_kwargs or {}
    save_markdown_report(evidence_list, f'{output_dir}/{prefix}_report.md',
                         title=title_prefix, **rpt_kw)
    save_reasoning_json(evidence_list, f'{output_dir}/{prefix}_evidence.json')

    # 4. Combined UMAP + Dotplot
    save_combined(adata, annotation_col, marker_dict, annotation_order,
                  title_prefix, f'{output_dir}/{prefix}_umap_dotplot')

    # 5. Full 4-panel report
    save_full_report(adata, annotation_col, marker_dict,
                     evidence_list, annotation_order,
                     title_prefix, f'{output_dir}/{prefix}_full_report')

    print(f"All visualizations saved to {output_dir}/{prefix}_*")


def save_figure(fig, filepath):
    """Save as PNG + SVG."""
    fig.savefig(f'{filepath}.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    fig.savefig(f'{filepath}.svg', bbox_inches='tight',
                facecolor='white', edgecolor='none')
```

---

## 3. Dotplot (No Overlap)

```python
import pandas as pd
import numpy as np

def save_dotplot(adata, marker_dict, annotation_col, annotation_order, filepath,
                 max_genes_per_row=15, swap_axes_threshold=20,
                 show_brackets=True):
    """
    Save dotplot with NO DOT OVERLAP and marker brackets.

    Parameters:
    -----------
    max_genes_per_row : int
        If more genes, split into multiple rows or swap axes
    swap_axes_threshold : int
        If more genes than this, swap axes (genes on y-axis)
    show_brackets : bool
        Show cell type brackets above gene names (default: True)

    Layout:
    - X-axis (bottom): Gene names
    - Y-axis: Cell type names
    - X-axis (top): Brackets showing marker ownership
    """
    # Collect markers in order (grouped by cell type)
    all_markers = []
    for ct in annotation_order:
        all_markers.extend(marker_dict.get(ct, []))
    seen = set()
    unique_markers = [m for m in all_markers if not (m in seen or seen.add(m))]

    n_genes = len(unique_markers)
    n_groups = len(annotation_order)

    # Set category order
    adata.obs[annotation_col] = pd.Categorical(
        adata.obs[annotation_col], categories=annotation_order, ordered=True
    )

    # Calculate optimal dot size to prevent overlap
    dot_config = calculate_dot_size(n_genes, n_groups)

    # Decide layout based on gene count
    if n_genes > swap_axes_threshold:
        # Many genes: swap axes (genes on y-axis, groups on x-axis)
        # Note: brackets not shown when swapped
        save_dotplot_swapped(adata, unique_markers, annotation_col,
                             annotation_order, dot_config, filepath,
                             marker_dict=marker_dict)
    elif n_genes > max_genes_per_row:
        # Moderate genes: split into panels with brackets
        save_dotplot_split(adata, unique_markers, annotation_col,
                           annotation_order, dot_config, filepath,
                           genes_per_panel=max_genes_per_row,
                           marker_dict=marker_dict if show_brackets else None)
    else:
        # Few genes: standard layout with brackets
        save_dotplot_standard(adata, unique_markers, annotation_col,
                              annotation_order, dot_config, filepath,
                              marker_dict=marker_dict if show_brackets else None)


def calculate_dot_size(n_genes, n_groups):
    """
    Calculate dot size parameters to prevent overlap.

    Key principle: dot_max=1.0 (always show full 0-100% fraction range),
    control physical dot size via largest_dot instead.

    Returns dict with:
    - dot_max: always 1.0 (full fraction range)
    - smallest_dot, largest_dot: physical dot size range (in points^2)
    """
    total_dots = n_genes * n_groups

    if total_dots <= 50:
        return {
            'dot_max': 1.0,
            'smallest_dot': 0,
            'largest_dot': 100,
            'figsize_per_gene': 0.35,
            'figsize_per_group': 0.4
        }
    elif total_dots <= 150:
        return {
            'dot_max': 1.0,
            'smallest_dot': 0,
            'largest_dot': 80,
            'figsize_per_gene': 0.3,
            'figsize_per_group': 0.35
        }
    elif total_dots <= 300:
        return {
            'dot_max': 1.0,
            'smallest_dot': 0,
            'largest_dot': 60,
            'figsize_per_gene': 0.25,
            'figsize_per_group': 0.3
        }
    else:
        return {
            'dot_max': 1.0,
            'smallest_dot': 0,
            'largest_dot': 40,
            'figsize_per_gene': 0.2,
            'figsize_per_group': 0.25
        }


def save_dotplot_standard(adata, markers, annotation_col, annotation_order,
                          dot_config, filepath, marker_dict=None):
    """
    Standard dotplot using sc.pl.DotPlot OOP API for precise dot size control.

    If marker_dict is provided, pass it as var_names (dict) for automatic
    bracket grouping by scanpy. Otherwise pass flat marker list.

    Layout:
    - X-axis (bottom): Gene names
    - Y-axis: Cell type names
    - X-axis (top): Brackets showing which cell type each gene belongs to
    """
    n_genes = len(markers)
    n_groups = len(annotation_order)

    # Dynamic figure size
    width = max(4.0, n_genes * dot_config['figsize_per_gene'] + 1.5)
    height = max(3.5, n_groups * dot_config['figsize_per_group'] + 1.5)

    # Use dict var_names for automatic bracket grouping
    var_names = marker_dict if marker_dict else markers

    dp = sc.pl.DotPlot(
        adata,
        var_names=var_names,
        groupby=annotation_col,
        standard_scale='var',
        cmap='RdYlBu_r',
        figsize=(width, height),
    )
    dp.style(
        dot_max=dot_config['dot_max'],
        smallest_dot=dot_config['smallest_dot'],
        largest_dot=dot_config['largest_dot'],
        dot_edge_color='none',
        dot_edge_lw=0,
    )
    dp.var_group_rotation = 90
    dp.savefig(f'{filepath}.png', dpi=300, bbox_inches='tight')
    plt.close('all')

    # Also save SVG
    dp = sc.pl.DotPlot(
        adata,
        var_names=var_names,
        groupby=annotation_col,
        standard_scale='var',
        cmap='RdYlBu_r',
        figsize=(width, height),
    )
    dp.style(
        dot_max=dot_config['dot_max'],
        smallest_dot=dot_config['smallest_dot'],
        largest_dot=dot_config['largest_dot'],
        dot_edge_color='none',
        dot_edge_lw=0,
    )
    dp.var_group_rotation = 90
    dp.savefig(f'{filepath}.svg', bbox_inches='tight')
    plt.close('all')


## NOTE: Manual bracket drawing is no longer needed.
## When using sc.pl.DotPlot OOP API with var_names as a dict,
## scanpy automatically draws brackets above the gene columns.
## Set dp.var_group_rotation = 90 for vertical bracket labels.


def save_dotplot_swapped(adata, markers, annotation_col, annotation_order,
                         dot_config, filepath, marker_dict=None):
    """
    Swapped dotplot (genes on y-axis) for many genes.
    Uses sc.pl.DotPlot OOP API for precise dot size control.
    """
    n_genes = len(markers)
    n_groups = len(annotation_order)

    # Swapped dimensions
    width = max(4.0, n_groups * dot_config['figsize_per_group'] + 2.0)
    height = max(4.0, n_genes * dot_config['figsize_per_gene'] + 1.0)

    var_names = marker_dict if marker_dict else markers

    dp = sc.pl.DotPlot(
        adata,
        var_names=var_names,
        groupby=annotation_col,
        standard_scale='var',
        cmap='RdYlBu_r',
        figsize=(width, height),
    )
    dp.style(
        dot_max=dot_config['dot_max'],
        smallest_dot=dot_config['smallest_dot'],
        largest_dot=dot_config['largest_dot'],
        dot_edge_color='none',
        dot_edge_lw=0,
    )
    dp.swap_axes()
    dp.var_group_rotation = 90
    dp.savefig(f'{filepath}.png', dpi=300, bbox_inches='tight')
    plt.close('all')

    # Also save SVG
    dp = sc.pl.DotPlot(
        adata,
        var_names=var_names,
        groupby=annotation_col,
        standard_scale='var',
        cmap='RdYlBu_r',
        figsize=(width, height),
    )
    dp.style(
        dot_max=dot_config['dot_max'],
        smallest_dot=dot_config['smallest_dot'],
        largest_dot=dot_config['largest_dot'],
        dot_edge_color='none',
        dot_edge_lw=0,
    )
    dp.swap_axes()
    dp.var_group_rotation = 90
    dp.savefig(f'{filepath}.svg', bbox_inches='tight')
    plt.close('all')


## NOTE: Manual y-axis bracket drawing also replaced by DotPlot OOP API.
## Use dp.swap_axes() for swapped layout with automatic bracket handling.


def save_dotplot_split(adata, markers, annotation_col, annotation_order,
                       dot_config, filepath, genes_per_panel=15, marker_dict=None):
    """
    Split dotplot into multiple panels using sc.pl.DotPlot OOP API.
    Each panel gets its own DotPlot with bracket grouping.
    """
    n_genes = len(markers)
    n_panels = int(np.ceil(n_genes / genes_per_panel))

    for i in range(n_panels):
        start_idx = i * genes_per_panel
        end_idx = min((i + 1) * genes_per_panel, n_genes)
        panel_markers = markers[start_idx:end_idx]

        # Build subset marker_dict for this panel
        panel_var_names = panel_markers
        if marker_dict:
            panel_marker_dict = {}
            for ct in annotation_order:
                ct_markers = marker_dict.get(ct, [])
                panel_ct_markers = [m for m in ct_markers if m in panel_markers]
                if panel_ct_markers:
                    panel_marker_dict[ct] = panel_ct_markers
            if panel_marker_dict:
                panel_var_names = panel_marker_dict

        suffix = f'_part{i+1}' if n_panels > 1 else ''

        dp = sc.pl.DotPlot(
            adata,
            var_names=panel_var_names,
            groupby=annotation_col,
            standard_scale='var',
            cmap='RdYlBu_r',
        )
        dp.style(
            dot_max=dot_config['dot_max'],
            smallest_dot=dot_config['smallest_dot'],
            largest_dot=dot_config['largest_dot'],
            dot_edge_color='none',
            dot_edge_lw=0,
        )
        dp.var_group_rotation = 90
        dp.savefig(f'{filepath}{suffix}.png', dpi=300, bbox_inches='tight')
        plt.close('all')
```

---

## 4. Markdown Report Output

Generate a consolidated verification/annotation report in the standard format.
This is the **primary deliverable** — a structured markdown with marker tables,
TF activity, and hyperlinked literature evidence per cell type.

### Report Structure

```
# Report Title
## Subtitle

**Species**: ...
**Date**: YYYY-MM-DD
**PMID Verification**: N/N VERIFIED

---

# I. Lineage Name (N cells, N original types)

**Source**: `file.h5ad`
**Reclassifications**: N (details)
**Removals**: N clusters (details)

---

## 1. original_name -> Verified_Name | N cells | STATUS

### Marker Genes and Functions

| Gene | pct | Enrichment | Biological Function |
|------|-----|-----------|-------------------|
| **GENE** (alias) | XX.X% | XX.Xx | Description... |

### Key TFs: TF1 (score), TF2 (score), ...

### Literature Evidence

| PMID | First Author | Year | Journal | Title |
|------|-------------|------|---------|-------|
| [PMID](https://pubmed.ncbi.nlm.nih.gov/PMID/) | Author | Year | *Journal* | Title |

---
```

### Report Generation Code

```python
from datetime import date

def save_markdown_report(evidence_list, output_path, title="Annotation Report",
                         species=None, source_file=None,
                         lineage_name=None,
                         n_cells_total=None, n_original_types=None,
                         reclassifications=None, removals=None,
                         n_cells_clean=None):
    """
    Generate consolidated markdown report from annotation evidence.

    Output format matches the reference structure:
      # Title
      ## Marker Genes, Functions, and Literature Evidence
      **Species** / **Date** / **PMID Verification**
      # I. Lineage (N cells, N original types)
      **Source** / **Clean output** / **Reclassifications** / **Removals**
      ## 1. original -> Verified | N cells | STATUS
        ### Marker Genes and Functions  (table)
        ### Key TFs: ...
        ### Literature Evidence  (table with hyperlinked PMIDs)

    Parameters
    ----------
    evidence_list : list[dict]
        List of evidence entries from annotation_evidence.json.
        Each entry has: annotation, markers[], tf_activity[], references[],
        confidence_level, n_cells, original_name, status, etc.
    output_path : str
        Output .md file path.
    title : str
        Report title.
    species : str, optional
        Species name (e.g., "Ferret (*Mustela putorius furo*)").
    source_file : str, optional
        Source h5ad filename.
    lineage_name : str, optional
        Lineage name for section header (e.g., "B Lineage", "T Lineage").
    n_cells_total : int, optional
        Total cell count before filtering.
    n_original_types : int, optional
        Number of original cell types.
    reclassifications : list[str], optional
        List of reclassification descriptions (e.g., ["NB3->Pre_B"]).
    removals : list[str], optional
        List of removed cluster descriptions (e.g., ["PB2=T/NK contamination"]).
    n_cells_clean : int, optional
        Cell count after filtering.
    """
    lines = []

    # --- Header ---
    lines.append(f"# {title}")
    lines.append("## Marker Genes, Functions, and Literature Evidence")
    lines.append("")
    if species:
        lines.append(f"**Species**: {species}")
    lines.append(f"**Date**: {date.today().isoformat()}")

    # Count verified PMIDs
    total_pmids = sum(len(e.get('references', [])) for e in evidence_list)
    verified_pmids = sum(
        1 for e in evidence_list
        for r in e.get('references', [])
        if r.get('status', '').startswith('VERIFIED') or
           r.get('status', '') == 'DOUBLE_VERIFIED'
    )
    lines.append(f"**PMID Verification**: {verified_pmids}/{total_pmids} VERIFIED via NCBI Entrez")
    lines.append("")
    lines.append("---")
    lines.append("")

    # --- Lineage section header ---
    lineage = lineage_name or "Lineage"
    cells_str = f"{n_cells_total:,} cells" if n_cells_total else ""
    types_str = f"{n_original_types} original types" if n_original_types else ""
    header_parts = [p for p in [cells_str, types_str] if p]
    if header_parts:
        lines.append(f"# I. {lineage} ({', '.join(header_parts)})")
    else:
        lines.append(f"# I. {lineage}")
    lines.append("")
    if source_file:
        lines.append(f"**Source**: `{source_file}`")
    if n_cells_clean and n_cells_total:
        removed = n_cells_total - n_cells_clean
        pct = removed / n_cells_total * 100
        lines.append(f"**Clean output**: {n_cells_clean:,} cells "
                     f"(removed {removed:,} = {pct:.0f}%)")
    if reclassifications:
        lines.append(f"**Reclassifications**: {len(reclassifications)} "
                     f"({', '.join(reclassifications)})")
    if removals:
        lines.append(f"**Removals**: {len(removals)} clusters "
                     f"({', '.join(removals)})")
    lines.append("")
    lines.append("---")
    lines.append("")

    # --- Per cell type sections ---
    for i, entry in enumerate(evidence_list, 1):
        annotation = entry.get('annotation', 'Unknown')
        original = entry.get('original_name', entry.get('cluster_id', ''))
        n_cells = entry.get('n_cells', '')
        status = entry.get('status', entry.get('confidence_level', 'CONFIRMED'))

        # Section header
        cell_count = f" | {n_cells:,} cells" if isinstance(n_cells, int) else ""
        if original and original != annotation:
            lines.append(f"## {i}. {original} -> {annotation}{cell_count} | {status}")
        else:
            lines.append(f"## {i}. {annotation}{cell_count} | {status}")
        lines.append("")

        # --- Marker table ---
        markers = entry.get('markers', [])
        if markers:
            lines.append("### Marker Genes and Functions")
            lines.append("")
            lines.append("| Gene | pct | Enrichment | Biological Function |")
            lines.append("|------|-----|-----------|-------------------|")

            for m in markers:
                gene = m.get('gene', '')
                alias = m.get('alias', '')
                pct = m.get('pct_in', m.get('pct', ''))
                enrichment = m.get('enrichment', m.get('log2fc', ''))
                function = m.get('function', m.get('biological_function', ''))

                # Format gene name
                gene_str = f"**{gene}**"
                if alias:
                    gene_str += f" ({alias})"

                # Format pct
                if isinstance(pct, (int, float)):
                    pct_str = f"{pct:.1f}%" if pct <= 1 else f"{pct:.1f}%"
                    if pct <= 1:  # fraction format
                        pct_str = f"{pct*100:.1f}%"
                else:
                    pct_str = str(pct)

                # Format enrichment
                if isinstance(enrichment, (int, float)):
                    enrich_str = f"{enrichment:.1f}x"
                else:
                    enrich_str = str(enrichment)

                lines.append(f"| {gene_str} | {pct_str} | {enrich_str} | {function} |")

            lines.append("")

        # --- Key TFs ---
        tf_activity = entry.get('tf_activity', [])
        if tf_activity:
            tf_parts = []
            for tf in tf_activity[:6]:  # Top 6 TFs
                name = tf.get('tf', tf.get('name', ''))
                score = tf.get('score', tf.get('activity', ''))
                if isinstance(score, (int, float)):
                    tf_parts.append(f"{name} ({score:.2f})")
                else:
                    tf_parts.append(f"{name} ({score})")
            lines.append(f"### Key TFs: {', '.join(tf_parts)}")
            lines.append("")

        # --- Notes ---
        notes = entry.get('notes', entry.get('note', ''))
        if notes:
            lines.append(f"### Note: {notes}")
            lines.append("")

        # --- Literature table ---
        references = entry.get('references', [])
        if references:
            lines.append("### Literature Evidence")
            lines.append("")
            lines.append("| PMID | First Author | Year | Journal | Title |")
            lines.append("|------|-------------|------|---------|-------|")

            for ref in references:
                pmid = ref.get('pmid', '')
                first_author = ref.get('first_author', ref.get('author', ''))
                year = ref.get('year', '')
                journal = ref.get('journal', '')
                ref_title = ref.get('title', ref.get('finding', ''))

                # Hyperlink PMID
                pmid_str = f"[{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
                # Italicize journal
                journal_str = f"*{journal}*" if journal else ""

                lines.append(f"| {pmid_str} | {first_author} | {year} "
                             f"| {journal_str} | {ref_title} |")

            lines.append("")

        lines.append("---")
        lines.append("")

    # Write
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f"Report saved: {output_path}")


def save_reasoning_json(evidence_list, output_path):
    """Save evidence as structured JSON (machine-readable companion to .md report)."""
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump({'annotations': evidence_list}, f, indent=2, ensure_ascii=False)
    print(f"Evidence JSON saved: {output_path}")
```

---

## 5. Combined Figures

```python
from matplotlib.gridspec import GridSpec

def save_combined(adata, annotation_col, marker_dict, annotation_order,
                  title_prefix, filepath, show_brackets=True):
    """
    Save UMAP + Dotplot (2-panel) with no dot overlap and marker brackets.

    Layout:
    - Left: UMAP
    - Right: Dotplot with brackets showing marker ownership
    """
    # Collect markers
    all_markers = []
    for ct in annotation_order:
        all_markers.extend(marker_dict.get(ct, []))
    seen = set()
    unique_markers = [m for m in all_markers if not (m in seen or seen.add(m))]

    n_genes = len(unique_markers)
    n_groups = len(annotation_order)

    # Calculate dot size
    dot_config = calculate_dot_size(n_genes, n_groups)

    # Dynamic figure size (add space for brackets)
    dot_width = max(4.0, n_genes * dot_config['figsize_per_gene'] + 1.5)
    total_width = 3.5 + dot_width + 0.5
    height = max(4.0, n_groups * dot_config['figsize_per_group'] + 1.5)

    fig = plt.figure(figsize=(total_width, height))
    gs = GridSpec(1, 2, width_ratios=[3.5, dot_width], wspace=0.3)

    # UMAP
    umap_config = calculate_umap_size(adata.n_obs)
    ax_umap = fig.add_subplot(gs[0])
    sc.pl.umap(adata, color=annotation_col, title='Cell Types',
               legend_loc='on data',
               legend_fontsize=umap_config['legend_fontsize'],
               legend_fontoutline=umap_config['legend_fontoutline'],
               size=umap_config['size'],
               frameon=True, ax=ax_umap, show=False)

    # Dotplot with no overlap (use DotPlot OOP API on separate figure, then embed)
    # Note: For combined figures, save dotplot separately then compose,
    # or use the functional API with controlled dot sizes.
    ax_dot = fig.add_subplot(gs[1])
    adata.obs[annotation_col] = pd.Categorical(
        adata.obs[annotation_col], categories=annotation_order, ordered=True
    )

    sc.pl.dotplot(
        adata,
        var_names=marker_dict if show_brackets else unique_markers,
        groupby=annotation_col,
        categories_order=annotation_order,
        ax=ax_dot,
        show=False,
        dendrogram=False,
        dot_max=dot_config['dot_max'],
        smallest_dot=dot_config['smallest_dot'],
        largest_dot=dot_config['largest_dot'],
        color_map='RdYlBu_r',
        standard_scale='var',
        title='Marker Expression'
    )

    # Rotate labels
    plt.sca(ax_dot)
    plt.xticks(rotation=90, ha='center', fontsize=7)

    fig.suptitle(title_prefix, fontsize=10, fontweight='bold', y=1.02)
    save_figure(fig, filepath)
    plt.close(fig)


def save_full_report(adata, annotation_col, marker_dict,
                     evidence_list, annotation_order,
                     title_prefix, filepath):
    """Save 4-panel figure report with no dot overlap."""
    all_markers = []
    for ct in annotation_order:
        all_markers.extend(marker_dict.get(ct, []))
    seen = set()
    unique_markers = [m for m in all_markers if not (m in seen or seen.add(m))]

    n_genes = len(unique_markers)
    n_groups = len(annotation_order)
    dot_config = calculate_dot_size(n_genes, n_groups)

    fig = plt.figure(figsize=FIGURE_SIZES['full_report'])
    gs = GridSpec(2, 2, height_ratios=[1.2, 1], width_ratios=[1, 1.2],
                  wspace=0.25, hspace=0.35)

    # TOP-LEFT: UMAP
    umap_config = calculate_umap_size(adata.n_obs)
    ax_umap = fig.add_subplot(gs[0, 0])
    sc.pl.umap(adata, color=annotation_col, title='Cell Types',
               legend_loc='on data',
               legend_fontsize=umap_config['legend_fontsize'],
               legend_fontoutline=umap_config['legend_fontoutline'],
               size=umap_config['size'],
               frameon=True, ax=ax_umap, show=False)

    # TOP-RIGHT: Dotplot
    ax_dot = fig.add_subplot(gs[0, 1])
    adata.obs[annotation_col] = pd.Categorical(
        adata.obs[annotation_col], categories=annotation_order, ordered=True
    )
    sc.pl.dotplot(
        adata,
        var_names=marker_dict,
        groupby=annotation_col,
        categories_order=annotation_order,
        ax=ax_dot,
        show=False,
        dendrogram=False,
        dot_max=dot_config['dot_max'],
        smallest_dot=dot_config['smallest_dot'],
        largest_dot=dot_config['largest_dot'],
        color_map='RdYlBu_r',
        standard_scale='var',
        title='Markers'
    )
    plt.sca(ax_dot)
    plt.xticks(rotation=90, ha='center', fontsize=6)

    # BOTTOM-LEFT: Summary from evidence
    ax_reason = fig.add_subplot(gs[1, 0])
    ax_reason.axis('off')
    reason_lines = []
    for e in evidence_list[:8]:
        ann = e.get('annotation', '?')
        conf = e.get('confidence_level', '?')
        reason_lines.append(f"* {ann}: {conf}")
    if len(evidence_list) > 8:
        reason_lines.append(f"  ... +{len(evidence_list)-8} more")
    ax_reason.text(0.02, 0.98, '\n'.join(reason_lines),
                   transform=ax_reason.transAxes, fontsize=7,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='#f0f0f0', alpha=0.5))
    ax_reason.set_title('Annotations', fontsize=9, fontweight='bold', loc='left')

    # BOTTOM-RIGHT: Key references
    ax_ref = fig.add_subplot(gs[1, 1])
    ax_ref.axis('off')
    ref_lines = []
    for e in evidence_list:
        for r in e.get('references', [])[:1]:
            status = 'VV' if r.get('status') == 'DOUBLE_VERIFIED' else 'V'
            ref_lines.append(f"[{status}] PMID:{r.get('pmid', 'N/A')}")
            if len(ref_lines) >= 10:
                break
        if len(ref_lines) >= 10:
            break
    ax_ref.text(0.02, 0.98, '\n'.join(ref_lines),
                transform=ax_ref.transAxes, fontsize=7,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='#fffacd', alpha=0.5))
    ax_ref.set_title('References', fontsize=9, fontweight='bold', loc='left')

    fig.suptitle(f'{title_prefix} - Report', fontsize=11, fontweight='bold', y=0.98)
    save_figure(fig, filepath)
    plt.close(fig)
```

---

## 6. PPTX Report

```python
from pptx import Presentation
from pptx.util import Inches, Pt

def create_pptx_report(tier_name, figure_paths, evidence_list,
                       output_path='annotation_report.pptx'):
    """Create PPTX report (no API calls)."""
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    layout = prs.slide_layouts[6]  # Blank

    # Title slide
    slide = prs.slides.add_slide(layout)
    title = slide.shapes.add_textbox(Inches(0.5), Inches(2.5), Inches(12), Inches(1))
    p = title.text_frame.paragraphs[0]
    p.text = f"Hierarchical Annotation: {tier_name}"
    p.font.size = Pt(44)
    p.font.bold = True
    p.font.name = 'Arial'

    # Figures slide
    slide = prs.slides.add_slide(layout)
    if 'umap' in figure_paths:
        slide.shapes.add_picture(f"{figure_paths['umap']}.png",
                                 Inches(0.5), Inches(1), Inches(5.5), Inches(5.5))
    if 'dotplot' in figure_paths:
        slide.shapes.add_picture(f"{figure_paths['dotplot']}.png",
                                 Inches(6.5), Inches(1), Inches(6), Inches(5.5))

    prs.save(output_path)
    print(f"PPTX saved: {output_path}")
```

---

## Checklist

```
Per tier completion:
- [ ] {prefix}_report.md           ← PRIMARY: Consolidated markdown report
- [ ] {prefix}_evidence.json       ← Machine-readable evidence
- [ ] {prefix}_umap.png/.svg
- [ ] {prefix}_dotplot.png/.svg
- [ ] {prefix}_umap_dotplot.png/.svg
- [ ] {prefix}_full_report.png/.svg
- [ ] {prefix}_report.pptx (optional)

Report format requirements:
- [ ] Per cell type: ## N. original -> Verified | N cells | STATUS
- [ ] Marker table: Gene (bold+alias) | pct | Enrichment | Biological Function
- [ ] Key TFs line with scores
- [ ] Literature table: PMID (hyperlinked) | First Author | Year | Journal | Title
- [ ] Header: Species, Date, PMID verification count

Figure requirements:
- [ ] Font: Arial
- [ ] DPI: 150
```
