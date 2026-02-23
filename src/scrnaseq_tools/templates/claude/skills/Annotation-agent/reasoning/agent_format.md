# Agent Reasoning Format

Compact 3-iteration reasoning template for cell type annotation.
Reasoning output maps to `evidence.json` → `save_markdown_report()` generates `report.md` (PRIMARY deliverable).

---

## Overview

```
Iteration 1: Evidence Collection + PRE-Reasoning Literature
Iteration 2: Reasoning + Verification
Iteration 3: POST-Reasoning Literature + Final Decision
```

---

## Iteration 1: Evidence Collection

### Bioinformatician Role

```markdown
**Bioinformatician**:
Cluster X - Top 50 Valid Markers (pct>=25%, LFC>=1, padj<0.05):

| Rank | Gene | pct_in | log2FC | padj |
|------|------|--------|--------|------|
| 1 | [gene] | [val] | [val] | [val] |
| 2 | [gene] | [val] | [val] | [val] |
...
[Show all 50]

Functional Groups Identified:
- Group A: [genes]
- Group B: [genes]
- Group C: [genes]

TF Activity: [TF1] ([score]), [TF2] ([score])
Pathway Activity: [Pathway1] ([score]), [Pathway2] ([score])
```

### Computational Biologist Role

```markdown
**Computational Biologist**:
Marker Combination Analysis:

Candidate 1: [Cell type name]
  -> Combination: [marker1] + [marker2] + [marker3]
  -> Rationale: [brief explanation]

Candidate 2: [Cell type name]
  -> Combination: [marker1] + [marker2] + [marker3]
  -> Rationale: [brief explanation]

Candidate 3: [Cell type name]
  -> Combination: [marker1] + [marker2] + [marker3]
  -> Rationale: [brief explanation]

PRE-REASONING LITERATURE SEARCH (PMID REQUIRED):

Query 1: "[marker combo] [cell type context]"
  - **PMID: [number]** - "[key finding]" ([journal], [year])
  - **PMID: [number]** - "[key finding]" ([journal], [year])

Query 2: "[marker combo] [cell type context]"
  - **PMID: [number]** - "[key finding]" ([journal], [year])
  - **PMID: [number]** - "[key finding]" ([journal], [year])

Query 3: "[marker combo] [cell type context]"
  - **PMID: [number]** - "[key finding]" ([journal], [year])
  - **PMID: [number]** - "[key finding]" ([journal], [year])

⚠️ WARNING: Literature search WITHOUT valid PMID is INVALID
```

---

## Iteration 2: Reasoning + Verification

### Scientific Critic Role

```markdown
**Scientific Critic**:
REASONING - Comparing Candidates:

Candidate 1 ([name]):
  [+] Evidence supporting
  [-] Evidence against
  -> Assessment

Candidate 2 ([name]):
  [+] Evidence supporting
  [-] Evidence against
  -> Assessment

Candidate 3 ([name]):
  [+] Evidence supporting
  [-] Evidence against
  -> Assessment

Reference Verification:
PMID:xxx ([description]):
  - Content: "[relevant quote/finding]"
  - Supports: [marker-cell type link]
  - Status: **VERIFIED**

PMID:xxx ([description]):
  - Content: "[relevant quote/finding]"
  - Supports: [marker-cell type link]
  - Status: **VERIFIED**

Preliminary decision: [Candidate X] seems best fit because...
Need additional verification for: [specific aspect]
```

---

## Iteration 3: Final Decision

### Domain Expert Role

```markdown
**Domain Expert**:
POST-REASONING LITERATURE SEARCH:

Query: "[chosen cell type] [specific characteristic] single cell"
  - PMID:xxx - "[confirms/supports]"
  - PMID:xxx - "[confirms/supports]"

Re-Verification (Double Check):

Marker Combination: [final markers]
  - PMID:xxx: [verification result]
  - PMID:xxx: [verification result]
  - **DOUBLE VERIFIED**
```

### PI Role (Final Decision)

```markdown
**PI**:
Final Evidence Summary:

| Marker Combo | Candidate | Literature | Decision |
|--------------|-----------|------------|----------|
| [combo 1] | [name] | [status] | [reject/consider] |
| [combo 2] | [name] | [status] | [reject/consider] |
| [combo 3] | [name] | [status] | **SELECTED** |

**Final Assignment**: [Cell type name]
**Marker Combination**: [marker1] + [marker2] + [marker3]
**Confidence**: [High/Medium/Low]
**Key References**: PMID:xxx, PMID:xxx
**Reasoning Path**: Literature -> Compare -> Literature -> Confirm
```

---

## Compact Format (Alternative)

For efficiency, use this compact single-block format:

```markdown
## Cluster X Annotation

### Evidence
- Top markers: [gene1], [gene2], [gene3], [gene4], [gene5]
- TF activity: [TF1] (high), [TF2] (moderate)
- Pathway: [Pathway1] (active)

### Candidates
1. [Type1]: [marker combo] - [brief rationale]
2. [Type2]: [marker combo] - [brief rationale]
3. [Type3]: [marker combo] - [brief rationale]

### Literature (PRE)
- Query: "[search terms]" -> PMID:xxx ([finding])
- Query: "[search terms]" -> PMID:xxx ([finding])

### Reasoning
- [Type1]: [+points] vs [-points]
- [Type2]: [+points] vs [-points]
- Best fit: [TypeX] because [reason]

### Literature (POST)
- Verification: PMID:xxx -> VERIFIED
- Double-check: PMID:xxx -> DOUBLE_VERIFIED

### Decision
- **Annotation**: [Final cell type]
- **Markers**: [combo]
- **Confidence**: [High/Medium/Low]
- **References**: [PMIDs]
```

---

## Key Rules

1. **Always generate 2-3 candidates** before deciding
2. **Search literature BEFORE reasoning** (PRE-search)
3. **Search literature AFTER reasoning** (POST-verification)
4. **Double-verify** all references before final assignment
5. **Use 2-5 markers** per candidate (never single marker)
6. **Document reasoning path** for reproducibility

---

## MANDATORY: Reasoning → evidence.json Mapping

```
╔══════════════════════════════════════════════════════════════════════╗
║  Every reasoning output MUST map to evidence.json fields below       ║
║  If a field cannot be filled, the reasoning is INCOMPLETE            ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Field Mapping Table

| Reasoning Element | JSON Field | Required | Used in report.md |
|-------------------|------------|----------|-------------------|
| Final Assignment | `annotation` | YES | Section header |
| Original name | `original_name` | YES | `## N. original -> Verified` |
| Cell count | `n_cells` | YES | `| N cells` |
| Verification status | `status` | YES | `| STATUS` |
| Marker Combination | `markers[]` | YES | Marker table |
| — gene alias | `markers[].alias` | NO | `**GENE** (alias)` |
| — enrichment | `markers[].enrichment` | YES | `Nx` column |
| — function | `markers[].function` | YES | Biological Function column |
| Confidence (High/Med/Low) | `confidence_level` | YES | — |
| Confidence scoring | `confidence_score` | YES (0-12) | — |
| Key References | `references[]` | YES | Literature table |
| — first author | `references[].first_author` | YES | First Author column |
| — year | `references[].year` | YES | Year column |
| — journal | `references[].journal` | YES | *Journal* column |
| — title | `references[].title` | YES | Title column |
| Candidate comparison | `candidates[]` | YES | — |
| Supporting evidence | `reasons[]` | YES | — |
| Conflicting evidence | `conflicts[]` | NO (if none) | — |
| TF Activity | `tf_activity[]` | if computed | `### Key TFs:` |
| Pathway Activity | `pathway_activity[]` | if computed | — |
| Novel flag | `is_novel` | YES | — |

### evidence.json Entry Schema

Fields marked with `→ report` are used by `save_markdown_report()` to generate the report.

```json
{
  "cluster_id": "tier3_{major_type}_{dev_state}_{cluster_num}",
  "tier": "tier3",
  "subset_id": "{major_type}_{dev_state}",

  "annotation": "{func_state}",                   // → report: section header
  "original_name": "{original_sbcl_ct}",           // → report: "## N. original -> Verified"
  "full_label": "{major_type}_{dev_state}_{func_state}",
  "n_cells": 3046,                                 // → report: "| N cells"
  "status": "CONFIRMED",                           // → report: "| STATUS" (CONFIRMED/RECLASSIFICATION/PARTIALLY CONFIRMED)

  "markers": [
    {
      "gene": "{marker1}",
      "alias": "{common_name}",                    // → report: "**GENE** (alias)"
      "pct_in": 95.4,                              // → report: pct column
      "enrichment": 211.8,                         // → report: Enrichment column (Nx)
      "function": "Biological function description",// → report: Biological Function column
      "log2fc": "{val}",
      "padj": "{val}"
    }
  ],

  "confidence_score": "{0-12}",
  "confidence_level": "{High/Medium/Low}",

  "candidates": [
    {
      "name": "{candidate1}",
      "markers": ["{marker1}", "{marker2}", "{marker3}"],
      "support": ["{supporting_evidence_1}", "{supporting_evidence_2}"],
      "against": [],
      "selected": true
    },
    {
      "name": "{candidate2}",
      "markers": ["{marker1}", "{marker2}", "{marker3}"],
      "support": ["{supporting_evidence}"],
      "against": ["{conflicting_evidence}"],
      "selected": false
    }
  ],

  "reasons": [
    "{marker_combo} combination matches {func_state} signature",
    "{TF} TF activity high ({score})",
    "No {negative_marker} expression rules out {alternative}"
  ],

  "conflicts": [],

  "references": [
    {
      "pmid": "{number}",                          // → report: [PMID](hyperlink)
      "first_author": "{Last Name} {Initial}",     // → report: First Author column
      "year": 2022,                                // → report: Year column
      "journal": "{Journal Name}",                 // → report: *Journal* column
      "title": "{Full article title}",             // → report: Title column
      "query": "{marker_combo} {cell_type_context}",
      "finding": "{key_finding}",
      "first_verify": "VERIFIED",
      "second_verify": "VERIFIED",
      "status": "DOUBLE_VERIFIED"
    }
  ],

  "tf_activity": [                                 // → report: "### Key TFs: TF1 (score), ..."
    {"tf": "{TF_name}", "score": "{val}", "pval": "{val}"}
  ],

  "pathway_activity": [
    {"pathway": "{pathway_name}", "score": "{val}", "pval": "{val}"}
  ],

  "is_novel": false,
  "novel_evidence": null
}
```

### Conversion Function

```python
def reasoning_to_evidence(cluster_id, tier, subset_id, reasoning_output):
    """Convert structured reasoning to evidence.json entry.

    This function MUST be called for every annotation decision.
    The output feeds directly into save_markdown_report() for report generation.
    """
    evidence = {
        "cluster_id": cluster_id,
        "tier": tier,
        "subset_id": subset_id,

        # From PI Final Decision → report section header
        "annotation": reasoning_output['final_assignment'],
        "original_name": reasoning_output.get('original_name', ''),
        "full_label": reasoning_output['full_label'],
        "n_cells": reasoning_output.get('n_cells', 0),
        "status": reasoning_output.get('status', 'CONFIRMED'),

        # Markers → report marker table
        # Each marker must have: gene, pct_in, enrichment, function
        # Optional: alias, log2fc, padj
        "markers": reasoning_output['marker_details'],

        # From decision/matrices.md scoring
        "confidence_score": reasoning_output['confidence_score'],
        "confidence_level": reasoning_output['confidence_level'],

        # From Computational Biologist candidates
        "candidates": reasoning_output['candidates'],

        # From Scientific Critic reasoning
        "reasons": reasoning_output['supporting_reasons'],
        "conflicts": reasoning_output.get('conflicts', []),

        # Literature → report literature table
        # Each ref must have: pmid, first_author, year, journal, title
        "references": reasoning_output['references'],

        # TF activity → report "### Key TFs: TF1 (score), ..."
        "tf_activity": reasoning_output.get('tf_activity', []),
        "pathway_activity": reasoning_output.get('pathway_activity', []),

        # Novel detection
        "is_novel": reasoning_output.get('is_novel', False),
        "novel_evidence": reasoning_output.get('novel_evidence', None)
    }

    # Validate required fields (including report-critical fields)
    required = ['annotation', 'markers', 'confidence_level', 'references',
                'candidates', 'reasons', 'n_cells', 'status']
    missing = [f for f in required if not evidence.get(f)]
    if missing:
        raise ValueError(f"INCOMPLETE reasoning: missing {missing}")

    # Validate marker fields needed for report
    for m in evidence['markers']:
        for key in ['gene', 'pct_in', 'enrichment', 'function']:
            if key not in m:
                raise ValueError(f"Marker {m.get('gene','?')} missing '{key}' (needed for report)")

    # Validate reference fields needed for report
    for r in evidence['references']:
        for key in ['pmid', 'first_author', 'year', 'journal', 'title']:
            if key not in r:
                raise ValueError(f"Reference PMID:{r.get('pmid','?')} missing '{key}' (needed for report)")

    return evidence
```
