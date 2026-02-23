# Integrated Reasoning Format (v2)

**Key Change**: Markers + TF/Pathway + Trajectory = Combined Evidence

---

## Overview

```
╔══════════════════════════════════════════════════════════════════════╗
║  v2 Reasoning: MUST include functional analysis in decision-making   ║
║                                                                      ║
║  Tier 2: Markers + TF Activity + Trajectory → Developmental State    ║
║  Tier 3: Markers + Pathway Activity → Functional State               ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Iteration Structure

```
Iteration 1: Evidence Collection (Markers + Functional)
Iteration 2: Candidate Reasoning + Literature
Iteration 3: Integrated Decision + Verification
```

---

## Tier 2 Reasoning Template

### Iteration 1: Evidence Collection

```markdown
**Bioinformatician**:

## Cluster {X} Evidence Summary

### DE Markers (Top 20 Valid)
| Rank | Gene | pct_in | log2FC | padj |
|------|------|--------|--------|------|
| 1 | {gene} | {val} | {val} | {val} |
| 2 | {gene} | {val} | {val} | {val} |
...

### TF Activity (Top 10)
| TF | Activity | p-value | Known Association |
|----|----------|---------|-------------------|
| {TF1} | {score} | {pval} | {known_function} |
| {TF2} | {score} | {pval} | {known_function} |
| {TF3} | {score} | {pval} | {known_function} |
...

### Trajectory Position
- Mean pseudotime: {value}
- Category: **{EARLY/MID/LATE}** (0-0.33 / 0.33-0.67 / 0.67-1.0)
- Interpretation: {interpretation based on pseudotime}

### Evidence Summary
- Marker signature: [{marker1}, {marker2}, {marker3}] → {candidate_type}
- TF activity: {TF1}+, {TF2}+ → {supported_state}
- Trajectory: {category} → {supported_state}
- **Convergent evidence**: {summary of convergence or conflict}
```

### Iteration 2: Candidate Reasoning + Literature

```markdown
**Computational Biologist**:

## Candidate Analysis

### Candidate 1: {type_name}
- **Marker support**: {marker1} (pct={val}%), {marker2} (pct={val}%), {marker3} (pct={val}%)
- **TF support**: {TF1} active ({score}), {TF2} active ({score})
- **Trajectory support**: {category} pseudotime ({value})
- **Confidence**: {HIGH/MEDIUM/LOW} ({rationale})

### Candidate 2: {type_name}
- **Marker support**: {description of marker evidence}
- **TF support**: {description of TF evidence}
- **Trajectory support**: {matches/conflicts} - {reason}
- **Confidence**: {HIGH/MEDIUM/LOW} ({rationale})

### Candidate 3: {type_name}
- **Marker support**: {description of marker evidence}
- **TF support**: {description of TF evidence}
- **Trajectory support**: {matches/conflicts}
- **Confidence**: {HIGH/MEDIUM/LOW} ({rationale})

## PRE-REASONING LITERATURE (PMID REQUIRED)

Query 1: "{marker_combo} {cell_type_context}"
- **PMID: {number}** - "{key finding}" ({journal}, {year})
- **PMID: {number}** - "{key finding}" ({journal}, {year})

Query 2: "{TF_combo} {cell_type_context}"
- **PMID: {number}** - "{key finding}" ({journal}, {year})
```

### Iteration 3: Integrated Decision

```markdown
**PI - Final Decision**:

## Evidence Integration Table

| Evidence Type | Candidate 1 ({name}) | Candidate 2 ({name}) | Candidate 3 ({name}) |
|---------------|---------------------|----------------------|-------------------------|
| Markers | ✅/⚠️/❌ {marker_combo} | ✅/⚠️/❌ {description} | ✅/⚠️/❌ {description} |
| TF Activity | ✅/⚠️/❌ {TF_support} | ✅/⚠️/❌ {TF_support} | ✅/⚠️/❌ {TF_support} |
| Trajectory | ✅/⚠️/❌ {category} ({value}) | ✅/⚠️/❌ {description} | ✅/⚠️/❌ {description} |
| Literature | ✅ {N} PMIDs | - | - |

## Final Assignment

**Annotation**: {selected_candidate}
**Confidence**: {HIGH/MEDIUM/LOW} ({score}/12 points)
- Markers: {0-3} pts ({N} markers meeting criteria)
- Refs: {0-3} pts ({N} double-verified)
- TF consistency: {0-3} pts ({TFs} match)
- Trajectory alignment: {0-3} pts ({category} matches)

**Key Evidence**:
1. Marker combination {marker1}+{marker2}+{marker3} supports {annotation}
2. {TF1} and {TF2} TF activity confirms {transcriptional_program}
3. {EARLY/MID/LATE} pseudotime ({value}) consistent with {state}
4. Literature verified (PMID: {number}, {number})

**References**:
- PMID:{number} - DOUBLE_VERIFIED
- PMID:{number} - DOUBLE_VERIFIED
```

---

## Tier 3 Reasoning Template

### Iteration 1: Evidence Collection

```markdown
**Bioinformatician**:

## Cluster {X} Evidence Summary (Parent: {Major}_{Dev})

### DE Markers (Top 20 Valid)
| Rank | Gene | pct_in | log2FC | padj |
|------|------|--------|--------|------|
| 1 | {gene} | {val} | {val} | {val} |
| 2 | {gene} | {val} | {val} | {val} |
| 3 | {gene} | {val} | {val} | {val} |
...

### Pathway Activity (Top 10)
| Pathway | Activity | p-value | Interpretation |
|---------|----------|---------|----------------|
| {pathway1} | {score} | {pval} | {interpretation} |
| {pathway2} | {score} | {pval} | {interpretation} |
| {pathway3} | {score} | {pval} | {interpretation} |
| {pathway4} | {score} | {pval} | {interpretation} |
...

### Parent Context
- Tier 1: {major_type}
- Tier 2: {dev_state}
- Tissue: {tissue_context}

### Evidence Summary
- Marker signature: [{marker1}, {marker2}, {marker3}] → {candidate_type}
- Pathway activity: {pathway1}+, {pathway2}+ → {functional_interpretation}
- **Convergent evidence**: {summary of convergence or conflict}
```

### Iteration 2: Candidate Reasoning + Literature

```markdown
**Computational Biologist**:

## Candidate Analysis

### Candidate 1: {func_state_name}
- **Marker support**: {marker1} ({pct}%), {marker2} ({pct}%), {marker3} ({pct}%)
- **Pathway support**: {pathway1} high, {pathway2} high → {interpretation}
- **Confidence**: {HIGH/MEDIUM/LOW}

### Candidate 2: {func_state_name}
- **Marker support**: {description of marker evidence}
- **Pathway support**: {description of pathway evidence}
- **Confidence**: {HIGH/MEDIUM/LOW}

### Candidate 3: {func_state_name}
- **Marker support**: {description of marker evidence}
- **Pathway support**: {description of pathway evidence}
- **Confidence**: {HIGH/MEDIUM/LOW}

## PRE-REASONING LITERATURE (PMID REQUIRED)

Query 1: "{marker_combo} {cell_type_context}"
- **PMID: {number}** - "{key finding}" ({journal}, {year})

Query 2: "{pathway} {func_state} {cell_type_context}"
- **PMID: {number}** - "{key finding}" ({journal}, {year})
```

### Iteration 3: Integrated Decision

```markdown
**PI - Final Decision**:

## Evidence Integration Table

| Evidence Type | Candidate 1 ({name}) | Candidate 2 ({name}) | Candidate 3 ({name}) |
|---------------|-------------------------|-------------------------|-------------------|
| Markers | ✅/⚠️/❌ {marker_combo} | ✅/⚠️/❌ {description} | ✅/⚠️/❌ {description} |
| Pathway | ✅/⚠️/❌ {pathway_support} | ✅/⚠️/❌ {description} | ✅/⚠️/❌ {description} |
| Literature | ✅ {N} PMIDs | - | - |

## Final Assignment

**Annotation**: {selected_candidate}
**Full Label**: {major_type}_{dev_state}_{func_state}
**Confidence**: {HIGH/MEDIUM/LOW} ({score}/12 points)
- Markers: {0-3} pts ({N} markers meeting criteria)
- Refs: {0-3} pts ({N} double-verified)
- Pathway consistency: {0-3} pts ({pathways} match)
- Literature: {0-3} pts

**Key Evidence**:
1. Marker combination {marker1}+{marker2}+{marker3} supports {annotation}
2. {pathway1} and {pathway2} pathway activity indicates {functional_state}
3. Absence of {negative_markers} rules out {alternative_candidate}
4. Literature verified (PMID: {number}, {number})

**References**:
- PMID:{number} - DOUBLE_VERIFIED
- PMID:{number} - VERIFIED
```

---

## Evidence Integration Scoring

### Scoring Matrix (v2)

| Criteria | High (3 pts) | Medium (2 pts) | Low (1 pt) | None (0 pts) |
|----------|--------------|----------------|------------|--------------|
| **Markers** | >= 3 meeting criteria | 2 markers | 1 marker | 0 |
| **References** | >= 2 DOUBLE_VERIFIED | 1 DOUBLE or 2 VERIFIED | 1 VERIFIED | None |
| **TF/Pathway** | Strong match (>= 2 TFs/pathways support) | Partial (1 supports) | Weak/conflicting | Not computed |
| **Trajectory** | Matches category | Borderline | Conflicts | Not computed |

### Confidence Thresholds

| Total Score | Confidence Level | Action |
|-------------|------------------|--------|
| **>= 10** | HIGH | Proceed with annotation |
| **7-9** | MEDIUM | Proceed with caution, note limitations |
| **4-6** | LOW | Additional evidence needed, consider flagging |
| **< 4** | INSUFFICIENT | Cannot annotate, flag for review |

---

## Evidence.json Schema (v2)

```json
{
  "cluster_id": "tier2_{major_type}_{cluster_num}",
  "tier": 2,
  "subset_id": "{major_type}",

  "annotation": "{dev_state}",
  "full_label": "{major_type}_{dev_state}",

  "markers": [
    {"gene": "{marker1}", "pct_in": "{val}", "log2fc": "{val}", "padj": "{val}"},
    {"gene": "{marker2}", "pct_in": "{val}", "log2fc": "{val}", "padj": "{val}"},
    {"gene": "{marker3}", "pct_in": "{val}", "log2fc": "{val}", "padj": "{val}"}
  ],

  "tf_activity": [
    {"tf": "{TF1}", "activity": "{score}", "pval": "{pval}", "supports": "{dev_state}"},
    {"tf": "{TF2}", "activity": "{score}", "pval": "{pval}", "supports": "{dev_state}"}
  ],

  "trajectory": {
    "mean_pseudotime": "{value}",
    "category": "{Early/Mid/Late}",
    "interpretation": "{interpretation}"
  },

  "pathway_activity": null,

  "confidence_score": "{0-12}",
  "confidence_level": "{High/Medium/Low}",

  "scoring_breakdown": {
    "markers": "{0-3}",
    "references": "{0-3}",
    "tf_pathway_consistency": "{0-3}",
    "trajectory_alignment": "{0-3}"
  },

  "candidates": [
    {
      "name": "{candidate1}",
      "markers": ["{marker1}", "{marker2}", "{marker3}"],
      "tf_support": ["{TF1}", "{TF2}"],
      "trajectory_support": true,
      "selected": true
    },
    {
      "name": "{candidate2}",
      "markers": ["{marker}"],
      "tf_support": ["{TF}"],
      "trajectory_support": false,
      "selected": false
    }
  ],

  "references": [
    {
      "pmid": "{number}",
      "query": "{marker_combo} {cell_type_context}",
      "status": "DOUBLE_VERIFIED"
    },
    {
      "pmid": "{number}",
      "query": "{TF_combo} {cell_type_context}",
      "status": "DOUBLE_VERIFIED"
    }
  ],

  "is_novel": false,
  "novel_evidence": null
}
```

---

## Quick Reference: Evidence Sources

### Tier 2 Evidence

| Source | Location | Key Columns |
|--------|----------|-------------|
| DE Markers | `rank_genes_groups` | names, pct_nz, logfoldchanges, pvals_adj |
| TF Activity | `obsm['ulm_estimate']` | TF names, activity scores |
| TF p-values | `obsm['ulm_pvals']` | TF names, p-values |
| Pseudotime | `obs['pseudotime']` | 0-1 scale |

### Tier 3 Evidence

| Source | Location | Key Columns |
|--------|----------|-------------|
| DE Markers | `rank_genes_groups` (re-computed!) | names, pct_nz, logfoldchanges, pvals_adj |
| Pathway Activity | `obsm['mlm_estimate']` | Pathway names, activity scores |
| Pathway p-values | `obsm['mlm_pvals']` | Pathway names, p-values |
