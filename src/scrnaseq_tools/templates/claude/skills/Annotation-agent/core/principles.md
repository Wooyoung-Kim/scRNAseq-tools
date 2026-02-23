# Core Principles (v3)

Critical rules with MANDATORY functional analysis enforcement + DATA-DRIVEN outlier detection.

```
╔══════════════════════════════════════════════════════════════════════╗
║  🚨 v3 CHANGES:                                                       ║
║  - NO hardcoded TF/Pathway associations                              ║
║  - Statistical outlier detection (z-score > 2.5)                     ║
║  - LLM-guided literature search for outliers                         ║
║  - Confidence penalties for unexplained outliers                     ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## 1. Marker Criteria (MANDATORY)

All markers must meet ALL criteria:

| Criterion | Threshold |
|-----------|-----------|
| `pct_in` | >= 25% |
| `log2FC` | >= 1 |
| `padj` | < 0.05 |

```python
valid_markers = de_results[
    (de_results['pct_in'] >= 0.25) &
    (de_results['log2fc'] >= 1.0) &
    (de_results['padj'] < 0.05)
].head(50)
```

---

## 2. DE Re-computation Rule (HARD FAIL)

```
╔══════════════════════════════════════════════════════════════════════╗
║  HARD FAIL: Using markers from previous tier's DE table              ║
║  Each tier MUST have unique subset_id and fresh DE computation       ║
╚══════════════════════════════════════════════════════════════════════╝
```

| Tier | Data Subset | subset_id format |
|------|-------------|------------------|
| Tier 1 | Full dataset | `full_dataset` |
| Tier 2 | Per major type | `{major_type}` |
| Tier 3 | Per dev state | `{major_type}_{dev_state}` |

---

## 3. Functional Analysis (MANDATORY - v3)

```
╔══════════════════════════════════════════════════════════════════════╗
║  ⚠️ v3 ENFORCEMENT: Cannot annotate without functional analysis       ║
║                                                                      ║
║  Tier 2: TF Activity (MANDATORY) + Trajectory (if >= 2000 cells)     ║
║         + Statistical Outlier Detection (z-score)                    ║
║  Tier 3: Pathway Activity (MANDATORY)                                ║
║         + Statistical Outlier Detection (z-score)                    ║
║                                                                      ║
║  🆕 v3: Outlier detection MANDATORY before reasoning                 ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Verification Code

```python
def verify_functional_analysis(adata, tier):
    """
    MUST call before annotation.
    Raises AssertionError if functional analysis missing.
    """
    errors = []

    if tier == 2:
        if 'ulm_estimate' not in adata.obsm:
            errors.append("❌ TF activity (ulm_estimate) not computed")
        if 'ulm_pvals' not in adata.obsm:
            errors.append("❌ TF p-values (ulm_pvals) not computed")
        if adata.n_obs >= 2000 and 'pseudotime' not in adata.obs.columns:
            errors.append("❌ Trajectory (pseudotime) not computed for >= 2000 cells")

    elif tier == 3:
        if 'mlm_estimate' not in adata.obsm:
            errors.append("❌ Pathway activity (mlm_estimate) not computed")
        if 'mlm_pvals' not in adata.obsm:
            errors.append("❌ Pathway p-values (mlm_pvals) not computed")

    if errors:
        raise AssertionError(
            f"\n{'='*60}\n"
            f"FUNCTIONAL ANALYSIS INCOMPLETE (Tier {tier})\n"
            f"{'='*60}\n" +
            "\n".join(errors) +
            f"\n\n⚠️ CANNOT PROCEED TO ANNOTATION!\n"
            f"Run functional analysis first.\n"
            f"{'='*60}"
        )

    print(f"✅ Functional analysis verified for Tier {tier}")
    return True
```

---

## 4. PMID Requirement (MANDATORY)

```
╔══════════════════════════════════════════════════════════════════════╗
║  EVERY literature reference MUST include PMID                        ║
║  References WITHOUT PMID are INVALID and will be REJECTED            ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Search Methods (in order)

1. MCP `pubmed_search` tool
2. Bio.Entrez (Biopython)
3. WebSearch with "PMID" query

### Reference Format

```json
{
  "pmid": "12345678",
  "query": "CCR7 TCF7 naive T cell",
  "title": "...",
  "journal": "...",
  "year": "...",
  "first_verify": "VERIFIED",
  "second_verify": "VERIFIED",
  "status": "DOUBLE_VERIFIED"
}
```

---

## 5. DATA-DRIVEN MARKER ENFORCEMENT (ABSOLUTE RULE)

```
╔══════════════════════════════════════════════════════════════════════╗
║  NEVER HARDCODE CANONICAL MARKERS                                     ║
║  ALL markers MUST be derived from DE results of the CURRENT dataset   ║
║                                                                      ║
║  Canonical markers from textbooks/papers may NOT match your data:     ║
║  - A "known" marker may be anti-enriched in your dataset              ║
║  - Enrichment patterns differ between species/tissues/conditions      ║
║  - Always verify before using                                        ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Anti-Enrichment Detection (MANDATORY)

Before using ANY marker for annotation or visualization:

```python
def validate_marker_enrichment(de_df, gene, target_group, groupby_col):
    """
    Validate that a marker is ACTUALLY enriched in the target group.
    MUST run before using any marker — even "canonical" ones.

    Returns:
    --------
    dict with pct, enrichment_ratio, is_valid
    """
    row = de_df[(de_df['gene'] == gene) & (de_df['group'] == target_group)]
    if row.empty:
        return {'gene': gene, 'pct': 0, 'enrichment': 0, 'is_valid': False,
                'reason': 'Gene not found in DE results'}

    pct_in = row['pct_in'].values[0]
    pct_out = row['pct_out'].values[0]
    enrichment = pct_in / max(pct_out, 0.001)

    is_valid = (pct_in >= 0.25) and (enrichment >= 1.5)

    return {
        'gene': gene,
        'pct': pct_in,
        'enrichment': round(enrichment, 2),
        'is_valid': is_valid,
        'reason': 'ANTI-ENRICHED' if enrichment < 1.0 else
                  'LOW_PCT' if pct_in < 0.25 else
                  'LOW_ENRICHMENT' if enrichment < 1.5 else 'VALID'
    }
```

### Rules

| Rule | Description |
|------|-------------|
| **No hardcoded marker lists** | Every marker must come from `rank_genes_groups` or equivalent DE |
| **No hardcoded TF/Pathway associations** | ❌ NEVER create dictionaries mapping TF/Pathway → cell type |
| **Anti-enrichment check** | If enrichment < 1.0, the marker is CONTRA-INDICATED for that cell type |
| **Minimum pct** | pct_in >= 25% in target group |
| **Minimum enrichment** | enrichment >= 1.5x vs other groups |
| **Dotplot markers** | Only include markers passing ALL criteria; omit bracket if none pass |
| **Cross-species caution** | Canonical human markers may not apply to ferret/mouse/other species |

### ❌ FORBIDDEN: Hardcoded TF/Pathway Associations

```python
# ❌ NEVER DO THIS
TF_ASSOCIATIONS = {
    'TBX21': 'Th1/Effector',  # Misses ABC (age-associated B cells)!
    'TCF7': 'Naive/Memory',
}

PATHWAY_ASSOCIATIONS = {
    'WNT': 'Stemness',  # Misses Tpex (precursor exhausted)!
    'TNFa': 'Inflammatory',
}

# ✅ INSTEAD: Use statistical z-scores + LLM literature search
z_score = (cluster_value - global_mean) / global_std
if abs(z_score) > 2.5:
    prompt = f"Search '{TF_name} {cell_lineage}' in literature"
```

---

## 5.5 Statistical Outlier Detection (NEW - v3)

```
╔══════════════════════════════════════════════════════════════════════╗
║  MANDATORY: Compute z-scores for TF/Pathway activity before reasoning║
║  Outliers (|z| > 2.5) trigger specialized literature searches        ║
╚══════════════════════════════════════════════════════════════════════╝
```

### Z-score Calculation

```python
def detect_outliers(adata, activity_key, cluster_col):
    """
    Detect statistical outliers in TF/Pathway activity.
    NO hardcoding - pure statistics.
    """
    scores = adata.obsm[activity_key]

    # Global statistics
    global_mean = scores.mean()
    global_std = scores.std()

    outliers = []
    for cluster in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == cluster
        cluster_mean = scores[mask].mean()

        for feature in scores.columns:
            # Z-score: how many standard deviations from mean?
            z = (cluster_mean[feature] - global_mean[feature]) / (global_std[feature] + 1e-6)

            # Flag if |z| > 2.5 (statistical outlier)
            if abs(z) > 2.5 and abs(cluster_mean[feature]) > 0.5:
                outliers.append({
                    'cluster': cluster,
                    'feature': feature,
                    'z_score': z,
                    'activity': cluster_mean[feature]
                })

    return outliers
```

### Outlier Interpretation Rules

| Outlier Count | Action |
|---------------|--------|
| **0** | Standard annotation workflow |
| **1-2** | Note outliers, investigate if marker pattern also unusual |
| **>= 3** | **HIGH PRIORITY** - specialized/rare cell type likely |

### Reasoning Workflow for Outliers

```
IF outliers detected:
  1. Generate reasoning prompt:
     "⚠️ OUTLIERS: {feature} (z={score:.2f})"
     "→ Search '{feature} {parent_context}'"

  2. LLM performs literature search:
     - "{outlier_feature} {cell_lineage}"
     - "{outlier_feature} {marker1} {marker2}"
     - Focus: "age-associated", "atypical", "specialized"

  3. Decision:
     IF literature explains outlier:
       → Use specialized annotation (e.g., ABC, Tpex)
       → FULL confidence (12 pts)
     ELSE IF standard annotation chosen:
       → CONFIDENCE PENALTY (max 6-7 pts)
       → Unexplained outliers = insufficient evidence
```

---

## 6. Marker Combination Strategy

```
ALWAYS use 2-5 markers per cell type (NEVER single marker)
```

| Good | Bad |
|------|-----|
| 2-5 markers | Single marker |
| Literature-supported | No support |
| Functionally coherent | Contradictory |
| TF/Pathway aligned | TF/Pathway conflicts |

---

## 7. Integrated Reasoning (v3)

```
Evidence Sources:
├── DE Markers (Top 50)
├── TF Activity (Tier 2)
├── Trajectory Position (Tier 2, if computed)
├── Pathway Activity (Tier 3)
├── 🆕 Statistical Outliers (TF/Pathway z-scores)
└── 🆕 Evidence Conflicts (automated detection)

All sources MUST be considered in reasoning.
Cannot rely on markers alone.
Outliers MUST trigger specialized literature searches.
```

### Evidence Weight

| Source | Weight in Decision | Notes |
|--------|-------------------|-------|
| Markers | 20% | Primary evidence |
| TF/Pathway | 20% | Functional context |
| **Outliers** | **20%** | 🆕 Signals specialized states |
| Trajectory | 20% (Tier 2) | Developmental context |
| Literature | 20% | Final validation |

### Reasoning Priority (v3)

```
1. Check for statistical outliers (|z| > 2.5)
   IF outliers detected → PRIORITY literature search

2. Check for evidence conflicts
   IF conflicts detected → Investigate specialized states

3. Standard annotation workflow
   IF no outliers/conflicts → Normal processing

4. Confidence scoring
   Unexplained outliers → PENALTY (max 6-7 pts)
   Literature-explained outliers → FULL confidence (12 pts)
```

---

## 8. Confidence Scoring (v3 - Outlier-Aware)

### Tier 2 - Standard (NO outliers)

| Criteria | Points |
|----------|--------|
| Markers (>= 3) | 3 |
| References (>= 2 DOUBLE_VERIFIED) | 3 |
| TF consistency | 3 |
| Trajectory alignment | 3 |
| **Max** | **12** |

### Tier 2 - WITH Outliers

| Criteria | Standard Annotation | Specialized Annotation |
|----------|---------------------|------------------------|
| Markers (>= 3) | 3 | 3 |
| References (>= 2 DOUBLE_VERIFIED) | 3 | **3 (REQUIRED)** |
| TF consistency | **0-1** ⚠️ (outlier unexplained) | **3** (if lit. explains outlier) |
| Trajectory alignment | 3 | 3 OR N/A |
| **Max** | **7-10** ⚠️ | **12** ✅ |

### Tier 3 - Standard (NO outliers)

| Criteria | Points |
|----------|--------|
| Markers (>= 3) | 3 |
| References (>= 2 DOUBLE_VERIFIED) | 3 |
| Pathway consistency | 3 |
| Literature support | 3 |
| **Max** | **12** |

### Tier 3 - WITH Outliers

| Criteria | Standard Annotation | Specialized Annotation |
|----------|---------------------|------------------------|
| Markers (>= 3) | 3 | 3 |
| References (>= 2 DOUBLE_VERIFIED) | 3 | **3 (REQUIRED)** |
| Pathway consistency | **0-1** ⚠️ (outlier unexplained) | **3** (if lit. explains outlier) |
| Conflict resolution | **0** ⚠️ | **3** (if lit. resolves conflict) |
| **Max** | **6-7** ⚠️ | **12** ✅ |

### Thresholds

| Score | Confidence | Notes |
|-------|------------|-------|
| >= 10 | HIGH | ✅ Well-supported |
| 7-9 | MEDIUM | ⚠️ Some uncertainty |
| 4-6 | LOW | ⚠️ **Typical for unexplained outliers** |
| < 4 | INSUFFICIENT | ❌ Do not annotate |

### Key Rules (v3)

```
RULE 1: Outliers present + standard annotation chosen
        → MANDATORY confidence penalty (max 6-7 pts)

RULE 2: Outliers present + literature explains outliers
        → FULL confidence possible (12 pts)

RULE 3: Multiple outliers (>= 3) unexplained
        → SEVERE penalty (max 6 pts = LOW confidence)

RULE 4: Forces LLM to search literature for outliers
        → Cannot ignore anomalies
```

---

## 9. NEVER DO (Safeguards - v3 Updated)

### 1. Skip functional analysis

```
BAD:  Annotate without TF/Pathway analysis
GOOD: Compute TF/Pathway → Verify → Then annotate
```

### 2. Use upstream markers

```
BAD:  CD3D (Tier 1) for T cell subset distinction
GOOD: Re-compute DE within T cell subset
```

### 3. Annotate with single marker

```
BAD:  "GZMB high = Cytotoxic"
GOOD: "GZMB + PRF1 + NKG7 + TNFa pathway = Cytotoxic"
```

### 4. Ignore TF/Pathway conflicts

```
BAD:  Markers say Effector, TF says Naive → Force Effector
GOOD: Note conflict, reduce confidence, investigate
```

### 5. Skip PMID verification

```
BAD:  Literature supports → Use without PMID
GOOD: Search PMID → Verify → Double-verify → DOUBLE_VERIFIED
```

### 6. Hardcode canonical markers

```
BAD:  markers = {'CellType': ['GENE_A', 'GENE_B']}  # Textbook markers, not verified in data
GOOD: markers = get_valid_markers_from_de(de_results, 'CellType')  # Data-driven from DE
```

### 7. Ignore anti-enrichment

```
BAD:  GENE_X is a "known" marker for CellType → Use it (even if pct < 5%, enrichment < 1.0x)
GOOD: GENE_X pct < 5%, enrichment < 1.0x → EXCLUDE, find data-supported alternatives from DE
```

### 8. Skip contamination check

```
BAD:  All clusters must be the expected lineage
GOOD: Check top DE markers for cross-lineage signatures (e.g., other lineage's defining genes)
      → Flag and remove contamination clusters before final annotation
```

### 🆕 9. Hardcode TF/Pathway associations (v3)

```
❌ ABSOLUTELY FORBIDDEN:
TF_ASSOCIATIONS = {
    'TBX21': 'Th1/Effector',  # Misses ABC in B cells!
    'TCF7': 'Naive/Memory',
}

PATHWAY_ASSOCIATIONS = {
    'WNT': 'Stemness',  # Misses Tpex in Effector context!
    'TNFa': 'Inflammatory',
}

✅ CORRECT APPROACH:
# Compute z-scores
z = (cluster_value - global_mean) / global_std

# Generate reasoning prompt
if abs(z) > 2.5:
    prompt = f"⚠️ OUTLIER: {feature} (z={z:.2f})\n"
    prompt += f"→ Search '{feature} {context}' in literature"

# Let LLM + literature determine meaning
```

**Why this is critical:**
- TBX21 in T cells ≠ TBX21 in B cells
- WNT in Naive ≠ WNT in Effector
- Context-dependent interpretation REQUIRES literature search

### 🆕 10. Ignore statistical outliers (v3)

```
BAD:  TBX21 z-score = 3.8 in B cells → Ignore, annotate as Memory_B
      (Outlier unexplained → Confidence penalty)

GOOD: TBX21 z-score = 3.8 → Search "TBX21 B cell"
      → Find PMID:25006127 (ABC)
      → Annotate as ABC (Outlier explained → Full confidence)
```

### 🆕 11. Use standard annotation with unexplained outliers (v3)

```
BAD:  Cluster has 3+ outliers
      → Annotate with standard name anyway
      → Confidence = 6-7 pts (LOW)

GOOD: Cluster has 3+ outliers
      → Prioritize specialized literature searches
      → Find specialized annotation OR flag as Novel
      → Confidence = 10-12 pts (HIGH)
```

---

## 10. Validation Checkpoints (v3 - Outlier-Aware)

### Before Tier 2 Annotation

```
[ ] DE computed within major type
[ ] TF activity computed (ulm_estimate)
[ ] Trajectory computed if >= 2000 cells
[ ] 🆕 Statistical outliers detected (tf_outliers in uns)
[ ] 🆕 Evidence conflicts detected (evidence_conflicts in uns)
[ ] verify_functional_analysis(adata, tier=2) PASSED
```

### Before Tier 3 Annotation

```
[ ] DE RE-computed within developmental state
[ ] Pathway activity computed (mlm_estimate)
[ ] 🆕 Statistical outliers detected (pathway_outliers in uns)
[ ] 🆕 Evidence conflicts detected (evidence_conflicts in uns)
[ ] verify_functional_analysis(adata, tier=3) PASSED
```

### Before Final Assignment

```
[ ] All evidence types considered (markers + TF/pathway + trajectory + 🆕 outliers)
[ ] 🆕 IF outliers present: Specialized literature searches performed
[ ] 🆕 IF outliers unexplained: Confidence penalty applied OR novel flagged
[ ] Confidence scored (max 12, penalties for unexplained outliers)
[ ] >= 2 references DOUBLE_VERIFIED
[ ] Novel populations flagged if applicable
[ ] Visualization saved with functional analysis
```

### 🆕 Outlier Handling Checklist

```
IF outliers detected (|z| > 2.5):
  [ ] Count outliers: N = {count}
  [ ] Generate reasoning prompts for each outlier
  [ ] Perform literature searches: "{outlier} {context}"
  [ ] Evaluate results:
      [ ] Literature explains outlier → Use specialized annotation
      [ ] No literature + coherent pattern → Flag as Novel
      [ ] No literature + incoherent → Reduce confidence, note for review
  [ ] Apply confidence scoring:
      [ ] Explained outliers → Full confidence (12 pts)
      [ ] Unexplained outliers → Penalty (6-10 pts)
```

---

## 11. Output Requirements (v3 - Outlier Information)

### Tier 2 Output

```python
# Required in saved .h5ad
assert 'tier2_annotation' in adata.obs.columns
assert 'ulm_estimate' in adata.obsm  # TF activity
assert 'ulm_pvals' in adata.obsm
if adata.n_obs >= 2000:
    assert 'pseudotime' in adata.obs.columns

# 🆕 v3: Outlier and conflict information
assert 'tf_outliers' in adata.uns  # Statistical outliers
assert 'evidence_conflicts' in adata.uns  # Conflict detection results

# Optional: Per-cluster metadata
# adata.uns['cluster_metadata'] = {
#     'cluster_0': {'has_outliers': True, 'outlier_count': 2, 'annotation_rationale': '...'},
# }
```

### Tier 3 Output

```python
assert 'tier3_annotation' in adata.obs.columns
assert 'final_annotation' in adata.obs.columns
assert 'mlm_estimate' in adata.obsm  # Pathway activity
assert 'mlm_pvals' in adata.obsm

# 🆕 v3: Outlier and conflict information
assert 'pathway_outliers' in adata.uns  # Statistical outliers
assert 'evidence_conflicts' in adata.uns  # Conflict detection results
```

### 🆕 Outlier Metadata Format

```python
# adata.uns['tf_outliers'] format
[
    {
        'cluster': '3',
        'tf': 'TBX21',
        'activity': 1.8,
        'z_score': 3.8,
        'pval': 0.0001,
        'direction': 'HIGH'
    },
    ...
]

# adata.uns['evidence_conflicts'] format
[
    {
        'cluster': '3',
        'type': 'TF_STATISTICAL_OUTLIER',
        'detail': 'TBX21 unusually HIGH',
        'severity': 'HIGH',
        'action': 'Search "TBX21 B lineage"'
    },
    ...
]
```

---

## 12. Summary: v2 → v3 Changes

### What Changed

| Aspect | v2 | v3 |
|--------|----|----|
| **TF/Pathway interpretation** | Hardcoded associations | Statistical z-scores + LLM search |
| **Outlier handling** | None | Mandatory detection + reasoning |
| **Confidence scoring** | Fixed (12 pts max) | Dynamic (6-12 pts based on outliers) |
| **Literature search** | Standard markers only | Priority for outliers |
| **Novel detection** | Manual flagging | Automatic via outlier detection |

### Migration Checklist

```
For existing workflows using v2:

[ ] Remove all TF_ASSOCIATIONS dictionaries
[ ] Remove all PATHWAY_ASSOCIATIONS dictionaries
[ ] Add outlier detection steps:
    - compute_tf_statistics() in Tier 2
    - compute_pathway_statistics() in Tier 3
[ ] Add outlier-aware evidence collection
[ ] Update confidence scoring logic
[ ] Update reasoning prompts to include outliers
[ ] Test with known specialized subsets (ABC, Tpex, etc.)
```
