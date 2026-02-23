# Decision Matrices (v2)

Confidence scoring with integrated functional analysis.

---

## 1. Confidence Level Assignment (v2 Updated)

### Scoring Table

| Criteria | High (3 pts) | Medium (2 pts) | Low (1 pt) | None (0 pts) |
|----------|--------------|----------------|------------|--------------|
| **Markers meeting criteria** | >= 3 | 2 | 1 | 0 |
| **Double-verified refs** | >= 2 | 1 | 0 | - |
| **TF/Pathway consistency** | Strong (>= 2 support) | Partial (1 supports) | Weak/Conflicting | Not computed |
| **Trajectory alignment** (Tier 2) | Matches | Borderline | Conflicts | Not computed |
| **Literature support** | Multiple PMIDs | 1 PMID | Indirect | None |

### Tier 2 Max Score: 12 points

| Component | Max Points |
|-----------|------------|
| Markers | 3 |
| References | 3 |
| TF consistency | 3 |
| Trajectory alignment | 3 |
| **Total** | **12** |

### Tier 3 Max Score: 12 points

| Component | Max Points |
|-----------|------------|
| Markers | 3 |
| References | 3 |
| Pathway consistency | 3 |
| Literature support | 3 |
| **Total** | **12** |

### Final Confidence

| Total Score | Confidence | Action |
|-------------|------------|--------|
| >= 10 | **High** | Proceed |
| 7-9 | **Medium** | Proceed with notes |
| 4-6 | **Low** | Additional evidence needed |
| < 4 | **Insufficient** | Cannot annotate |

---

## 2. TF/Pathway Consistency Scoring

### TF Consistency (Tier 2)

| Score | Criteria |
|-------|----------|
| **3 (Strong)** | >= 2 TFs active that support the candidate cell type |
| **2 (Partial)** | 1 TF supports, or 2 TFs with one conflicting |
| **1 (Weak)** | TFs present but don't clearly support candidate |
| **0 (None)** | No TF activity computed OR all TFs conflict |

**Example:**
- {TF1} active ({score}) + {TF2} active ({score}), both support {candidate} → 3 pts (Strong)
- {TF1} active only → 2 pts (Partial)
- {TF3} active (suggests {alternative}, not {candidate}) → 1 pt (Weak/Conflicting)

### Pathway Consistency (Tier 3)

| Score | Criteria |
|-------|----------|
| **3 (Strong)** | >= 2 pathways active that match functional state |
| **2 (Partial)** | 1 pathway supports, or mixed signals |
| **1 (Weak)** | Pathways present but don't match candidate |
| **0 (None)** | No pathway activity computed OR all pathways conflict |

**Example:**
- {pathway1} high + {pathway2} high + {pathway3} high, all support {func_state} → 3 pts (Strong)
- {pathway1} high only → 2 pts (Partial)
- {pathway4} high (suggests {alternative}, not {func_state}) → 1 pt (Weak)

---

## 3. Trajectory Alignment Scoring (Tier 2 Only)

> **⚠️ Pseudotime-to-cell type mapping is LINEAGE-SPECIFIC.**
> Do NOT assume a universal mapping. Determine expected positions from the data.

### Scoring Criteria

| Score | Criteria |
|-------|----------|
| **3 (Matches)** | Candidate's expected developmental stage matches observed pseudotime position |
| **2 (Borderline)** | Near category cutoff, or candidate could span two categories |
| **1 (Conflicts)** | Candidate's expected stage directly contradicts observed pseudotime |
| **0 (N/A)** | Pseudotime not computed |

### How to Score (Data-Driven)

1. **Determine candidate** from markers + TF evidence (NOT from pseudotime)
2. **Ask**: Is this candidate expected to be early, mid, or late in its lineage's differentiation?
   - Use literature / prior knowledge about the candidate's developmental position
   - Do NOT use a fixed lookup table — differentiation order varies by lineage and dataset
3. **Compare**: Does the observed mean pseudotime match that expectation?

### Pseudotime Categories

| Category | Range | Interpretation |
|----------|-------|----------------|
| Early | < 0.33 | Progenitor / undifferentiated states |
| Mid | 0.33 - 0.67 | Intermediate / transitional states |
| Late | > 0.67 | Terminally differentiated / effector states |

### Example
- Candidate: {cell_type} (from marker/TF evidence)
- Literature says {cell_type} is a {early/mid/late}-stage cell in this lineage
- Observed pseudotime: {value} ({category})
- Match → 3 pts / Borderline → 2 pts / Conflict → 1 pt

---

## 4. Quick Decision Trees

### Tier 2 Decision Tree

```
Are >= 3 markers meeting criteria?
├── NO → LOW confidence, need more markers
└── YES → Continue
        │
        v
    Is TF activity computed?
    ├── NO → MEDIUM confidence max, note missing TF
    └── YES → Continue
            │
            v
        Do >= 2 TFs support candidate?
        ├── NO → MEDIUM confidence
        └── YES → Continue
                │
                v
            Is trajectory computed (>= 2000 cells)?
            ├── NO (< 2000 cells) → Score without trajectory
            └── YES → Continue
                    │
                    v
                Does trajectory match candidate?
                ├── NO → MEDIUM confidence (TF vs Trajectory conflict)
                └── YES → Continue
                        │
                        v
                    Are >= 2 references DOUBLE_VERIFIED?
                    ├── NO → MEDIUM confidence
                    └── YES → HIGH confidence
```

### Tier 3 Decision Tree

```
Are >= 3 markers meeting criteria?
├── NO → LOW confidence
└── YES → Continue
        │
        v
    Is pathway activity computed?
    ├── NO → ❌ CANNOT PROCEED (v2 rule)
    └── YES → Continue
            │
            v
        Do >= 2 pathways support functional state?
        ├── NO → MEDIUM confidence
        └── YES → Continue
                │
                v
            Are >= 2 references DOUBLE_VERIFIED?
            ├── NO → MEDIUM confidence
            └── YES → HIGH confidence
```

---

## 5. Novel Population Criteria (Updated)

### Trigger Conditions

Flag as **NOVEL** if ANY:

| Condition | Check |
|-----------|-------|
| **No literature match** | 3+ PubMed queries return no relevant results |
| **Opposing markers** | Co-expression of contradictory markers |
| **Unusual TF/Pathway** | TF or pathway combination not reported together |
| **Trajectory outlier** | Falls between canonical pseudotime positions |
| **Condition-specific** | Only in specific samples |

### Novel Confirmation Checklist

```
REQUIRED (ALL must be met):
[ ] >= 3 distinguishing markers (pct >= 25%, LFC >= 1)
[ ] Unique TF or pathway signature
[ ] >= 50 cells
[ ] Document 3+ failed PubMed queries
[ ] Not explained by doublet contamination
[ ] Not technical artifact (check QC)

EVIDENCE DOCUMENTATION:
[ ] List failed queries with results
[ ] Note TF/pathway anomalies
[ ] Record trajectory position if applicable
```

### Novel Naming

```
Format: Novel_{DistinctiveMarkers}

Examples:
- Novel_{marker1}+{marker2}+{marker3}
  (Use the top 2-3 most distinctive DE markers for the cluster)
```

---

## 6. Evidence Integration Template

### Tier 2 Template

```markdown
## Evidence Summary

| Type | Finding | Score |
|------|---------|-------|
| Markers | {marker1}+, {marker2}+, {marker3}+ (3 markers) | {0-3}/3 |
| TF Activity | {TF1} ({score}), {TF2} ({score}) active | {0-3}/3 |
| Trajectory | Pseudotime {value} ({category}) | {0-3}/3 |
| References | PMID:{id1}, PMID:{id2} ({n} DOUBLE_VERIFIED) | {0-3}/3 |
| **Total** | | **{total}/12** |

## Final Decision
- **Annotation**: {dev_state}
- **Confidence**: {High/Medium/Low} ({total}/12)
- **Key Evidence**:
  1. Marker combination supports {dev_state}
  2. TF activity consistent with {dev_state} program
  3. Pseudotime position matches lineage expectation
  4. Literature verified
```

### Tier 3 Template

```markdown
## Evidence Summary

| Type | Finding | Score |
|------|---------|-------|
| Markers | {marker1}+, {marker2}+, {marker3}+ (3 markers) | {0-3}/3 |
| Pathway | {pathway1} ({score}), {pathway2} ({score}) | {0-3}/3 |
| References | PMID:{id} ({n} DOUBLE_VERIFIED) | {0-3}/3 |
| Literature | {literature_summary} | {0-3}/3 |
| **Total** | | **{total}/12** |

## Final Decision
- **Annotation**: {func_state}
- **Confidence**: {High/Medium/Low} ({total}/12)
- **Key Evidence**:
  1. Marker combination supports {func_state}
  2. Pathway activity consistent with {func_state}
  3. Literature verified
```

---

## 7. When to Flag for Review

Flag cluster for manual review if:

```
[ ] Confidence < 7 after reasoning
[ ] TF/Pathway conflicts with markers
[ ] Novel criteria partially met (2-3 of 5)
[ ] < 50 cells in cluster
[ ] High doublet score
[ ] Ambiguous between 2+ canonical types
[ ] Functional analysis not computed (v2 violation)
```
