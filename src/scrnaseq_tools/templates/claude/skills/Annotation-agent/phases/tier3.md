# Phase 3: Tier 3 - Functional States (v2)

```
╔══════════════════════════════════════════════════════════════════════╗
║  ⚠️ v2 CHANGE: Pathway Activity is MANDATORY                         ║
║  Cannot proceed to annotation without pathway analysis               ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Workflow Overview

```
Step 1: Subset Data (from Tier 2)
Step 2: Re-cluster
Step 3: RE-COMPUTE DE (Critical)
Step 4: Pathway Activity Analysis (MANDATORY) ← NEW ENFORCEMENT
Step 4.5: Statistical Outlier Detection (DATA-DRIVEN) ← NEW
Step 4.6: Evidence Conflict Detection (DATA-DRIVEN) ← NEW
Step 5: VERIFY Functional Analysis ← NEW
Step 6: Filter Valid Markers
Step 7: Integrated Pre-Analysis ← UPDATED with Outliers/Conflicts
Step 8: MCP Verification (PMID Required)
Step 9: 3-Iteration Reasoning with Integrated Evidence ← UPDATED
Step 10: Detect Novel Populations
Step 11: Assign Annotations
Step 12: Save Results
```

---

## Step-by-Step Workflow

### Step 1: Subset Data

```python
major_type = 'T cells'
dev_state = 'Effector'

# Load Tier 2 subset
tier2_subset = sc.read_h5ad(f'annotation_output/subsets/tier2/{major_type.replace(" ", "_")}.h5ad')

# Further subset to one developmental state
subset = tier2_subset[tier2_subset.obs['tier2_annotation'] == dev_state].copy()
print(f"Subsetting {major_type} -> {dev_state}: {subset.n_obs} cells")
```

### Step 2: Re-cluster

```python
# Higher resolution for functional states
sc.pp.neighbors(subset, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=50)
sc.tl.leiden(subset, resolution=1.2, key_added='tier3_cluster')
print(f"Clusters: {subset.obs['tier3_cluster'].nunique()}")
```

### Step 3: RE-COMPUTE DE (Critical)

```python
# MUST re-compute - Tier 2 markers no longer differentially expressed
sc.tl.rank_genes_groups(subset, groupby='tier3_cluster', method='wilcoxon')
de_df = sc.get.rank_genes_groups_df(subset, group=None)
print(f"⚠️ DE computed for subset: {major_type}_{dev_state} ({subset.n_obs} cells)")
```

### Step 4: Pathway Activity Analysis (MANDATORY)

```python
import decoupler as dc

print("🔬 Step 4: Pathway Activity Analysis (MANDATORY)")

# Get PROGENy network
net = dc.op.progeny(organism='human', top=500)
print(f"   PROGENy network: {len(net['source'].unique())} pathways")

# Run MLM
dc.run_mlm(
    mat=subset,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
    use_raw=False
)

# VERIFY
assert 'mlm_estimate' in subset.obsm, "❌ Pathway activity not computed!"
assert 'mlm_pvals' in subset.obsm, "❌ Pathway p-values not computed!"
print(f"✅ Pathway activity computed: {subset.obsm['mlm_estimate'].shape[1]} pathways")
```

### Step 4.5: Statistical Outlier Detection (NEW - DATA-DRIVEN)

```python
print("🔍 Step 4.5: Statistical Outlier Detection (data-driven, NO hardcoding)")

# Detect Pathway activity outliers using z-scores (completely data-driven)
cluster_col = 'tier3_cluster'
pw_key = 'score_mlm' if 'score_mlm' in subset.obsm else 'mlm_estimate'
pval_key = 'padj_mlm' if 'padj_mlm' in subset.obsm else 'mlm_pvals'

pw_scores = subset.obsm[pw_key]
pw_pvals = subset.obsm[pval_key]

outliers = []

for cluster in subset.obs[cluster_col].unique():
    mask = subset.obs[cluster_col] == cluster
    cluster_pw_mean = pw_scores[mask].mean()
    cluster_pw_pval = pw_pvals[mask].mean()

    for pathway in pw_scores.columns:
        # Global statistics across all clusters
        global_mean = pw_scores[pathway].mean()
        global_std = pw_scores[pathway].std()
        cluster_val = cluster_pw_mean[pathway]
        cluster_pval = cluster_pw_pval[pathway]

        # Z-score (statistical measure of how unusual this value is)
        if global_std > 1e-6:  # Avoid division by zero
            z_score = (cluster_val - global_mean) / global_std
        else:
            z_score = 0.0

        # Flag as outlier if:
        # 1. |z-score| > 2.5 (statistically significant deviation)
        # 2. |Activity| > 0.5 (meaningful activity, not just noise)
        # 3. p-value < 0.05 (significant in this cluster)
        if abs(z_score) > 2.5 and abs(cluster_val) > 0.5 and cluster_pval < 0.05:
            outliers.append({
                'cluster': str(cluster),
                'pathway': pathway,
                'activity': float(cluster_val),
                'z_score': float(z_score),
                'pval': float(cluster_pval),
                'direction': 'HIGH' if z_score > 0 else 'LOW'
            })

# Store outliers in AnnData for reasoning
subset.uns['pathway_outliers'] = outliers

if outliers:
    print(f"   ⚠️ FOUND {len(outliers)} Pathway activity outliers (statistical anomalies)")
    print(f"   → These may indicate specialized/rare functional states")
    print(f"   → Will be investigated in reasoning phase")

    # Show summary
    from collections import Counter
    cluster_counts = Counter(o['cluster'] for o in outliers)
    for cluster, count in cluster_counts.most_common(5):
        pathways = [o['pathway'] for o in outliers if o['cluster'] == cluster]
        print(f"      Cluster {cluster}: {count} outlier pathways ({', '.join(pathways[:3])}...)")
else:
    print(f"   ✅ No statistical outliers detected")
    print(f"   → Standard functional states expected")
```

### Step 4.6: Evidence Conflict Detection (NEW - DATA-DRIVEN)

```python
print("🔍 Step 4.6: Evidence Conflict Detection (data-driven)")

# Detect conflicts between markers and pathway activity (NO hardcoding)
conflicts = []

for cluster in subset.obs[cluster_col].unique():
    mask = subset.obs[cluster_col] == cluster

    # Get top markers for this cluster
    cluster_de = de_df[de_df['group'] == cluster].nlargest(20, 'logfoldchanges')
    top_markers = cluster_de['names'].tolist() if 'names' in cluster_de.columns else []

    # Get top pathways
    cluster_pw_mean = pw_scores[mask].mean()
    top_pw_indices = cluster_pw_mean.nlargest(10).index.tolist()
    low_pw_indices = cluster_pw_mean.nsmallest(10).index.tolist()

    # CONFLICT 1: Pathway activity outlier (unusual functional program)
    cluster_outliers = [o for o in outliers if o['cluster'] == str(cluster)]
    if cluster_outliers:
        high_outlier_pws = [o['pathway'] for o in cluster_outliers if o['direction'] == 'HIGH']
        low_outlier_pws = [o['pathway'] for o in cluster_outliers if o['direction'] == 'LOW']

        if high_outlier_pws:
            conflicts.append({
                'cluster': str(cluster),
                'type': 'PATHWAY_STATISTICAL_OUTLIER',
                'detail': f"Pathways with unusually HIGH activity: {', '.join(high_outlier_pws[:5])}",
                'severity': 'HIGH',
                'action': 'Search literature for these pathway combinations + parent context'
            })

        if low_outlier_pws:
            conflicts.append({
                'cluster': str(cluster),
                'type': 'PATHWAY_SUPPRESSION_OUTLIER',
                'detail': f"Pathways with unusually LOW activity: {', '.join(low_outlier_pws[:5])}",
                'severity': 'MEDIUM',
                'action': 'May indicate exhausted/anergic/quiescent state'
            })

    # CONFLICT 2: Opposing marker patterns
    # Example: Both activation markers and exhaustion markers
    activation_markers = ['CD69', 'CD38', 'IL2RA', 'ICOS']
    exhaustion_markers = ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'TOX']

    activ_count = sum(1 for m in activation_markers if m in top_markers)
    exhaust_count = sum(1 for m in exhaustion_markers if m in top_markers)

    if activ_count >= 2 and exhaust_count >= 2:
        conflicts.append({
            'cluster': str(cluster),
            'type': 'OPPOSING_MARKER_PATTERN',
            'detail': f'Both activation ({activ_count}) and exhaustion ({exhaust_count}) markers present',
            'severity': 'HIGH',
            'action': 'Investigate: chronic stimulation, progenitor exhausted, or transitional state'
        })

    # CONFLICT 3: Cytotoxic markers without inflammatory pathways
    cytotox_markers = ['GZMB', 'PRF1', 'NKG7', 'GNLY']
    cytotox_count = sum(1 for m in cytotox_markers if m in top_markers)

    if cytotox_count >= 2:
        # Check for inflammatory pathways
        inflammatory_pws = ['TNFa', 'NFkB', 'Trail', 'JAK-STAT']
        inflam_pw_count = sum(1 for pw in inflammatory_pws if pw in top_pw_indices)

        if inflam_pw_count == 0:
            conflicts.append({
                'cluster': str(cluster),
                'type': 'CYTOTOXIC_WITHOUT_INFLAMMATION',
                'detail': f'Cytotoxic markers ({cytotox_count}) but no inflammatory pathway activity',
                'severity': 'MEDIUM',
                'action': 'May indicate: exhausted cytotoxic, pre-activated, or tissue-resident'
            })

    # CONFLICT 4: Stem-like markers in terminal state context
    if dev_state in ['Effector', 'Memory', 'Plasma_cells']:
        stem_markers = ['TCF7', 'LEF1', 'ID3', 'SELL']
        stem_count = sum(1 for m in stem_markers if m in top_markers)

        if stem_count >= 2:
            conflicts.append({
                'cluster': str(cluster),
                'type': 'STEM_IN_TERMINAL_CONTEXT',
                'detail': f'Stem-like markers ({stem_count}) in {dev_state} context',
                'severity': 'HIGH',
                'action': 'May indicate: precursor exhausted (Tpex), stem-like memory, or plasticity'
            })

# Store conflicts
subset.uns['evidence_conflicts'] = conflicts

if conflicts:
    high_priority = [c for c in conflicts if c['severity'] == 'HIGH']
    print(f"   ⚠️ FOUND {len(conflicts)} evidence conflicts ({len(high_priority)} high priority)")
    print(f"   → These require careful investigation in reasoning")

    # Show high priority conflicts
    for conflict in high_priority[:5]:
        print(f"      Cluster {conflict['cluster']}: {conflict['type']}")
        print(f"         → {conflict['detail']}")
else:
    print(f"   ✅ No major evidence conflicts detected")
```

### Step 5: VERIFY Functional Analysis (NEW)

```python
print("🔍 Step 5: Verifying functional analysis...")

errors = []

if 'mlm_estimate' not in subset.obsm:
    errors.append("Pathway activity (mlm_estimate) not computed")
if 'mlm_pvals' not in subset.obsm:
    errors.append("Pathway p-values (mlm_pvals) not computed")

if errors:
    raise AssertionError(
        "❌ FUNCTIONAL ANALYSIS INCOMPLETE:\n" +
        "\n".join(f"  - {e}" for e in errors) +
        "\n\n⚠️ Cannot proceed to annotation without pathway analysis!"
    )

print("✅ Functional analysis verified - ready for annotation")
```

### Step 6: Filter Valid Markers

```python
valid_markers_dict = {}
for cluster in subset.obs['tier3_cluster'].unique():
    valid_markers_dict[cluster] = filter_valid_markers(de_df, cluster, top_n=50)
```

### Step 7: Integrated Pre-Analysis (UPDATED with Outliers/Conflicts)

**Present ALL evidence to Claude:**

```markdown
=== {Major Type} → {Dev State} → Cluster {X} Analysis ===

## DE Markers (Top 50 Valid)
| Rank | Gene | pct_in | log2FC | padj |
|------|------|--------|--------|------|
| 1 | GZMB | 0.92 | 3.5 | 1e-50 |
| 2 | PRF1 | 0.85 | 2.8 | 1e-45 |
[show all 50]

## Pathway Activity (Top 10)
| Pathway | Activity | p-value | Z-score | Outlier? |
|---------|----------|---------|---------|----------|
| {PW1} | {val} | {pval} | {z:.2f} | {✅ NORMAL / ⚠️ HIGH OUTLIER / ⚠️ LOW OUTLIER} |
| {PW2} | {val} | {pval} | {z:.2f} | {✅ NORMAL / ⚠️ HIGH OUTLIER / ⚠️ LOW OUTLIER} |
[show top 10 + ALL outliers even if not in top 10]

## 🔍 STATISTICAL OUTLIERS (if detected)

**Pathway Activity Outliers for Cluster {X}:**
{for each outlier in this cluster}
- **{Pathway}** (activity={val:.2f}, z-score={z:.2f}, p={pval:.2e})
  → This pathway activity is **statistically unusual** (|z| > 2.5)
  → **Action**: Search literature for "{Pathway} {parent_context}" combination
  → **Interpretation**: May indicate specialized functional state

**Total outliers in this cluster: {N}**
**Interpretation**:
- IF N = 0: Standard functional state expected
- IF N = 1-2: Minor variation, investigate if marker pattern also unusual
- IF N >= 3: **HIGH PRIORITY** - likely specialized functional state (e.g., Tpex, TRM, anergic)

## ⚠️ EVIDENCE CONFLICTS (if detected)

{for each conflict in this cluster}
**{conflict_type}**:
- {detail}
- Severity: {HIGH/MEDIUM/INFO}
- Action: {action_instruction}

## Parent Context
- Tier 1: {major_type}
- Tier 2: {dev_state}
- Tissue: {tissue}

## 🧠 REASONING INSTRUCTION

**Standard Workflow (NO outliers/conflicts)**:
1. Identify functional state based on markers + pathways + parent context
2. Search literature for standard marker+pathway combinations
3. Assign annotation with standard confidence scoring

**SPECIAL WORKFLOW (outliers or HIGH severity conflicts detected)**:
1. **DO NOT assume standard functional states**
2. **Priority literature searches**:
   - For EACH outlier pathway: "{pathway_name} {parent_context}"
   - For unusual marker combinations: "{marker1} {marker2} {major_type} {dev_state}"
   - For conflicts: "{conflict_keywords} {parent_context}"
   - Focus on: "precursor", "progenitor", "tissue-resident", "exhausted", "anergic", "plastic"
3. **Candidate generation**:
   - Candidate A: Literature-supported specialized state (HIGH confidence if found)
   - Candidate B: Standard functional state (LOW confidence, explain outlier conflict)
   - Candidate C: Novel/uncharacterized if no literature match (document failed queries)
4. **Decision criteria**:
   - IF literature explains outlier/conflict → Use specialized annotation
   - IF no literature but pattern coherent → Flag as NOVEL with descriptive name
   - IF pattern incoherent → Reduce confidence, manual review required

## Evidence Summary (for Claude)
- Marker signature suggests: [candidates from markers]
- Pathway activity suggests: [candidates from pathways]
- **Pathway outliers suggest**: {interpretation OR "N/A - no outliers"}
- **Evidence conflicts**: {summary OR "N/A - no conflicts"}
- Parent context supports: [candidates fitting dev_state]
- **Evidence convergence**: {CONVERGENT / CONFLICTED - explain}
```

**Claude Output Format:**

```yaml
cluster: {X}
parent_context: "{major_type} → {dev_state}"
has_outliers: {true/false}
outlier_count: {N}
has_conflicts: {true/false}

# IF OUTLIERS DETECTED:
outlier_investigation:
  outlier_pathways: ["{PW1}", "{PW2}"]  # List of outlier pathways
  literature_search_required: true
  priority: "HIGH"  # If >= 3 outliers

candidates:
  # IF OUTLIERS/CONFLICTS: Prioritize specialized functional state candidates
  - func_state: "{Specialized_State}"  # e.g., "Tpex", "TRM", "Anergic"
    confidence: "High/Medium/Low"
    key_markers: ["{marker1}", "{marker2}", "{marker3}"]
    key_pathways: ["{PW1}", "{PW2}"]
    outlier_pathways: ["{PW_outlier1}"]  # NEW - Pathways flagged as statistical outliers
    pathway_support: "{description}"
    pathway_outlier_explanation: "{how outlier supports this annotation}"  # NEW
    conflict_resolution: "{how this explains any conflicts}"  # NEW
    reasoning: "{integrated reasoning addressing outliers/conflicts}"

  # Standard functional state (lower confidence if outliers unexplained)
  - func_state: "{Standard_State}"  # e.g., "Cytotoxic"
    confidence: "Low"  # Reduced if outliers present
    key_markers: ["{marker1}", "{marker2}"]
    key_pathways: ["{PW1}"]
    outlier_pathways: []
    pathway_support: "{description}"
    pathway_outlier_conflict: "{cannot explain outlier pathways}"  # NEW
    reasoning: "{why this is less likely given outliers/conflicts}"

is_novel_candidate: {true/false}
novel_reason: "{if true, explain why no literature support AND pattern is coherent}"

mcp_queries:
  # IF OUTLIERS: Add specific outlier queries
  - query: "{outlier_pathway} {parent_context}"  # e.g., "WNT Effector T cell"
    purpose: "Investigate statistical pathway outlier"
    priority: "HIGH"

  - query: "{conflict_keywords} {parent_context}"  # e.g., "GZMB TCF7 exhausted T cell"
    purpose: "Resolve evidence conflict"
    priority: "HIGH"

  - query: "{marker1} {marker2} {pathway1} {parent_context}"
    purpose: "Standard marker+pathway validation"
    priority: "NORMAL"
```

### Step 8: MCP Verification (PMID Required)

```python
for candidate in claude_candidates:
    for query in candidate['mcp_queries']:
        result = mcp.pubmed_search(query['query'], max_results=3)

        if result['count'] > 0:
            pmid = result['results'][0]['pmid']
            verification = mcp.verify_reference(
                pmid=pmid,
                markers=candidate['key_markers'],
                cell_type=f"{candidate['func_state']} {parent_context}"
            )
            print(f"PMID:{pmid} - {verification['status']}")
```

### Step 9: 3-Iteration Reasoning with Integrated Evidence (UPDATED - GENERATES reasoning_results)

```python
print(f"\n🧠 Step 9: Performing 3-Iteration Reasoning for {major_type} -> {dev_state}...")

# This step generates reasoning_results which is used in:
# - Step 11.5: Create annotation_evidence
# - Step 12: Visualization

reasoning_results = {}

for cluster in subset.obs['tier3_cluster'].unique():
    cluster_str = str(cluster)
    print(f"\n--- Cluster {cluster_str} Reasoning ---")

    # Get cluster-specific data
    cluster_mask = subset.obs['tier3_cluster'] == cluster
    cluster_de = de_df[de_df['group'] == cluster]

    # Filter valid markers
    valid_markers = cluster_de[
        (cluster_de['pct_nz_group'] >= 0.25) &
        (cluster_de['logfoldchanges'] >= 1.0) &
        (cluster_de['pvals_adj'] < 0.05)
    ].nlargest(50, 'logfoldchanges')

    # Extract top markers
    top_genes = valid_markers.nlargest(10, 'logfoldchanges')['names'].tolist()

    # Get Pathway activity for this cluster
    pathway_scores = subset.obsm['mlm_estimate']
    cluster_pathway_mean = pathway_scores[cluster_mask].mean()
    top_pathways = cluster_pathway_mean.nlargest(5).index.tolist()
    top_pathway_scores = cluster_pathway_mean.nlargest(5).values.tolist()

    # Get outliers and conflicts for this cluster
    cluster_outliers = [o for o in subset.uns.get('pathway_outliers', [])
                       if o['cluster'] == cluster_str]
    cluster_conflicts = [c for c in subset.uns.get('evidence_conflicts', [])
                        if c['cluster'] == cluster_str]

    # Create structured reasoning output
    cluster_reasoning = {
        'cluster_id': cluster_str,
        'final_assignment': None,  # Will be filled by reasoning
        'full_label': None,

        # Markers
        'marker_details': [
            {
                'gene': row['names'],
                'pct_in': row['pct_nz_group'],
                'log2fc': row['logfoldchanges'],
                'padj': row['pvals_adj']
            }
            for _, row in valid_markers.head(5).iterrows()
        ],

        # Candidates
        'candidates': [],

        # References (from MCP verification)
        'references': [],

        # Supporting reasons
        'supporting_reasons': [
            f"Top markers: {', '.join(top_genes[:5])}",
            f"Top pathways: {', '.join([f'{pw}({score:.2f})' for pw, score in zip(top_pathways[:3], top_pathway_scores[:3])])}"
        ],

        # Conflicts
        'conflicts': cluster_conflicts,

        # Pathway activity
        'pathway_activity': [
            {'pathway': pw, 'score': float(score), 'pval': 0.01}
            for pw, score in zip(top_pathways, top_pathway_scores)
        ],

        # Outliers
        'pathway_outliers': cluster_outliers,

        # Confidence scoring
        'confidence_score': 0,
        'confidence_level': 'Low',

        # Novel flag
        'is_novel': False,
        'novel_evidence': None
    }

    # Check for conflict severity
    high_severity_conflicts = [c for c in cluster_conflicts if c.get('severity') == 'HIGH']
    if high_severity_conflicts:
        cluster_reasoning['supporting_reasons'].append(
            f"⚠️ HIGH severity conflicts: {len(high_severity_conflicts)}"
        )

    # Check for pathway outliers
    if cluster_outliers:
        outlier_pathways = [o['pathway'] for o in cluster_outliers]
        cluster_reasoning['supporting_reasons'].append(
            f"⚠️ Pathway outliers: {', '.join(outlier_pathways[:3])}"
        )

    reasoning_results[cluster_str] = cluster_reasoning
    print(f"  ✓ Reasoning structure created for cluster {cluster_str}")

print(f"\n✅ Reasoning completed for {len(reasoning_results)} clusters")
print("⚠️  NOTE: In production, complete the reasoning with:")
print("   → LLM 3-iteration reasoning (reasoning/integrated_format.md)")
print("   → MCP literature verification")
print("   → Pathway-based confidence scoring")
print("   → Outlier/Conflict investigation")
```

**Template Reference**: `reasoning/integrated_format.md`

**Reasoning MUST include:**
- **Iteration 1**: Markers + Pathway Activity + Outliers + Conflicts
- **Iteration 2**: Candidate comparison + outlier/conflict investigation
- **Iteration 3**: Integrated scoring and final decision

**Confidence Scoring:**

| Criteria | Standard (No outliers) | With Outliers (unexplained) | With Outliers (literature-explained) |
|----------|------------------------|----------------------------|-------------------------------------|
| Markers | 3 pts | 3 pts | 3 pts |
| References | 3 pts | 3 pts (REQUIRED) | 3 pts (REQUIRED) |
| Pathway consistency | 3 pts | 0-1 pts | 3 pts |
| Literature/Conflict | 3 pts | 0 pts (unresolved) | 3 pts (resolved) |
| **Total** | **12 pts** | **6-7 pts** | **12 pts** |

### Step 10: Detect Novel Populations (ENHANCED)

**Flag as NOVEL if:**

```
[ ] No literature match after 3+ queries (including outlier/conflict-specific queries)
[ ] Opposing markers present (e.g., GZMB + TCF7 without clear Tpex signature)
[ ] Unusual pathway combination not reported in literature (outlier pathways)
[ ] Statistical outliers (>= 3) with no literature explanation
[ ] Evidence conflicts (HIGH severity) unresolved by literature
[ ] Only appears in specific samples/conditions
```

**Required for Novel confirmation:**
```
[ ] >= 3 distinguishing markers
[ ] Unique pathway signature (with statistical support - outliers or specific combination)
[ ] >= 50 cells
[ ] Document 3+ failed queries including:
    - Standard markers: "{marker1} {marker2} {parent_context}"
    - Outlier pathways: "{outlier_pathway} {parent_context}"
    - Conflicts: "{conflict_keywords} {parent_context}"
[ ] Evidence is COHERENT (not just incoherent noise):
    - Markers form functional signature
    - Pathway pattern is interpretable
    - NOT due to doublets or contamination
```

**Descriptive Naming for Novel Populations:**
```
Format: {ParentContext}_{DistinguishingFeature}

Examples:
- T_Effector_GZMB+TCF7+  (unusual marker combo)
- T_Effector_WNT_high    (outlier pathway)
- B_GC_NFkB_low          (suppressed pathway)
- T_Memory_Prolif_anergi (conflict: proliferation + low activity)
```

### Step 11: Assign Annotations

```python
tier3_mapping = {
    '0': 'Cytotoxic',
    '1': 'Exhausted',
    '2': 'TRM',
    '3': 'Novel_GZMB+TCF7',  # Novel population
}

subset.obs['tier3_annotation'] = subset.obs['tier3_cluster'].map(tier3_mapping)

# Final hierarchical label
subset.obs['final_annotation'] = f"{major_type}_{dev_state}_" + subset.obs['tier3_annotation'].astype(str)
```

### Step 11.5: Create annotation_evidence (NEW - REQUIRED FOR VISUALIZATION)

```python
# Convert reasoning results to annotation_evidence format
# This is REQUIRED for visualization in Step 12

from reasoning.agent_format import reasoning_to_evidence

annotation_evidence = []
safe_name = f"{major_type}_{dev_state}".replace(' ', '_')

for cluster in subset.obs['tier3_cluster'].unique():
    # Get annotation for this cluster
    func_state = tier3_mapping.get(str(cluster))
    if not func_state:
        continue

    # Get cluster-specific data from reasoning process
    # NOTE: This assumes reasoning was performed for each cluster in Step 9
    cluster_reasoning = reasoning_results.get(str(cluster), {})

    # Create evidence entry
    evidence_entry = reasoning_to_evidence(
        cluster_id=f"tier3_{safe_name}_cluster_{cluster}",
        tier="tier3",
        subset_id=safe_name,
        reasoning_output=cluster_reasoning
    )

    annotation_evidence.append(evidence_entry)

print(f"✅ Created annotation_evidence for {major_type}_{dev_state}: {len(annotation_evidence)} entries")

# Save to JSON for persistence
import json
import os
os.makedirs('annotation_output/references', exist_ok=True)

# Load existing evidence if present
evidence_file = 'annotation_output/references/annotation_evidence.json'
all_evidence = []
if os.path.exists(evidence_file):
    with open(evidence_file, 'r') as f:
        all_evidence = json.load(f)

# Add new Tier 3 evidence
all_evidence.extend(annotation_evidence)

# Save combined evidence
with open(evidence_file, 'w') as f:
    json.dump(all_evidence, f, indent=2, ensure_ascii=False)

print(f"✅ Updated: {evidence_file} (total entries: {len(all_evidence)})")
```

### Step 12: Save Results

```python
safe_name = f"{major_type}_{dev_state}".replace(' ', '_')
subset.write(f'annotation_output/subsets/tier3/{safe_name}.h5ad')

# Verify saved data has pathway analysis
saved = sc.read_h5ad(f'annotation_output/subsets/tier3/{safe_name}.h5ad')
assert 'mlm_estimate' in saved.obsm, "Pathway activity not saved!"
print("✅ Subset saved with pathway analysis")

# ✅ Generate Visualizations (NEW)
print(f"\n🎨 Generating Tier 3 Visualizations for {major_type} -> {dev_state}...")
import os
os.makedirs('annotation_output/figures', exist_ok=True)

# Prepare marker dict from annotation_evidence (VALIDATED markers only)
# annotation_evidence는 Step 9에서 reasoning 결과로 생성됨
marker_dict = {}
annotation_order = sorted(subset.obs['tier3_annotation'].unique())

for entry in annotation_evidence:
    if entry['tier'] == 'tier3' and entry['subset_id'] == safe_name:
        func_state = entry['annotation']
        # Extract validated marker gene names from reasoning results
        validated_markers = [m['gene'] for m in entry['markers']]
        marker_dict[func_state] = validated_markers

# Import visualization functions
import sys
sys.path.append('annotation_output')
from tools.visualization import save_all_visualizations
from tools.functional_plots import (
    plot_pathway_activity_heatmap,
    plot_pathway_dotplot
)

# 1. Generate standard visualizations (UMAP, Dotplot, Combined)
save_all_visualizations(
    adata=subset,
    annotation_col='tier3_annotation',
    marker_dict=marker_dict,
    reasoning_dict={e['annotation']: e for e in annotation_evidence
                    if e['tier'] == 'tier3' and e['subset_id'] == safe_name},
    reference_dict={e['annotation']: e['references'] for e in annotation_evidence
                    if e['tier'] == 'tier3' and e['subset_id'] == safe_name},
    annotation_order=annotation_order,
    title_prefix=f'Tier 3 - {major_type} {dev_state}',
    output_dir='annotation_output/figures/'
)

# 2. Generate Pathway activity heatmap (Tier 3 MANDATORY)
plot_pathway_activity_heatmap(
    subset,
    cluster_col='tier3_annotation',
    figsize=(8, 6),
    save_path=f'annotation_output/figures/{safe_name}_tier3_pathway_heatmap'
)
print(f"   - {safe_name}_tier3_pathway_heatmap.png/.svg")

# 3. Generate Pathway dotplot (optional, comprehensive view)
plot_pathway_dotplot(
    subset,
    cluster_col='tier3_annotation',
    figsize=(10, 6),
    save_path=f'annotation_output/figures/{safe_name}_tier3_pathway_dotplot'
)
print(f"   - {safe_name}_tier3_pathway_dotplot.png/.svg")

print(f"✅ Tier 3 Visualizations for {major_type} -> {dev_state} completed:")
print(f"   - {safe_name}_tier3_umap.png/.svg")
print(f"   - {safe_name}_tier3_dotplot.png/.svg (with marker brackets, VALIDATED markers)")
print(f"   - {safe_name}_tier3_pathway_heatmap.png/.svg (MANDATORY)")
print(f"   - {safe_name}_tier3_combined.png/.svg (multi-panel)")
print(f"   - {safe_name}_tier3_reasoning.json/.txt")
```

---

## Pathway Interpretation Reference

### Key Pathway-Functional State Associations

| Pathway | High Activity Suggests | Low Activity Suggests | Outlier Context |
|---------|------------------------|----------------------|-----------------|
| **TNFa** | Pro-inflammatory, activated effector | Quiescent, anergic | **HIGH in Exhausted** → Chronic inflammation |
| **NFkB** | Activation, survival, inflammation | Apoptotic, suppressed | **LOW in Effector** → Dysfunction |
| **JAK-STAT** | Cytokine response, proliferation | Cytokine unresponsive | **HIGH without proliferation** → Chronic stimulation |
| **Trail** | Apoptosis signaling, cytotoxic function | Non-cytotoxic | **HIGH in non-Effector** → Atypical killing |
| **MAPK** | Proliferation, differentiation | Quiescent | **LOW in activated** → Cell cycle arrest |
| **PI3K** | Survival, metabolism | Metabolic stress | **HIGH in Exhausted** → Compensatory metabolism |
| **Hypoxia** | Tissue adaptation, tumor microenvironment | Normoxic | **HIGH outlier** → Tissue-resident or tumor |
| **TGFb** | Immunosuppression, regulatory | Pro-inflammatory | **HIGH in Effector** → Regulatory dysfunction |
| **WNT** | Stemness, self-renewal | Terminal differentiation | **HIGH in Effector** → Precursor/plastic state |
| **p53** | DNA damage, stress response | Unstressed | **HIGH outlier** → Stress/senescence |
| **VEGF** | Angiogenesis | - | **HIGH in lymphocytes** → Tissue remodeling |

### Standard Pathway Combinations

| Combination | Likely Functional State |
|-------------|------------------------|
| TNFa + NFkB + Trail | Cytotoxic/Effector |
| JAK-STAT + MAPK | Proliferating |
| TGFb high | Regulatory/Suppressive |
| Hypoxia + TNFa | Tumor-infiltrating effector |
| WNT + PI3K | Stem-like/Memory |

### **Outlier Pathway Combinations** (Specialized States)

| Unusual Combination | Potential Interpretation | Literature Keywords |
|---------------------|-------------------------|---------------------|
| **WNT high in Effector** | Precursor exhausted (Tpex), stem-like effector | "TCF7 effector", "Tpex", "stem-like CD8" |
| **TNFa low in Cytotoxic** | Exhausted cytotoxic, cold tumor infiltrating | "exhausted cytotoxic", "dysfunctional" |
| **TGFb high + GZMB high** | Regulatory cytotoxic (rare), suppressed TIL | "TGFb cytotoxic", "suppressed effector" |
| **Hypoxia high (outlier)** | Tissue-resident memory, tumor-infiltrating | "TRM", "hypoxic T cell", "tumor microenvironment" |
| **Multiple pathways LOW** | Anergic, exhausted terminal, quiescent | "anergy", "terminal exhaustion", "quiescence" |
| **p53 + low proliferation** | Senescent, DNA damage arrest | "senescent T cell", "DNA damage lymphocyte" |

---

## Checklist (Per Developmental State)

```
MANDATORY:
- [ ] Subset extracted from Tier 2
- [ ] Re-clustered (resolution 1.0-1.5)
- [ ] DE RE-COMPUTED within subset
- [ ] ⚠️ Pathway activity computed (mlm_estimate in obsm)
- [ ] ⚠️ Statistical outliers detected (pathway_outliers in uns) ← NEW
- [ ] ⚠️ Evidence conflicts detected (evidence_conflicts in uns) ← NEW
- [ ] ⚠️ Functional analysis VERIFIED before annotation

ANNOTATION:
- [ ] Valid markers filtered (Top 50)
- [ ] Integrated Pre-Analysis (markers + pathways + outliers + conflicts) ← UPDATED
- [ ] MCP verification with PMID
- [ ] 3-iteration reasoning with integrated evidence
- [ ] IF OUTLIERS: Specialized functional state literature searches performed ← NEW
- [ ] IF CONFLICTS: Conflict resolution literature searches performed ← NEW
- [ ] IF OUTLIERS/CONFLICTS: Outlier-aware confidence scoring ← NEW
- [ ] Novel populations checked (enhanced with outlier detection) ← UPDATED
- [ ] References double-verified
- [ ] Confidence scored (max 12 pts, penalties for unexplained outliers)

OUTPUT:
- [ ] Tier 3 annotations assigned (standard OR specialized based on outliers/conflicts)
- [ ] Final hierarchical labels constructed
- [ ] Subset saved with pathway analysis + outlier/conflict info
- [ ] Visualizations saved
```

---

## Output Format

```
{Major Type}_{Dev State}_{Functional State}

Examples:
- T cells_Effector_Cytotoxic
- T cells_Effector_Exhausted
- T cells_Effector_TRM
- T cells_Effector_Novel_GZMB+TCF7
- B cells_GC_Light_Zone
- B cells_GC_Dark_Zone
```
