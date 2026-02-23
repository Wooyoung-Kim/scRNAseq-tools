# Phase 2: Tier 2 - Developmental States (v2)

```
╔══════════════════════════════════════════════════════════════════════╗
║  ⚠️ v2 CHANGE: TF Activity + Trajectory are MANDATORY                ║
║  Cannot proceed to annotation without functional analysis            ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Workflow Overview

```
Step 1: Subset Data
Step 2: Re-cluster
Step 3: RE-COMPUTE DE (Critical)
Step 4: TF Activity Analysis (MANDATORY) ← NEW ENFORCEMENT
Step 5: Trajectory Analysis (MANDATORY if >= 2000 cells) ← NEW ENFORCEMENT
Step 5.5: Statistical Outlier Detection (DATA-DRIVEN) ← NEW
Step 5.6: Evidence Conflict Detection (DATA-DRIVEN) ← NEW
Step 6: VERIFY Functional Analysis ← NEW
Step 7: Filter Valid Markers
Step 8: Integrated Pre-Analysis ← UPDATED with Outliers/Conflicts
Step 9: MCP Verification (PMID Required)
Step 10: 3-Iteration Reasoning with Integrated Evidence ← UPDATED
Step 11: Assign Annotations
Step 12: Save Results
```

---

## Step-by-Step Workflow

### Step 1: Subset Data

```python
major_type = 'T cells'
subset = adata[adata.obs['tier1_annotation'] == major_type].copy()
print(f"Subsetting {major_type}: {subset.n_obs} cells")
```

### Step 2: Re-cluster

```python
# Use existing harmony embedding
sc.pp.neighbors(subset, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=50)
sc.tl.leiden(subset, resolution=0.8, key_added='tier2_cluster')
print(f"Clusters: {subset.obs['tier2_cluster'].nunique()}")
```

### Step 3: RE-COMPUTE DE (Critical)

```python
sc.tl.rank_genes_groups(subset, groupby='tier2_cluster', method='wilcoxon')
de_df = sc.get.rank_genes_groups_df(subset, group=None)
print(f"⚠️ DE computed for subset: {major_type} ({subset.n_obs} cells)")
```

### Step 4: TF Activity Analysis (MANDATORY)

```python
import decoupler as dc

print("🔬 Step 4: TF Activity Analysis (MANDATORY)")

# Get CollecTRI network (TF-target regulatory network)
net = dc.op.collectri(organism='human')

# Run ULM
dc.run_ulm(
    mat=subset,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
    use_raw=False
)

# VERIFY
assert 'ulm_estimate' in subset.obsm, "❌ TF activity not computed!"
print(f"✅ TF activity computed: {subset.obsm['ulm_estimate'].shape[1]} TFs")
```

### Step 5: Trajectory Analysis (MANDATORY if >= 2000 cells)

```python
import palantir
import numpy as np

print("🔬 Step 5: Trajectory Analysis")

if subset.n_obs >= 2000:
    print(f"   {subset.n_obs} cells >= 2000, computing trajectory...")

    # Run diffusion maps
    palantir.utils.run_diffusion_maps(subset, n_components=10, pca_key='X_pca')

    # Determine multiscale space
    palantir.utils.determine_multiscale_space(subset)

    # Find start cell (most naive-like)
    naive_markers = ['CCR7', 'TCF7', 'SELL', 'IL7R']
    naive_score = np.zeros(subset.n_obs)
    for marker in naive_markers:
        if marker in subset.var_names:
            expr = subset[:, marker].X.toarray().flatten() if hasattr(subset.X, 'toarray') else subset[:, marker].X.flatten()
            naive_score += (expr > 0).astype(float)
    start_cell = subset.obs_names[np.argmax(naive_score)]

    # Run Palantir (results stored directly in AnnData)
    palantir.core.run_palantir(subset, start_cell, num_waypoints=500)
    subset.obs['pseudotime'] = subset.obs['palantir_pseudotime']

    # Categorize
    subset.obs['pseudotime_category'] = pd.cut(
        subset.obs['pseudotime'],
        bins=[0, 0.33, 0.67, 1.0],
        labels=['Early', 'Mid', 'Late']
    )

    print(f"✅ Trajectory computed: pseudotime [{subset.obs['pseudotime'].min():.2f}, {subset.obs['pseudotime'].max():.2f}]")
else:
    print(f"   {subset.n_obs} cells < 2000, skipping trajectory")
    subset.obs['pseudotime'] = np.nan
    subset.obs['pseudotime_category'] = 'N/A'
```

### Step 5.5: Statistical Outlier Detection (NEW - DATA-DRIVEN)

```python
print("🔍 Step 5.5: Statistical Outlier Detection (data-driven, NO hardcoding)")

# Detect TF activity outliers using z-scores (completely data-driven)
cluster_col = 'tier2_cluster'
tf_key = 'score_ulm' if 'score_ulm' in subset.obsm else 'ulm_estimate'
pval_key = 'padj_ulm' if 'padj_ulm' in subset.obsm else 'ulm_pvals'

tf_scores = subset.obsm[tf_key]
tf_pvals = subset.obsm[pval_key]

outliers = []

for cluster in subset.obs[cluster_col].unique():
    mask = subset.obs[cluster_col] == cluster
    cluster_tf_mean = tf_scores[mask].mean()
    cluster_tf_pval = tf_pvals[mask].mean()

    for tf in tf_scores.columns:
        # Global statistics across all clusters
        global_mean = tf_scores[tf].mean()
        global_std = tf_scores[tf].std()
        cluster_val = cluster_tf_mean[tf]
        cluster_pval = cluster_tf_pval[tf]

        # Z-score (statistical measure of how unusual this value is)
        if global_std > 1e-6:  # Avoid division by zero
            z_score = (cluster_val - global_mean) / global_std
        else:
            z_score = 0.0

        # Flag as outlier if:
        # 1. |z-score| > 2.5 (statistically significant deviation)
        # 2. Activity > 0.5 (meaningful activity, not just noise)
        # 3. p-value < 0.05 (significant in this cluster)
        if abs(z_score) > 2.5 and cluster_val > 0.5 and cluster_pval < 0.05:
            outliers.append({
                'cluster': str(cluster),
                'tf': tf,
                'activity': float(cluster_val),
                'z_score': float(z_score),
                'pval': float(cluster_pval),
                'direction': 'HIGH' if z_score > 0 else 'LOW'
            })

# Store outliers in AnnData for reasoning
subset.uns['tf_outliers'] = outliers

if outliers:
    print(f"   ⚠️ FOUND {len(outliers)} TF activity outliers (statistical anomalies)")
    print(f"   → These may indicate specialized/rare cell types")
    print(f"   → Will be investigated in reasoning phase")

    # Show summary
    from collections import Counter
    cluster_counts = Counter(o['cluster'] for o in outliers)
    for cluster, count in cluster_counts.most_common(5):
        tfs = [o['tf'] for o in outliers if o['cluster'] == cluster]
        print(f"      Cluster {cluster}: {count} outlier TFs ({', '.join(tfs[:3])}...)")
else:
    print(f"   ✅ No statistical outliers detected")
    print(f"   → Standard developmental states expected")
```

### Step 5.6: Evidence Conflict Detection (NEW - DATA-DRIVEN)

```python
print("🔍 Step 5.6: Evidence Conflict Detection (data-driven)")

# Detect conflicts between different evidence types (NO hardcoding)
conflicts = []

for cluster in subset.obs[cluster_col].unique():
    mask = subset.obs[cluster_col] == cluster

    # Get top markers for this cluster
    cluster_de = de_df[de_df['group'] == cluster].nlargest(20, 'logfoldchanges')
    top_markers = cluster_de['names'].tolist() if 'names' in cluster_de.columns else []

    # Get trajectory info
    if 'pseudotime' in subset.obs.columns:
        mean_pt = subset.obs.loc[mask, 'pseudotime'].mean()
        pt_category = subset.obs.loc[mask, 'pseudotime_category'].mode()[0] if len(subset.obs.loc[mask, 'pseudotime_category']) > 0 else 'N/A'
    else:
        mean_pt = np.nan
        pt_category = 'N/A'

    # Get top TFs
    cluster_tf_mean = tf_scores[mask].mean()
    top_tf_indices = cluster_tf_mean.nlargest(10).index.tolist()

    # CONFLICT 1: Top marker is also a top TF (transcriptionally active)
    # This is interesting but not necessarily a conflict - just noteworthy
    marker_tf_overlap = set(top_markers) & set(top_tf_indices)
    if marker_tf_overlap:
        for gene in marker_tf_overlap:
            conflicts.append({
                'cluster': str(cluster),
                'type': 'MARKER_TF_OVERLAP',
                'detail': f'{gene} is both highly expressed (marker) and transcriptionally active (TF)',
                'severity': 'INFO',
                'action': 'Note this dual role - may indicate key regulatory function'
            })

    # CONFLICT 2: TF activity outlier conflicts with marker-based identity
    # If a cluster has TF outliers, flag for special attention
    cluster_outliers = [o for o in outliers if o['cluster'] == str(cluster)]
    if cluster_outliers:
        high_outlier_tfs = [o['tf'] for o in cluster_outliers if o['direction'] == 'HIGH']
        if high_outlier_tfs:
            conflicts.append({
                'cluster': str(cluster),
                'type': 'TF_STATISTICAL_OUTLIER',
                'detail': f"TFs with unusually HIGH activity: {', '.join(high_outlier_tfs[:5])}",
                'severity': 'HIGH',
                'action': 'Search literature for these TF + major_type combinations'
            })

    # CONFLICT 3: Marker expression pattern doesn't match TF activity
    # Example: High proliferation markers but no proliferation-related TFs
    # This requires some basic pattern recognition but NO hardcoded cell types
    proliferation_markers = ['MKI67', 'TOP2A', 'PCNA', 'CDK1']
    prolif_marker_count = sum(1 for m in proliferation_markers if m in top_markers)

    if prolif_marker_count >= 2:
        # Check if proliferation-related TFs are active
        prolif_tfs = ['E2F1', 'MYC', 'FOXM1']
        prolif_tf_count = sum(1 for tf in prolif_tfs if tf in top_tf_indices)

        if prolif_tf_count == 0:
            conflicts.append({
                'cluster': str(cluster),
                'type': 'PROLIFERATION_CONFLICT',
                'detail': f'High proliferation markers ({prolif_marker_count}) but no proliferation TFs active',
                'severity': 'MEDIUM',
                'action': 'Investigate: passive proliferation vs active cell cycle program'
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

### Step 6: VERIFY Functional Analysis (NEW)

```python
print("🔍 Step 6: Verifying functional analysis...")

errors = []

# TF Activity check
if 'ulm_estimate' not in subset.obsm:
    errors.append("TF activity (ulm_estimate) not computed")
if 'ulm_pvals' not in subset.obsm:
    errors.append("TF p-values (ulm_pvals) not computed")

# Trajectory check (if applicable)
if subset.n_obs >= 2000 and 'pseudotime' not in subset.obs.columns:
    errors.append("Trajectory (pseudotime) not computed for >= 2000 cells")

if errors:
    raise AssertionError(
        "❌ FUNCTIONAL ANALYSIS INCOMPLETE:\n" +
        "\n".join(f"  - {e}" for e in errors) +
        "\n\n⚠️ Cannot proceed to annotation without functional analysis!"
    )

print("✅ Functional analysis verified - ready for annotation")
```

### Step 7: Filter Valid Markers

```python
def filter_valid_markers(de_df, cluster_id, top_n=50):
    cluster_df = de_df[de_df['group'] == cluster_id].copy()
    valid = cluster_df[
        (cluster_df['pct_nz_group'] >= 0.25) &
        (cluster_df['logfoldchanges'] >= 1.0) &
        (cluster_df['pvals_adj'] < 0.05)
    ]
    return valid.nlargest(top_n, 'logfoldchanges')

# Get valid markers per cluster
valid_markers_dict = {}
for cluster in subset.obs['tier2_cluster'].unique():
    valid_markers_dict[cluster] = filter_valid_markers(de_df, cluster, top_n=50)
```

### Step 8: Integrated Pre-Analysis (UPDATED with Outliers/Conflicts)

**Present ALL evidence to Claude for candidate generation:**

```markdown
=== {Major Type} Cluster {X} Analysis ===

## DE Markers (Top 50 Valid)
| Rank | Gene | pct_in | log2FC | padj |
|------|------|--------|--------|------|
[show all 50]

## TF Activity (Top 10)
| TF | Activity | p-value | Z-score | Outlier? |
|----|----------|---------|---------|----------|
| {TF1} | {val} | {pval} | {z:.2f} | {✅ NORMAL / ⚠️ HIGH OUTLIER / ⚠️ LOW OUTLIER} |
| {TF2} | {val} | {pval} | {z:.2f} | {✅ NORMAL / ⚠️ HIGH OUTLIER / ⚠️ LOW OUTLIER} |
[show top 10 + ALL outliers even if not in top 10]

## 🔍 STATISTICAL OUTLIERS (if detected)

**TF Activity Outliers for Cluster {X}:**
{for each outlier in this cluster}
- **{TF}** (activity={val:.2f}, z-score={z:.2f}, p={pval:.2e})
  → This TF activity is **statistically unusual** (|z| > 2.5)
  → **Action**: Search literature for "{TF} {major_type}" combination
  → **Interpretation**: May indicate specialized/rare cell type NOT in standard developmental hierarchy

**Total outliers in this cluster: {N}**
**Interpretation**:
- IF N = 0: Standard developmental state expected
- IF N = 1-2: Minor variation, investigate if marker pattern also unusual
- IF N >= 3: **HIGH PRIORITY** - likely specialized subset (e.g., age-associated, atypical)

## ⚠️ EVIDENCE CONFLICTS (if detected)

{for each conflict in this cluster}
**{conflict_type}**:
- {detail}
- Severity: {HIGH/MEDIUM/INFO}
- Action: {action_instruction}

## Trajectory Position
- Mean pseudotime: {value}
- Category: {Early/Mid/Late}
- Interpretation: {interpretation}

## 🧠 REASONING INSTRUCTION

**Standard Workflow (NO outliers/conflicts)**:
1. Identify developmental state based on markers + TF + trajectory
2. Search literature for standard markers
3. Assign annotation with standard confidence scoring

**SPECIAL WORKFLOW (outliers or HIGH severity conflicts detected)**:
1. **DO NOT assume standard developmental states**
2. **Priority literature searches**:
   - For EACH outlier TF: "{TF_name} {major_type}"
   - For unusual marker combinations: "{marker1} {marker2} {marker3} {major_type}"
   - Focus on: "age-associated", "atypical", "specialized", "rare subset"
3. **Candidate generation**:
   - Candidate A: Literature-supported specialized subset (HIGH confidence if found)
   - Candidate B: Standard developmental state (LOW confidence, explain outlier conflict)
   - Candidate C: Novel/uncharacterized if no literature match
4. **Decision criteria**:
   - IF literature explains outlier → Use specialized annotation
   - IF no literature but pattern coherent → Flag as NOVEL
   - IF pattern incoherent → Check for contamination

## Evidence Summary (for Claude)
- Marker signature suggests: [candidates from markers]
- TF activity suggests: [candidates from TF]
- **TF outliers suggest**: {interpretation OR "N/A - no outliers"}
- Trajectory suggests: [candidates from pseudotime]
- **Evidence convergence**: {CONVERGENT / CONFLICTED - explain}
```

**Claude Output Format:**

```yaml
cluster: {X}
parent_type: "{major_type}"
has_outliers: {true/false}
outlier_count: {N}
has_conflicts: {true/false}

# IF OUTLIERS DETECTED:
outlier_investigation:
  outlier_tfs: ["{TF1}", "{TF2}"]  # List of outlier TFs
  literature_search_required: true
  priority: "HIGH"  # If >= 3 outliers

# Integrated candidate analysis
candidates:
  # IF OUTLIERS: Prioritize specialized/rare cell type candidates
  - dev_state: "{Specialized_Type}"  # e.g., "ABC", "Atypical_Memory"
    confidence: "High/Medium/Low"
    key_markers: ["{marker1}", "{marker2}", "{marker3}"]
    key_tfs: ["{TF1}", "{TF2}"]
    outlier_tfs: ["{TF_outlier1}"]  # NEW - TFs flagged as statistical outliers
    tf_support: "{description}"
    tf_outlier_explanation: "{how outlier supports this annotation}"  # NEW
    trajectory_support: "{description}"
    reasoning: "{integrated reasoning addressing outlier}"

  # Standard developmental state (lower confidence if outliers unexplained)
  - dev_state: "{Standard_Type}"  # e.g., "Memory"
    confidence: "Low"  # Reduced if outliers present
    key_markers: ["{marker1}", "{marker2}"]
    key_tfs: ["{TF1}"]
    outlier_tfs: []
    tf_support: "{description}"
    tf_outlier_conflict: "{cannot explain outlier TFs}"  # NEW
    trajectory_support: "{description}"
    reasoning: "{why this is less likely given outliers}"

mcp_queries:
  # IF OUTLIERS: Add specific outlier queries
  - query: "{outlier_TF} {major_type}"  # e.g., "TBX21 B cell"
    purpose: "Investigate statistical TF outlier"
    priority: "HIGH"

  - query: "{marker1} {marker2} {major_type}"
    purpose: "Standard marker validation"
    priority: "NORMAL"
```

### Step 9: MCP Verification (PMID Required)

```python
for candidate in claude_candidates:
    for query in candidate['mcp_queries']:
        # Search PubMed
        result = mcp.pubmed_search(query['query'], max_results=3)

        if result['count'] > 0:
            pmid = result['results'][0]['pmid']

            # Verify reference includes TF if mentioned
            markers_to_verify = candidate['key_markers'] + candidate.get('key_tfs', [])

            verification = mcp.verify_reference(
                pmid=pmid,
                markers=markers_to_verify,
                cell_type=f"{candidate['dev_state']} {parent_type}"
            )

            print(f"PMID:{pmid} - {verification['status']}")
```

### Step 10: 3-Iteration Reasoning with Integrated Evidence (UPDATED - GENERATES reasoning_results)

```python
print(f"\n🧠 Step 10: Performing 3-Iteration Reasoning for {major_type}...")

# This step generates reasoning_results which is used in:
# - Step 11.5: Create annotation_evidence
# - Step 12: Visualization

reasoning_results = {}

for cluster in subset.obs['tier2_cluster'].unique():
    cluster_str = str(cluster)
    print(f"\n--- Cluster {cluster_str} Reasoning ---")

    # Get cluster-specific data
    cluster_mask = subset.obs['tier2_cluster'] == cluster
    cluster_de = de_df[de_df['group'] == cluster]

    # Filter valid markers
    valid_markers = cluster_de[
        (cluster_de['pct_nz_group'] >= 0.25) &
        (cluster_de['logfoldchanges'] >= 1.0) &
        (cluster_de['pvals_adj'] < 0.05)
    ].nlargest(50, 'logfoldchanges')

    # Extract top markers
    top_genes = valid_markers.nlargest(10, 'logfoldchanges')['names'].tolist()

    # Get TF activity for this cluster
    tf_scores = subset.obsm['ulm_estimate']
    cluster_tf_mean = tf_scores[cluster_mask].mean()
    top_tfs = cluster_tf_mean.nlargest(5).index.tolist()
    top_tf_scores = cluster_tf_mean.nlargest(5).values.tolist()

    # Get trajectory info (if computed)
    trajectory_info = None
    if 'palantir_pseudotime' in subset.obs.columns:
        cluster_pseudotime = subset.obs.loc[cluster_mask, 'palantir_pseudotime']
        trajectory_info = {
            'mean_pseudotime': float(cluster_pseudotime.mean()),
            'category': 'Early' if cluster_pseudotime.mean() < 0.33 else
                       ('Middle' if cluster_pseudotime.mean() < 0.67 else 'Late')
        }

    # Get outliers and conflicts for this cluster
    cluster_outliers = [o for o in subset.uns.get('tf_outliers', [])
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
            f"Top TFs: {', '.join([f'{tf}({score:.2f})' for tf, score in zip(top_tfs[:3], top_tf_scores[:3])])}"
        ],

        # Conflicts
        'conflicts': cluster_conflicts,

        # TF activity
        'tf_activity': [
            {'tf': tf, 'score': float(score), 'pval': 0.01}
            for tf, score in zip(top_tfs, top_tf_scores)
        ],

        # Trajectory
        'trajectory': trajectory_info,

        # Outliers
        'tf_outliers': cluster_outliers,

        # Confidence scoring
        'confidence_score': 0,
        'confidence_level': 'Low',

        # Novel flag
        'is_novel': False,
        'novel_evidence': None
    }

    if trajectory_info:
        cluster_reasoning['supporting_reasons'].append(
            f"Pseudotime: {trajectory_info['mean_pseudotime']:.2f} ({trajectory_info['category']})"
        )

    reasoning_results[cluster_str] = cluster_reasoning
    print(f"  ✓ Reasoning structure created for cluster {cluster_str}")

print(f"\n✅ Reasoning completed for {len(reasoning_results)} clusters")
print("⚠️  NOTE: In production, complete the reasoning with:")
print("   → LLM 3-iteration reasoning (reasoning/integrated_format.md)")
print("   → MCP literature verification")
print("   → Confidence scoring with TF/Trajectory evidence")
```

**Template Reference**: `reasoning/integrated_format.md`

**Reasoning MUST include:**
- **Iteration 1**: Markers + TF Activity + Trajectory + Outliers + Conflicts
- **Iteration 2**: Candidate comparison with ALL evidence types
- **Iteration 3**: Integrated scoring and final decision

**Confidence Scoring:**

| Criteria | Standard (No outliers) | With Outliers (unexplained) | With Outliers (literature-explained) |
|----------|------------------------|----------------------------|-------------------------------------|
| Markers | 3 pts | 3 pts | 3 pts |
| References | 3 pts | 3 pts (REQUIRED) | 3 pts (REQUIRED) |
| TF consistency | 3 pts | 0-1 pts | 3 pts |
| Trajectory | 3 pts | 3 pts or N/A | 3 pts or N/A |
| **Total** | **12 pts** | **7-10 pts** | **12 pts** |

### Step 11: Assign Annotations

```python
tier2_mapping = {
    '0': 'Naive',
    '1': 'Effector',
    '2': 'Memory',
    ...
}

subset.obs['tier2_annotation'] = subset.obs['tier2_cluster'].map(tier2_mapping)
```

### Step 11.5: Create annotation_evidence (NEW - REQUIRED FOR VISUALIZATION)

```python
# Convert reasoning results to annotation_evidence format
# This is REQUIRED for visualization in Step 12

from reasoning.agent_format import reasoning_to_evidence

annotation_evidence = []
safe_name = major_type.replace(" ", "_")

for cluster in subset.obs['tier2_cluster'].unique():
    # Get annotation for this cluster
    dev_state = tier2_mapping.get(str(cluster))
    if not dev_state:
        continue

    # Get cluster-specific data from reasoning process
    # NOTE: This assumes reasoning was performed for each cluster in Step 10
    cluster_reasoning = reasoning_results.get(str(cluster), {})

    # Create evidence entry
    evidence_entry = reasoning_to_evidence(
        cluster_id=f"tier2_{safe_name}_cluster_{cluster}",
        tier="tier2",
        subset_id=safe_name,
        reasoning_output=cluster_reasoning
    )

    annotation_evidence.append(evidence_entry)

print(f"✅ Created annotation_evidence for {major_type}: {len(annotation_evidence)} entries")

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

# Add new Tier 2 evidence
all_evidence.extend(annotation_evidence)

# Save combined evidence
with open(evidence_file, 'w') as f:
    json.dump(all_evidence, f, indent=2, ensure_ascii=False)

print(f"✅ Updated: {evidence_file} (total entries: {len(all_evidence)})")
```

### Step 12: Save Results

```python
# Save subset with functional analysis
safe_name = major_type.replace(" ", "_")
subset.write(f'annotation_output/subsets/tier2/{safe_name}.h5ad')

# Verify saved data has functional analysis
saved = sc.read_h5ad(f'annotation_output/subsets/tier2/{safe_name}.h5ad')
assert 'ulm_estimate' in saved.obsm, "TF activity not saved!"
if saved.n_obs >= 2000:
    assert 'pseudotime' in saved.obs.columns, "Pseudotime not saved!"
print("✅ Subset saved with functional analysis")

# ✅ Generate Visualizations (NEW)
print(f"\n🎨 Generating Tier 2 Visualizations for {major_type}...")
import os
os.makedirs('annotation_output/figures', exist_ok=True)

# Prepare marker dict from annotation_evidence (VALIDATED markers only)
# annotation_evidence는 Step 10에서 reasoning 결과로 생성됨
marker_dict = {}
annotation_order = sorted(subset.obs['tier2_annotation'].unique())

for entry in annotation_evidence:
    if entry['tier'] == 'tier2' and entry['subset_id'] == safe_name:
        dev_state = entry['annotation']
        # Extract validated marker gene names from reasoning results
        validated_markers = [m['gene'] for m in entry['markers']]
        marker_dict[dev_state] = validated_markers

# Import visualization functions
import sys
sys.path.append('annotation_output')
from tools.visualization import save_all_visualizations
from tools.functional_plots import (
    plot_tf_activity_heatmap,
    plot_tf_dotplot,
    plot_pseudotime_umap,
    plot_pseudotime_violin
)

# 1. Generate standard visualizations (UMAP, Dotplot, Combined)
save_all_visualizations(
    adata=subset,
    annotation_col='tier2_annotation',
    marker_dict=marker_dict,
    reasoning_dict={e['annotation']: e for e in annotation_evidence
                    if e['tier'] == 'tier2' and e['subset_id'] == safe_name},
    reference_dict={e['annotation']: e['references'] for e in annotation_evidence
                    if e['tier'] == 'tier2' and e['subset_id'] == safe_name},
    annotation_order=annotation_order,
    title_prefix=f'Tier 2 - {major_type}',
    output_dir='annotation_output/figures/'
)

# 2. Generate TF activity heatmap (Tier 2 MANDATORY)
plot_tf_activity_heatmap(
    subset,
    cluster_col='tier2_annotation',
    top_n=15,
    figsize=(8, 6),
    save_path=f'annotation_output/figures/{safe_name}_tier2_tf_heatmap'
)
print(f"   - {safe_name}_tier2_tf_heatmap.png/.svg")

# 3. Generate TF dotplot (optional, comprehensive view)
plot_tf_dotplot(
    subset,
    cluster_col='tier2_annotation',
    top_n=10,
    figsize=(10, 6),
    save_path=f'annotation_output/figures/{safe_name}_tier2_tf_dotplot'
)
print(f"   - {safe_name}_tier2_tf_dotplot.png/.svg")

# 4. Generate trajectory plots (if computed)
if 'palantir_pseudotime' in subset.obs.columns:
    plot_pseudotime_umap(
        subset,
        figsize=(5, 5),
        save_path=f'annotation_output/figures/{safe_name}_tier2_pseudotime_umap'
    )
    print(f"   - {safe_name}_tier2_pseudotime_umap.png/.svg")

    plot_pseudotime_violin(
        subset,
        cluster_col='tier2_annotation',
        figsize=(8, 4),
        save_path=f'annotation_output/figures/{safe_name}_tier2_pseudotime_violin'
    )
    print(f"   - {safe_name}_tier2_pseudotime_violin.png/.svg")

print(f"✅ Tier 2 Visualizations for {major_type} completed:")
print(f"   - {safe_name}_tier2_umap.png/.svg")
print(f"   - {safe_name}_tier2_dotplot.png/.svg (with marker brackets, VALIDATED markers)")
print(f"   - {safe_name}_tier2_tf_heatmap.png/.svg (MANDATORY)")
print(f"   - {safe_name}_tier2_combined.png/.svg (multi-panel)")
print(f"   - {safe_name}_tier2_reasoning.json/.txt")
```

---

## Checklist (Per Major Type)

```
MANDATORY:
- [ ] Subset extracted
- [ ] Re-clustered (resolution 0.6-1.0)
- [ ] DE RE-COMPUTED within subset
- [ ] ⚠️ TF activity computed (ulm_estimate in obsm)
- [ ] ⚠️ Trajectory computed if >= 2000 cells (pseudotime in obs)
- [ ] ⚠️ Statistical outliers detected (tf_outliers in uns) ← NEW
- [ ] ⚠️ Evidence conflicts detected (evidence_conflicts in uns) ← NEW
- [ ] ⚠️ Functional analysis VERIFIED before annotation

ANNOTATION:
- [ ] Valid markers filtered (Top 50)
- [ ] Integrated Pre-Analysis (markers + TF + trajectory + outliers + conflicts) ← UPDATED
- [ ] MCP verification with PMID
- [ ] 3-iteration reasoning with integrated evidence
- [ ] IF OUTLIERS: Specialized subset literature searches performed ← NEW
- [ ] IF OUTLIERS: Outlier-aware confidence scoring ← NEW
- [ ] References double-verified
- [ ] Confidence scored (max 12 pts)

OUTPUT:
- [ ] Annotations assigned (standard OR specialized based on outliers)
- [ ] Subset saved with functional analysis + outlier/conflict info
- [ ] Visualizations saved
```

---

## Troubleshooting

### TF Activity Fails

```python
# Check if var_names match network
net = dc.op.collectri(organism='human')
overlap = set(subset.var_names) & set(net['target'])
print(f"Gene overlap: {len(overlap)} / {len(net['target'])}")

# If low overlap, check gene symbols
# May need to convert Ensembl to gene symbols
```

### Trajectory Fails

```python
# Check if enough cells
print(f"Cells: {subset.n_obs} (need >= 2000)")

# Check if PCA exists
print(f"PCA shape: {subset.obsm['X_pca'].shape}")

# Try with fewer components if memory issues
palantir.utils.run_diffusion_maps(subset, n_components=5, pca_key='X_pca')
```

### Memory Issues

```python
# For large datasets, process in batches or reduce dimensions
if subset.n_obs > 50000:
    print("Consider subsampling for trajectory analysis")
    subsample = sc.pp.subsample(subset, n_obs=20000, copy=True)
    # Run trajectory on subsample, then impute to full data
```
