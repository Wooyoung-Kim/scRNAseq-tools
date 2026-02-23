# Phase 1: Tier 1 - Major Cell Types (v3 - DATA-DRIVEN)

```
╔══════════════════════════════════════════════════════════════════════╗
║  GOAL: "이 세포가 어떤 계통(lineage)에 속하는가?"                      ║
║                                                                      ║
║  분류 대상: T cells, B cells, NK cells, Myeloid, Epithelial, etc.   ║
║  분류 기준: Data-driven DE markers + Statistical outlier detection   ║
║                                                                      ║
║  🆕 v3 CHANGES:                                                       ║
║  - NO hardcoded lineage marker lists                                ║
║  - Statistical outlier detection (marker expression z-scores)       ║
║  - Automated cross-lineage contamination detection                  ║
║  - Doublet detection via mixed lineage signatures                   ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## 목표 명확화

| 이 Tier에서 결정 | 이 Tier에서 결정 안함 |
|------------------|----------------------|
| T cell vs B cell vs Myeloid | Naive vs Effector (Tier 2) |
| 세포 계통 정체성 | 발달 단계 (Tier 2) |
| Lineage markers | 기능 상태 (Tier 3) |

---

## Input

```python
# QC 완료된 full dataset
adata = sc.read_h5ad('preprocessed_data.h5ad')

# 필요 조건
assert 'X_pca' in adata.obsm or 'X_pca_harmony' in adata.obsm
assert adata.n_obs > 0
```

## Output

```python
adata.obs['tier1_annotation']  # Major cell type labels
# 예: "T cells", "B lineage", "Myeloid", "NK cells"
```

---

## Workflow

### Step 1: Clustering

```python
import scanpy as sc

# Use harmony-corrected embedding if available
use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'

sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, n_pcs=50)
sc.tl.leiden(adata, resolution=0.8, key_added='tier1_cluster')

print(f"Clusters: {adata.obs['tier1_cluster'].nunique()}")
```

### Step 2: Compute DE

```python
sc.tl.rank_genes_groups(adata, groupby='tier1_cluster', method='wilcoxon')
de_df = sc.get.rank_genes_groups_df(adata, group=None)

print(f"⚠️ DE computed for subset: full_dataset ({adata.n_obs} cells)")
```

### Step 2.5: Statistical Outlier Detection (NEW - DATA-DRIVEN)

```python
print("🔍 Step 2.5: Statistical Outlier Detection (data-driven, NO hardcoding)")

# Detect marker expression outliers using z-scores
import numpy as np
import pandas as pd

cluster_col = 'tier1_cluster'

# Compute marker expression statistics per cluster
marker_stats = {}
for cluster in adata.obs[cluster_col].unique():
    mask = adata.obs[cluster_col] == cluster
    cluster_cells = adata[mask]

    # Get expression for all genes
    if hasattr(cluster_cells.X, 'toarray'):
        expr = cluster_cells.X.toarray()
    else:
        expr = cluster_cells.X

    # Mean expression per gene
    mean_expr = expr.mean(axis=0)
    marker_stats[cluster] = pd.Series(mean_expr, index=adata.var_names)

# Global statistics
global_means = pd.concat(marker_stats.values(), axis=1).mean(axis=1)
global_stds = pd.concat(marker_stats.values(), axis=1).std(axis=1)

# Detect outliers
outliers = []
for cluster, cluster_expr in marker_stats.items():
    for gene in cluster_expr.index:
        if global_stds[gene] > 1e-6:  # Avoid division by zero
            z_score = (cluster_expr[gene] - global_means[gene]) / global_stds[gene]

            # Flag if expression is unusually high (z > 2.5) and meaningful (> 0.1)
            if z_score > 2.5 and cluster_expr[gene] > 0.1:
                outliers.append({
                    'cluster': str(cluster),
                    'gene': gene,
                    'expression': float(cluster_expr[gene]),
                    'z_score': float(z_score),
                    'direction': 'HIGH'
                })

# Store outliers
adata.uns['marker_outliers'] = outliers

if outliers:
    print(f"   ⚠️ FOUND {len(outliers)} marker expression outliers")
    print(f"   → These may indicate specialized cell types or contamination")

    # Show summary
    from collections import Counter
    cluster_counts = Counter(o['cluster'] for o in outliers)
    for cluster, count in cluster_counts.most_common(5):
        genes = [o['gene'] for o in outliers if o['cluster'] == cluster][:5]
        print(f"      Cluster {cluster}: {count} outlier markers ({', '.join(genes)}...)")
else:
    print(f"   ✅ No statistical outliers detected")
```

### Step 2.6: Cross-Lineage Contamination Detection (NEW - DATA-DRIVEN)

```python
print("🔍 Step 2.6: Cross-Lineage Contamination Detection (data-driven)")

# Detect cross-lineage contamination without hardcoded marker lists
# Strategy: Use top DE markers from each cluster as "lineage signatures"

contamination_issues = []

for cluster in adata.obs[cluster_col].unique():
    mask = adata.obs[cluster_col] == cluster
    cluster_cells = adata[mask]

    # Get top 10 DE markers for this cluster
    cluster_de = de_df[de_df['group'] == cluster].nlargest(10, 'logfoldchanges')
    top_markers = cluster_de['names'].tolist()

    # For each OTHER cluster, check if those markers are suspiciously high
    for other_cluster in adata.obs[cluster_col].unique():
        if other_cluster == cluster:
            continue

        other_mask = adata.obs[cluster_col] == other_cluster
        other_cells = adata[other_mask]

        # Check how many top markers are expressed in the other cluster
        contamination_score = 0
        contaminating_markers = []

        for marker in top_markers:
            if marker in adata.var_names:
                # Expression in other cluster
                if hasattr(other_cells.X, 'toarray'):
                    expr = other_cells[:, marker].X.toarray().flatten()
                else:
                    expr = other_cells[:, marker].X.flatten()

                pct_expr = (expr > 0).mean()

                # If > 10% of other cluster expresses this marker, potential contamination
                if pct_expr > 0.10:
                    contamination_score += pct_expr
                    contaminating_markers.append(f"{marker} ({pct_expr:.1%})")

        # Flag if contamination score is high
        if contamination_score > 0.5:  # More than 50% contamination signal
            contamination_issues.append({
                'source_cluster': str(cluster),
                'target_cluster': str(other_cluster),
                'contamination_score': float(contamination_score),
                'markers': contaminating_markers[:3],  # Top 3
                'severity': 'HIGH' if contamination_score > 1.0 else 'MEDIUM'
            })

# Store contamination info
adata.uns['contamination_issues'] = contamination_issues

if contamination_issues:
    high_priority = [c for c in contamination_issues if c['severity'] == 'HIGH']
    print(f"   ⚠️ FOUND {len(contamination_issues)} potential contamination issues "
          f"({len(high_priority)} high priority)")

    # Show high priority
    for issue in high_priority[:5]:
        print(f"      Cluster {issue['target_cluster']} contaminated by Cluster {issue['source_cluster']}")
        print(f"         Markers: {', '.join(issue['markers'])}")
else:
    print(f"   ✅ No cross-lineage contamination detected")
```

### Step 3: Filter Valid Markers

```python
def filter_valid_markers(de_df, cluster_id, top_n=50):
    cluster_df = de_df[de_df['group'] == cluster_id].copy()
    valid = cluster_df[
        (cluster_df['pct_nz_group'] >= 0.25) &
        (cluster_df['logfoldchanges'] >= 1.0) &
        (cluster_df['pvals_adj'] < 0.05)
    ]
    return valid.nlargest(top_n, 'logfoldchanges')
```

### Step 4: Identify Major Cell Types (DATA-DRIVEN with Dynamic Knowledge)

**⚠️ NO HARDCODED MARKER LISTS - USE DYNAMIC KNOWLEDGE BASE**

Strategy:
1. Extract top DE markers from current dataset
2. Use dynamic knowledge base to search literature
3. Get cell type candidates from PubMed (data-driven)
4. Rank candidates by confidence
5. Flag novel populations if no match

```python
from tools.dynamic_knowledge import (
    search_cell_type_from_markers,
    verify_cell_type_assignment,
    detect_novel_population
)

# For each cluster, identify cell type via dynamic literature search
def identify_cluster_cell_type(adata, cluster_id, de_df, outliers, contamination):
    """
    Identify cell type using dynamic knowledge base.
    NO hardcoded cell types or markers.
    """
    # Get valid markers
    cluster_de = de_df[de_df['group'] == cluster_id]
    valid_markers = cluster_de[
        (cluster_de['pct_nz_group'] >= 0.25) &
        (cluster_de['logfoldchanges'] >= 1.0) &
        (cluster_de['pvals_adj'] < 0.05)
    ].nlargest(50, 'logfoldchanges')

    top_markers = valid_markers['names'].tolist()

    # Get outliers for this cluster
    cluster_outliers = [o for o in outliers if o['cluster'] == str(cluster_id)]

    # Get contamination info
    cluster_contamination = [
        c for c in contamination
        if c['target_cluster'] == str(cluster_id) or c['source_cluster'] == str(cluster_id)
    ]

    print(f"\n🔍 Cluster {cluster_id}: Dynamic cell type search...")
    print(f"   Top markers: {top_markers[:5]}")

    # DYNAMIC KNOWLEDGE SEARCH (NO hardcoding)
    candidates = search_cell_type_from_markers(
        markers=top_markers[:10],
        tissue=adata.uns.get('tissue', None),
        species=adata.uns.get('species', 'human'),
        max_combinations=10
    )

    if not candidates or candidates[0].is_novel:
        print(f"   ⚠️ Novel population detected!")
        print(f"   → No strong literature match found")
        print(f"   → Flag for manual review and deeper literature search")

        return {
            'cluster_id': cluster_id,
            'n_cells': (adata.obs['tier1_cluster'] == cluster_id).sum(),
            'top_markers': valid_markers.to_dict('records')[:10],
            'cell_type_candidates': [],
            'is_novel': True,
            'outliers': cluster_outliers,
            'contamination': cluster_contamination
        }

    # Verify top candidate
    top_candidate = candidates[0]
    is_valid, reason = verify_cell_type_assignment(
        cluster_markers=top_markers,
        candidate=top_candidate,
        min_marker_overlap=2
    )

    print(f"   ✓ Top candidate: {top_candidate.name}")
    print(f"   → Confidence: {top_candidate.confidence:.2f}")
    print(f"   → Supporting markers: {top_candidate.supporting_markers}")
    print(f"   → PMIDs: {top_candidate.pmids}")
    print(f"   → Valid: {is_valid} ({reason})")

    return {
        'cluster_id': cluster_id,
        'n_cells': (adata.obs['tier1_cluster'] == cluster_id).sum(),
        'top_markers': valid_markers.to_dict('records')[:10],
        'cell_type_candidates': [
            {
                'name': c.name,
                'confidence': c.confidence,
                'supporting_markers': c.supporting_markers,
                'pmids': c.pmids,
                'evidence': c.evidence
            }
            for c in candidates[:3]  # Top 3 candidates
        ],
        'is_novel': False,
        'outliers': cluster_outliers,
        'contamination': cluster_contamination
    }
```

**Example Output:**
```
🔍 Cluster 0: Dynamic cell type search...
   Top markers: ['CD3D', 'CD3E', 'TRAC', 'CD8A', 'GZMB']
   ✓ Top candidate: T cells
   → Confidence: 0.92
   → Supporting markers: ['CD3D', 'CD3E', 'TRAC', 'CD8A']
   → PMIDs: ['12345678', '23456789', '34567890', '45678901']
   → Valid: True (Valid)

🔍 Cluster 5: Dynamic cell type search...
   Top markers: ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5']
   ⚠️ Novel population detected!
   → No strong literature match found
   → Flag for manual review and deeper literature search
```

### Step 5: LLM Pre-Analysis (UPDATED - Outlier-Aware)

```markdown
=== Cluster {X} Analysis (Tier 1) ===

## Top 50 DE Markers
| Rank | Gene | pct_in | log2FC | padj | Z-score | Outlier? |
|------|------|--------|--------|------|---------|----------|
| 1 | CD3D | 0.95 | 4.5 | 1e-100 | 3.2 | ⚠️ HIGH OUTLIER |
| 2 | CD3E | 0.92 | 4.2 | 1e-95 | 2.8 | ⚠️ HIGH OUTLIER |
| 3 | TRAC | 0.88 | 3.8 | 1e-80 | 2.6 | ⚠️ HIGH OUTLIER |
...

## 🔍 STATISTICAL OUTLIERS (if detected)

**Marker Expression Outliers for Cluster {X}:**
- **CD3D** (expr=0.95, z=3.2) - Unusually HIGH expression
- **CD3E** (expr=0.92, z=2.8) - Unusually HIGH expression
- **TRAC** (expr=0.88, z=2.6) - Unusually HIGH expression

**Total outliers: 3**
**Interpretation**:
→ Multiple outliers suggest this is a distinct major lineage
→ Search literature: "CD3D CD3E TRAC" to identify lineage

## ⚠️ CONTAMINATION ALERTS (if detected)

**Potential Contamination:**
- Source Cluster 2 → This Cluster: CD79A (15%), MS4A1 (12%)
  → Severity: MEDIUM
  → Action: Check if this is doublet contamination or mixed population

## 🧠 REASONING INSTRUCTION

**Standard Workflow (NO outliers/contamination)**:
1. Identify lineage based on top DE markers
2. Search literature for marker combinations
3. Assign major cell type annotation

**SPECIAL WORKFLOW (outliers or contamination detected)**:
1. **Priority investigation**:
   - For EACH outlier marker combination: Search "{marker1} {marker2} {marker3}"
   - Focus on lineage-defining marker combinations
2. **Contamination handling**:
   - IF contamination detected: Check if doublets (multiple lineage markers)
   - Search: "doublet {marker1} {marker2}" or specific lineages
3. **Decision**:
   - IF single lineage signature → Assign major type
   - IF mixed lineage signatures + high contamination → Flag as doublet/contamination
   - IF novel combination → Search literature, may be rare cell type

## Candidate Analysis

Based on outliers and DE:
- Marker combination: CD3D + CD3E + TRAC (all HIGH outliers)
- Literature search required: "CD3D CD3E TRAC lineage"
- Contamination status: NONE / {details}
```

### Step 6: MCP Verification (PMID Required)

```python
# Search for lineage marker verification
result = mcp.pubmed_search("CD3D CD3E T cell lineage marker", max_results=3)
pmid = result['results'][0]['pmid']

verification = mcp.verify_reference(
    pmid=pmid,
    markers=['CD3D', 'CD3E', 'TRAC'],
    cell_type='T cells'
)
# Status: VERIFIED / DOUBLE_VERIFIED
```

### Step 6.5: 3-Iteration Reasoning (NEW - GENERATES reasoning_results)

```python
print("\n🧠 Step 6.5: Performing 3-Iteration Reasoning for each cluster...")

# This step generates reasoning_results which is used in:
# - Step 7.5: Create annotation_evidence
# - Step 9: Visualization

reasoning_results = {}

for cluster in adata.obs['tier1_cluster'].unique():
    cluster_str = str(cluster)
    print(f"\n--- Cluster {cluster_str} Reasoning ---")

    # Get cluster-specific data
    cluster_mask = adata.obs['tier1_cluster'] == cluster
    cluster_de = de_df[de_df['group'] == cluster]

    # Filter valid markers
    valid_markers = cluster_de[
        (cluster_de['pct_nz_group'] >= 0.25) &
        (cluster_de['logfoldchanges'] >= 1.0) &
        (cluster_de['pvals_adj'] < 0.05)
    ].nlargest(50, 'logfoldchanges')

    # Extract top markers for reasoning
    top_genes = valid_markers.nlargest(10, 'logfoldchanges')['names'].tolist()

    # Get outliers for this cluster (if any)
    cluster_outliers = [o for o in adata.uns.get('marker_outliers', [])
                       if o['cluster'] == cluster_str]

    # ITERATION 1: Evidence Collection + PRE-Reasoning Literature
    # (This would involve MCP pubmed_search for each marker combination)

    # ITERATION 2: Reasoning + Verification
    # (Compare candidates, verify references)

    # ITERATION 3: POST-Reasoning Literature + Final Decision
    # (Double-check references, make final assignment)

    # For now, create a structured output based on the reasoning template
    # In practice, this would be generated through LLM reasoning

    # Example reasoning output structure:
    cluster_reasoning = {
        'cluster_id': cluster_str,
        'final_assignment': None,  # Will be filled by reasoning
        'full_label': None,

        # Markers (from valid_markers)
        'marker_details': [
            {
                'gene': row['names'],
                'pct_in': row['pct_nz_group'],
                'log2fc': row['logfoldchanges'],
                'padj': row['pvals_adj']
            }
            for _, row in valid_markers.head(5).iterrows()
        ],

        # Candidates (would be generated by reasoning)
        'candidates': [
            {
                'name': 'Candidate_1',
                'markers': top_genes[:3],
                'support': ['Evidence 1', 'Evidence 2'],
                'against': [],
                'selected': True
            }
        ],

        # References (from MCP verification)
        'references': [
            # Format: {pmid, genes, title, journal, year, authors_short, finding,
            #          first_verify, second_verify, status}
        ],

        # Supporting reasons
        'supporting_reasons': [
            f"Top markers: {', '.join(top_genes[:5])}",
            f"Cluster size: {cluster_mask.sum()} cells"
        ],

        # Conflicts (if any)
        'conflicts': [],

        # Confidence scoring (0-12)
        'confidence_score': 0,  # Will be calculated
        'confidence_level': 'Low',  # Low/Medium/High

        # Outliers (if any)
        'marker_outliers': cluster_outliers,

        # Novel flag
        'is_novel': False,
        'novel_evidence': None
    }

    # Store reasoning result
    reasoning_results[cluster_str] = cluster_reasoning

    print(f"  ✓ Reasoning structure created for cluster {cluster_str}")

print(f"\n✅ Reasoning completed for {len(reasoning_results)} clusters")
print("⚠️  NOTE: In production, this step would involve:")
print("   1. LLM-based 3-iteration reasoning (Bioinformatician → Computational Biologist → Scientific Critic → Domain Expert → PI)")
print("   2. MCP pubmed_search for literature verification")
print("   3. Candidate comparison and final decision")
print("   4. Confidence scoring based on evidence quality")
print("\n   For now, reasoning_results contains structured templates.")
print("   → Complete the reasoning manually or integrate LLM reasoning pipeline.")
```

**IMPORTANT**: This step creates the `reasoning_results` dictionary that is used by:
- Step 7.5: `annotation_evidence` generation
- Step 9: Visualization

The actual reasoning process should follow the template in `reasoning/integrated_format.md`:
1. **Iteration 1**: Evidence collection + PRE-reasoning literature search
2. **Iteration 2**: Reasoning + verification (compare candidates)
3. **Iteration 3**: POST-reasoning literature + final decision

### Step 7: Assign Annotations

```python
tier1_mapping = {
    '0': 'T cells',
    '1': 'T cells',
    '2': 'B lineage',
    '3': 'Myeloid',
    '4': 'NK cells',
    ...
}

adata.obs['tier1_annotation'] = adata.obs['tier1_cluster'].map(tier1_mapping)
```

### Step 7.5: Create annotation_evidence (NEW - REQUIRED FOR VISUALIZATION)

```python
# Convert reasoning results to annotation_evidence format
# This is REQUIRED for visualization in Step 9

from reasoning.agent_format import reasoning_to_evidence

annotation_evidence = []

for cluster in adata.obs['tier1_cluster'].unique():
    # Get annotation for this cluster
    cell_type = tier1_mapping.get(str(cluster))
    if not cell_type:
        continue

    # Get cluster-specific data from reasoning process
    # NOTE: This assumes reasoning was performed for each cluster
    # reasoning_results should be created during Step 6 (MCP verification + reasoning)

    cluster_reasoning = reasoning_results.get(str(cluster), {})

    # Create evidence entry
    evidence_entry = reasoning_to_evidence(
        cluster_id=f"tier1_cluster_{cluster}",
        tier="tier1",
        subset_id="full_dataset",
        reasoning_output=cluster_reasoning
    )

    annotation_evidence.append(evidence_entry)

print(f"✅ Created annotation_evidence: {len(annotation_evidence)} entries")

# Save to JSON for persistence
import json
import os
os.makedirs('annotation_output/references', exist_ok=True)

with open('annotation_output/references/tier1_annotation_evidence.json', 'w') as f:
    json.dump(annotation_evidence, f, indent=2, ensure_ascii=False)

print("✅ Saved: annotation_output/references/tier1_annotation_evidence.json")
```

### Step 8: Cross-Contamination Verification (DATA-DRIVEN)

```python
# Verify contamination after annotation (NO hardcoded markers)
def verify_contamination_post_annotation(adata, annotation_col, de_df):
    """
    Verify cross-contamination using data-driven approach.
    For each annotation, check if cells express TOP markers from OTHER annotations.
    NO hardcoded marker lists.
    """
    issues = []

    # Build signature for each annotation from top DE markers
    annotation_signatures = {}
    for ann in adata.obs[annotation_col].unique():
        # Get cells with this annotation
        ann_mask = adata.obs[annotation_col] == ann

        # Get clusters that belong to this annotation
        clusters_in_ann = adata.obs[ann_mask]['tier1_cluster'].unique()

        # Get top markers for these clusters
        top_markers = []
        for cluster in clusters_in_ann:
            cluster_de = de_df[de_df['group'] == cluster]
            cluster_top = cluster_de[
                (cluster_de['pct_nz_group'] >= 0.50) &  # High expression
                (cluster_de['logfoldchanges'] >= 2.0) &  # Strong enrichment
                (cluster_de['pvals_adj'] < 0.001)  # Highly significant
            ].nlargest(5, 'logfoldchanges')['names'].tolist()
            top_markers.extend(cluster_top)

        # Keep unique top markers as signature
        annotation_signatures[ann] = list(set(top_markers))[:10]

    # Check each annotation for contamination from others
    for ann in adata.obs[annotation_col].unique():
        mask = adata.obs[annotation_col] == ann
        subset = adata[mask]

        # Check for markers from OTHER annotations
        for other_ann, other_markers in annotation_signatures.items():
            if other_ann == ann:
                continue

            contamination_count = 0
            contaminating_genes = []

            for marker in other_markers:
                if marker in subset.var_names:
                    if hasattr(subset.X, 'toarray'):
                        expr = subset[:, marker].X.toarray().flatten()
                    else:
                        expr = subset[:, marker].X.flatten()

                    pct = (expr > 0).mean()
                    if pct > 0.05:  # > 5% expressing
                        contamination_count += 1
                        contaminating_genes.append(f"{marker} ({pct:.1%})")

            # Flag if multiple markers from another lineage are present
            if contamination_count >= 2:
                issues.append({
                    'annotation': ann,
                    'contaminated_by': other_ann,
                    'n_markers': contamination_count,
                    'markers': contaminating_genes[:3],
                    'severity': 'HIGH' if contamination_count >= 3 else 'MEDIUM'
                })

    return issues

# Run verification
contamination = verify_contamination_post_annotation(adata, 'tier1_annotation', de_df)

if contamination:
    print("⚠️ Post-annotation contamination check:")
    for issue in contamination:
        print(f"  {issue['annotation']} contaminated by {issue['contaminated_by']}")
        print(f"    → {issue['n_markers']} markers: {', '.join(issue['markers'])}")
        print(f"    → Severity: {issue['severity']}")
        print(f"    → Action: Review these cells, may be doublets or mis-annotated")
else:
    print("✅ No significant cross-contamination detected")
```

### Step 9: Save Results

```python
# Save full dataset with Tier 1 annotations
adata.write('annotation_output/tier1_annotated.h5ad')

# Save subsets per major type
import os
os.makedirs('annotation_output/subsets/tier1', exist_ok=True)

for major_type in adata.obs['tier1_annotation'].unique():
    subset = adata[adata.obs['tier1_annotation'] == major_type].copy()
    safe_name = major_type.replace(' ', '_')
    subset.write(f'annotation_output/subsets/tier1/{safe_name}.h5ad')
    print(f"Saved: {safe_name} ({subset.n_obs} cells)")

# ✅ Generate Visualizations (NEW)
print("\n🎨 Generating Tier 1 Visualizations...")
os.makedirs('annotation_output/figures', exist_ok=True)

# Prepare marker dict from annotation_evidence (VALIDATED markers only)
# annotation_evidence는 Step 7에서 reasoning 결과로 생성됨
marker_dict = {}
annotation_order = sorted(adata.obs['tier1_annotation'].unique())

for entry in annotation_evidence:
    if entry['tier'] == 'tier1':
        cell_type = entry['annotation']
        # Extract validated marker gene names from reasoning results
        validated_markers = [m['gene'] for m in entry['markers']]
        marker_dict[cell_type] = validated_markers

# Import visualization functions
import sys
sys.path.append('annotation_output')  # Add to path if needed
from tools.visualization import save_all_visualizations

# Generate all visualizations
save_all_visualizations(
    adata=adata,
    annotation_col='tier1_annotation',
    marker_dict=marker_dict,
    reasoning_dict={entry['annotation']: entry for entry in annotation_evidence if entry['tier'] == 'tier1'},
    reference_dict={entry['annotation']: entry['references'] for entry in annotation_evidence if entry['tier'] == 'tier1'},
    annotation_order=annotation_order,
    title_prefix='Tier 1 - Major Cell Types',
    output_dir='annotation_output/figures/'
)

print("✅ Tier 1 Visualizations saved:")
print("   - tier1_umap.png/.svg")
print("   - tier1_dotplot.png/.svg (with marker brackets)")
print("   - tier1_umap_dotplot.png/.svg (2-panel combined)")
print("   - tier1_full_report.png/.svg (4-panel)")
print("   - tier1_reasoning.json/.txt")
```

---

## Checklist (v3 - Outlier-Aware)

```
ANALYSIS:
- [ ] Leiden clustering 완료 (resolution 0.5-1.0)
- [ ] DE 계산 완료 (rank_genes_groups)
- [ ] 🆕 Statistical outliers detected (marker_outliers in uns)
- [ ] 🆕 Cross-lineage contamination detected (contamination_issues in uns)
- [ ] Valid markers 필터링 (pct >= 25%, LFC >= 1)

ANNOTATION:
- [ ] 각 클러스터에 2-4개 lineage marker 확인 (from DE, NOT hardcoded)
- [ ] 🆕 IF outliers: Marker combination literature search performed
- [ ] 🆕 IF contamination: Doublet check performed
- [ ] MCP로 PMID 검증 완료
- [ ] Cross-contamination < 5% (data-driven verification)
- [ ] 모든 클러스터 annotation 완료

OUTPUT:
- [ ] tier1_annotated.h5ad 저장 (with outlier/contamination info)
- [ ] 각 major type subset 저장 (tier1/{type}.h5ad)
- [ ] UMAP + Dotplot 시각화 저장
```

---

## 성공 기준

| 기준 | 설명 | 통과 |
|------|------|------|
| Marker count | 각 annotation에 2-4개 lineage marker | >= 2 |
| PMID | 모든 annotation에 PMID 참조 | 100% |
| Contamination | Cross-type contamination | < 5% |
| Coverage | 모든 클러스터 annotation | 100% |

---

## 주의사항 (v3 Updated)

```
❌ 하지 말 것:
- 발달 상태 포함 (예: "Naive T cells" → Tier 2에서)
- 기능 상태 포함 (예: "Cytotoxic T cells" → Tier 3에서)
- 단일 마커로 annotation (예: "CD3D+ cells")
- 🆕 Hardcoded lineage marker lists 사용
  ❌ t_markers = ['CD3D', 'CD3E', 'TRAC']  # 하드코딩
  ❌ b_markers = ['CD79A', 'CD79B', 'MS4A1']  # 하드코딩
- 🆕 Statistical outliers 무시
- 🆕 Contamination warnings 무시

✅ 해야 할 것:
- 계통만 결정 (예: "T cells", "B lineage")
- 2-4개 마커 조합 사용 (from current dataset DE)
- 🆕 Outlier 기반 우선순위 설정
  - Outlier markers (z > 2.5) → 문헌 검색 우선
  - Marker combination → "{marker1} {marker2} {marker3}" 검색
- 🆕 Contamination 체크
  - Mixed lineage signatures → Doublet 가능성 검토
  - High contamination score → 재검토 또는 제외
- PMID 검증
```

---

## 🆕 Doublet Detection Strategy (v3)

```python
def detect_potential_doublets(adata, annotation_col, de_df):
    """
    Detect potential doublets based on mixed lineage signatures.
    NO hardcoded markers - uses top DE markers from each annotation.
    """
    doublet_candidates = []

    # Get signature for each major type
    signatures = {}
    for ann in adata.obs[annotation_col].unique():
        mask = adata.obs[annotation_col] == ann
        clusters = adata.obs[mask]['tier1_cluster'].unique()

        top_genes = []
        for cluster in clusters:
            cluster_top = de_df[
                (de_df['group'] == cluster) &
                (de_df['pct_nz_group'] >= 0.5) &
                (de_df['logfoldchanges'] >= 2.0)
            ].nlargest(5, 'logfoldchanges')['names'].tolist()
            top_genes.extend(cluster_top)

        signatures[ann] = list(set(top_genes))[:10]

    # Check each cluster for mixed signatures
    for cluster in adata.obs['tier1_cluster'].unique():
        mask = adata.obs['tier1_cluster'] == cluster
        cluster_cells = adata[mask]

        # Count how many signatures are present
        signature_scores = {}
        for ann, genes in signatures.items():
            score = 0
            for gene in genes:
                if gene in adata.var_names:
                    if hasattr(cluster_cells.X, 'toarray'):
                        expr = cluster_cells[:, gene].X.toarray().flatten()
                    else:
                        expr = cluster_cells[:, gene].X.flatten()
                    pct = (expr > 0).mean()
                    if pct > 0.3:  # 30% expressing
                        score += pct

            signature_scores[ann] = score

        # If 2+ signatures are strong, likely doublet
        strong_signatures = [ann for ann, score in signature_scores.items() if score > 1.0]
        if len(strong_signatures) >= 2:
            doublet_candidates.append({
                'cluster': str(cluster),
                'signatures': strong_signatures,
                'scores': {k: v for k, v in signature_scores.items() if k in strong_signatures},
                'recommendation': 'Review for doublets - consider filtering'
            })

    return doublet_candidates

# Run doublet detection
doublets = detect_potential_doublets(adata, 'tier1_annotation', de_df)

if doublets:
    print("⚠️ Potential doublet clusters detected:")
    for d in doublets:
        print(f"  Cluster {d['cluster']}: Mixed signatures {d['signatures']}")
        print(f"    → {d['recommendation']}")
else:
    print("✅ No obvious doublet clusters detected")
```
