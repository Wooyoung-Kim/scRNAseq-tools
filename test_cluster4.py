"""Test Cluster 4 (expected: ABC / Age-associated B Cells) with TF activity integration."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
import csv

from dynamic_knowledge import (
    search_subtype_with_tf_activity,
    search_subtype_from_markers,
    filter_informative_markers,
    extract_cell_subtypes_from_text,
    _get_cached_aliases,
    _build_tf_subtype_queries,
    pubmed_search,
)

# Load Cluster 4 markers
markers_4 = []
with open('/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/B_cells_cluster_4_markers.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        markers_4.append(row['names'])

print(f"Raw markers (top 20): {markers_4[:20]}")

filtered = filter_informative_markers(markers_4, max_markers=10)
print(f"Filtered markers: {filtered}")

# Simulated TF activities for Cluster 4 (from decoupler-like output)
tf_activities_4 = {
    'TBX21': 3.5,    # T-bet: strongly active → ABCs
    'SOX5': 2.8,     # High DE marker
    'POU2F2': 1.8,   # Oct-2: B cell TF
    'BATF': 1.5,     # Active in ABCs
    'BHLHE41': 1.2,  # DEC2: active in ABCs
    'ZEB2': 1.0,     # EMT TF, also in ABCs
}

# First, check what aliases we get for key TFs
print("\n--- TF Aliases ---")
for tf in ['TBX21', 'BATF', 'ZEB2']:
    aliases = _get_cached_aliases(tf)
    print(f"  {tf} → {aliases}")

# Check what queries are built for TBX21
print("\n--- TBX21 Queries ---")
queries = _build_tf_subtype_queries('TBX21', 'B cells', tier=2, species='human')
for i, q in enumerate(queries):
    print(f"  Q{i}: {q}")

# Run a quick TBX21 T-bet query manually to see what papers we find
print("\n--- Manual T-bet query test ---")
results = pubmed_search('T-bet[Title/Abstract] AND B[Title/Abstract] AND (differentiation OR development OR subset OR subtype)', max_results=5)
for r in results[:3]:
    title = r.get('title', '')
    pmid = r.get('pmid', '')
    print(f"  PMID {pmid}: {title[:100]}")
    # Extract subtypes from this paper
    text = (title or '') + ' ' + (r.get('abstract') or '')
    subtypes = extract_cell_subtypes_from_text(text, 'B cells')
    print(f"    Subtypes found: {subtypes[:10]}")

# Run full TF-integrated search
print("\n" + "=" * 60)
print("=== Cluster 4: Full TF-integrated search ===")
print("=" * 60)
candidates = search_subtype_with_tf_activity(
    markers_4, tf_activities_4,
    parent_cell_type="B cells", tier=2, species="human"
)
for i, c in enumerate(candidates[:10]):
    print(f"  {i+1}. {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}, markers={c.supporting_markers[:5]}")
    # Show evidence titles
    for ev in c.evidence[:2]:
        print(f"       → {ev.get('title', 'N/A')[:80]}")
