"""Quick test: Cluster 0 only, to verify validation fix."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
import csv

from dynamic_knowledge import search_subtype_with_tf_activity, filter_informative_markers

# Load Cluster 0 markers
markers = []
with open('/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/B_cells_cluster_0_markers.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        markers.append(row['names'])

tf_acts = {'KLF2': 2.5, 'FOXO1': 2.0, 'BACH2': 1.5, 'TCF3': 1.2}

print("Cluster 0 (expected: Naive/Follicular B)")
print(f"  Top markers: {filter_informative_markers(markers, max_markers=10)[:5]}")
print(f"  TFs: {list(tf_acts.items())}")

candidates = search_subtype_with_tf_activity(
    markers, tf_acts,
    parent_cell_type="B cells", tier=2, species="human"
)

print(f"\n  Top 10 candidates:")
for i, c in enumerate(candidates[:10]):
    print(f"    {i+1}. {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}")
