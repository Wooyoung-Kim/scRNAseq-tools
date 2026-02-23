"""Debug: trace marker-only vs TF-only vs merged for Cluster 0."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
import csv

from dynamic_knowledge import (
    search_subtype_from_markers,
    search_subtype_with_tf_activity,
    filter_informative_markers,
)

# Load Cluster 0 markers
markers = []
with open('/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/B_cells_cluster_0_markers.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        markers.append(row['names'])

tf_acts = {'KLF2': 2.5, 'FOXO1': 2.0, 'BACH2': 1.5, 'TCF3': 1.2}

# Step 1: Marker-only results
print("=" * 60)
print("MARKER-ONLY (search_subtype_from_markers)")
print("=" * 60)
marker_cands = search_subtype_from_markers(
    markers, parent_cell_type="B cells", tier=2, species="human"
)
for i, c in enumerate(marker_cands[:10]):
    print(f"  {i+1}. {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}, markers={c.supporting_markers[:3]}")

# Step 2: Full pipeline
print(f"\n{'=' * 60}")
print("FULL PIPELINE (with TF)")
print("=" * 60)
full_cands = search_subtype_with_tf_activity(
    markers, tf_acts,
    parent_cell_type="B cells", tier=2, species="human"
)
for i, c in enumerate(full_cands[:10]):
    print(f"  {i+1}. {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}, markers={c.supporting_markers[:5]}")
