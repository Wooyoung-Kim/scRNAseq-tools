"""Debug: check PMID overlap between ABC and CD11c+ B Cells."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
import csv

from dynamic_knowledge import search_subtype_with_tf_activity, filter_informative_markers

markers_4 = []
with open('/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/B_cells_cluster_4_markers.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        markers_4.append(row['names'])

tf_activities_4 = {
    'TBX21': 3.5, 'SOX5': 2.8, 'POU2F2': 1.8,
    'BATF': 1.5, 'BHLHE41': 1.2, 'ZEB2': 1.0,
}

candidates = search_subtype_with_tf_activity(
    markers_4, tf_activities_4,
    parent_cell_type="B cells", tier=2, species="human"
)

# Show PMIDs for ABC-related candidates
for c in candidates:
    c_lower = c.name.lower()
    if any(k in c_lower for k in ['age', 'cd11c', 'effector', 'atypical', 'dn ']):
        print(f"\n{c.name} (conf={c.confidence:.3f})")
        print(f"  PMIDs: {c.pmids}")
        for ev in c.evidence:
            print(f"  Evidence: PMID {ev.get('pmid','?')} - {ev.get('title','?')[:80]}")
