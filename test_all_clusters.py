"""Comprehensive test of all 6 B cell clusters using per-marker voting system."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
import csv

from dynamic_knowledge import (
    search_subtype_by_voting,
    filter_informative_markers,
    _zscore_tf_activities,
    _select_top_tfs,
)

# Load all cluster markers
clusters = {}
for i in range(6):
    markers = []
    path = f'/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/B_cells_cluster_{i}_markers.csv'
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            markers.append(row['names'])
    clusters[i] = markers

# Load REAL TF activities per cluster (from decoupler output, 745 TFs each)
tf_activities_per_cluster = {}
for i in range(6):
    tf_path = f'/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/tf_activity/cluster_{i}_tf_activity.csv'
    tf_acts = {}
    with open(tf_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf_acts[row['TF']] = float(row['activity_score'])
    tf_activities_per_cluster[i] = tf_acts

# Load REAL pseudotime data
pseudotime_per_cluster = {}
pt_path = '/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/trajectory/cluster_pseudotime_summary.csv'
with open(pt_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        cid = int(row['tier2_cluster'])
        pseudotime_per_cluster[cid] = {
            'mean': float(row['mean']),
            'std': float(row['std']),
        }

expected = {
    0: "Naive/Follicular B",
    1: "Plasma cell precursors / Plasmablasts",
    2: "Transitional/Immature B",
    3: "Pro-B / Pre-B",
    4: "ABC (Age-associated B cells)",
    5: "Pro-B / Pre-B (early)",
}

print("=" * 70)
print("PER-MARKER VOTING B CELL CLUSTER TEST")
print("=" * 70)

for cluster_id in sorted(clusters.keys()):
    markers = clusters[cluster_id]
    tf_acts = tf_activities_per_cluster.get(cluster_id, {})
    pt_stats = pseudotime_per_cluster.get(cluster_id)
    exp = expected.get(cluster_id, "Unknown")

    filtered = filter_informative_markers(markers, max_markers=10)
    z_scores = _zscore_tf_activities(tf_acts, tf_activities_per_cluster) if tf_acts else {}
    z_sorted = sorted(z_scores.items(), key=lambda x: x[1], reverse=True)

    # Select top 3 TFs for display (matching the voting system's selection logic)
    import math
    active_tfs = {k: v for k, v in tf_acts.items() if v > 0.5}
    tf_sel_w = {}
    for k, v in active_tfs.items():
        z = max(z_scores.get(k, 0), 0)
        raw_weight = math.log1p(v)
        tf_sel_w[k] = z * raw_weight + 0.01 * raw_weight
    top_tfs = _select_top_tfs(tf_sel_w, max_tfs=3)

    print(f"\n{'='*70}")
    print(f"Cluster {cluster_id} (expected: {exp})")
    print(f"  Top markers: {filtered[:5]}")
    print(f"  Top TFs (z): {[(n, f'{z_scores.get(n,0):.2f}') for n, _ in top_tfs]}")
    if pt_stats:
        print(f"  Pseudotime: mean={pt_stats['mean']:.3f}, std={pt_stats['std']:.3f}")

    candidates = search_subtype_by_voting(
        markers, tf_acts,
        parent_cell_type="B cells", tier=2, species="human",
        n_voting_markers=8, n_voting_tfs=3,
        pseudotime_stats=pt_stats,
        all_cluster_tf_activities=tf_activities_per_cluster,
        all_cluster_pseudotime=pseudotime_per_cluster,
    )

    print(f"  Top 5 candidates:")
    for i, c in enumerate(candidates[:5]):
        # Show voting details if available
        vote_info = ""
        if c.evidence and isinstance(c.evidence[0], dict):
            ev = c.evidence[0]
            vm = ev.get('voting_markers', [])
            vt = ev.get('voting_tfs', [])
            cross = ev.get('cross_validated', False)
            vote_info = f" [markers:{len(vm)}, TFs:{len(vt)}, cross:{cross}]"
        print(f"    {i+1}. {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}{vote_info}")

print(f"\n{'='*70}")
print("TEST COMPLETE")
print("=" * 70)
