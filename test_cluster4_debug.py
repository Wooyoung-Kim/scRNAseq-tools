"""Debug: trace where "Human and Mouse Memory B Cells" name comes from."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
import csv

from dynamic_knowledge import (
    search_subtype_from_markers,
    filter_informative_markers,
    extract_cell_subtypes_from_text,
    extract_cell_types_from_text,
    normalize_cell_type_name,
    aggregate_candidates,
    _build_subtype_queries,
    _score_cell_types_by_marker_context,
    pubmed_search,
    generate_marker_combinations,
)

# Load Cluster 4 markers
markers_4 = []
with open('/home/kwy7605/data_61/Vaccine_V3/annotation_output/markers/tier2/B_cells_cluster_4_markers.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        markers_4.append(row['names'])

filtered = filter_informative_markers(markers_4, max_markers=10)
print(f"Filtered markers: {filtered}")

# Trace marker-based search manually
parent_cell_type = "B cells"
marker_combinations = generate_marker_combinations(filtered[:8], max_size=3)

# Also individual markers
for m in filtered[:3]:
    if [m] not in marker_combinations:
        marker_combinations.append([m])

all_raw_candidates = []

for combo_idx, combo in enumerate(marker_combinations[:5]):
    queries = _build_subtype_queries(combo, parent_cell_type, tier=2, species="human")

    results = []
    for query in queries:
        results = pubmed_search(query, max_results=5)
        if results:
            break

    for result in results:
        title = result.get('title') or ''
        abstract = result.get('abstract') or ''
        text = title + " " + abstract
        pmid = result['pmid']

        # Extract subtypes
        cell_types = extract_cell_subtypes_from_text(text, parent_cell_type)

        # Check if "Human and Mouse" appears in any extracted type
        for ct in cell_types:
            if 'human' in ct.lower() or 'mouse' in ct.lower():
                print(f"\n!!! FOUND species-contaminated cell type: '{ct}'")
                print(f"    From PMID {pmid}: {title[:120]}")
                print(f"    Combo: {combo}")
                # Check which pattern extracted it
                # Re-run extractions separately
                print(f"    --- Subtype extraction yields ---")
                for st in cell_types:
                    print(f"      - {st}")
                print(f"    --- General extraction from same text ---")
                general = extract_cell_types_from_text(text)
                for g in general:
                    print(f"      - {g}")
                print(f"    --- Normalized: '{normalize_cell_type_name(ct, preserve_subtypes=True)}'")

        scored_types = _score_cell_types_by_marker_context(cell_types, combo, text)
        for cell_type, context_score in scored_types:
            norm = normalize_cell_type_name(cell_type, preserve_subtypes=True)
            if 'human' in norm.lower() or 'mouse' in norm.lower():
                print(f"  AFTER NORM still has species: '{norm}' (from '{cell_type}')")
            all_raw_candidates.append({
                'cell_type': cell_type,
                'markers': combo,
                'pmid': pmid,
                'title': title,
                'journal': result.get('journal', 'N/A'),
                'year': result.get('year', 'N/A'),
                'context_score': context_score
            })

# Aggregate and show
aggregated = aggregate_candidates(all_raw_candidates, preserve_subtypes=True)
print(f"\n--- Aggregated marker candidates ---")
for c in aggregated[:10]:
    norm_lower = c.name.lower()
    flag = " *** SPECIES" if ('human' in norm_lower or 'mouse' in norm_lower) else ""
    print(f"  {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}{flag}")
