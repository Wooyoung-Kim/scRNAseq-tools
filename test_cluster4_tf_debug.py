"""Debug: trace TF channel to find where "Human and Mouse Memory B Cells" comes from."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')

from dynamic_knowledge import (
    extract_cell_subtypes_from_text,
    extract_cell_types_from_text,
    normalize_cell_type_name,
    aggregate_candidates,
    _get_cached_aliases,
    _build_tf_subtype_queries,
    _score_cell_types_by_marker_context,
    _select_top_tfs,
    pubmed_search,
)

parent_cell_type = "B cells"
tf_activities = {
    'TBX21': 3.5, 'SOX5': 2.8, 'POU2F2': 1.8,
    'BATF': 1.5, 'BHLHE41': 1.2, 'ZEB2': 1.0,
}

top_tfs = _select_top_tfs(tf_activities, max_tfs=3)
print(f"Top TFs: {top_tfs}")

tf_candidates_raw = []

for tf_name, tf_score in top_tfs:
    print(f"\n{'='*60}")
    print(f"TF: {tf_name} (score={tf_score})")

    aliases = _get_cached_aliases(tf_name)
    print(f"  Aliases: {aliases}")

    queries = _build_tf_subtype_queries(tf_name, parent_cell_type, tier=2, species='human')

    all_results = []
    seen_pmids = set()
    for qi, query in enumerate(queries):
        results = pubmed_search(query, max_results=3)
        for r in results:
            if r['pmid'] not in seen_pmids:
                all_results.append(r)
                seen_pmids.add(r['pmid'])
        if len(all_results) >= 8:
            break

    print(f"  Total unique papers: {len(all_results)}")

    for result in all_results:
        title = result.get('title') or ''
        abstract = result.get('abstract') or ''
        text = title + " " + abstract
        pmid = result['pmid']

        cell_types = extract_cell_subtypes_from_text(text, parent_cell_type)

        # Check each extracted type for species contamination
        for ct in cell_types:
            ct_lower = ct.lower()
            if 'human' in ct_lower or 'mouse' in ct_lower:
                print(f"\n  !!! SPECIES HIT: '{ct}' from PMID {pmid}")
                print(f"      Title: {title[:120]}")
                # Show first 500 chars of abstract to find the source text
                print(f"      Abstract (first 300): {abstract[:300]}")
                norm = normalize_cell_type_name(ct, preserve_subtypes=True)
                print(f"      Normalized: '{norm}'")

        scored_types = _score_cell_types_by_marker_context(cell_types, [tf_name], text)
        for cell_type, context_score in scored_types:
            tf_candidates_raw.append({
                'cell_type': cell_type,
                'markers': [tf_name],
                'pmid': pmid,
                'title': title,
                'journal': result.get('journal', 'N/A'),
                'year': result.get('year', 'N/A'),
                'context_score': context_score,
                'tf_name': tf_name,
                'tf_score': tf_score
            })

# Aggregate
tf_aggregated = aggregate_candidates(tf_candidates_raw, preserve_subtypes=True)
print(f"\n--- Aggregated TF candidates ---")
for c in tf_aggregated[:10]:
    flag = " *** SPECIES" if ('human' in c.name.lower() or 'mouse' in c.name.lower()) else ""
    print(f"  {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}{flag}")
