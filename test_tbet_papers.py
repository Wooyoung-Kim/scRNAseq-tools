"""Check what subtypes we extract from all TBX21/T-bet papers."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')

from dynamic_knowledge import (
    extract_cell_subtypes_from_text,
    _get_cached_aliases,
    _build_tf_subtype_queries,
    pubmed_search,
)

queries = _build_tf_subtype_queries('TBX21', 'B cells', tier=2, species='human')

all_results = []
seen_pmids = set()
for q in queries:
    results = pubmed_search(q, max_results=3)
    for r in results:
        if r['pmid'] not in seen_pmids:
            all_results.append(r)
            seen_pmids.add(r['pmid'])
    if len(all_results) >= 8:
        break

print(f"Total papers from TBX21 queries: {len(all_results)}\n")

for r in all_results:
    title = r.get('title') or ''
    abstract = r.get('abstract') or ''
    text = title + " " + abstract
    pmid = r['pmid']

    subtypes = extract_cell_subtypes_from_text(text, 'B cells')

    print(f"PMID {pmid}: {title[:100]}")
    print(f"  Extracted subtypes: {subtypes}")

    # Check if "age" or "atypical" appears anywhere in text
    text_lower = text.lower()
    for keyword in ['age-associated', 'age associated', 'atypical', 'abc ', 'abcs', 'double-negative', 'dn2', 'cd11c']:
        if keyword in text_lower:
            # Find the context
            idx = text_lower.find(keyword)
            context = text[max(0,idx-50):idx+80]
            print(f"  ** Found '{keyword}' in text: ...{context}...")
    print()
