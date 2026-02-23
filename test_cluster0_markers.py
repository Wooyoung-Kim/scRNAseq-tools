"""Quick check: what do individual marker searches find for Cluster 0?"""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')

from dynamic_knowledge import (
    pubmed_search, extract_cell_subtypes_from_text,
    _build_subtype_queries, filter_informative_markers,
)

markers = ['BTLA', 'FCER2', 'BANK1', 'CD200', 'LTB']

# Check what queries are built for individual markers
print("=== Individual marker query check ===")
for m in markers[:3]:
    queries = _build_subtype_queries([m], "B cells", tier=2, species="human")
    print(f"\n--- {m} ---")
    for q in queries[:4]:
        print(f"  Query: {q}")
        results = pubmed_search(q, max_results=3)
        print(f"  Results: {len(results)}")
        for r in results[:2]:
            title = (r.get('title') or '')[:100]
            text = (r.get('title') or '') + " " + (r.get('abstract') or '')
            subtypes = extract_cell_subtypes_from_text(text, "B cells")
            naive_related = [s for s in subtypes if any(k in s.lower()
                for k in ['naive', 'follicular', 'mature', 'resting'])]
            print(f"    PMID {r['pmid']}: {title}")
            print(f"      Naive-related: {naive_related}")
            print(f"      All: {subtypes[:5]}")
        if results:
            break  # Same break logic as actual code

# Also check specific targeted queries
print("\n\n=== Targeted naive B queries ===")
targeted_queries = [
    'FCER2[Title/Abstract] AND naive[Title/Abstract] AND B cell',
    'CD23[Title/Abstract] AND follicular[Title/Abstract] AND B cell',
    'BTLA[Title/Abstract] AND naive[Title/Abstract] AND B cell',
    'FCER2 follicular B cell',
    'CD23 naive B cell human',
]
for q in targeted_queries:
    results = pubmed_search(q, max_results=3)
    print(f"\nQuery: {q}")
    print(f"  Results: {len(results)}")
    for r in results[:2]:
        title = (r.get('title') or '')[:100]
        text = (r.get('title') or '') + " " + (r.get('abstract') or '')
        subtypes = extract_cell_subtypes_from_text(text, "B cells")
        print(f"    PMID {r['pmid']}: {title}")
        naive_sub = [s for s in subtypes if any(k in s.lower()
            for k in ['naive', 'follicular', 'mature', 'resting'])]
        print(f"      Naive-related: {naive_sub}")
