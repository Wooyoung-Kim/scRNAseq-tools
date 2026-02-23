"""Check what "T-bet + CD11c + B cell" queries find."""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
from dynamic_knowledge import pubmed_search, extract_cell_subtypes_from_text

# Direct T-bet + CD11c searches
queries = [
    'T-bet[Title/Abstract] AND CD11c[Title/Abstract] AND B[Title/Abstract]',
    '"T-bet" "CD11c" "B cell" subset',
    '"age-associated B cells" T-bet',
    '"T-bet" "atypical" B cell',
    'T-bet[Title/Abstract] AND (age-associated OR atypical OR ABC) AND B[Title/Abstract]',
]

for q in queries:
    print(f"\nQuery: {q}")
    results = pubmed_search(q, max_results=5)
    print(f"  Results: {len(results)}")
    for r in results[:3]:
        title = r.get('title', '')[:120]
        pmid = r['pmid']
        text = (r.get('title') or '') + ' ' + (r.get('abstract') or '')
        subtypes = extract_cell_subtypes_from_text(text, 'B cells')
        print(f"  PMID {pmid}: {title}")
        abc_related = [s for s in subtypes if any(k in s.lower() for k in ['age', 'atypical', 'cd11c', 'abc', 'effector', 'double'])]
        print(f"    ABC-related: {abc_related}")
        print(f"    All subtypes: {subtypes}")
