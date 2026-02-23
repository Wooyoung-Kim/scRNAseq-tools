# Dynamic Knowledge Base System

**완전한 데이터 드라이븐 세포 타입 발견 시스템**

NO hardcoded cell types. NO hardcoded markers.
Cell type identification is purely literature-driven from actual DE results.

---

## Philosophy

```
╔══════════════════════════════════════════════════════════════════════╗
║  Traditional Approach (REJECTED):                                    ║
║    1. Check if markers match hardcoded list                          ║
║    2. Assign predefined cell type                                    ║
║    → Misses novel populations                                        ║
║    → Biased towards common cell types                                ║
║                                                                      ║
║  Dynamic Knowledge Approach (THIS SYSTEM):                           ║
║    1. Take top DE markers from data                                  ║
║    2. Search literature for marker combinations                      ║
║    3. Extract cell types mentioned in papers                         ║
║    4. Aggregate evidence → Ranked candidates                         ║
║    5. If no match → Flag as novel                                    ║
║    → Discovers rare/novel cell types                                 ║
║    → Unbiased, data-driven                                           ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## Workflow

### Input: DE Markers from Data

```python
# Example: Cluster 5 from tier1_cluster
cluster_markers = [
    'CD3D',    # LFC=4.5, padj=1e-100
    'CD3E',    # LFC=4.2, padj=1e-95
    'TRAC',    # LFC=3.8, padj=1e-80
    'CD8A',    # LFC=3.2, padj=1e-70
    'GZMB'     # LFC=2.8, padj=1e-60
]
```

### Step 1: Generate Marker Combinations

```python
from tools.dynamic_knowledge import generate_marker_combinations

combinations = generate_marker_combinations(cluster_markers, max_size=5)

# Results:
# ['CD3D', 'CD3E']                          # Top 2
# ['CD3D', 'CD3E', 'TRAC']                  # Top 3
# ['CD3D', 'CD3E', 'TRAC', 'CD8A']          # Top 4
# ['CD3D', 'CD3E', 'TRAC', 'CD8A', 'GZMB']  # Top 5
# ['CD3D', 'TRAC', 'CD8A']                  # Skip pattern
# ...
```

### Step 2: Search Literature for Each Combination

```python
from tools.dynamic_knowledge import search_cell_type_from_markers

candidates = search_cell_type_from_markers(
    markers=cluster_markers,
    tissue="peripheral blood",  # Optional
    species="human",
    max_combinations=10
)

# For each combination, searches PubMed:
# Query 1: "CD3D CD3E human cell type marker"
# Query 2: "CD3D CD3E TRAC human cell type marker"
# Query 3: "CD3D CD3E TRAC CD8A peripheral blood human cell type marker"
# ...
```

### Step 3: Extract Cell Types from Results

```python
# From paper titles/abstracts:
# "CD3+ T cells express CD8A in peripheral blood..."
# "Cytotoxic T lymphocytes express GZMB..."
# "Memory T cells show high CD3D expression..."

cell_types_found = [
    "T cells",           # Found in 5 papers
    "T lymphocytes",     # Found in 3 papers (normalized to "T cells")
    "CD8+ T cells",      # Found in 4 papers (normalized to "T cells")
    "B cells"            # Found in 1 paper (low confidence)
]
```

### Step 4: Aggregate and Rank Candidates

```python
# Scoring criteria:
# - Frequency: How many papers mention this cell type
# - Marker overlap: How many markers support it
# - Recency: Newer papers weighted higher

candidates = [
    CellTypeCandidate(
        name="T cells",
        confidence=0.92,
        supporting_markers=['CD3D', 'CD3E', 'TRAC', 'CD8A'],
        pmids=['12345678', '23456789', '34567890', '45678901', '56789012'],
        evidence=[...],
        is_novel=False
    ),
    CellTypeCandidate(
        name="B cells",
        confidence=0.12,
        supporting_markers=['CD3D'],  # Weak
        pmids=['99999999'],
        evidence=[...],
        is_novel=False
    )
]

# Top candidate: "T cells" (92% confidence)
```

### Step 5: Verify Assignment

```python
from tools.dynamic_knowledge import verify_cell_type_assignment

is_valid, reason = verify_cell_type_assignment(
    cluster_markers=cluster_markers,
    candidate=candidates[0],
    min_marker_overlap=2
)

# is_valid=True, reason="Valid"
# → Assign cluster as "T cells"
```

---

## Novel Population Detection

If no strong candidate is found, flag as novel:

```python
from tools.dynamic_knowledge import detect_novel_population

is_novel = detect_novel_population(
    cluster_markers=['GENE1', 'GENE2', 'GENE3'],
    candidates=[...],
    confidence_threshold=0.3
)

if is_novel:
    print("⚠️ Novel population detected!")
    print("Markers: GENE1, GENE2, GENE3")
    print("→ Perform deeper literature review")
    print("→ Consider tissue-specific or disease-specific cell types")
```

**Novel criteria:**
1. Top candidate confidence < 0.3
2. No papers found for marker combination
3. Mixed signals (multiple weak conflicting candidates)

---

## Integration with Annotation Workflow

### Tier 1 Example

```python
import scanpy as sc
from tools.dynamic_knowledge import search_cell_type_from_markers

# Load data
adata = sc.read_h5ad('preprocessed_data.h5ad')

# Cluster
sc.tl.leiden(adata, resolution=0.8, key_added='tier1_cluster')

# DE
sc.tl.rank_genes_groups(adata, groupby='tier1_cluster', method='wilcoxon')
de_df = sc.get.rank_genes_groups_df(adata, group=None)

# For each cluster
for cluster in adata.obs['tier1_cluster'].unique():
    # Get top markers
    cluster_de = de_df[
        (de_df['group'] == cluster) &
        (de_df['pct_nz_group'] >= 0.25) &
        (de_df['logfoldchanges'] >= 1.0) &
        (de_df['pvals_adj'] < 0.05)
    ].nlargest(20, 'logfoldchanges')

    top_markers = cluster_de['names'].tolist()

    # Search literature (data-driven, NO hardcoding)
    candidates = search_cell_type_from_markers(
        markers=top_markers,
        tissue=adata.uns.get('tissue', None),
        species=adata.uns.get('species', 'human')
    )

    # Check if novel
    if candidates[0].is_novel:
        print(f"Cluster {cluster}: Novel population")
        print(f"  Top markers: {top_markers[:5]}")
        print(f"  → Flag for manual review")
    else:
        # Verify and assign
        is_valid, reason = verify_cell_type_assignment(
            cluster_markers=top_markers,
            candidate=candidates[0]
        )

        if is_valid:
            print(f"Cluster {cluster}: {candidates[0].name}")
            print(f"  Confidence: {candidates[0].confidence:.2f}")
            print(f"  PMIDs: {candidates[0].pmids}")
        else:
            print(f"Cluster {cluster}: Verification failed")
            print(f"  Reason: {reason}")
```

---

## Advantages Over Hardcoded Systems

| Aspect | Hardcoded System | Dynamic Knowledge System |
|--------|------------------|--------------------------|
| **Discovery** | Only finds known types | Finds novel/rare types |
| **Bias** | Biased to common types | Unbiased, data-driven |
| **Species** | Human-centric | Works for any species |
| **Tissue** | Generic markers | Tissue-specific discovery |
| **Updates** | Manual updates | Auto-updates (new papers) |
| **Flexibility** | Rigid categories | Flexible, context-aware |

---

## Implementation Notes

### PubMed API Integration

The `pubmed_search()` function in `dynamic_knowledge.py` is a **stub**.
You must implement it using:

**Option 1: MCP Tool (Recommended)**
```python
def pubmed_search(query: str, max_results: int = 5) -> List[Dict]:
    """Use MCP pubmed_search tool."""
    # Call MCP tool via agent
    result = mcp_call('pubmed_search', {
        'query': query,
        'max_results': max_results
    })
    return result['results']
```

**Option 2: Biopython (Fallback)**
```python
from Bio import Entrez

Entrez.email = "your_email@example.com"
Entrez.api_key = "your_api_key"

def pubmed_search(query: str, max_results: int = 5) -> List[Dict]:
    """Use Bio.Entrez for PubMed search."""
    # Search
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    pmids = record['IdList']

    # Fetch details
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="xml")
    abstracts = Entrez.read(handle)

    # Parse and return
    results = []
    for article in abstracts['PubmedArticle']:
        results.append({
            'pmid': str(article['MedlineCitation']['PMID']),
            'title': article['MedlineCitation']['Article']['ArticleTitle'],
            'abstract': article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0],
            'journal': article['MedlineCitation']['Article']['Journal']['Title'],
            'year': int(article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'])
        })

    return results
```

See `tools/literature_verification.md` for full Bio.Entrez implementation.

---

## Testing

```python
# Test marker combination generation
from tools.dynamic_knowledge import generate_marker_combinations

markers = ['A', 'B', 'C', 'D', 'E']
combos = generate_marker_combinations(markers, max_size=3)
assert ['A', 'B'] in combos
assert ['A', 'B', 'C'] in combos
print("✓ Combination generation works")

# Test cell type extraction
from tools.dynamic_knowledge import extract_cell_types_from_text

text = "CD3+ T cells and CD19+ B cells were identified."
types = extract_cell_types_from_text(text)
assert 'T cells' in types or 'T' in types
assert 'B cells' in types or 'B' in types
print("✓ Cell type extraction works")

# Test normalization
from tools.dynamic_knowledge import normalize_cell_type_name

assert normalize_cell_type_name("T cell") == "T Cells"
assert normalize_cell_type_name("CD8+ T cells") == "T Cells"
assert normalize_cell_type_name("memory B cell") == "B Cells"
print("✓ Normalization works")
```

---

## Troubleshooting

### Issue: No candidates found

**Cause**: Markers too specific or not in literature
**Solution**:
1. Try fewer markers (2-3 instead of 5)
2. Remove tissue-specific modifiers
3. Try synonyms (e.g., "TdT" → "DNTT")
4. Check if markers are human gene symbols

### Issue: Multiple weak candidates

**Cause**: Mixed cell population or doublets
**Solution**:
1. Check cross-contamination (Step 2.6 in tier1.md)
2. Run doublet detection
3. Increase clustering resolution
4. Flag as mixed population

### Issue: Wrong cell type assigned

**Cause**: Marker overlap with wrong type
**Solution**:
1. Use more markers (increase `min_marker_overlap`)
2. Add tissue context
3. Check PMID evidence manually
4. Verify with TF/Pathway analysis (Tier 2/3)

---

## Future Enhancements

1. **Species-specific databases**: Mouse, rat, zebrafish markers
2. **Tissue-specific priors**: Weight tissue-relevant papers higher
3. **Disease context**: Cancer vs. healthy cell types
4. **Graph-based reasoning**: Build cell type ontology from papers
5. **LLM integration**: Use LLM to read abstracts and extract relationships

---

## Key Takeaways

```
✓ NO hardcoded cell types or markers
✓ Purely literature-driven from DE results
✓ Discovers novel/rare populations automatically
✓ Unbiased, data-driven approach
✓ Works for any species, tissue, disease context
✓ Integrates with existing annotation workflow
```

**This is the future of cell type annotation: data-driven, literature-verified, and discovery-oriented.**
