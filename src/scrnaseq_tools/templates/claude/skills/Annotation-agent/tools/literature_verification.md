# Literature Verification Tool

PubMed API를 사용한 실시간 마커 문헌 검증.

---

## API 설정

```python
# NCBI E-utilities API (server.py에서 사용)
NCBI_API_KEY = "40b96e1094387e03e7f9133ec6e33e881108"
NCBI_EMAIL = "kwy7605@gmail.com"
# Rate limit with API key: 10 requests/second
```

---

## 검색 방법 우선순위

```
╔══════════════════════════════════════════════════════════════════════╗
║  1. [PRIMARY] MCP pubmed_search tool                                 ║
║     - Rate limit: 10 req/sec with API key                           ║
║     - Returns: PMID, title, authors, journal, year                  ║
║                                                                      ║
║  2. [FALLBACK] Bio.Entrez (Biopython)                               ║
║     - Direct Python interface to NCBI                               ║
║     - Same API, different wrapper                                   ║
║                                                                      ║
║  3. [FALLBACK] WebSearch with "PMID" query                          ║
║     - For when API is unavailable                                   ║
╚══════════════════════════════════════════════════════════════════════╝
```

---

## MCP Tool Usage

### pubmed_search

```json
{
  "name": "pubmed_search",
  "arguments": {
    "query": "DNTT Pro-B cell marker",
    "max_results": 3
  }
}
```

**Response:**
```json
{
  "query": "DNTT Pro-B cell marker",
  "count": 3,
  "total_found": 245,
  "results": [
    {
      "pmid": "2785044",
      "title": "Terminal deoxynucleotidyl transferase...",
      "authors": ["Choudhury BA", "Ebtehaj K", "..."],
      "journal": "J Immunol",
      "year": "1989"
    }
  ]
}
```

### verify_reference

```json
{
  "name": "verify_reference",
  "arguments": {
    "pmid": "11007475",
    "markers": ["AICDA", "AID"],
    "cell_type": "germinal center B cell"
  }
}
```

**Response:**
```json
{
  "pmid": "11007475",
  "title": "Class switch recombination and hypermutation...",
  "status": "VERIFIED",
  "found_markers": ["AICDA", "AID"],
  "missing_markers": [],
  "marker_ratio": 1.0,
  "cell_type_found": true
}
```

---

## Python Implementation

```python
from Bio import Entrez
import time

# NCBI credentials
Entrez.email = "kwy7605@gmail.com"
Entrez.api_key = "40b96e1094387e03e7f9133ec6e33e881108"


def search_marker_reference(marker: str, cell_type: str, max_results: int = 3) -> list:
    """
    Search PubMed for marker-cell type associations.

    Args:
        marker: Gene symbol (e.g., 'DNTT')
        cell_type: Cell type context (e.g., 'Pro-B cell')
        max_results: Maximum results to return

    Returns:
        List of dicts with pmid, title, journal, year
    """
    query = f"{marker} {cell_type}"

    try:
        # Search
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()

        pmids = record.get("IdList", [])
        if not pmids:
            return []

        # Get summaries
        handle = Entrez.esummary(db="pubmed", id=",".join(pmids))
        summaries = Entrez.read(handle)
        handle.close()

        results = []
        for doc in summaries:
            results.append({
                'pmid': doc.get('Id', 'N/A'),
                'title': doc.get('Title', 'N/A'),
                'journal': doc.get('Source', 'N/A'),
                'year': doc.get('PubDate', 'N/A')[:4]
            })

        return results

    except Exception as e:
        print(f"PubMed search error: {e}")
        return []


def verify_markers_with_pmid(markers: list, cell_type: str) -> dict:
    """
    Verify a list of markers and get PMIDs.

    Args:
        markers: List of gene symbols
        cell_type: Cell type context

    Returns:
        dict with 'verified' list and 'not_found' list
    """
    verified = []
    not_found = []

    for marker in markers:
        results = search_marker_reference(marker, cell_type)
        time.sleep(0.15)  # Rate limiting

        if results:
            verified.append({
                'marker': marker,
                'pmid': results[0]['pmid'],
                'title': results[0]['title'],
                'status': 'VERIFIED'
            })
        else:
            not_found.append(marker)

    return {
        'verified': verified,
        'not_found': not_found,
        'rate': len(verified) / len(markers) if markers else 0
    }
```

---

## Reasoning Output Format

### Literature Search (PMID Required)

```markdown
**Literature Verification**:

Query: "DNTT Pro-B cell"
Results:
| PMID | Title | Journal | Year | Status |
|------|-------|---------|------|--------|
| 2785044 | TdT+ B cell precursors... | J Immunol | 1989 | ✓ VERIFIED |
| 34530764 | RAG in lymphocyte development | Front Immunol | 2023 | ✓ VERIFIED |

Query: "AICDA germinal center"
Results:
| PMID | Title | Journal | Year | Status |
|------|-------|---------|------|--------|
| 11007475 | AID class switch recombination | Cell | 2000 | ✓ VERIFIED |
```

---

## Error Handling

```python
def search_with_retry(query: str, max_retries: int = 3) -> list:
    """Search with automatic retry on failure."""
    for attempt in range(max_retries):
        try:
            return search_marker_reference(query.split()[0], " ".join(query.split()[1:]))
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(1)  # Wait before retry
                continue
            raise e
    return []
```

---

## Integration with Annotation Agent

### In Reasoning Template

```markdown
### Literature (PRE)

**Marker Verification (API Search)**:

| Marker | Cell Type | PMID | Title | Status |
|--------|-----------|------|-------|--------|
| DNTT | Pro-B | 2785044 | TdT+ B cell precursors | ✓ |
| RAG1 | Pro-B | 34530764 | RAG in lymphocyte development | ✓ |
| ERG | Pro-B | 19471024 | ERG in hematopoiesis | ✓ |

Search method: NCBI E-utilities API
API status: OK (10 req/sec)
```

---

## Reference Output Format (MANDATORY — Per-Marker Dict)

Every `pmid_verified` field MUST be a **per-marker dictionary**, NOT a flat list.
Each marker gene is a key, and its value is a list of verified references.

```json
"pmid_verified": {
  "DNTT": [
    {
      "pmid": "2785044",
      "author": "Choudhury BA et al.",
      "year": 1989,
      "journal": "J Immunol",
      "title": "TdT+ B cell precursors..."
    }
  ],
  "RAG1": [
    {
      "pmid": "34530764",
      "author": "Song Y et al.",
      "year": 2023,
      "journal": "Front Immunol",
      "title": "RAG in lymphocyte development"
    }
  ],
  "ERG": []
}
```

```
❌ FORBIDDEN: "pmid_verified": [{"pmid": "...", ...}]        ← flat list
✅ REQUIRED:  "pmid_verified": {"Gene1": [...], "Gene2": [...]}  ← per-marker dict
```

| Field | Required | Description |
|-------|----------|-------------|
| `pmid` | YES | PubMed ID (numeric string) |
| `author` | YES | "FirstAuthor et al." format |
| `year` | YES | Publication year (integer) |
| `journal` | YES | Journal name |
| `title` | YES | Full paper title |

---

## Key Points

1. **모든 참조에 PMID 필수** - API 검색으로 실시간 확보
2. **Per-Marker 검증 필수** - 각 marker gene마다 개별 검증
3. **Rate limiting 준수** - 0.1-0.15초 간격
4. **Retry 로직** - API 오류 시 자동 재시도
5. **미검증 마커 = 빈 리스트** - reference 못 찾은 마커도 key로 포함 (빈 `[]`)
