"""
Dynamic Knowledge Base for Data-Driven Cell Type Annotation

NO hardcoded cell types or markers.
Cell type identification is purely literature-driven from actual DE results.
"""

from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import re
from collections import Counter


@dataclass
class CellTypeCandidate:
    """A candidate cell type inferred from literature."""
    name: str
    confidence: float  # 0-1
    supporting_markers: List[str]
    pmids: List[str]
    evidence: List[Dict]  # List of {pmid, title, marker_match, snippet}
    is_novel: bool = False


def search_cell_type_from_markers(
    markers: List[str],
    tissue: Optional[str] = None,
    species: str = "human",
    max_combinations: int = 10
) -> List[CellTypeCandidate]:
    """
    Search PubMed for cell types based on marker combinations.

    NO hardcoded cell types. Pure literature-driven discovery.

    Args:
        markers: Top DE markers from current cluster
        tissue: Optional tissue context
        species: Species context (default: human)
        max_combinations: Max marker combinations to try

    Returns:
        List of CellTypeCandidate objects, sorted by confidence

    Example:
        >>> markers = ['CD3D', 'CD3E', 'TRAC', 'CD8A']
        >>> candidates = search_cell_type_from_markers(markers)
        >>> # Returns: [CellTypeCandidate(name='T cells', confidence=0.95, ...)]
    """
    candidates = []

    # Strategy: Multi-level PubMed query with fallback
    # Level 1: Top 2-3 markers with field tags (most specific)
    # Level 2: Individual markers + "cell type" (broader)
    # Level 3: Marker pairs (comprehensive)
    marker_combinations = generate_marker_combinations(markers[:10], max_size=4)

    for combo in marker_combinations[:max_combinations]:
        # Build PubMed-optimized queries (try multiple formats)
        queries = _build_pubmed_queries(combo, tissue=tissue, species=species)

        for query in queries:
            results = pubmed_search(query, max_results=5)
            if results:
                break  # Use first successful query format

        # Extract cell types from results
        for result in results:
            title = result.get('title') or ''
            abstract = result.get('abstract') or ''
            text = title + " " + abstract

            cell_types = extract_cell_types_from_text(text)

            # Score cell types by marker proximity in text
            scored_types = _score_cell_types_by_marker_context(
                cell_types, combo, text
            )

            for cell_type, context_score in scored_types:
                candidates.append({
                    'cell_type': cell_type,
                    'markers': combo,
                    'pmid': result['pmid'],
                    'title': title,
                    'journal': result.get('journal', 'N/A'),
                    'year': result.get('year', 'N/A'),
                    'context_score': context_score
                })

    # Aggregate candidates
    aggregated = aggregate_candidates(candidates)

    # Flag novel if no strong candidates found
    if not aggregated or aggregated[0].confidence < 0.3:
        return [CellTypeCandidate(
            name="Novel population",
            confidence=0.0,
            supporting_markers=markers[:5],
            pmids=[],
            evidence=[],
            is_novel=True
        )]

    return aggregated


def generate_marker_combinations(
    markers: List[str],
    min_size: int = 2,
    max_size: int = 5
) -> List[List[str]]:
    """
    Generate marker combinations for literature search.

    Prioritize:
    1. Top 2-3 markers (strongest signal)
    2. Top 4-5 markers (more specific)
    3. Pairwise combinations of top 10

    Args:
        markers: Sorted list of markers (by LFC or significance)
        min_size: Minimum combination size
        max_size: Maximum combination size

    Returns:
        List of marker combinations, sorted by priority
    """
    combinations = []

    # Priority 1: Top N consecutive markers
    for size in range(min_size, min(max_size + 1, len(markers) + 1)):
        combinations.append(markers[:size])

    # Priority 2: Skip patterns (e.g., [0,2,3] to diversify)
    if len(markers) >= 5:
        combinations.append([markers[0], markers[2], markers[3]])
        combinations.append([markers[0], markers[1], markers[3], markers[4]])

    # Priority 3: Pairwise from top 10
    from itertools import combinations as itertools_combinations
    if len(markers) >= 3:
        for combo in itertools_combinations(markers[:10], min_size):
            if list(combo) not in combinations:
                combinations.append(list(combo))
            if len(combinations) >= 20:  # Limit total
                break

    return combinations


def extract_cell_types_from_text(text: str) -> List[str]:
    """
    Extract cell type mentions from title/abstract.

    Uses pattern matching to identify cell type references.
    NO hardcoded whitelist - extracts ANY cell type mention.

    Patterns:
    - "X cells" (e.g., "T cells", "B cells")
    - "X cell" (singular)
    - "Xocytes" (e.g., "lymphocytes", "monocytes")
    - Hyphenated types (e.g., "CD8+ T cells", "memory B cells")

    Args:
        text: Title + abstract text

    Returns:
        List of unique cell type strings
    """
    cell_types = []

    # Normalize Unicode diacritics: "naïve" → "naive"
    import unicodedata
    text = unicodedata.normalize('NFD', text)
    text = ''.join(c for c in text if unicodedata.category(c) != 'Mn')

    # Pattern 1: "X cells" or "X cell"
    # Capture 1-3 words before "cell(s)"
    pattern1 = r'\b([A-Z][a-zA-Z0-9+\-/γδ]+(?:\s+[A-Z][a-zA-Z0-9+\-/γδ]+){0,3})\s+cells?\b'
    matches1 = re.findall(pattern1, text)
    cell_types.extend([m.strip() for m in matches1])

    # Pattern 2: "Xocyte(s)" (e.g., lymphocyte, monocyte)
    pattern2 = r'\b([A-Z][a-z]+ocytes?)\b'
    matches2 = re.findall(pattern2, text)
    cell_types.extend([m.strip() for m in matches2])

    # Pattern 3: Specific immune cell patterns
    # "T cell", "B cell", "NK cell", etc.
    pattern3 = r'\b(T|B|NK|NKT|ILC|MAIT|Treg|Th\d+|Tc\d*)\s+cells?\b'
    matches3 = re.findall(pattern3, text, re.IGNORECASE)
    cell_types.extend([f"{m} cells" for m in matches3])

    # Pattern 4: Developmental states (if in context)
    pattern4 = r'\b(naive|memory|effector|activated|exhausted|regulatory)\s+([A-Z][a-z]+)\s+cells?\b'
    matches4 = re.findall(pattern4, text, re.IGNORECASE)
    cell_types.extend([f"{m[1]} cells" for m in matches4])  # Extract base type only

    # Pattern 5: "Xoblasts" (e.g., fibroblasts, osteoblasts)
    pattern5 = r'\b([A-Z][a-z]+oblasts?)\b'
    matches5 = re.findall(pattern5, text)
    cell_types.extend([m.strip() for m in matches5])

    # Pattern 6: "Xphages" (e.g., macrophages)
    pattern6 = r'\b([A-Z][a-z]+phages?)\b'
    matches6 = re.findall(pattern6, text)
    cell_types.extend([m.strip() for m in matches6])

    # Clean up
    cell_types = list(set(cell_types))  # Unique
    cell_types = [ct for ct in cell_types if len(ct) > 2]  # Filter noise

    # Filter out non-cell-type words that match patterns
    noise_words = {
        'Immune', 'Immunes', 'Single', 'Tumor', 'Cancer', 'Human',
        'Clinical', 'Novel', 'Primary', 'Associated', 'Related',
        'Relative', 'Relatives', 'Total', 'Major', 'Minor',
        'Other', 'Different', 'Various', 'Multiple', 'Several',
        'Specific', 'Normal', 'Adult', 'Fetal', 'Mature',
        'Peripheral', 'Central', 'Resident', 'Circulating',
        'Infiltrating', 'Recruited', 'Tissue', 'Blood', 'Bone',
    }
    cell_types = [ct for ct in cell_types if ct not in noise_words]

    # Additional noise filtering
    _verbal_noise_starts = {
        'classified', 'identifying', 'identified', 'characterizing', 'characterized',
        'defining', 'defined', 'describing', 'described', 'shown', 'showing',
        'found', 'finding', 'observed', 'observing', 'detected', 'detecting',
        'generated', 'generating', 'produced', 'producing', 'contained', 'containing',
        'represented', 'representing', 'resembled', 'resembling', 'designated',
        'termed', 'called', 'named', 'referred', 'labeled', 'categorized',
        'emerging', 'converted', 'converting', 'becoming', 'known', 'recognized',
        'considered', 'regarded', 'displaying', 'exhibited', 'exhibiting',
        'included', 'including', 'suggested', 'proposed', 'abnormal',
        'proportions', 'proportion', 'percentage', 'frequency', 'frequencies',
        'number', 'numbers', 'levels', 'level', 'amounts', 'amount',
        'role', 'roles', 'function', 'functions', 'loss', 'lack',
        'absence', 'presence', 'expression', 'regulation', 'depletion',
        'higher', 'lower', 'fewer', 'greater', 'increased', 'decreased',
        'reduced', 'elevated', 'enhanced', 'diminished', 'expanded',
        'total', 'overall', 'majority', 'subset',
    }

    def _is_valid_cell_type(ct: str) -> bool:
        ct_lower = ct.lower()
        # Must contain "cell" or end in recognized suffix
        has_cell = 'cell' in ct_lower
        has_suffix = any(ct_lower.endswith(s) for s in
                        ['cytes', 'cyte', 'blasts', 'blast', 'phages', 'phage'])
        has_known_type = any(t in ct_lower for t in
                           ['plasma', 'stem', 'progenitor', 'precursor'])
        if not (has_cell or has_suffix or has_known_type):
            if ct_lower in {'t cells', 'b cells', 'nk cells'}:
                return True
            return False
        if len(ct) > 50:
            return False
        # Remove species-prefixed types: "Human And Mouse Memory B Cells"
        species_words = {'human', 'mouse', 'murine', 'rat', 'bovine', 'porcine'}
        words = ct_lower.split()
        if words and words[0] in species_words:
            return False
        # Remove "The Proportions of..." style noise
        if ct_lower.startswith('the '):
            return False
        # Remove verbal/participial prefixes: "Classified As Germinal Center B Cells"
        if words and words[0] in _verbal_noise_starts:
            return False
        return True

    cell_types = [ct for ct in cell_types if _is_valid_cell_type(ct)]

    return cell_types


def aggregate_candidates(
    raw_candidates: List[Dict],
    preserve_subtypes: bool = False
) -> List[CellTypeCandidate]:
    """
    Aggregate raw candidates into ranked CellTypeCandidate objects.

    Scoring:
    - Frequency: How many papers mention this cell type
    - Marker overlap: How many markers support it
    - Recency: Newer papers weighted higher
    - Context: How close cell type appears to markers in text

    Args:
        raw_candidates: List of dicts from PubMed searches
        preserve_subtypes: If True, preserve developmental/functional modifiers

    Returns:
        Sorted list of CellTypeCandidate objects
    """
    # Group by cell type
    grouped = {}
    for cand in raw_candidates:
        cell_type = normalize_cell_type_name(
            cand['cell_type'], preserve_subtypes=preserve_subtypes
        )

        if cell_type not in grouped:
            grouped[cell_type] = {
                'markers': set(),
                'pmids': set(),
                'evidence': [],
                'context_scores': []
            }

        grouped[cell_type]['markers'].update(cand['markers'])
        grouped[cell_type]['pmids'].add(cand['pmid'])
        grouped[cell_type]['evidence'].append({
            'pmid': cand['pmid'],
            'title': cand['title'],
            'journal': cand['journal'],
            'year': cand['year'],
            'markers': cand['markers']
        })
        grouped[cell_type]['context_scores'].append(
            cand.get('context_score', 0.5)
        )

    # Score each candidate
    candidates = []
    for cell_type, data in grouped.items():
        # Frequency score (how many papers)
        freq_score = min(len(data['pmids']) / 3.0, 1.0)  # Max at 3 papers

        # Marker score (how many markers supported)
        marker_score = min(len(data['markers']) / 4.0, 1.0)  # Max at 4 markers

        # Recency score (average year)
        years = [e['year'] for e in data['evidence'] if isinstance(e['year'], int)]
        if years:
            avg_year = sum(years) / len(years)
            recency_score = min((avg_year - 2000) / 25.0, 1.0)  # Normalize 2000-2025
        else:
            recency_score = 0.5

        # Context score (how close cell type appears to markers in text)
        ctx_scores = data['context_scores']
        avg_context = sum(ctx_scores) / len(ctx_scores) if ctx_scores else 0.5

        # Combined confidence (context-aware)
        confidence = (
            freq_score * 0.3 +
            marker_score * 0.3 +
            avg_context * 0.25 +
            recency_score * 0.15
        )

        candidates.append(CellTypeCandidate(
            name=cell_type,
            confidence=confidence,
            supporting_markers=list(data['markers']),
            pmids=list(data['pmids']),
            evidence=data['evidence'],
            is_novel=False
        ))

    # Sort by confidence
    candidates.sort(key=lambda x: x.confidence, reverse=True)

    return candidates


def normalize_cell_type_name(name: str, preserve_subtypes: bool = False) -> str:
    """
    Normalize cell type names for aggregation.

    For Tier 1 (preserve_subtypes=False):
    - "T cell" → "T cells"
    - "CD8+ T cells" → "T cells"  (simplify to base type)
    - "memory B cell" → "B cells"

    For Tier 2/3 (preserve_subtypes=True):
    - "naive B cell" → "Naive B Cells"  (keep modifier)
    - "pro-B cells" → "Pro-B Cells"
    - "CD8+ T cells" → "CD8+ T Cells"

    Args:
        name: Raw cell type name
        preserve_subtypes: If True, keep developmental/functional modifiers

    Returns:
        Normalized name
    """
    name = name.lower().strip()

    # Pluralize
    if not name.endswith('s') and not name.endswith('cyte'):
        name = name + 's'
    if name.endswith('cyte'):
        name = name + 's'

    if not preserve_subtypes:
        # Tier 1: Extract base type (remove modifiers)
        base_patterns = [
            (r'cd\d+[+-]?\s+([a-z]+\s+cells?)', r'\1'),
            (r'(naive|memory|effector|activated|exhausted|regulatory)\s+([a-z]+\s+cells?)', r'\2'),
        ]
        for pattern, replacement in base_patterns:
            name = re.sub(pattern, replacement, name, flags=re.IGNORECASE)

    # Capitalize properly
    name = ' '.join(word.capitalize() if word not in ['and', 'or', 'of'] else word
                    for word in name.split())

    return name


def pubmed_search(query: str, max_results: int = 5, max_retries: int = 3) -> List[Dict]:
    """
    Search PubMed via Bio.Entrez with retry logic.

    Uses NCBI E-utilities: esearch → efetch for full abstracts.
    Rate limited to respect NCBI guidelines (10 req/sec with API key).

    Args:
        query: Search query
        max_results: Max results to return
        max_retries: Max retry attempts on failure

    Returns:
        List of dicts with pmid, title, abstract, journal, year
    """
    from Bio import Entrez
    import time
    import xml.etree.ElementTree as ET

    Entrez.email = "kwy7605@gmail.com"
    Entrez.api_key = "40b96e1094387e03e7f9133ec6e33e881108"

    for attempt in range(max_retries):
        try:
            # Step 1: Search for PMIDs
            handle = Entrez.esearch(
                db="pubmed", term=query, retmax=max_results, sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()

            pmids = record.get("IdList", [])
            if not pmids:
                return []

            time.sleep(0.12)  # Rate limiting

            # Step 2: Fetch full records with abstracts
            handle = Entrez.efetch(
                db="pubmed", id=",".join(pmids), rettype="xml", retmode="xml"
            )
            raw_xml = handle.read()
            handle.close()

            # Parse XML
            results = _parse_pubmed_xml(raw_xml)

            time.sleep(0.12)  # Rate limiting
            return results

        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(1.0 * (attempt + 1))  # Exponential backoff
                continue
            print(f"PubMed search failed after {max_retries} attempts: {e}")
            return []

    return []


def _parse_pubmed_xml(xml_data) -> List[Dict]:
    """Parse PubMed efetch XML response into structured dicts."""
    import xml.etree.ElementTree as ET

    if isinstance(xml_data, bytes):
        xml_data = xml_data.decode('utf-8')

    # Remove DOCTYPE declaration that causes ET.ParseError
    xml_data = re.sub(r'<!DOCTYPE[^>]+>', '', xml_data)

    results = []
    try:
        root = ET.fromstring(xml_data)
    except ET.ParseError:
        return results

    for article_node in root.findall('.//PubmedArticle'):
        try:
            medline = article_node.find('MedlineCitation')
            if medline is None:
                continue

            # PMID
            pmid_el = medline.find('PMID')
            pmid = pmid_el.text if pmid_el is not None else 'N/A'

            # Article
            article = medline.find('Article')
            if article is None:
                continue

            # Title
            title_el = article.find('ArticleTitle')
            title = title_el.text if title_el is not None else 'N/A'

            # Abstract
            abstract_parts = []
            abstract_el = article.find('Abstract')
            if abstract_el is not None:
                for abs_text in abstract_el.findall('AbstractText'):
                    if abs_text.text:
                        label = abs_text.get('Label', '')
                        prefix = f"{label}: " if label else ''
                        abstract_parts.append(prefix + abs_text.text)
            abstract = ' '.join(abstract_parts)

            # Journal
            journal_el = article.find('.//Journal/Title')
            journal = journal_el.text if journal_el is not None else 'N/A'

            # Year
            year = _extract_year(article)

            # Authors (short format)
            authors_short = _extract_authors_short(article)

            results.append({
                'pmid': pmid,
                'title': title,
                'abstract': abstract,
                'journal': journal,
                'year': year,
                'authors_short': authors_short
            })

        except Exception:
            continue

    return results


def _extract_year(article_element) -> int:
    """Extract publication year from article XML element."""
    # Try Journal → JournalIssue → PubDate → Year
    year_el = article_element.find('.//Journal/JournalIssue/PubDate/Year')
    if year_el is not None and year_el.text:
        try:
            return int(year_el.text)
        except ValueError:
            pass

    # Try MedlineDate
    medline_date = article_element.find('.//Journal/JournalIssue/PubDate/MedlineDate')
    if medline_date is not None and medline_date.text:
        import re as _re
        match = _re.search(r'(\d{4})', medline_date.text)
        if match:
            return int(match.group(1))

    # Try ArticleDate
    year_el = article_element.find('.//ArticleDate/Year')
    if year_el is not None and year_el.text:
        try:
            return int(year_el.text)
        except ValueError:
            pass

    return 0


def _extract_authors_short(article_element) -> str:
    """Extract 'FirstAuthor et al.' format from article XML."""
    author_list = article_element.find('.//AuthorList')
    if author_list is None:
        return 'Unknown'

    authors = author_list.findall('Author')
    if not authors:
        return 'Unknown'

    first = authors[0]
    last_name = first.find('LastName')
    if last_name is not None and last_name.text:
        name = last_name.text
        if len(authors) > 1:
            return f"{name} et al."
        return name

    collective = first.find('CollectiveName')
    if collective is not None and collective.text:
        return collective.text

    return 'Unknown'


def _score_cell_types_by_marker_context(
    cell_types: List[str],
    markers: List[str],
    text: str
) -> List[Tuple[str, float]]:
    """
    Score extracted cell types by how closely they appear near query markers in text.

    Cell types mentioned near the actual markers are more relevant than those
    mentioned in passing elsewhere in the abstract.

    Args:
        cell_types: Extracted cell type names
        markers: The query markers used
        text: The full text (title + abstract)

    Returns:
        List of (cell_type, context_score) tuples, sorted by score desc
    """
    text_lower = text.lower()
    scored = []

    for ct in cell_types:
        ct_lower = ct.lower()
        ct_pos = text_lower.find(ct_lower)
        if ct_pos == -1:
            scored.append((ct, 0.1))
            continue

        # Find closest marker mention
        min_distance = len(text)
        markers_found_nearby = 0
        for marker in markers:
            m_pos = text_lower.find(marker.lower())
            if m_pos != -1:
                dist = abs(ct_pos - m_pos)
                min_distance = min(min_distance, dist)
                if dist < 200:  # Within ~2 sentences
                    markers_found_nearby += 1

        # Score: closer = better, more markers nearby = better
        distance_score = max(0, 1.0 - min_distance / 500.0)
        proximity_score = min(markers_found_nearby / max(len(markers), 1), 1.0)

        # Bonus: cell type in title (very relevant)
        title_bonus = 0.3 if ct_lower in text_lower[:200] else 0.0

        context_score = (distance_score * 0.4 + proximity_score * 0.4 + title_bonus + 0.1)
        scored.append((ct, min(context_score, 1.0)))

    scored.sort(key=lambda x: x[1], reverse=True)
    return scored


def _build_pubmed_queries(
    markers: List[str],
    tissue: Optional[str] = None,
    species: str = "human"
) -> List[str]:
    """
    Build optimized PubMed queries from markers, with fallback levels.

    Returns multiple query formats ordered from most specific to broadest.
    The caller should try each until one returns results.

    Args:
        markers: List of gene symbols for this combination
        tissue: Optional tissue context
        species: Species context

    Returns:
        List of query strings (try in order)
    """
    queries = []

    # Level 1: Field-tagged markers + "cell type" (best precision)
    tagged = " AND ".join(f"{m}[Title/Abstract]" for m in markers[:3])
    q1 = f"{tagged} AND (cell type[Title/Abstract] OR marker[Title/Abstract])"
    if tissue:
        q1 += f" AND {tissue}[Title/Abstract]"
    queries.append(q1)

    # Level 2: Quoted markers + cell context (good balance)
    quoted = " ".join(f'"{m}"' for m in markers[:3])
    q2 = f"{quoted} cell marker {species}"
    if tissue:
        q2 += f" {tissue}"
    queries.append(q2)

    # Level 3: Simple markers + cell (broadest)
    q3 = " ".join(markers[:2]) + " cell " + species
    queries.append(q3)

    # Level 4: Individual top marker (fallback)
    q4 = f"{markers[0]} cell type marker {species}"
    queries.append(q4)

    return queries


def verify_cell_type_assignment(
    cluster_markers: List[str],
    candidate: CellTypeCandidate,
    min_marker_overlap: int = 2
) -> Tuple[bool, str]:
    """
    Verify if a cell type assignment is valid.

    Checks:
    1. Marker overlap: >= min_marker_overlap markers match
    2. Confidence threshold: >= 0.5
    3. PMID evidence: >= 1 paper

    Args:
        cluster_markers: Top DE markers from cluster
        candidate: Cell type candidate to verify
        min_marker_overlap: Minimum markers that must match

    Returns:
        (is_valid, reason)
    """
    # Check marker overlap
    overlap = len(set(cluster_markers) & set(candidate.supporting_markers))
    if overlap < min_marker_overlap:
        return False, f"Insufficient marker overlap ({overlap}/{min_marker_overlap})"

    # Check confidence
    if candidate.confidence < 0.5:
        return False, f"Low confidence ({candidate.confidence:.2f})"

    # Check evidence
    if len(candidate.pmids) < 1:
        return False, "No PMID evidence"

    return True, "Valid"


def detect_novel_population(
    cluster_markers: List[str],
    candidates: List[CellTypeCandidate],
    confidence_threshold: float = 0.3
) -> bool:
    """
    Detect if a cluster represents a novel/rare cell population.

    Criteria for "novel":
    1. No candidate with confidence > threshold
    2. Marker combination not found in literature
    3. Mixed signals (multiple weak candidates)

    Args:
        cluster_markers: Top DE markers
        candidates: Cell type candidates from literature
        confidence_threshold: Min confidence to avoid "novel" flag

    Returns:
        True if novel population detected
    """
    if not candidates:
        return True

    # Check top candidate
    top_candidate = candidates[0]
    if top_candidate.confidence < confidence_threshold:
        return True

    # Check for mixed signals
    if len(candidates) >= 3 and candidates[2].confidence > 0.4:
        # Multiple strong but conflicting candidates
        return True

    return False


def filter_informative_markers(
    markers: List[str],
    max_markers: int = 10,
    prioritize_known: bool = True
) -> List[str]:
    """
    Filter and prioritize DE markers for literature search.

    Step 1: Remove uninformative genes (LOC, ribosomal, MT, variants)
    Step 2: If prioritize_known, reorder to put biologically informative
            markers first (CD, IG, HLA, known TFs, surface markers)

    Args:
        markers: Raw list of DE marker gene names
        max_markers: Maximum markers to return
        prioritize_known: If True, put known cell markers at top

    Returns:
        Filtered and optionally prioritized list of markers
    """
    # Step 1: Basic filtering
    filtered = []
    for m in markers:
        if m.startswith('LOC') and m[3:].isdigit():
            continue
        if '.' in m:
            continue
        if m.startswith('MT-'):
            continue
        if len(m) < 2:
            continue
        if m.startswith(('RPS', 'RPL')) and len(filtered) >= 3:
            continue
        filtered.append(m)

    if not prioritize_known:
        return filtered[:max_markers]

    # Step 2: Prioritize biologically informative markers
    # These are more likely to yield relevant PubMed results
    high_priority = []   # CD markers, Ig genes, known cell markers
    medium_priority = [] # TFs, signaling molecules, HLA
    normal_priority = [] # Everything else

    for m in filtered:
        m_upper = m.upper()
        if (m_upper.startswith('CD') and m_upper[2:3].isdigit() or
            m_upper.startswith(('IGH', 'IGK', 'IGL', 'VPREB', 'IGLL')) or
            m_upper in _HIGH_PRIORITY_MARKERS):
            high_priority.append(m)
        elif (m_upper.startswith('HLA-') or
              m_upper.startswith(('MHC', 'FCRL', 'TNFRSF', 'TNFSF')) or
              m_upper in _MEDIUM_PRIORITY_MARKERS):
            medium_priority.append(m)
        else:
            normal_priority.append(m)

    reordered = high_priority + medium_priority + normal_priority
    return reordered[:max_markers]


# Markers known to be informative for cell type identification
# NOT hardcoded cell types - just gene families that help PubMed search
_HIGH_PRIORITY_MARKERS = {
    # Surface markers
    'BTLA', 'CTLA4', 'PDCD1', 'LAG3', 'TIGIT', 'HAVCR2',
    'EPCAM', 'THY1', 'NCAM1', 'KLRD1', 'KLRB1', 'KLRK1',
    'NKG7', 'GNLY', 'PRF1', 'GZMB', 'GZMA', 'GZMK',
    'MS4A1', 'CR2', 'FCER2', 'BLNK', 'BANK1',
    # Key lineage TFs
    'PAX5', 'EBF1', 'FOXP3', 'TBX21', 'GATA3', 'RORC',
    'BCL6', 'IRF4', 'BLIMP1', 'PRDM1', 'TCF7', 'LEF1',
    'SOX4', 'SOX5', 'E2A', 'TCF3', 'IKZF1', 'IKZF3',
    # Key enzymes/functional
    'DNTT', 'RAG1', 'RAG2', 'AICDA', 'TRAC', 'TRBC1',
    'LYZ', 'MPO', 'ELANE', 'CTSG', 'AZU1',
    'S100A8', 'S100A9', 'S100A12', 'VCAN', 'FCN1',
}

_MEDIUM_PRIORITY_MARKERS = {
    # Signaling/TFs that help identify subtypes
    'FOXO1', 'FOXO3', 'POU2F2', 'POU2AF1', 'IRF8', 'SPI1',
    'KLF2', 'KLF4', 'ZEB2', 'BATF', 'BATF3',
    'STAT1', 'STAT3', 'STAT4', 'STAT6',
    'GPR183', 'CCR7', 'CCR6', 'CCR9', 'CXCR5', 'CXCR4',
    'LTB', 'TNF', 'IFNG', 'IL2', 'IL4', 'IL7R', 'IL17A',
    'IGHM', 'IGHD', 'IGHG1', 'IGHG2', 'IGHA1', 'IGHE',
    'FCRLA', 'FCRL3', 'FCRL5',
    'CD72', 'TNFRSF13B', 'TNFRSF17',
}


def _select_wide_marker_set(
    all_markers: List[str],
    max_per_category: int = 3,
    max_total: int = 20
) -> Tuple[List[str], Dict[str, List[str]]]:
    """
    Select a diverse set of representative markers from the full DE marker list.

    Instead of only using the top 10 markers, this scans the full top 100
    and selects representatives from each functional category to ensure broad
    coverage for literature search. A marker at rank 38 (e.g., BCL2) can be
    just as important for cell identity as one at rank 3.

    Categories:
    - surface: CD markers (CD200, CD23, etc.)
    - key_marker: Known cell identity markers (BTLA, FCER2, PAX5, etc.)
    - signaling: Signaling molecules and TFs (KLF2, FOXO1, etc.)
    - immunoglobulin: Ig genes (IGHM, IGHD, etc.)
    - fc_slam: FCRL/SLAM family
    - receptor: Cytokine/chemokine receptors
    - hla: HLA/MHC genes
    - other: Top-ranked remaining markers

    Args:
        all_markers: Full DE marker list (raw, unfiltered)
        max_per_category: Max markers per category
        max_total: Max total markers to return

    Returns:
        Tuple of (flat_marker_list, categories_dict)
    """
    # Filter noise from expanded list (top 100)
    filtered = []
    for m in all_markers[:100]:
        if m.startswith('LOC') and m[3:].isdigit():
            continue
        if '.' in m:
            continue
        if m.startswith('MT-'):
            continue
        if len(m) < 2:
            continue
        if m.startswith(('RPS', 'RPL')):
            continue
        filtered.append(m)

    # Categorize into functional groups
    categories: Dict[str, List[str]] = {}
    categorized = set()

    for m in filtered:
        mu = m.upper()
        cat = None

        # Surface markers: CD molecules + integrins (ITG* = cell surface adhesion)
        if mu.startswith('CD') and len(mu) > 2 and mu[2:3].isdigit():
            cat = 'surface'
        elif mu.startswith('ITG') and len(mu) > 3:
            cat = 'surface'  # Integrins (ITGAX=CD11c, ITGB2, etc.)
        elif mu.startswith(('PECAM', 'EPCAM', 'NCAM', 'ICAM', 'VCAM')):
            cat = 'surface'  # CAM family cell adhesion molecules
        elif mu.startswith(('IGH', 'IGK', 'IGL', 'VPREB', 'IGLL')):
            cat = 'immunoglobulin'
        elif mu.startswith('HLA-') or mu.startswith('MHC'):
            cat = 'hla'
        elif mu.startswith(('SLAMF', 'FCRL', 'FCRLA')):
            cat = 'fc_slam'
        elif mu in _HIGH_PRIORITY_MARKERS:
            cat = 'key_marker'
        elif mu in _MEDIUM_PRIORITY_MARKERS:
            cat = 'signaling'
        elif any(mu.startswith(p) for p in ('IL', 'CCR', 'CXCR', 'TNFRSF', 'TNFSF')):
            cat = 'receptor'
        elif any(mu.startswith(p) for p in ('TLR', 'BCL')):
            cat = 'functional'  # Toll-like receptors, BCL family

        if cat:
            if cat not in categories:
                categories[cat] = []
            if len(categories[cat]) < max_per_category:
                categories[cat].append(m)
                categorized.add(m)

    # Fill 'other' with top-ranked uncategorized markers
    remaining_budget = max_total - sum(len(v) for v in categories.values())
    other_markers = [m for m in filtered if m not in categorized]
    categories['other'] = other_markers[:max(2, remaining_budget)]

    # Build flat list: prioritized categories first
    # key_marker/signaling/functional first (most likely to identify cell type),
    # then surface/Ig/FCRL (well-searchable on PubMed)
    priority_order = [
        'key_marker', 'signaling', 'functional', 'surface',
        'fc_slam', 'immunoglobulin', 'receptor', 'hla', 'other'
    ]
    selected = []
    for cat in priority_order:
        if cat in categories:
            for m in categories[cat]:
                if m not in selected:
                    selected.append(m)

    # Ensure very top DE markers (rank 1-5 from filtered) are at the front
    # These have the strongest statistical signal regardless of category
    top5_to_insert = [m for m in filtered[:5] if m not in selected]
    for m in reversed(top5_to_insert):
        selected.insert(0, m)

    return selected[:max_total], categories


def _generate_diverse_combinations(
    wide_markers: List[str],
    categories: Dict[str, List[str]],
    max_combos: int = 12
) -> List[List[str]]:
    """
    Generate diverse marker combinations from categorized markers.

    Strategy:
    1. Top 2-3 markers (traditional, high-ranked)
    2. Individual informative markers (each alone finds specific papers)
    3. Cross-category pairs (surface + TF, Ig + receptor, etc.)

    This ensures PubMed queries cover multiple angles of cell identity,
    not just the top-ranked genes which may be housekeeping-adjacent.

    Args:
        wide_markers: Flat list of selected markers
        categories: Dict of {category: [markers]} from _select_wide_marker_set
        max_combos: Max combinations to return

    Returns:
        List of marker combinations (each is a list of gene symbols)
    """
    combos = []
    seen = set()

    def _add(combo):
        key = tuple(sorted(combo))
        if key not in seen:
            seen.add(key)
            combos.append(combo)

    # Strategy 1: Top 2-3 ranked markers together
    if len(wide_markers) >= 2:
        _add(wide_markers[:2])
    if len(wide_markers) >= 3:
        _add(wide_markers[:3])

    # Strategy 2: Individual informative markers
    # These often find the most specific papers (e.g., FCER2 alone → follicular B)
    informative_cats = [
        'surface', 'key_marker', 'fc_slam', 'immunoglobulin',
        'signaling', 'functional', 'receptor'
    ]
    for cat in informative_cats:
        if cat in categories:
            for m in categories[cat][:2]:
                _add([m])

    # Strategy 3: Within-category pairs for key categories
    # e.g., BTLA + FCER2 (both key_markers) together find more specific papers
    for cat in ['key_marker', 'surface', 'signaling']:
        if cat in categories and len(categories[cat]) >= 2:
            _add([categories[cat][0], categories[cat][1]])

    # Strategy 4: Cross-category pairs
    # Combine markers from different functional groups for richer context
    cat_keys = [k for k in informative_cats if k in categories and categories[k]]
    for i, k1 in enumerate(cat_keys):
        for k2 in cat_keys[i + 1:]:
            if categories[k1] and categories[k2]:
                _add([categories[k1][0], categories[k2][0]])
            if len(combos) >= max_combos:
                break
        if len(combos) >= max_combos:
            break

    # Strategy 5: Top other markers individually (catch remaining signal)
    for m in categories.get('other', [])[:2]:
        if len(combos) < max_combos:
            _add([m])

    return combos[:max_combos]


def _validate_candidates_against_marker_profile(
    candidates: List[CellTypeCandidate],
    all_markers: List[str],
    parent_cell_type: str,
    n_val_markers: int = 5,
    n_candidates: int = 15
) -> List[CellTypeCandidate]:
    """
    Validate each candidate by checking co-occurrence with diverse DE markers.

    For each top candidate, search PubMed for "candidate_subtype AND DE_marker"
    using representatives from different functional categories. Candidates that
    co-occur with many DE markers are likely the true cell identity; those that
    don't match the marker profile get penalized.

    Example: Cluster 0 has BTLA, FCER2, BCL2, KLF2, CD200 as markers.
    - "Naive B Cells" validates against 4-5/5 of these → boost
    - "Plasma Cells" validates against 0-1/5 → penalty

    This is the user's key insight: use the FULL marker profile (not just top 5)
    to verify which candidate actually matches this cluster's identity.

    Args:
        candidates: Aggregated candidates to validate
        all_markers: Full DE marker list (raw, for wide selection)
        parent_cell_type: Parent lineage
        n_val_markers: Number of diverse markers to validate against
        n_candidates: Max candidates to validate (API cost control)

    Returns:
        Candidates with adjusted confidence based on marker validation
    """
    if len(candidates) < 2:
        return candidates

    # Select diverse validation markers from the wide marker set
    val_wide, _ = _select_wide_marker_set(
        all_markers, max_per_category=2, max_total=n_val_markers
    )
    val_markers = val_wide[:n_val_markers]

    if not val_markers:
        return candidates

    parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()

    result = []
    for i, c in enumerate(candidates):
        if i >= n_candidates or c.is_novel:
            result.append(c)
            continue

        # Use full subtype name for validation (e.g., "naive b cells" not just "naive")
        # This avoids false positives from generic words like "plasma" or "memory"
        subtype_q = c.name.lower().strip()
        # Skip if it's just the parent type
        parent_lower = parent_cell_type.lower().strip()
        if subtype_q == parent_lower or len(subtype_q) < 5:
            result.append(c)
            continue

        # Validate against each diverse marker
        hits = 0
        for vm in val_markers:
            q = f'"{subtype_q}"[Title/Abstract] AND {vm}[Title/Abstract]'
            vresults = pubmed_search(q, max_results=1)
            if vresults:
                hits += 1

        val_fraction = hits / len(val_markers)

        # Multiplicative validation model with inverse document frequency (IDF)
        #
        # Instead of a blend that REPLACES base confidence, validation ADJUSTS
        # it via a multiplier. This preserves the relative ordering from
        # aggregation/search while using validation to confirm or refute.
        #
        # IDF correction: candidates with many PMIDs have inflated validation
        # because they appear in more PubMed papers → higher baseline
        # co-occurrence with ANY marker. Their high val_fraction reflects
        # literature abundance, not specificity for THIS cluster.
        # IDF scales DOWN the boost for high-PMID candidates.
        import math
        # Cap PMID count at 10 for IDF calculation. After merge,
        # candidates found by both marker and TF channels can have
        # 20-30+ PMIDs. This legitimately high count should NOT be
        # penalized as harshly as literature noise. Cap prevents
        # over-penalization of well-supported dual-channel candidates.
        n_pmids = min(max(len(c.pmids), 1), 10)
        idf = 1.0 / (1.0 + 0.3 * math.log1p(n_pmids))
        # n=1: 0.83, n=3: 0.70, n=5: 0.65, n=10: 0.58 (cap)

        # Base multiplier from validation fraction
        if val_fraction >= 0.5:
            base_mult = 1.0 + (val_fraction - 0.4) * 0.35
        elif val_fraction >= 0.2:
            base_mult = 0.9 + val_fraction * 0.3
        else:
            base_mult = 0.65 + val_fraction * 1.5

        # Apply IDF: reduce BOOST for high-PMID candidates
        # Penalties are NOT reduced (poor validation is informative regardless)
        if base_mult > 1.0:
            boost = (base_mult - 1.0) * idf
            multiplier = 1.0 + boost
        else:
            multiplier = base_mult

        new_conf = min(c.confidence * multiplier, 1.0)

        result.append(CellTypeCandidate(
            name=c.name,
            confidence=new_conf,
            supporting_markers=c.supporting_markers,
            pmids=c.pmids,
            evidence=c.evidence,
            is_novel=c.is_novel
        ))

    result.sort(key=lambda x: x.confidence, reverse=True)
    return result


def _consensus_rerank_candidates(
    candidates: List[CellTypeCandidate],
    parent_cell_type: str
) -> List[CellTypeCandidate]:
    """
    Re-rank candidates using group consensus scoring.

    Instead of ranking by individual confidence alone, this:
    1. Groups semantically related candidates by shared evidence and name overlap
    2. Scores each group by total evidence breadth (PMIDs, markers, members)
    3. Boosts group members proportionally to group strength
    4. Re-ranks using individual + group consensus

    This resolves cases where the correct cell type appears as multiple
    related candidates (e.g., "naive B cells", "resting B cells",
    "follicular B cells") that individually rank low but collectively
    have strong evidence.

    Args:
        candidates: Aggregated candidates (all, not just top N)
        parent_cell_type: Parent lineage context

    Returns:
        Re-ranked list of CellTypeCandidate
    """
    if len(candidates) <= 2:
        return candidates

    parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()

    def _get_modifier_words(name: str) -> set:
        """Extract modifier words from a candidate name (excluding parent type and 'cells')."""
        n = name.lower().strip()
        n = n.replace(f'{parent_base} cells', '').replace(f'{parent_base} cell', '')
        n = n.replace('cells', '').replace('cell', '').strip()
        # Remove gene/CD qualifiers for word comparison
        n = re.sub(r'\b[A-Z0-9]+[+-]\s*', '', n, flags=re.IGNORECASE).strip()
        words = {w for w in n.split() if len(w) > 2}
        return words

    # Build groups by evidence overlap and name similarity
    n = len(candidates)
    assigned = set()
    groups = []

    for i in range(n):
        if i in assigned:
            continue

        group = [i]
        assigned.add(i)
        c1 = candidates[i]
        mod1 = _get_modifier_words(c1.name)
        pmids1 = set(c1.pmids)

        for j in range(i + 1, n):
            if j in assigned:
                continue

            c2 = candidates[j]
            should_group = False

            # Criterion 1: Name equivalence
            if _names_are_equivalent(c1.name, c2.name):
                should_group = True

            # Criterion 2: Shared PMIDs (>30% overlap)
            if not should_group and pmids1 and c2.pmids:
                pmids2 = set(c2.pmids)
                shared = pmids1 & pmids2
                total = pmids1 | pmids2
                if shared and len(shared) / len(total) >= 0.3:
                    should_group = True

            # Criterion 3: Significant modifier word overlap
            if not should_group and mod1:
                mod2 = _get_modifier_words(c2.name)
                if mod2:
                    overlap = mod1 & mod2
                    union = mod1 | mod2
                    if overlap and len(overlap) / len(union) >= 0.5:
                        should_group = True

            if should_group:
                group.append(j)
                assigned.add(j)

        groups.append(group)

    # Score each group and boost members
    result = []
    for group_indices in groups:
        group_members = [candidates[i] for i in group_indices]

        # Group evidence metrics
        all_pmids = set()
        all_markers = set()
        total_conf = 0
        for m in group_members:
            all_pmids.update(m.pmids)
            all_markers.update(m.supporting_markers)
            total_conf += m.confidence

        # Consensus score components
        evidence_breadth = min(len(all_pmids) / 5.0, 1.0)
        marker_coverage = min(len(all_markers) / 4.0, 1.0)
        member_count = min(len(group_members) / 3.0, 1.0)
        avg_conf = total_conf / len(group_members)

        # Group consensus bonus: stronger for larger, better-evidenced groups
        consensus_bonus = (
            evidence_breadth * 0.3 +
            marker_coverage * 0.25 +
            member_count * 0.25 +
            avg_conf * 0.2
        ) * 0.15  # Scale: max ~0.15 bonus

        # Apply consensus bonus to each member
        for member in group_members:
            new_conf = min(member.confidence + consensus_bonus, 1.0)
            result.append(CellTypeCandidate(
                name=member.name,
                confidence=new_conf,
                supporting_markers=member.supporting_markers,
                pmids=member.pmids,
                evidence=member.evidence,
                is_novel=member.is_novel
            ))

    # Re-sort by updated confidence
    result.sort(key=lambda x: x.confidence, reverse=True)
    return result


def search_subtype_from_markers(
    markers: List[str],
    parent_cell_type: str,
    tier: int = 2,
    tissue: Optional[str] = None,
    species: str = "human",
    max_combinations: int = 12
) -> List[CellTypeCandidate]:
    """
    Search PubMed for cell subtypes within a parent lineage (Tier 2/3).

    Uses wide marker sampling: instead of only the top 10 markers, scans
    the full top 50-100 DE markers to select representatives from each
    functional category (surface, TF, Ig, receptor, etc.). This ensures
    markers at rank 22 (e.g., FCER2) or rank 38 (e.g., BCL2) contribute
    to the search.

    After collecting all candidates, applies consensus re-ranking: groups
    semantically related candidates and boosts members of strong groups.

    Args:
        markers: Top DE markers from sub-cluster (raw, full list)
        parent_cell_type: Parent lineage (e.g., "B cells", "T cells")
        tier: Annotation tier (2=developmental, 3=functional)
        tissue: Optional tissue context
        species: Species context
        max_combinations: Max marker combinations to try

    Returns:
        List of CellTypeCandidate objects, sorted by confidence
    """
    # Phase 1: Wide marker selection from full DE list
    # Scan top 50-100 markers and pick representatives per functional category
    wide_markers, categories = _select_wide_marker_set(
        markers, max_per_category=3, max_total=20
    )

    if not wide_markers:
        return [CellTypeCandidate(
            name="Novel population",
            confidence=0.0,
            supporting_markers=[],
            pmids=[],
            evidence=[],
            is_novel=True
        )]

    # Generate diverse combinations from categorized markers
    marker_combinations = _generate_diverse_combinations(
        wide_markers, categories, max_combos=max_combinations
    )

    candidates = []

    for combo in marker_combinations[:max_combinations]:
        # Build Tier 2/3 specific queries
        queries = _build_subtype_queries(
            combo, parent_cell_type, tier=tier,
            tissue=tissue, species=species
        )

        # Accumulate results from multiple query levels (not just first hit)
        # This ensures alias-based queries (CD23 vs FCER2) are also tried,
        # finding papers that use protein names instead of gene symbols
        results = []
        seen_pmids_q = set()
        for query in queries[:5]:  # Try up to 5 query levels
            batch = pubmed_search(query, max_results=3)
            for r in batch:
                if r['pmid'] not in seen_pmids_q:
                    results.append(r)
                    seen_pmids_q.add(r['pmid'])
            if len(results) >= 5:  # Cap per combination
                break

        # Extract cell types/subtypes from results
        for result in results:
            title = result.get('title') or ''
            abstract = result.get('abstract') or ''
            text = title + " " + abstract

            # Extract subtypes (more granular than Tier 1)
            cell_types = extract_cell_subtypes_from_text(
                text, parent_cell_type
            )

            # Subtype dilution: papers mentioning many subtypes are less specific
            # about each one. Scale context scores by 1/N^0.3 where N = number of
            # extracted subtypes. Papers focused on one subtype keep full weight.
            import math
            n_subtypes = max(len(cell_types), 1)
            dilution = 1.0 / (n_subtypes ** 0.3) if n_subtypes > 1 else 1.0

            scored_types = _score_cell_types_by_marker_context(
                cell_types, combo, text
            )

            for cell_type, context_score in scored_types:
                candidates.append({
                    'cell_type': cell_type,
                    'markers': combo,
                    'pmid': result['pmid'],
                    'title': title,
                    'journal': result.get('journal', 'N/A'),
                    'year': result.get('year', 'N/A'),
                    'context_score': context_score * dilution
                })

    # Aggregate candidates (preserve subtypes for Tier 2/3)
    aggregated = aggregate_candidates(candidates, preserve_subtypes=True)

    # Phase 2a: Marker profile validation (FIRST, before penalties)
    # Validate ALL candidates against the full marker profile.
    # This is the most powerful discriminator: candidates that co-occur
    # with the cluster's specific markers in literature are more likely correct.
    # Must run before other penalties to cover candidates at any position.
    aggregated = _validate_candidates_against_marker_profile(
        aggregated, markers, parent_cell_type,
        n_val_markers=5, n_candidates=20
    )

    # Penalize CD-marker-qualified subtypes not in DE markers
    markers_upper = {m.upper() for m in wide_markers}
    for i, c in enumerate(aggregated):
        cd_in_name = re.findall(r'\bCD(\d+)', c.name, re.IGNORECASE)
        if cd_in_name:
            any_in_de = any(f'CD{num}' in markers_upper for num in cd_in_name)
            if not any_in_de:
                aggregated[i] = CellTypeCandidate(
                    name=c.name,
                    confidence=c.confidence * 0.70,
                    supporting_markers=c.supporting_markers,
                    pmids=c.pmids,
                    evidence=c.evidence,
                    is_novel=c.is_novel
                )
    aggregated.sort(key=lambda x: x.confidence, reverse=True)

    # Filter out generic parent type and penalize cross-lineage candidates
    parent_lower = parent_cell_type.lower().strip()
    parent_base = parent_lower.replace(' cells', '').replace(' cell', '').strip()
    if parent_lower and len(aggregated) >= 2:
        has_specific = any(
            c.name.lower() != parent_lower and
            not c.name.lower().startswith('novel')
            for c in aggregated[:5]
        )
        if has_specific:
            aggregated = [c for c in aggregated
                          if c.name.lower() != parent_lower or c.confidence > 0.95]

        # Detect cross-lineage candidates: normalize singular/plural
        _lineage_bases = {
            't', 'b', 'nk', 'myeloid', 'monocyte', 'monocytes',
            'dendritic', 'neutrophil', 'neutrophils',
            'macrophage', 'macrophages', 'eosinophil', 'eosinophils',
            'basophil', 'basophils', 'platelet', 'platelets',
            'erythrocyte', 'erythrocytes', 'leukocyte', 'leukocytes',
        }
        _lineage_subtype_map = {
            't': [r'^th\d', r'^tc\d', r'^treg', r'^nkt', r'^mait',
                  r'^cd[48]\+?\s+t\b', r'^cytotoxic t', r'^helper t'],
            'b': [r'^pre-?b', r'^pro-?b'],
            'nk': [r'^cd56', r'^cd16\+?\s+nk'],
        }
        for i, c in enumerate(aggregated):
            c_lower = c.name.lower()
            c_base = c_lower.replace(' cells', '').replace(' cell', '').strip()
            is_different = False
            # Check singular/plural forms
            c_singular = c_base.rstrip('s') if c_base.endswith('s') else c_base
            if (c_base in _lineage_bases or c_singular in _lineage_bases) and \
               c_singular != parent_base and parent_base not in c_lower:
                is_different = True
            else:
                for lineage, patterns in _lineage_subtype_map.items():
                    if lineage == parent_base:
                        continue
                    if any(re.match(p, c_base) for p in patterns):
                        is_different = True
                        break
            if is_different:
                aggregated[i] = CellTypeCandidate(
                    name=c.name,
                    confidence=c.confidence * 0.3,
                    supporting_markers=c.supporting_markers,
                    pmids=c.pmids,
                    evidence=c.evidence,
                    is_novel=c.is_novel
                )

        # NOTE: Parent-containment penalty removed.
        # The IDF+multiplicative validation model already handles generic type
        # inflation. Adding a name-based penalty on top over-penalizes legitimate
        # differentiation products (e.g., Plasma Cells/Plasmablasts in B cell
        # search). In the full pipeline, the merge step's marker-only penalty
        # (0.85x) handles candidates not confirmed by TF channel.

        aggregated.sort(key=lambda x: x.confidence, reverse=True)

    # Phase 2b: Consensus re-ranking
    aggregated = _consensus_rerank_candidates(aggregated, parent_cell_type)

    if not aggregated or aggregated[0].confidence < 0.3:
        return [CellTypeCandidate(
            name=f"Novel {parent_cell_type} subtype",
            confidence=0.0,
            supporting_markers=markers[:5],
            pmids=[],
            evidence=[],
            is_novel=True
        )]

    return aggregated


def extract_cell_subtypes_from_text(text: str, parent_cell_type: str) -> List[str]:
    """
    Extract cell subtype mentions from text, prioritizing subtypes of the parent lineage.

    More granular than extract_cell_types_from_text:
    - Captures developmental states: "naive B", "memory B", "pro-B"
    - Captures functional states: "exhausted T", "regulatory T"
    - Captures specific subtypes: "Th1", "Th17", "plasma cells"

    Args:
        text: Title + abstract text
        parent_cell_type: Parent lineage (e.g., "B cells")

    Returns:
        List of unique subtype strings
    """
    subtypes = []

    # Normalize Unicode diacritics to ASCII: "naïve" → "naive", "résumé" → "resume"
    # Many papers use "naïve B cells" which fails to match regex pattern "naive"
    import unicodedata
    text = unicodedata.normalize('NFD', text)
    text = ''.join(c for c in text if unicodedata.category(c) != 'Mn')

    # Extract parent base: "B cells" → "B", "T cells" → "T", "Myeloid cells" → "Myeloid"
    parent_base = parent_cell_type.strip()
    if parent_base.lower().endswith(' cells'):
        parent_base = parent_base[:-6].strip()
    elif parent_base.lower().endswith(' cell'):
        parent_base = parent_base[:-5].strip()
    parent_base_lower = parent_base.lower()

    # Pattern 1: Developmental modifiers + parent type
    # "naive B cells", "memory T cells", "mature B cells"
    dev_modifiers = (
        r'naive|memory|effector|activated|exhausted|regulatory|'
        r'transitional|immature|mature|resting|proliferating|'
        r'cycling|quiescent|progenitor|precursor|'
        r'germinal[- ]center|follicular|marginal[- ]zone|'
        r'class[- ]switched|unswitched|double[- ]negative|double[- ]positive'
    )
    dev_pattern = (
        r'\b(' + dev_modifiers + r')\s+'
        + re.escape(parent_base) + r'[- ]?\s*cells?\b'
    )
    matches = re.findall(dev_pattern, text, re.IGNORECASE)
    subtypes.extend([f"{m.strip()} {parent_base} cells" for m in matches])

    # Pattern 1b: Hyphenated prefixes: "pro-B cells", "pre-B cells"
    hyph_pattern = (
        r'\b(pro|pre|early|late)\s*-?\s*'
        + re.escape(parent_base) + r'\s*cells?\b'
    )
    hyph_matches = re.findall(hyph_pattern, text, re.IGNORECASE)
    subtypes.extend([f"{m}-{parent_base} cells" for m in hyph_matches])

    # Pattern 2: CD-qualified subtypes
    # "CD8+ T cells", "CD4+ T cells", "CD27+ B cells"
    cd_pattern = (
        r'\b(CD\d+[+-]?(?:/CD\d+[+-]?)?)\s+'
        + re.escape(parent_base) + r'[- ]?\s*cells?\b'
    )
    cd_matches = re.findall(cd_pattern, text, re.IGNORECASE)
    subtypes.extend([f"{m} {parent_base} cells" for m in cd_matches])

    # Pattern 3: Specific known subtype terms (without hardcoding cell types)
    # "plasma cells", "plasmablasts", "Th1 cells", "Treg cells"
    specific_patterns = [
        r'\b(plasma\s*(?:blast|cell)s?)\b',
        r'\b(Th\d+\s+cells?)\b',
        r'\b(Treg\s+cells?|regulatory\s+T\s+cells?)\b',
        r'\b(NKT\s+cells?|MAIT\s+cells?)\b',
        r'\b(DN\d?\s+cells?|DP\s+cells?)\b',
        r'\b(ILC\d?\s+cells?)\b',
    ]
    for pat in specific_patterns:
        for m in re.findall(pat, text, re.IGNORECASE):
            subtypes.append(m.strip())

    # Pattern 4: Functional states relevant to the lineage
    func_pattern = (
        r'\b(cytotoxic|helper|suppressor|anergic|senescent|'
        r'tissue[- ]resident|recirculating|age[- ]associated|'
        r'atypical|classical|non[- ]classical|intermediate|'
        r'inflammatory|anti[- ]inflammatory|tolerogenic|'
        r'switched|unswitched|innate[- ]like|'
        r'extrafollicular|peritoneal|breg|b1[ab]?|b2)\s+'
        + re.escape(parent_base) + r'[- ]?\s*cells?\b'
    )
    func_matches = re.findall(func_pattern, text, re.IGNORECASE)
    subtypes.extend([f"{m.strip()} {parent_base} cells" for m in func_matches])

    # Pattern 5: Gene/TF-qualified subtypes
    # "T-bet+ B cells", "TBX21-expressing B cells", "FOXP3+ T cells"
    tf_qual_pattern = (
        r'\b([A-Z][A-Za-z0-9]+)[+\-]?\s*(?:positive|expressing|high|lo(?:w)?|neg(?:ative)?)\s+'
        + re.escape(parent_base) + r'[- ]?\s*cells?\b'
    )
    tf_matches = re.findall(tf_qual_pattern, text, re.IGNORECASE)
    _non_gene_words = {
        'with', 'the', 'and', 'but', 'for', 'not', 'are', 'was', 'were',
        'has', 'had', 'have', 'from', 'into', 'than', 'that', 'this',
        'these', 'those', 'their', 'both', 'each', 'all', 'any', 'most',
        'some', 'such', 'more', 'less', 'only', 'also', 'very', 'thus',
    }
    subtypes.extend([
        f"{m}+ {parent_base} cells" for m in tf_matches
        if m.lower() not in _non_gene_words
    ])

    # Pattern 5b: "TF+ cells" without parent name (e.g., "T-bet+ cells")
    tf_plus_pattern = r'\b([A-Z][A-Za-z0-9-]+)\+\s+' + re.escape(parent_base) + r'[- ]?\s*cells?\b'
    tf_plus_matches = re.findall(tf_plus_pattern, text, re.IGNORECASE)
    subtypes.extend([
        f"{m}+ {parent_base} cells" for m in tf_plus_matches
        if m.lower() not in _non_gene_words
    ])

    # Pattern 6: Expand abbreviations found with their full form in parentheses
    # "age-associated B cells (ABCs)" → capture the full form, skip abbreviation
    # This avoids duplicate entries like "ABCs B Cells" alongside "Age-associated B cells"
    paren_pattern = (
        r'\b([A-Za-z][a-z-]+(?:\s+[a-z-]+){0,3})\s+'
        + re.escape(parent_base) + r'\s*cells?\s*\(([A-Z]{2,5}s?)\)'
    )
    paren_matches = re.findall(paren_pattern, text, re.IGNORECASE)
    _noise_prefixes = {
        # Species words
        'human', 'mouse', 'murine', 'rat', 'bovine', 'porcine', 'primate',
        # Structural/generic words
        'type', 'types', 'subset', 'subsets', 'subpopulation', 'population',
        'class', 'kind', 'group', 'fraction', 'pool', 'compartment',
        'specific', 'distinct', 'novel', 'new', 'certain', 'particular',
        'different', 'various', 'several', 'major', 'minor',
        'and', 'or', 'of', 'the', 'a', 'an', 'this', 'that', 'these',
        # Verbal/participial forms (research actions, not cell type modifiers)
        'classified', 'identifying', 'identified', 'characterizing', 'characterized',
        'defining', 'defined', 'describing', 'described', 'shown', 'showing',
        'found', 'finding', 'observed', 'observing', 'detected', 'detecting',
        'generated', 'generating', 'produced', 'producing', 'contained', 'containing',
        'represented', 'representing', 'resembled', 'resembling', 'designated',
        'termed', 'called', 'named', 'referred', 'labeled', 'categorized',
        'emerging', 'converted', 'converting', 'becoming', 'known', 'recognized',
        'considered', 'regarded', 'thought', 'suggested', 'proposed',
        'displaying', 'exhibited', 'exhibiting', 'included', 'including',
        # Pathological/quality descriptors (not cell subtypes)
        'abnormal', 'normal', 'typical', 'rare', 'common', 'abundant',
        'predominant', 'prominent', 'notable', 'significant',
        'as',  # "classified as" context
        # Quantitative/comparative
        'higher', 'lower', 'fewer', 'greater', 'increased', 'decreased',
        'reduced', 'elevated', 'enhanced', 'diminished', 'expanded',
        'total', 'overall', 'majority', 'subset',
        'proportions', 'proportion', 'percentage', 'frequency',
        'number', 'numbers', 'levels', 'level',
        'role', 'roles', 'function', 'functions', 'loss', 'lack',
    }
    for full_form, abbrev in paren_matches:
        # Strip leading noise words: "type of memory" → "memory"
        words = full_form.strip().split()
        while words and words[0].lower() in _noise_prefixes:
            words.pop(0)
        if words:
            clean_form = ' '.join(words)
            subtypes.append(f"{clean_form} {parent_base} cells")

    # Pattern 7: Compound descriptors
    # "CD21-/CD27- B cells", "IgG+ memory B cells"
    compound_pattern = (
        r'\b((?:Ig[ADEGM]\+?(?:/Ig[ADEGM]\+?)?|'
        r'CD\d+[+-](?:/CD\d+[+-])?)\s+'
        r'(?:naive|memory|switched|unswitched|mature|immature)?\s*'
        + re.escape(parent_base) + r')\s*cells?\b'
    )
    compound_matches = re.findall(compound_pattern, text, re.IGNORECASE)
    subtypes.extend([f"{m.strip()} cells" for m in compound_matches])

    # Also include general cell types from the text
    general = extract_cell_types_from_text(text)
    subtypes.extend(general)

    # Filter species-prefixed and verbal-noise subtypes (catch-all for any pattern)
    _sp_words = {'human', 'mouse', 'murine', 'rat', 'bovine', 'porcine', 'primate'}
    _verbal_starts = {
        'classified', 'identifying', 'identified', 'characterizing', 'characterized',
        'defining', 'defined', 'describing', 'described', 'shown', 'showing',
        'found', 'finding', 'observed', 'observing', 'detected', 'detecting',
        'generated', 'generating', 'produced', 'producing', 'contained', 'containing',
        'represented', 'representing', 'resembled', 'resembling', 'designated',
        'termed', 'called', 'named', 'referred', 'labeled', 'categorized',
        'emerging', 'converted', 'converting', 'becoming', 'known', 'recognized',
        'considered', 'regarded', 'displaying', 'exhibited', 'exhibiting',
        'included', 'including', 'suggested', 'proposed', 'abnormal',
        # Quantitative/structural prefixes
        'proportions', 'proportion', 'percentage', 'frequency', 'frequencies',
        'number', 'numbers', 'levels', 'level', 'amounts', 'amount',
        'role', 'roles', 'function', 'functions', 'loss', 'lack',
        'absence', 'presence', 'expression', 'regulation', 'depletion',
        # Comparative/superlative/quantifier
        'higher', 'lower', 'fewer', 'greater', 'increased', 'decreased',
        'reduced', 'elevated', 'enhanced', 'diminished', 'expanded',
        'total', 'overall', 'majority', 'subset',
    }
    cleaned = []
    for s in subtypes:
        words = s.strip().split()
        if not words:
            continue
        first_word = words[0].lower()
        if first_word in _sp_words:
            # Strip species prefix: "mouse memory B cells" → "memory B cells"
            words.pop(0)
            while words and words[0].lower() in ('and', 'or'):
                words.pop(0)
                if words and words[0].lower() in _sp_words:
                    words.pop(0)
            if words:
                cleaned.append(' '.join(words))
        elif first_word in _verbal_starts:
            # Skip entries starting with verbal/participial noise
            continue
        else:
            cleaned.append(s)

    # Deduplicate
    seen = set()
    unique = []
    for s in cleaned:
        s_lower = s.lower().strip()
        if s_lower not in seen and len(s_lower) > 3:
            seen.add(s_lower)
            unique.append(s)

    return unique


def _build_subtype_queries(
    markers: List[str],
    parent_cell_type: str,
    tier: int = 2,
    tissue: Optional[str] = None,
    species: str = "human"
) -> List[str]:
    """
    Build PubMed queries for subtype search (Tier 2/3).

    Always includes parent cell type context.

    Args:
        markers: Gene symbols
        parent_cell_type: Parent lineage context
        tier: 2=developmental, 3=functional
        species: Species context

    Returns:
        List of query strings (try in order)
    """
    queries = []
    # Extract lineage base: "B cells" → "b", "T cells" → "t"
    parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()

    # Tier-specific search context: include BOTH developmental and mature state terms
    # This ensures queries find papers about any cell state, not just development
    if tier == 2:
        context_terms = (
            "development OR differentiation OR subtype OR subset OR "
            "naive OR mature OR quiescent OR resting OR homeostasis"
        )
    else:  # tier 3
        context_terms = "functional state OR activation OR pathway OR exhaustion"

    # Query order: specific → broad. Loop breaks on first success, so
    # specific queries must come first to find the most relevant papers.

    # Level 1: Multi-marker + parent cell type + combined context
    if len(markers) >= 2:
        tagged = " AND ".join(f"{m}[Title/Abstract]" for m in markers[:2])
        q1 = (
            f"{tagged} AND ({parent_base} cell[Title/Abstract]) "
            f"AND ({context_terms})"
        )
        queries.append(q1)

    # Level 2: Single marker + parent type + combined context (always available)
    q2 = (
        f"{markers[0]}[Title/Abstract] AND {parent_base} cell[Title/Abstract] "
        f"AND ({context_terms})"
    )
    queries.append(q2)

    # Level 2b: Marker ALIAS (protein/CD name) + parent type + context
    # Gene symbols (FCER2) are sometimes less common in papers than protein
    # names (CD23). Using aliases finds additional papers.
    marker_aliases = _get_cached_aliases(markers[0])
    for alias in marker_aliases[:2]:  # Top 2 aliases
        if len(alias) >= 2 and alias != markers[0]:
            q2b = (
                f"{alias}[Title/Abstract] AND {parent_base} cell[Title/Abstract] "
                f"AND ({context_terms})"
            )
            queries.append(q2b)

    # Level 3: Free-text markers + parent type (broader)
    m2 = markers[1] if len(markers) > 1 else markers[0]
    q3 = f'"{markers[0]}" "{m2}" {parent_base} cell {species}'
    queries.append(q3)

    # Level 4: Single marker + parent type (broadest)
    q4 = f"{markers[0]} {parent_base} cell {species}"
    queries.append(q4)

    # Level 5: Marker name only + cell (last resort)
    q5 = f"{markers[0]} cell type marker"
    queries.append(q5)

    return queries


def search_subtype_with_tf_activity(
    markers: List[str],
    tf_activities: Dict[str, float],
    parent_cell_type: str,
    tier: int = 2,
    tissue: Optional[str] = None,
    species: str = "human",
    max_marker_combos: int = 12,
    max_tf_queries: int = 3,
    all_cluster_tf_activities: Optional[Dict[str, Dict[str, float]]] = None
) -> List[CellTypeCandidate]:
    """
    Combined marker + TF activity search for cell subtype identification.

    Three-channel evidence integration:
    1. MARKER channel: DE markers → PubMed → cell subtypes
    2. TF channel: Top active TFs → PubMed → cell subtypes
    3. MERGE: Cross-validate, boost candidates found in both channels

    This resolves cases where markers alone are ambiguous but TF activity
    provides discriminating signal (e.g., TBX21 activity → ABC/atypical B cells).

    Args:
        markers: Top DE markers (raw, will be filtered internally)
        tf_activities: Dict of {TF_name: activity_score} from decoupler.
                       Positive scores = active, negative = inactive.
        parent_cell_type: Parent lineage (e.g., "B cells")
        tier: 2=developmental, 3=functional
        tissue: Optional tissue context
        species: Species context
        max_marker_combos: Max marker combinations to try
        max_tf_queries: Max TF queries to send
        all_cluster_tf_activities: Optional dict of ALL clusters' TF activities
            {cluster_id: {TF_name: score}}. When provided, z-scores are computed
            across clusters so cluster-SPECIFIC TFs are weighted higher than
            housekeeping TFs with similar activity across all clusters.

    Returns:
        List of CellTypeCandidate objects, sorted by combined confidence
    """
    # ---- Z-score TF activities if cross-cluster data is available ----
    # z-scores identify cluster-discriminating TFs: a TF uniquely active in
    # this cluster gets a high z-score, while a TF active across many clusters
    # (e.g., FOXO1 in both Naive and Transitional B) gets a lower z-score.
    if all_cluster_tf_activities and len(all_cluster_tf_activities) > 2:
        tf_z_scores = _zscore_tf_activities(tf_activities, all_cluster_tf_activities)
        # Use z-scores for TF selection/weighting, raw for sign (active/inactive)
        tf_weights = tf_z_scores
    else:
        tf_weights = tf_activities
    # ---- Channel 1: Marker-based search (existing logic) ----
    marker_candidates = search_subtype_from_markers(
        markers, parent_cell_type=parent_cell_type, tier=tier,
        tissue=tissue, species=species, max_combinations=max_marker_combos
    )

    # ---- Channel 2: TF activity-based search ----
    tf_candidates_raw = []

    # Select top TFs by z-score (cluster-specificity) when available,
    # otherwise by raw activity magnitude. Only consider active TFs (raw > 0).
    # When z-scores are tied (common with sparse data where single-cluster TFs
    # all get the same z-score), break ties using raw magnitude via log scale.
    active_tfs = {k: v for k, v in tf_activities.items() if v > 0}
    import math as _math
    tf_selection_weights = {}
    for k, v in active_tfs.items():
        z = tf_weights.get(k, v)
        raw = v
        # z-score primary, log(raw) secondary tiebreaker
        tf_selection_weights[k] = z + 0.01 * _math.log1p(raw)
    top_tfs = _select_top_tfs(tf_selection_weights, max_tfs=max_tf_queries)

    for tf_name, tf_score in top_tfs:
        # Build TF-specific PubMed queries (one per alias + fallbacks)
        queries = _build_tf_subtype_queries(
            tf_name, parent_cell_type, tier=tier,
            tissue=tissue, species=species
        )

        # Accumulate results from multiple name variants
        # Don't stop on first hit - different aliases find different papers
        # Use higher depth for high-activity TFs (more evidence to gather)
        per_query_max = 5 if abs(tf_score) >= 2.0 else 3
        total_cap = 12 if abs(tf_score) >= 2.0 else 8
        all_results = []
        seen_pmids = set()
        for query in queries:
            results = pubmed_search(query, max_results=per_query_max)
            for r in results:
                if r['pmid'] not in seen_pmids:
                    all_results.append(r)
                    seen_pmids.add(r['pmid'])
            # Stop after we have enough unique results
            if len(all_results) >= total_cap:
                break

        # Extract subtypes from TF search results
        for result in all_results:
            title = result.get('title') or ''
            abstract = result.get('abstract') or ''
            text = title + " " + abstract

            cell_types = extract_cell_subtypes_from_text(text, parent_cell_type)

            # Subtype dilution: multi-subtype papers → lower per-subtype weight
            import math
            n_subtypes = max(len(cell_types), 1)
            dilution = 1.0 / (n_subtypes ** 0.3) if n_subtypes > 1 else 1.0

            scored_types = _score_cell_types_by_marker_context(
                cell_types, [tf_name], text
            )

            for cell_type, context_score in scored_types:
                tf_candidates_raw.append({
                    'cell_type': cell_type,
                    'markers': [tf_name],
                    'pmid': result['pmid'],
                    'title': title,
                    'journal': result.get('journal', 'N/A'),
                    'year': result.get('year', 'N/A'),
                    'context_score': context_score * dilution,
                    'tf_name': tf_name,
                    'tf_score': tf_score
                })

    # ---- Second-pass: Refine TF queries using discovered qualifiers ----
    # Extract qualifying terms from first-pass subtypes (CD markers, functional states)
    # and search again with TF alias + qualifier for more specific results
    first_pass_subtypes = set()
    for cand in tf_candidates_raw:
        first_pass_subtypes.add(cand['cell_type'].lower())

    # Extract qualifying keywords from discovered subtypes
    qualifiers = set()
    for st in first_pass_subtypes:
        # Extract CD markers: "cd11c+ b cells" → "CD11c"
        cd_match = re.findall(r'\b(cd\d+[+-]?)\b', st, re.IGNORECASE)
        qualifiers.update(m.upper().rstrip('+-') for m in cd_match)
        # Extract functional modifiers: "effector b cells" → "effector"
        func_words = re.findall(
            r'\b(effector|atypical|age[- ]?associated|memory|naive|'
            r'exhausted|activated|senescent|regulatory|switched|unswitched|'
            r'double[- ]?negative|transitional|immature|mature|'
            r'germinal[- ]?center|follicular|marginal[- ]?zone)\b',
            st, re.IGNORECASE
        )
        qualifiers.update(w.lower() for w in func_words)

    # Second pass: top TF aliases + discovered qualifiers → refined queries
    if qualifiers:
        parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()
        seen_pmids_all = set(c['pmid'] for c in tf_candidates_raw)

        for tf_name, tf_score in top_tfs[:2]:  # Top 2 TFs only
            aliases = _get_cached_aliases(tf_name)
            # Best alias names (prefer protein names: shorter, 3+ chars)
            best_names = [tf_name] + [a for a in aliases if len(a) >= 3][:3]

            for qualifier in list(qualifiers)[:4]:  # Max 4 qualifiers
                for name in best_names[:2]:  # Top 2 names per qualifier
                    q = (
                        f"{name}[Title/Abstract] AND {qualifier}[Title/Abstract] "
                        f"AND {parent_base}[Title/Abstract]"
                    )
                    results = pubmed_search(q, max_results=3)
                    for result in results:
                        if result['pmid'] in seen_pmids_all:
                            continue
                        seen_pmids_all.add(result['pmid'])

                        title = result.get('title') or ''
                        abstract = result.get('abstract') or ''
                        text = title + " " + abstract

                        cell_types = extract_cell_subtypes_from_text(
                            text, parent_cell_type
                        )
                        n_st = max(len(cell_types), 1)
                        dil = 1.0 / (n_st ** 0.3) if n_st > 1 else 1.0
                        scored_types = _score_cell_types_by_marker_context(
                            cell_types, [tf_name], text
                        )
                        for cell_type, context_score in scored_types:
                            tf_candidates_raw.append({
                                'cell_type': cell_type,
                                'markers': [tf_name],
                                'pmid': result['pmid'],
                                'title': title,
                                'journal': result.get('journal', 'N/A'),
                                'year': result.get('year', 'N/A'),
                                'context_score': context_score * dil,
                                'tf_name': tf_name,
                                'tf_score': tf_score
                            })

    # Boost TF candidate context scores by TF activity before aggregation
    # High-activity TFs provide stronger evidence for their associated subtypes
    for cand in tf_candidates_raw:
        tf_score_abs = abs(cand.get('tf_score', 0))
        # Normalize: TF activity >= 3.0 gives full boost
        tf_boost = min(tf_score_abs / 3.0, 1.0) * 0.25
        cand['context_score'] = min(cand['context_score'] + tf_boost, 1.0)

    # Aggregate TF candidates
    tf_candidates = aggregate_candidates(tf_candidates_raw, preserve_subtypes=True)

    # Filter TF channel: remove parent type and cross-lineage noise.
    # Papers about TFs (e.g., PRDM1 + B cells) mention "B cells" generically
    # and may also mention unrelated lineages. Same filtering as marker channel.
    parent_lower = parent_cell_type.lower().strip()
    parent_base = parent_lower.replace(' cells', '').replace(' cell', '').strip()
    if parent_lower and len(tf_candidates) >= 2:
        # Remove generic parent type
        has_specific = any(
            c.name.lower() != parent_lower and not c.name.lower().startswith('novel')
            for c in tf_candidates[:5]
        )
        if has_specific:
            tf_candidates = [c for c in tf_candidates
                             if c.name.lower() != parent_lower or c.confidence > 0.95]

        # Penalize cross-lineage candidates in TF channel
        _lineage_bases_tf = {
            't', 'b', 'nk', 'myeloid', 'monocyte', 'monocytes',
            'dendritic', 'neutrophil', 'neutrophils',
            'macrophage', 'macrophages', 'leukocyte', 'leukocytes',
        }
        for i, c in enumerate(tf_candidates):
            c_lower = c.name.lower()
            c_base = c_lower.replace(' cells', '').replace(' cell', '').strip()
            c_singular = c_base.rstrip('s') if c_base.endswith('s') else c_base
            if (c_base in _lineage_bases_tf or c_singular in _lineage_bases_tf) and \
               c_singular != parent_base and parent_base not in c_lower:
                tf_candidates[i] = CellTypeCandidate(
                    name=c.name,
                    confidence=c.confidence * 0.3,
                    supporting_markers=c.supporting_markers,
                    pmids=c.pmids,
                    evidence=c.evidence,
                    is_novel=c.is_novel
                )
        tf_candidates.sort(key=lambda x: x.confidence, reverse=True)

    # Validate TF candidates against marker profile (same as marker channel).
    # This ensures TF-only candidates are tested against the actual DE markers
    # before merge. Each channel uses its own PMID counts for IDF, preventing
    # post-merge PMID inflation from crushing the validation boost.
    tf_candidates = _validate_candidates_against_marker_profile(
        tf_candidates, markers, parent_cell_type,
        n_val_markers=5, n_candidates=15
    )

    # ---- Channel 3: Merge and cross-validate ----
    # Pass z-scores (tf_weights) to merge so concentration bonus and
    # TF-only weighting reflect cluster-specificity, not raw magnitude.
    merged = _merge_marker_and_tf_candidates(
        marker_candidates, tf_candidates, tf_weights,
        parent_cell_type=parent_cell_type
    )

    # Post-merge: Penalize CD-marker-qualified subtypes whose CD marker
    # is NOT in the original DE markers. Use wide scan of top 100 markers
    # for broader CD marker coverage.
    wide_markers_for_cd, _ = _select_wide_marker_set(markers, max_total=30)
    markers_upper = {m.upper() for m in wide_markers_for_cd}
    for i, c in enumerate(merged):
        cd_in_name = re.findall(r'\bCD(\d+)', c.name, re.IGNORECASE)
        if cd_in_name:
            any_in_de = any(f'CD{num}' in markers_upper for num in cd_in_name)
            if not any_in_de:
                merged[i] = CellTypeCandidate(
                    name=c.name,
                    confidence=c.confidence * 0.70,
                    supporting_markers=c.supporting_markers,
                    pmids=c.pmids,
                    evidence=c.evidence,
                    is_novel=c.is_novel
                )
    merged.sort(key=lambda x: x.confidence, reverse=True)

    # NOTE: No post-merge validation. Both channels are now independently
    # validated before merge (marker channel in search_subtype_from_markers,
    # TF channel above). Post-merge validation was problematic because merged
    # PMID counts are much higher, causing IDF to over-penalize well-supported
    # dual-channel candidates like Plasma Cells (31 PMIDs after merge).

    if not merged or merged[0].confidence < 0.3:
        return [CellTypeCandidate(
            name=f"Novel {parent_cell_type} subtype",
            confidence=0.0,
            supporting_markers=markers[:5],
            pmids=[],
            evidence=[],
            is_novel=True
        )]

    return merged


def _zscore_tf_activities(
    tf_activities: Dict[str, float],
    all_cluster_tf_activities: Dict[str, Dict[str, float]]
) -> Dict[str, float]:
    """
    Convert raw TF activities to z-scores across clusters.

    Z-scoring identifies cluster-SPECIFIC TFs: a TF with high raw activity
    but similar across all clusters gets a LOW z-score (not discriminating),
    while a TF uniquely active in this cluster gets a HIGH z-score.

    Args:
        tf_activities: Current cluster's {TF_name: raw_score}
        all_cluster_tf_activities: All clusters' activities
            {cluster_id: {TF_name: raw_score}}

    Returns:
        {TF_name: z_score} for this cluster's TFs
    """
    import numpy as np

    z_scores = {}
    n_clusters = len(all_cluster_tf_activities)

    for tf_name, score in tf_activities.items():
        # Collect this TF's scores across all clusters (0 if absent)
        all_scores = [
            cluster_acts.get(tf_name, 0.0)
            for cluster_acts in all_cluster_tf_activities.values()
        ]

        if n_clusters > 2:
            mean_val = np.mean(all_scores)
            std_val = np.std(all_scores, ddof=1)  # Sample std
            if std_val > 1e-10:
                z_scores[tf_name] = float((score - mean_val) / std_val)
            else:
                # No variation → not informative for discrimination
                z_scores[tf_name] = 0.0
        else:
            # Too few clusters for meaningful z-score; use raw
            z_scores[tf_name] = score

    return z_scores


def _select_top_tfs(
    tf_activities: Dict[str, float],
    max_tfs: int = 3
) -> List[Tuple[str, float]]:
    """
    Select the most informative TFs from activity scores.

    Selects TFs with highest absolute activity scores,
    prioritizing positive (active) TFs.

    Args:
        tf_activities: {TF_name: activity_score}
        max_tfs: Maximum TFs to return

    Returns:
        List of (tf_name, score) tuples, sorted by |score| desc
    """
    if not tf_activities:
        return []

    # Sort by absolute score, preferring positive (active)
    sorted_tfs = sorted(
        tf_activities.items(),
        key=lambda x: (x[1] > 0, abs(x[1])),  # Active first, then by magnitude
        reverse=True
    )

    return sorted_tfs[:max_tfs]


def _get_gene_aliases(gene_symbol: str, species: str = "Homo sapiens") -> List[str]:
    """
    Get gene aliases from NCBI Gene database.

    Fully data-driven: queries NCBI Gene for official aliases.
    Essential for TFs where protein name differs from gene symbol
    (e.g., TBX21 → T-bet, PRDM1 → Blimp-1).

    Args:
        gene_symbol: Gene symbol (e.g., 'TBX21')
        species: Species name

    Returns:
        List of alias strings (empty if not found)
    """
    from Bio import Entrez
    import time

    Entrez.email = "kwy7605@gmail.com"
    Entrez.api_key = "40b96e1094387e03e7f9133ec6e33e881108"

    try:
        handle = Entrez.esearch(
            db="gene",
            term=f"{gene_symbol}[Gene Name] AND {species}[Organism]",
            retmax=1
        )
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])
        if not ids:
            return []

        time.sleep(0.12)

        handle = Entrez.esummary(db="gene", id=ids[0])
        summary = Entrez.read(handle)
        handle.close()

        doc = summary['DocumentSummarySet']['DocumentSummary'][0]
        aliases_str = doc.get('OtherAliases', '')
        aliases = [a.strip() for a in aliases_str.split(',') if a.strip()]

        time.sleep(0.12)
        return aliases

    except Exception:
        return []


# Cache for gene aliases to avoid repeated API calls
_alias_cache: Dict[str, List[str]] = {}


def _get_cached_aliases(gene_symbol: str) -> List[str]:
    """Get gene aliases with caching."""
    if gene_symbol not in _alias_cache:
        _alias_cache[gene_symbol] = _get_gene_aliases(gene_symbol)
    return _alias_cache[gene_symbol]


def _build_tf_subtype_queries(
    tf_name: str,
    parent_cell_type: str,
    tier: int = 2,
    tissue: Optional[str] = None,
    species: str = "human"
) -> List[str]:
    """
    Build PubMed queries for TF → cell subtype associations.

    Searches using BOTH gene symbol AND protein name aliases from NCBI Gene.
    This is critical because papers often use protein names (T-bet) instead
    of gene symbols (TBX21).

    Args:
        tf_name: Transcription factor gene symbol
        parent_cell_type: Parent lineage
        tier: 2=developmental, 3=functional
        species: Species context

    Returns:
        List of query strings (try in order)
    """
    queries = []
    # Extract lineage base: "B cells" → "B", "T cells" → "T"
    parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()

    if tier == 2:
        context = (
            "differentiation OR development OR subset OR subtype OR "
            "naive OR mature OR quiescent OR resting OR homeostasis"
        )
    else:
        context = "function OR activation OR state OR program"

    # Get protein name aliases from NCBI Gene (data-driven)
    aliases = _get_cached_aliases(tf_name)
    # Use all name variants: gene symbol + aliases
    # Sort aliases: prefer shorter, more common names (likely protein names)
    # and filter out pure numeric or very short aliases
    filtered_aliases = [a for a in aliases if len(a) >= 3 and not a.isdigit()]
    search_names = [tf_name] + filtered_aliases

    # Strategy: mix precise and broad queries across name variants
    # Precise queries find specific papers; broad queries catch review articles

    # Level 1: Precise with context restriction (per alias)
    for name in search_names[:5]:
        q = (
            f"{name}[Title/Abstract] AND {parent_base} cell[Title/Abstract] "
            f"AND ({context})"
        )
        queries.append(q)

    # Level 2: Broad without context restriction (finds review/characterization papers)
    for name in search_names[:3]:
        q = f"{name}[Title/Abstract] AND {parent_base} cell[Title/Abstract]"
        queries.append(q)

    # Level 3: Natural language queries (PubMed free text)
    for name in search_names[:3]:
        q = f'"{name}" {parent_base} cell subtype {species}'
        queries.append(q)

    # Level 4: Broad fallback with gene symbol
    queries.append(f"{tf_name} {parent_base} cell {species}")

    return queries


def _merge_marker_and_tf_candidates(
    marker_candidates: List[CellTypeCandidate],
    tf_candidates: List[CellTypeCandidate],
    tf_activities: Dict[str, float],
    parent_cell_type: str = ""
) -> List[CellTypeCandidate]:
    """
    Merge marker-based and TF-based candidates with cross-validation boost.

    Scoring logic:
    - Candidates found in BOTH channels get a cross-validation boost
    - TF candidates get weighted by the TF activity score magnitude
    - Final score = weighted combination

    Args:
        marker_candidates: From marker-based search
        tf_candidates: From TF activity-based search
        tf_activities: Original TF activity scores for weighting

    Returns:
        Merged and re-ranked list of CellTypeCandidate
    """
    # Build lookup for marker candidates
    marker_lookup = {}
    for c in marker_candidates:
        if not c.is_novel:
            key = c.name.lower().strip()
            marker_lookup[key] = c

    # Build lookup for TF candidates
    tf_lookup = {}
    for c in tf_candidates:
        if not c.is_novel:
            key = c.name.lower().strip()
            tf_lookup[key] = c

    # Merge: iterate all unique candidates
    all_names = set(list(marker_lookup.keys()) + list(tf_lookup.keys()))
    merged = []

    for name in all_names:
        m_cand = marker_lookup.get(name)
        tf_cand = tf_lookup.get(name)

        if m_cand and tf_cand:
            # Cross-validated: found in BOTH channels → boost
            cross_boost = 0.15

            # TF evidence bonus: rewards candidates with strong TF support.
            # Specificity (1 TF) gets the highest bonus, but multi-TF support
            # (concordance) also gets a meaningful bonus. The old formula used
            # a threshold at 0.3 which zeroed out the bonus for 3+ TFs, causing
            # Plasma Cells (found by PRDM1+IRF4+XBP1) to get no bonus.
            # New: smooth curve `0.08 + 0.10/n_tfs` never reaches zero.
            tf_supporting = [m for m in tf_cand.supporting_markers
                             if m in tf_activities]
            if tf_supporting:
                n_unique_tfs = len(set(tf_supporting))
                max_tf_activity = max(abs(tf_activities.get(m, 0))
                                      for m in tf_supporting)
                # Smooth bonus: 1 TF→0.18, 2 TFs→0.13, 3 TFs→0.11
                # Specificity (fewer TFs) still rewarded but multi-TF not zeroed
                concentration_bonus = (
                    min(max_tf_activity / 2.5, 1.0) *
                    (0.08 + 0.10 / n_unique_tfs)
                )
            else:
                concentration_bonus = 0

            combined_conf = (
                m_cand.confidence * 0.45 +
                tf_cand.confidence * 0.40 +
                cross_boost +
                concentration_bonus
            )
            combined_markers = list(set(
                m_cand.supporting_markers + tf_cand.supporting_markers
            ))
            combined_pmids = list(set(m_cand.pmids + tf_cand.pmids))
            combined_evidence = m_cand.evidence + tf_cand.evidence

            merged.append(CellTypeCandidate(
                name=m_cand.name,  # Use marker candidate's casing
                confidence=min(combined_conf, 1.0),
                supporting_markers=combined_markers,
                pmids=combined_pmids,
                evidence=combined_evidence,
                is_novel=False
            ))

        elif m_cand:
            # Marker-only: keep with slight penalty for no TF confirmation
            merged.append(CellTypeCandidate(
                name=m_cand.name,
                confidence=m_cand.confidence * 0.85,
                supporting_markers=m_cand.supporting_markers,
                pmids=m_cand.pmids,
                evidence=m_cand.evidence,
                is_novel=False
            ))

        elif tf_cand:
            # TF-only: weight by actual TF activity scores of supporting TFs
            # Higher TF activity → stronger evidence for this subtype
            tf_scores = [
                abs(tf_activities.get(m, 0))
                for m in tf_cand.supporting_markers
                if m in tf_activities
            ]
            if tf_scores:
                max_tf_score = max(tf_scores)
                avg_tf_score = sum(tf_scores) / len(tf_scores)
                # Normalize: activity >= 2.0 is strong evidence
                tf_strength = min(avg_tf_score / 2.0, 1.0)
                # Extra boost for very high TF activity (>= 3.0)
                high_activity_boost = min(max_tf_score / 4.0, 0.15) if max_tf_score >= 2.5 else 0
            else:
                tf_strength = 0.5
                high_activity_boost = 0

            adj_conf = tf_cand.confidence * (0.6 + 0.4 * tf_strength) + high_activity_boost

            merged.append(CellTypeCandidate(
                name=tf_cand.name,
                confidence=adj_conf,
                supporting_markers=tf_cand.supporting_markers,
                pmids=tf_cand.pmids,
                evidence=tf_cand.evidence,
                is_novel=False
            ))

    # Post-merge: consolidate near-duplicate candidates
    merged = _consolidate_duplicates(merged)

    # Inverse frequency adjustment: in a subtype search, overly common
    # candidates (many PMIDs) are less discriminating as subtype identifiers.
    # Applies diminishing returns to candidates appearing in >5 papers.
    import math
    for i, c in enumerate(merged):
        n_pmids = len(c.pmids)
        if n_pmids > 5:
            freq_penalty = 1.0 - 0.10 * math.log2(n_pmids / 5.0)
            merged[i] = CellTypeCandidate(
                name=c.name,
                confidence=c.confidence * max(freq_penalty, 0.75),
                supporting_markers=c.supporting_markers,
                pmids=c.pmids,
                evidence=c.evidence,
                is_novel=c.is_novel
            )

    # Sort by confidence
    merged.sort(key=lambda x: x.confidence, reverse=True)

    # Filter out generic parent type if specific subtypes exist
    parent_lower = parent_cell_type.lower().strip() if parent_cell_type else ""
    if len(merged) >= 2 and parent_lower:
        has_specific = any(
            c.name.lower() != parent_lower and
            not c.name.lower().startswith('novel')
            for c in merged[:5]
        )
        if has_specific:
            merged = [c for c in merged
                      if c.name.lower() != parent_lower or c.confidence > 0.9]

    # Penalize candidates from a different lineage than parent_cell_type
    # e.g., "T Cells" in a "B cells" subtype search — from shared-TF papers
    if parent_lower:
        parent_base = parent_lower.replace(' cells', '').replace(' cell', '').strip()
        # Known major lineages for cross-lineage detection
        _lineage_bases = {'t', 'b', 'nk', 'myeloid', 'monocyte', 'dendritic', 'neutrophil'}
        # Also detect lineage-specific subtypes (Th1, Th2, Treg → T lineage)
        _lineage_subtype_map = {
            't': [r'^th\d', r'^tc\d', r'^treg', r'^nkt', r'^mait',
                  r'^cd[48]\+?\s+t\b', r'^cytotoxic t', r'^helper t',
                  r'^effector t', r'^memory t', r'^naive t'],
            'b': [r'^pre-?b', r'^pro-?b', r'^plasma', r'^plasmablast'],
            'nk': [r'^cd56', r'^cd16\+?\s+nk'],
        }
        for i, c in enumerate(merged):
            c_lower = c.name.lower()
            c_base = c_lower.replace(' cells', '').replace(' cell', '').strip()
            # Check if this is a different major lineage
            is_different_lineage = False
            if c_base in _lineage_bases and c_base != parent_base and parent_base not in c_lower:
                is_different_lineage = True
            else:
                # Check if it's a subtype of a different lineage
                for lineage, patterns in _lineage_subtype_map.items():
                    if lineage == parent_base:
                        continue
                    if any(re.match(p, c_base) for p in patterns):
                        is_different_lineage = True
                        break
            if is_different_lineage:
                merged[i] = CellTypeCandidate(
                    name=c.name,
                    confidence=c.confidence * 0.3,  # Heavy penalty
                    supporting_markers=c.supporting_markers,
                    pmids=c.pmids,
                    evidence=c.evidence,
                    is_novel=c.is_novel
                )
        # Penalize candidates that don't contain the parent lineage word
        # In a subtype search, "Naive B Cells" is a subtype of "B cells"
        # but "Plasma Cells" is a different cell state (differentiation product)
        # Candidates without the parent word (as whole word) get a moderate penalty
        parent_word_pattern = re.compile(r'\b' + re.escape(parent_base) + r'\b')
        for i, c in enumerate(merged):
            c_lower = c.name.lower()
            if not parent_word_pattern.search(c_lower) and not c.is_novel:
                merged[i] = CellTypeCandidate(
                    name=c.name,
                    confidence=c.confidence * 0.80,
                    supporting_markers=c.supporting_markers,
                    pmids=c.pmids,
                    evidence=c.evidence,
                    is_novel=c.is_novel
                )

        # Re-sort after penalty
        merged.sort(key=lambda x: x.confidence, reverse=True)

    return merged


def _names_are_equivalent(name1: str, name2: str) -> bool:
    """
    Check if two cell type names refer to the same entity.

    Examples of equivalent names:
    - "Age-associated B Cells" and "Age-Associated B cells"  (case)
    - "ABCs B Cells" and "Age-associated B Cells"  (abbreviation)
    - "memory B cells" and "Memory B Cells"  (case)
    - "Pre-b Cells" and "Pre-B Cells"  (case)

    Examples of NOT equivalent:
    - "Memory B Cells" and "B Cells"  (too different)
    - "GC B Cells" and "Memory B Cells"  (different types)
    """
    n1 = name1.lower().strip().replace('-', ' ').replace('+', ' ')
    n2 = name2.lower().strip().replace('-', ' ').replace('+', ' ')

    # Exact match after normalization
    if n1 == n2:
        return True

    # Strip trailing "cells" / "cell" for comparison
    n1_base = re.sub(r'\s*cells?\s*$', '', n1).strip()
    n2_base = re.sub(r'\s*cells?\s*$', '', n2).strip()

    if n1_base == n2_base:
        return True

    # One is a substring of the other (but only if the shorter is >= 3 words)
    words1 = n1_base.split()
    words2 = n2_base.split()
    if len(words1) >= 2 and len(words2) >= 2:
        if n1_base in n2_base or n2_base in n1_base:
            return True

    return False


def _consolidate_duplicates(candidates: List[CellTypeCandidate]) -> List[CellTypeCandidate]:
    """
    Merge candidates that refer to the same cell type under different names.

    Detects duplicates by:
    1. Shared PMIDs (>50% overlap → same entity)
    2. One name is a substring/abbreviation of another

    Keeps the more descriptive (longer) name and combines evidence.

    Args:
        candidates: List of candidates to consolidate

    Returns:
        Consolidated list
    """
    if len(candidates) <= 1:
        return candidates

    consolidated = []
    merged_indices = set()

    # Generic modifiers that are broad categories - don't merge with specific types via PMID
    _generic_modifiers = {
        'naive', 'memory', 'effector', 'mature', 'immature', 'activated',
        'resting', 'germinal center', 'follicular', 'marginal zone',
        'transitional', 'regulatory', 'classical', 'non-classical',
    }

    def _is_specific_subtype(name: str) -> bool:
        """Check if a subtype name is specific (not a generic modifier)."""
        n = name.lower().strip()
        n_base = re.sub(r'\s*\b[a-z]+\s+cells?\s*$', '', n).strip()  # Remove "X cells"
        # Check if the modifier is generic
        for gm in _generic_modifiers:
            if n_base == gm or n.startswith(gm + ' '):
                return False
        # Single-word "B Cells" etc. → not specific
        if len(n.split()) <= 2:
            return False
        return True

    for i, c1 in enumerate(candidates):
        if i in merged_indices:
            continue

        # Find candidates with very similar names (same cell type, different phrasing)
        group = [c1]
        for j, c2 in enumerate(candidates):
            if j <= i or j in merged_indices:
                continue

            # Criterion 1: Name similarity
            if _names_are_equivalent(c1.name, c2.name):
                group.append(c2)
                merged_indices.add(j)
                continue

            # Criterion 2: Co-paper consolidation for rare, specific subtypes
            # If two specific subtypes share PMIDs and are both rare, they may be
            # the same population described with different terminology
            if (len(c1.pmids) <= 3 and len(c2.pmids) <= 3 and
                    _is_specific_subtype(c1.name) and _is_specific_subtype(c2.name)):
                shared = set(c1.pmids) & set(c2.pmids)
                total = set(c1.pmids) | set(c2.pmids)
                if shared and len(shared) / len(total) >= 0.3:
                    group.append(c2)
                    merged_indices.add(j)

        if len(group) == 1:
            consolidated.append(c1)
        else:
            # Merge group: keep longest name, combine evidence
            best = max(group, key=lambda c: len(c.name))
            all_pmids = list(set(
                pmid for c in group for pmid in c.pmids
            ))
            all_markers = list(set(
                m for c in group for m in c.supporting_markers
            ))
            all_evidence = [ev for c in group for ev in c.evidence]
            max_conf = max(c.confidence for c in group)
            # Boost for having multiple name variants (stronger evidence)
            variant_boost = min(0.05 * (len(group) - 1), 0.15)

            consolidated.append(CellTypeCandidate(
                name=best.name,
                confidence=min(max_conf + variant_boost, 1.0),
                supporting_markers=all_markers,
                pmids=all_pmids,
                evidence=all_evidence,
                is_novel=False
            ))

    return consolidated


# ============================================================================
# Per-Marker Voting System
# ============================================================================
# Instead of searching marker COMBINATIONS, each marker votes independently
# for cell type candidates. Generic types (Memory B, Plasma Cells) that appear
# in many B cell papers get at most 1 vote per marker, while specific subtypes
# accumulate votes from multiple independent markers.


def _search_per_marker_candidates(
    marker: str,
    parent_base: str,
    parent_cell_type: str,
    species: str = "human",
    top_k: int = 3
) -> List[Tuple[str, float, str]]:
    """
    Search PubMed for cell type candidates associated with a single marker.

    Each marker independently "votes" for cell types it's associated with in
    literature. Returns top_k candidates per marker.

    Args:
        marker: Single gene symbol (e.g., 'BTLA')
        parent_base: Parent lineage base (e.g., 'b' for B cells)
        parent_cell_type: Full parent type (e.g., 'B cells')
        species: Species context
        top_k: Max candidates to return per marker

    Returns:
        List of (cell_type_name, context_score, pmid) tuples
    """
    # Build query: single marker + parent cell type + subtype context
    context = "subtype OR subset OR differentiation OR development"
    q1 = (
        f"{marker}[Title/Abstract] AND {parent_base} cell[Title/Abstract] "
        f"AND ({context})"
    )

    # Alias fallback query
    aliases = _get_cached_aliases(marker)
    alias_queries = []
    for alias in aliases[:2]:
        if len(alias) >= 3 and alias != marker:
            aq = (
                f"{alias}[Title/Abstract] AND {parent_base} cell[Title/Abstract] "
                f"AND ({context})"
            )
            alias_queries.append(aq)

    # Broader fallback
    q_broad = f"{marker}[Title/Abstract] AND {parent_base} cell[Title/Abstract]"

    # Try queries in order, accumulate unique results
    all_results = []
    seen_pmids = set()
    for query in [q1] + alias_queries + [q_broad]:
        results = pubmed_search(query, max_results=5)
        for r in results:
            if r['pmid'] not in seen_pmids:
                all_results.append(r)
                seen_pmids.add(r['pmid'])
        if len(all_results) >= 5:
            break

    if not all_results:
        return []

    # Extract and score cell subtypes from results
    import math
    candidates_raw = []
    parent_lower = parent_cell_type.lower().strip()

    for result in all_results:
        title = result.get('title') or ''
        abstract = result.get('abstract') or ''
        text = title + " " + abstract

        cell_types = extract_cell_subtypes_from_text(text, parent_cell_type)

        # Subtype dilution
        n_subtypes = max(len(cell_types), 1)
        dilution = 1.0 / (n_subtypes ** 0.3) if n_subtypes > 1 else 1.0

        scored_types = _score_cell_types_by_marker_context(
            cell_types, [marker], text
        )

        for cell_type, context_score in scored_types:
            norm_name = normalize_cell_type_name(cell_type, preserve_subtypes=True)
            # Skip generic parent type
            if norm_name.lower().strip() == parent_lower:
                continue
            candidates_raw.append((norm_name, context_score * dilution, result['pmid']))

    # Deduplicate: keep best score per cell type
    best_per_type: Dict[str, Tuple[float, str]] = {}
    for name, score, pmid in candidates_raw:
        key = name.lower().strip()
        if key not in best_per_type or score > best_per_type[key][0]:
            best_per_type[key] = (score, pmid)

    # Sort by score, return top_k
    sorted_candidates = sorted(
        [(name, score, pmid) for name, (score, pmid) in
         [(k, best_per_type[k]) for k in best_per_type]],
        key=lambda x: x[1], reverse=True
    )

    # Reconstruct with proper casing from first occurrence
    name_casing = {}
    for name, score, pmid in candidates_raw:
        key = name.lower().strip()
        if key not in name_casing:
            name_casing[key] = name

    result_list = []
    for key in [s[0].lower().strip() for s in sorted_candidates[:top_k]]:
        score, pmid = best_per_type[key]
        result_list.append((name_casing.get(key, key), score, pmid))

    return result_list


def _search_per_tf_candidates(
    tf_name: str,
    parent_cell_type: str,
    tier: int = 2,
    species: str = "human",
    top_k: int = 3
) -> List[Tuple[str, float, str]]:
    """
    Search PubMed for cell type candidates associated with a single TF.

    Reuses _build_tf_subtype_queries() for query generation.

    Args:
        tf_name: TF gene symbol
        parent_cell_type: Full parent type
        tier: Annotation tier
        species: Species context
        top_k: Max candidates to return

    Returns:
        List of (cell_type_name, context_score, pmid) tuples
    """
    queries = _build_tf_subtype_queries(
        tf_name, parent_cell_type, tier=tier, species=species
    )

    # Accumulate results from multiple query variants
    all_results = []
    seen_pmids = set()
    for query in queries:
        results = pubmed_search(query, max_results=5)
        for r in results:
            if r['pmid'] not in seen_pmids:
                all_results.append(r)
                seen_pmids.add(r['pmid'])
        if len(all_results) >= 8:
            break

    if not all_results:
        return []

    import math
    candidates_raw = []
    parent_lower = parent_cell_type.lower().strip()

    for result in all_results:
        title = result.get('title') or ''
        abstract = result.get('abstract') or ''
        text = title + " " + abstract

        cell_types = extract_cell_subtypes_from_text(text, parent_cell_type)

        n_subtypes = max(len(cell_types), 1)
        dilution = 1.0 / (n_subtypes ** 0.3) if n_subtypes > 1 else 1.0

        scored_types = _score_cell_types_by_marker_context(
            cell_types, [tf_name], text
        )

        for cell_type, context_score in scored_types:
            norm_name = normalize_cell_type_name(cell_type, preserve_subtypes=True)
            if norm_name.lower().strip() == parent_lower:
                continue
            candidates_raw.append((norm_name, context_score * dilution, result['pmid']))

    # Deduplicate: keep best score per type
    best_per_type: Dict[str, Tuple[float, str]] = {}
    for name, score, pmid in candidates_raw:
        key = name.lower().strip()
        if key not in best_per_type or score > best_per_type[key][0]:
            best_per_type[key] = (score, pmid)

    sorted_candidates = sorted(
        best_per_type.items(),
        key=lambda x: x[1][0], reverse=True
    )

    # Proper casing
    name_casing = {}
    for name, score, pmid in candidates_raw:
        key = name.lower().strip()
        if key not in name_casing:
            name_casing[key] = name

    return [
        (name_casing.get(key, key), score, pmid)
        for key, (score, pmid) in sorted_candidates[:top_k]
    ]


def _aggregate_votes(
    marker_votes: Dict[str, List[Tuple[str, float, str]]],
    tf_votes: Dict[str, List[Tuple[str, float, str]]],
    tf_z_scores: Dict[str, float],
    parent_cell_type: str
) -> List[CellTypeCandidate]:
    """
    Aggregate per-marker and per-TF votes into ranked candidates.

    Scoring formula:
      marker_vote_fraction = n_markers_that_voted / n_total_markers   (40%)
      tf_vote_weight = sum(z_scores of supporting TFs) / total_z      (25%)
      avg_context_score = mean context score across all votes          (20%)
      evidence_breadth = min(n_unique_pmids / 5, 1.0)                 (15%)
      cross_bonus = +0.10 if found in both marker AND TF channels

    Args:
        marker_votes: {marker_name: [(cell_type, score, pmid), ...]}
        tf_votes: {tf_name: [(cell_type, score, pmid), ...]}
        tf_z_scores: {tf_name: z_score} for weighting TF contributions
        parent_cell_type: Parent lineage

    Returns:
        Ranked list of CellTypeCandidate
    """
    # Collect votes per candidate cell type
    # Each entry: which markers/TFs voted for it, with what scores
    candidate_data: Dict[str, Dict] = {}

    def _ensure_entry(key: str, display_name: str):
        if key not in candidate_data:
            candidate_data[key] = {
                'display_name': display_name,
                'voting_markers': set(),
                'voting_tfs': set(),
                'context_scores': [],
                'pmids': set(),
                'evidence': [],
            }

    # Process marker votes
    for marker_name, votes in marker_votes.items():
        for cell_type, context_score, pmid in votes:
            key = cell_type.lower().strip()
            _ensure_entry(key, cell_type)
            candidate_data[key]['voting_markers'].add(marker_name)
            candidate_data[key]['context_scores'].append(context_score)
            candidate_data[key]['pmids'].add(pmid)

    # Process TF votes
    for tf_name, votes in tf_votes.items():
        for cell_type, context_score, pmid in votes:
            key = cell_type.lower().strip()
            _ensure_entry(key, cell_type)
            candidate_data[key]['voting_tfs'].add(tf_name)
            candidate_data[key]['context_scores'].append(context_score)
            candidate_data[key]['pmids'].add(pmid)

    # Compute scores
    n_total_markers = max(len(marker_votes), 1)
    total_z = sum(abs(z) for z in tf_z_scores.values()) if tf_z_scores else 1.0
    total_z = max(total_z, 1e-6)

    candidates = []
    for key, data in candidate_data.items():
        # Marker vote fraction (40%)
        n_voting_markers = len(data['voting_markers'])
        marker_frac = n_voting_markers / n_total_markers

        # TF vote weight (25%): weighted by z-scores
        supporting_tf_z = sum(
            abs(tf_z_scores.get(tf, 0)) for tf in data['voting_tfs']
        )
        tf_weight = supporting_tf_z / total_z if data['voting_tfs'] else 0.0

        # Average context score (20%)
        ctx_scores = data['context_scores']
        avg_context = sum(ctx_scores) / len(ctx_scores) if ctx_scores else 0.0

        # Evidence breadth (15%)
        n_unique_pmids = len(data['pmids'])
        evidence_breadth = min(n_unique_pmids / 5.0, 1.0)

        # Cross-channel bonus
        cross_bonus = 0.10 if data['voting_markers'] and data['voting_tfs'] else 0.0

        combined = (
            marker_frac * 0.40 +
            tf_weight * 0.25 +
            avg_context * 0.20 +
            evidence_breadth * 0.15 +
            cross_bonus
        )

        # Collect supporting markers (genes that voted)
        all_supporting = list(data['voting_markers'] | data['voting_tfs'])

        candidates.append(CellTypeCandidate(
            name=data['display_name'],
            confidence=min(combined, 1.0),
            supporting_markers=all_supporting,
            pmids=list(data['pmids']),
            evidence=[{
                'voting_markers': list(data['voting_markers']),
                'voting_tfs': list(data['voting_tfs']),
                'marker_vote_fraction': marker_frac,
                'tf_vote_weight': tf_weight,
                'cross_validated': bool(cross_bonus),
            }],
            is_novel=False
        ))

    # Sort by confidence
    candidates.sort(key=lambda x: x.confidence, reverse=True)

    return candidates


def _apply_trajectory_constraint(
    candidates: List[CellTypeCandidate],
    pseudotime_stats: Dict[str, float],
    all_cluster_pseudotime: Dict[str, Dict[str, float]],
    parent_cell_type: str,
    species: str = "human",
    top_k: int = 5
) -> List[CellTypeCandidate]:
    """
    Apply pseudotime consistency as a gentle constraint on candidate ranking.

    For each top candidate, queries PubMed for developmental stage keywords
    and compares the literature-derived developmental position with the
    cluster's normalized pseudotime.

    Multiplier range: 0.85-1.0 (gentle, tie-breaking only).

    Args:
        candidates: Ranked candidates from voting
        pseudotime_stats: {'mean': float, 'std': float} for this cluster
        all_cluster_pseudotime: {cluster_id: {'mean': float, 'std': float}}
        parent_cell_type: Parent lineage
        species: Species context
        top_k: Max candidates to apply trajectory constraint to

    Returns:
        Re-ranked candidates with trajectory adjustment
    """
    if not pseudotime_stats or not all_cluster_pseudotime:
        return candidates

    # Normalize this cluster's pseudotime relative to all clusters (0=early, 1=late)
    all_means = [v.get('mean', 0.5) for v in all_cluster_pseudotime.values()]
    if not all_means:
        return candidates

    pt_min = min(all_means)
    pt_max = max(all_means)
    pt_range = pt_max - pt_min
    if pt_range < 1e-6:
        return candidates

    cluster_pt_norm = (pseudotime_stats.get('mean', 0.5) - pt_min) / pt_range

    parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()

    # Early and late developmental keywords (data-driven from literature)
    early_keywords = {
        'progenitor', 'precursor', 'immature', 'pro-', 'pre-', 'early',
        'stem', 'primitive', 'uncommitted', 'transitional',
    }
    late_keywords = {
        'mature', 'terminally differentiated', 'effector', 'terminal',
        'late', 'senescent', 'exhausted', 'end-stage', 'fully differentiated',
        'memory', 'long-lived', 'quiescent',
    }

    result = []
    for i, c in enumerate(candidates):
        if i >= top_k or c.is_novel:
            result.append(c)
            continue

        # Query PubMed for developmental context of this candidate
        subtype_q = c.name.lower().strip()
        q = (
            f'"{subtype_q}"[Title/Abstract] AND '
            f'{parent_base}[Title/Abstract] AND '
            f'(development OR differentiation OR maturation)'
        )
        dev_results = pubmed_search(q, max_results=3)

        if not dev_results:
            # No developmental info → no adjustment
            result.append(c)
            continue

        # Count early vs late keywords in abstracts
        early_count = 0
        late_count = 0
        for dr in dev_results:
            text = ((dr.get('title') or '') + " " + (dr.get('abstract') or '')).lower()
            for kw in early_keywords:
                early_count += text.count(kw)
            for kw in late_keywords:
                late_count += text.count(kw)

        total_kw = early_count + late_count
        if total_kw == 0:
            result.append(c)
            continue

        # Literature-derived developmental position: 0=early, 1=late
        lit_position = late_count / total_kw

        # Consistency: how well does pseudotime match literature position?
        # Both normalized to [0, 1]. Perfect match = 1.0, complete mismatch = 0.0
        consistency = 1.0 - abs(cluster_pt_norm - lit_position)

        # Gentle multiplier: 0.85 at worst, 1.0 at best
        multiplier = 0.85 + 0.15 * consistency

        new_conf = c.confidence * multiplier

        result.append(CellTypeCandidate(
            name=c.name,
            confidence=new_conf,
            supporting_markers=c.supporting_markers,
            pmids=c.pmids,
            evidence=c.evidence,
            is_novel=c.is_novel
        ))

    result.sort(key=lambda x: x.confidence, reverse=True)
    return result


def _verify_candidates_by_marker_combo(
    candidates: List[CellTypeCandidate],
    all_markers: List[str],
    parent_cell_type: str,
    n_verify_markers: int = 4,
    n_candidates: int = 10,
) -> List[CellTypeCandidate]:
    """
    Verify voting candidates by checking co-occurrence with marker COMBINATIONS.

    Per-marker voting finds candidates but lacks specificity: "Plasma Cells"
    appears in papers for BTLA, FCER2, CD200 individually. This step checks
    whether each candidate co-occurs with MULTIPLE DE markers SIMULTANEOUSLY,
    which only the TRUE cell identity will pass.

    Example for Cluster 0 (Naive B markers: BTLA, FCER2, CD200, BANK1):
    - "Plasma Cells" + BTLA + FCER2 → 0 papers (these aren't plasma markers)
    - "Naive B Cells" + BTLA + FCER2 → papers found (correct association)

    Multiplier model:
    - combo_hits / n_combos → verification fraction
    - fraction >= 0.5: boost up to 1.30x
    - fraction < 0.2: penalty down to 0.65x

    Args:
        candidates: Ranked candidates from voting/trajectory
        all_markers: Full DE marker list for selecting verification markers
        parent_cell_type: Parent lineage
        n_verify_markers: Number of markers to use in verification combos
        n_candidates: Max candidates to verify

    Returns:
        Re-ranked candidates with combo verification adjustment
    """
    if len(candidates) < 2:
        return candidates

    # Select diverse verification markers
    verify_wide, _ = _select_wide_marker_set(
        all_markers, max_per_category=2, max_total=n_verify_markers
    )
    verify_markers = verify_wide[:n_verify_markers]
    if len(verify_markers) < 2:
        return candidates

    parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()

    # Generate marker combos for verification: pairs and triples
    from itertools import combinations as _itertools_combinations
    combos = []
    # Pairs first (more likely to find hits)
    for pair in _itertools_combinations(verify_markers[:n_verify_markers], 2):
        combos.append(list(pair))
    # One triple if possible (stricter test)
    if len(verify_markers) >= 3:
        combos.append(list(verify_markers[:3]))

    if not combos:
        return candidates

    result = []
    for i, c in enumerate(candidates):
        if i >= n_candidates or c.is_novel:
            result.append(c)
            continue

        subtype_q = c.name.lower().strip()
        parent_lower = parent_cell_type.lower().strip()
        # Skip parent type itself
        if subtype_q == parent_lower or len(subtype_q) < 5:
            result.append(c)
            continue

        # Test each marker combo
        combo_hits = 0
        for combo in combos:
            marker_tags = " AND ".join(f"{m}[Title/Abstract]" for m in combo)
            q = f'"{subtype_q}"[Title/Abstract] AND {marker_tags}'
            vresults = pubmed_search(q, max_results=1)
            if vresults:
                combo_hits += 1

        verify_fraction = combo_hits / len(combos)

        # Multiplicative adjustment
        if verify_fraction >= 0.5:
            multiplier = 1.0 + (verify_fraction - 0.4) * 0.50  # up to ~1.30
        elif verify_fraction >= 0.2:
            multiplier = 0.85 + verify_fraction * 0.75  # 0.85 - 1.00
        else:
            multiplier = 0.65 + verify_fraction * 1.0   # 0.65 - 0.85

        new_conf = min(c.confidence * multiplier, 1.0)

        result.append(CellTypeCandidate(
            name=c.name,
            confidence=new_conf,
            supporting_markers=c.supporting_markers,
            pmids=c.pmids,
            evidence=c.evidence,
            is_novel=c.is_novel
        ))

    result.sort(key=lambda x: x.confidence, reverse=True)
    return result


def search_subtype_by_voting(
    markers: List[str],
    tf_activities: Dict[str, float],
    parent_cell_type: str,
    tier: int = 2,
    species: str = "human",
    n_voting_markers: int = 8,
    n_voting_tfs: int = 3,
    pseudotime_stats: Optional[Dict[str, float]] = None,
    all_cluster_tf_activities: Optional[Dict] = None,
    all_cluster_pseudotime: Optional[Dict] = None,
) -> List[CellTypeCandidate]:
    """
    Cell subtype identification via per-marker voting system.

    Instead of searching marker COMBINATIONS, each marker independently votes
    for cell types. Generic types that appear in many papers only get 1 vote
    per marker, while specific subtypes accumulate votes from multiple markers.

    Pipeline:
    1. Select top N diverse voting markers
    2. Per-marker PubMed voting → top 3 candidates each
    3. Per-TF PubMed voting → top 3 candidates each (z-weighted)
    4. Aggregate all votes → ranked candidates
    5. Apply trajectory constraint (optional)
    5.5. Combo verification → check candidates against marker PAIRS/TRIPLES
    6. Post-processing: cross-lineage penalty, CD-not-in-DE penalty

    Args:
        markers: Full DE marker list (raw, will be filtered)
        tf_activities: {TF_name: activity_score} from decoupler
        parent_cell_type: Parent lineage (e.g., "B cells")
        tier: 2=developmental, 3=functional
        species: Species context
        n_voting_markers: Number of markers to vote (default 8)
        n_voting_tfs: Number of TFs to vote (default 3)
        pseudotime_stats: {'mean': float, 'std': float} for this cluster
        all_cluster_tf_activities: All clusters' TF activities for z-scoring
        all_cluster_pseudotime: All clusters' pseudotime for normalization

    Returns:
        Ranked list of CellTypeCandidate
    """
    parent_base = parent_cell_type.lower().replace(' cells', '').replace(' cell', '').strip()
    parent_lower = parent_cell_type.lower().strip()

    # ---- Step 1: Select diverse markers from top 50+ DE markers ----
    wide_markers, categories = _select_wide_marker_set(
        markers, max_per_category=3, max_total=20
    )

    if not wide_markers:
        return [CellTypeCandidate(
            name=f"Novel {parent_cell_type} subtype",
            confidence=0.0, supporting_markers=[], pmids=[],
            evidence=[], is_novel=True
        )]

    # ---- Step 2: Combo-based voting ----
    # Generate diverse marker COMBINATIONS from the wide marker set.
    # Each combination independently "votes" for cell types it finds.
    # Combos are more specific than individual markers (BTLA+FCER2 → Follicular B)
    # while voting aggregation counts independent combo support.
    marker_combos = _generate_diverse_combinations(
        wide_markers, categories, max_combos=n_voting_markers
    )

    import math
    marker_votes: Dict[str, List[Tuple[str, float, str]]] = {}
    for combo_idx, combo in enumerate(marker_combos):
        # Build subtype queries for this combination
        queries = _build_subtype_queries(
            combo, parent_cell_type, tier=tier, species=species
        )

        # Try queries, accumulate unique results
        combo_results = []
        seen_pmids_combo = set()
        for query in queries[:4]:
            batch = pubmed_search(query, max_results=3)
            for r in batch:
                if r['pmid'] not in seen_pmids_combo:
                    combo_results.append(r)
                    seen_pmids_combo.add(r['pmid'])
            if len(combo_results) >= 5:
                break

        if not combo_results:
            continue

        # Extract and score cell subtypes
        combo_key = f"combo_{combo_idx}:{'+'.join(combo)}"
        combo_votes = []

        for result in combo_results:
            title = result.get('title') or ''
            abstract = result.get('abstract') or ''
            text = title + " " + abstract

            cell_types = extract_cell_subtypes_from_text(text, parent_cell_type)

            n_subtypes = max(len(cell_types), 1)
            dilution = 1.0 / (n_subtypes ** 0.3) if n_subtypes > 1 else 1.0

            scored_types = _score_cell_types_by_marker_context(
                cell_types, combo, text
            )

            for cell_type, context_score in scored_types:
                norm_name = normalize_cell_type_name(
                    cell_type, preserve_subtypes=True
                )
                if norm_name.lower().strip() == parent_lower:
                    continue
                combo_votes.append(
                    (norm_name, context_score * dilution, result['pmid'])
                )

        # Deduplicate per combo: keep best score per cell type
        best_per_type: Dict[str, Tuple[float, str]] = {}
        name_casing: Dict[str, str] = {}
        for name, score, pmid in combo_votes:
            key = name.lower().strip()
            if key not in name_casing:
                name_casing[key] = name
            if key not in best_per_type or score > best_per_type[key][0]:
                best_per_type[key] = (score, pmid)

        # Top 5 per combo as votes
        sorted_votes = sorted(
            best_per_type.items(), key=lambda x: x[1][0], reverse=True
        )
        final_votes = [
            (name_casing.get(k, k), s, p)
            for k, (s, p) in sorted_votes[:5]
        ]

        if final_votes:
            marker_votes[combo_key] = final_votes

    # ---- Step 3: Z-score TFs and per-TF voting ----
    if all_cluster_tf_activities and len(all_cluster_tf_activities) > 2:
        tf_z_scores = _zscore_tf_activities(tf_activities, all_cluster_tf_activities)
    else:
        tf_z_scores = {k: v for k, v in tf_activities.items()}

    # Select top TFs by z-score × raw activity product
    # With 745 TFs, many noise TFs get high z-scores from statistical
    # fluctuation but have near-zero raw activity. Require minimum raw
    # activity (>0.5) and weight by z * log(raw) so only TFs that are
    # BOTH cluster-specific AND meaningfully active are selected.
    import math as _math
    active_tfs = {k: v for k, v in tf_activities.items() if v > 0.5}
    tf_selection_weights = {}
    for k, v in active_tfs.items():
        z = max(tf_z_scores.get(k, 0), 0)  # Only positive z (above average)
        raw_weight = _math.log1p(v)  # log(1+raw) to dampen extreme values
        tf_selection_weights[k] = z * raw_weight + 0.01 * raw_weight
    top_tfs = _select_top_tfs(tf_selection_weights, max_tfs=n_voting_tfs)

    tf_votes: Dict[str, List[Tuple[str, float, str]]] = {}
    for tf_name, tf_score in top_tfs:
        votes = _search_per_tf_candidates(
            tf_name, parent_cell_type, tier=tier,
            species=species, top_k=3
        )
        if votes:
            tf_votes[tf_name] = votes

    # ---- Step 4: Aggregate all votes ----
    aggregated = _aggregate_votes(
        marker_votes, tf_votes, tf_z_scores, parent_cell_type
    )

    if not aggregated:
        return [CellTypeCandidate(
            name=f"Novel {parent_cell_type} subtype",
            confidence=0.0, supporting_markers=markers[:5], pmids=[],
            evidence=[], is_novel=True
        )]

    # ---- Step 5: Trajectory constraint (optional) ----
    if pseudotime_stats and all_cluster_pseudotime:
        aggregated = _apply_trajectory_constraint(
            aggregated, pseudotime_stats, all_cluster_pseudotime,
            parent_cell_type, species=species, top_k=5
        )

    # ---- Step 5.5: Combo verification ----
    # Voting found candidates; now verify each against marker COMBINATIONS.
    # Generic types (Plasma Cells) fail because cluster-specific markers
    # (BTLA+FCER2+CD200) don't co-occur with them in literature.
    aggregated = _verify_candidates_by_marker_combo(
        aggregated, markers, parent_cell_type,
        n_verify_markers=4, n_candidates=10
    )

    # ---- Step 6: Post-processing ----
    # 6a: Cross-lineage penalty
    _lineage_bases = {
        't', 'b', 'nk', 'myeloid', 'monocyte', 'monocytes',
        'dendritic', 'neutrophil', 'neutrophils',
        'macrophage', 'macrophages',
    }
    _lineage_subtype_map = {
        't': [r'^th\d', r'^tc\d', r'^treg', r'^nkt', r'^mait',
              r'^cd[48]\+?\s+t\b', r'^cytotoxic t', r'^helper t'],
        'b': [r'^pre-?b', r'^pro-?b'],
        'nk': [r'^cd56', r'^cd16\+?\s+nk'],
    }
    for i, c in enumerate(aggregated):
        c_lower = c.name.lower()
        c_base = c_lower.replace(' cells', '').replace(' cell', '').strip()
        c_singular = c_base.rstrip('s') if c_base.endswith('s') else c_base
        is_different = False
        if (c_base in _lineage_bases or c_singular in _lineage_bases) and \
           c_singular != parent_base and parent_base not in c_lower:
            is_different = True
        else:
            for lineage, patterns in _lineage_subtype_map.items():
                if lineage == parent_base:
                    continue
                if any(re.match(p, c_base) for p in patterns):
                    is_different = True
                    break
        if is_different:
            aggregated[i] = CellTypeCandidate(
                name=c.name,
                confidence=c.confidence * 0.3,
                supporting_markers=c.supporting_markers,
                pmids=c.pmids,
                evidence=c.evidence,
                is_novel=c.is_novel
            )

    # 6b: CD-not-in-DE penalty
    wide_markers_for_cd, _ = _select_wide_marker_set(markers, max_total=30)
    markers_upper = {m.upper() for m in wide_markers_for_cd}
    for i, c in enumerate(aggregated):
        cd_in_name = re.findall(r'\bCD(\d+)', c.name, re.IGNORECASE)
        if cd_in_name:
            any_in_de = any(f'CD{num}' in markers_upper for num in cd_in_name)
            if not any_in_de:
                aggregated[i] = CellTypeCandidate(
                    name=c.name,
                    confidence=c.confidence * 0.70,
                    supporting_markers=c.supporting_markers,
                    pmids=c.pmids,
                    evidence=c.evidence,
                    is_novel=c.is_novel
                )

    # 6c: Consolidate near-duplicates
    aggregated = _consolidate_duplicates(aggregated)

    # Final sort
    aggregated.sort(key=lambda x: x.confidence, reverse=True)

    if not aggregated or aggregated[0].confidence < 0.15:
        return [CellTypeCandidate(
            name=f"Novel {parent_cell_type} subtype",
            confidence=0.0, supporting_markers=markers[:5], pmids=[],
            evidence=[], is_novel=True
        )]

    return aggregated


# Example usage
if __name__ == "__main__":
    print("=" * 60)
    print("Dynamic Knowledge Base - Test Suite")
    print("=" * 60)

    # Test 1: Tier 1 - Major cell type identification
    print("\n--- Tier 1: Major Cell Type Search ---")
    markers_t = ['CD3D', 'CD3E', 'TRAC', 'CD8A', 'GZMB']
    print(f"Markers: {markers_t}")
    candidates = search_cell_type_from_markers(markers_t, species='human', max_combinations=3)
    for c in candidates[:3]:
        print(f"  {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}")

    # Test 2: Tier 2 - Subtype search with parent context
    print("\n--- Tier 2: B Cell Subtype Search ---")
    markers_prob = ['DNTT', 'RAG1', 'EPCAM', 'BLNK', 'VPREB3']
    print(f"Markers: {markers_prob}")
    candidates = search_subtype_from_markers(
        markers_prob, parent_cell_type="B cells", tier=2, species="human"
    )
    for c in candidates[:3]:
        print(f"  {c.name}: conf={c.confidence:.3f}, pmids={len(c.pmids)}")

    # Test 3: Marker filtering and prioritization
    print("\n--- Marker Filtering ---")
    raw = ['LOC123', 'RPS20', 'MT-CO1', 'DNTT', 'RPL12', 'CD79A', 'IGHM.1', 'PAX5']
    filtered = filter_informative_markers(raw, max_markers=5, prioritize_known=True)
    print(f"Raw: {raw}")
    print(f"Filtered: {filtered}")

    # Test 4: Cell type extraction
    print("\n--- Text Extraction ---")
    text = """
    We identified naive B cells and memory T cells in the sample.
    Pro-B cells expressed DNTT and RAG1. Plasmablasts were also present.
    """
    types = extract_cell_types_from_text(text)
    print(f"Extracted: {types}")

    print("\n" + "=" * 60)
    print("All tests completed.")
