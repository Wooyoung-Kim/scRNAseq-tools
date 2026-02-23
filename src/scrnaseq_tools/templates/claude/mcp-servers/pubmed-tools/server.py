#!/usr/bin/env python3
"""
PubMed MCP Tools for Hierarchical Annotation

Tools:
- pubmed_search: Search PubMed for marker/cell type references
- verify_reference: Verify if PMID supports marker-cell type link
- fetch_abstract: Get full abstract for a PMID
"""

import asyncio
import json
from typing import Optional
from xml.etree import ElementTree as ET
import httpx
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

# NCBI E-utilities base URL
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI API credentials (increases rate limit from 3 to 10 requests/second)
# Register at: https://www.ncbi.nlm.nih.gov/account/settings/
NCBI_API_KEY = "40b96e1094387e03e7f9133ec6e33e881108"
NCBI_EMAIL = "kwy7605@gmail.com"

server = Server("pubmed-tools")


@server.list_tools()
async def list_tools() -> list[Tool]:
    """List available tools."""
    return [
        Tool(
            name="pubmed_search",
            description="Search PubMed for marker/cell type references. Returns PMIDs with titles.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query (e.g., 'GZMB PRF1 cytotoxic T cell')"
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Maximum number of results (default: 5)",
                        "default": 5
                    }
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="verify_reference",
            description="Verify if a PMID mentions specific markers in context of a cell type.",
            inputSchema={
                "type": "object",
                "properties": {
                    "pmid": {
                        "type": "string",
                        "description": "PubMed ID to verify"
                    },
                    "markers": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of markers to check for"
                    },
                    "cell_type": {
                        "type": "string",
                        "description": "Expected cell type context"
                    }
                },
                "required": ["pmid", "markers"]
            }
        ),
        Tool(
            name="fetch_abstract",
            description="Fetch the full abstract for a PMID.",
            inputSchema={
                "type": "object",
                "properties": {
                    "pmid": {
                        "type": "string",
                        "description": "PubMed ID"
                    }
                },
                "required": ["pmid"]
            }
        )
    ]


@server.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    """Handle tool calls."""

    if name == "pubmed_search":
        result = await pubmed_search(
            arguments["query"],
            arguments.get("max_results", 5)
        )
    elif name == "verify_reference":
        result = await verify_reference(
            arguments["pmid"],
            arguments["markers"],
            arguments.get("cell_type")
        )
    elif name == "fetch_abstract":
        result = await fetch_abstract(arguments["pmid"])
    else:
        result = {"error": f"Unknown tool: {name}"}

    return [TextContent(type="text", text=json.dumps(result, indent=2))]


async def pubmed_search(query: str, max_results: int = 5) -> dict:
    """
    Search PubMed and return results with titles.

    Args:
        query: Search query (e.g., "GZMB PRF1 cytotoxic T cell")
        max_results: Maximum number of results

    Returns:
        dict with PMIDs, titles, and search metadata
    """
    async with httpx.AsyncClient(timeout=30.0) as client:
        # Step 1: Search for PMIDs
        search_params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "sort": "relevance",
            "email": NCBI_EMAIL
        }
        # Only add API key if valid (prevents 400 errors from invalid keys)
        if NCBI_API_KEY:
            search_params["api_key"] = NCBI_API_KEY

        search_resp = await client.get(f"{EUTILS_BASE}/esearch.fcgi", params=search_params)
        search_data = search_resp.json()

        pmids = search_data.get("esearchresult", {}).get("idlist", [])

        if not pmids:
            return {
                "query": query,
                "count": 0,
                "results": [],
                "message": "No results found"
            }

        # Step 2: Fetch summaries for PMIDs
        summary_params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json",
            "email": NCBI_EMAIL
        }
        if NCBI_API_KEY:
            summary_params["api_key"] = NCBI_API_KEY

        summary_resp = await client.get(f"{EUTILS_BASE}/esummary.fcgi", params=summary_params)
        summary_data = summary_resp.json()

        results = []
        for pmid in pmids:
            if pmid in summary_data.get("result", {}):
                article = summary_data["result"][pmid]
                results.append({
                    "pmid": pmid,
                    "title": article.get("title", "N/A"),
                    "authors": article.get("authors", [])[:3],  # First 3 authors
                    "journal": article.get("source", "N/A"),
                    "year": article.get("pubdate", "N/A")[:4],
                    "doi": article.get("elocationid", "N/A")
                })

        return {
            "query": query,
            "count": len(results),
            "total_found": search_data.get("esearchresult", {}).get("count", 0),
            "results": results
        }


async def verify_reference(pmid: str, markers: list[str], cell_type: Optional[str] = None) -> dict:
    """
    Verify if a PMID mentions specific markers.

    Args:
        pmid: PubMed ID
        markers: List of gene/marker names to check
        cell_type: Optional cell type context to check

    Returns:
        Verification result with found/missing markers
    """
    # Fetch abstract
    abstract_data = await fetch_abstract(pmid)

    if "error" in abstract_data:
        return {
            "pmid": pmid,
            "verified": False,
            "error": abstract_data["error"]
        }

    # Combine title and abstract for searching
    text = f"{abstract_data.get('title', '')} {abstract_data.get('abstract', '')}".lower()

    # Check markers
    found_markers = []
    missing_markers = []

    for marker in markers:
        # Check various forms: GZMB, Granzyme B, granzyme-b
        marker_lower = marker.lower()
        if marker_lower in text or marker in text:
            found_markers.append(marker)
        else:
            missing_markers.append(marker)

    # Check cell type if provided
    cell_type_found = True
    if cell_type:
        cell_type_lower = cell_type.lower()
        cell_type_found = cell_type_lower in text

    # Determine verification status
    marker_ratio = len(found_markers) / len(markers) if markers else 0

    if marker_ratio >= 0.5 and cell_type_found:
        status = "VERIFIED"
    elif marker_ratio >= 0.3:
        status = "PARTIALLY_VERIFIED"
    else:
        status = "NOT_VERIFIED"

    return {
        "pmid": pmid,
        "title": abstract_data.get("title", "N/A"),
        "status": status,
        "markers_checked": markers,
        "found_markers": found_markers,
        "missing_markers": missing_markers,
        "marker_ratio": round(marker_ratio, 2),
        "cell_type_checked": cell_type,
        "cell_type_found": cell_type_found if cell_type else None
    }


def parse_pubmed_xml(xml_text: str) -> dict:
    """
    Parse PubMed XML safely using ElementTree.

    Args:
        xml_text: Raw XML response from PubMed efetch

    Returns:
        dict with title, abstract, journal, year or error
    """
    try:
        root = ET.fromstring(xml_text)

        # Extract title
        title_elem = root.find(".//ArticleTitle")
        title = ""
        if title_elem is not None:
            # Handle mixed content (text + child elements)
            title = "".join(title_elem.itertext()).strip()

        # Extract abstract (handle multiple AbstractText elements)
        abstract_parts = []
        for abstract_text in root.findall(".//AbstractText"):
            # Get label attribute if exists (e.g., "BACKGROUND", "METHODS")
            label = abstract_text.get("Label", "")
            # Handle mixed content within AbstractText
            text = "".join(abstract_text.itertext()).strip()
            if text:
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
        abstract = " ".join(abstract_parts)

        # Extract journal title
        journal_elem = root.find(".//Journal/Title")
        journal = ""
        if journal_elem is not None and journal_elem.text:
            journal = journal_elem.text.strip()

        # Extract year (try multiple paths)
        year = ""
        year_elem = root.find(".//PubDate/Year")
        if year_elem is not None and year_elem.text:
            year = year_elem.text.strip()
        else:
            # Try MedlineDate as fallback
            medline_date = root.find(".//PubDate/MedlineDate")
            if medline_date is not None and medline_date.text:
                # Extract year from formats like "2023 Jan-Feb"
                year = medline_date.text.strip()[:4]

        return {
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year
        }

    except ET.ParseError as e:
        return {"error": f"XML parsing failed: {e}"}


async def fetch_abstract(pmid: str) -> dict:
    """
    Fetch full abstract for a PMID.

    Args:
        pmid: PubMed ID

    Returns:
        dict with title, abstract, and metadata
    """
    async with httpx.AsyncClient(timeout=30.0) as client:
        params = {
            "db": "pubmed",
            "id": pmid,
            "retmode": "xml",
            "email": NCBI_EMAIL
        }
        # Only add API key if valid
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        resp = await client.get(f"{EUTILS_BASE}/efetch.fcgi", params=params)

        if resp.status_code != 200:
            return {"error": f"Failed to fetch PMID {pmid}: HTTP {resp.status_code}"}

        # Parse XML using ElementTree
        parsed = parse_pubmed_xml(resp.text)

        if "error" in parsed:
            return {"pmid": pmid, **parsed}

        return {
            "pmid": pmid,
            "title": parsed["title"],
            "abstract": parsed["abstract"],
            "journal": parsed["journal"],
            "year": parsed["year"]
        }


async def main():
    """Run the MCP server."""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, server.create_initialization_options())


if __name__ == "__main__":
    asyncio.run(main())
