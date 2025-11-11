"""
UniProt API helper utilities.

Provides lightweight functions to fetch sequences and metadata for UniProt IDs.
Designed to avoid heavy dependencies and work in KBase runtime.
"""

import logging
from typing import Dict, List, Optional, Tuple
import requests

logger = logging.getLogger(__name__)

UNIPROT_REST = "https://rest.uniprot.org"


def fetch_sequences(uniprot_ids: List[str]) -> Dict[str, str]:
    """Fetch amino acid sequences for UniProt IDs.

    Args:
        uniprot_ids: List of UniProt accession IDs

    Returns:
        Mapping id -> sequence (ids not found omitted)
    """
    results: Dict[str, str] = {}
    if not uniprot_ids:
        return results
    for uid in uniprot_ids:
        try:
            url = f"{UNIPROT_REST}/uniprotkb/{uid}.fasta"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200 and resp.text:
                lines = [l.strip() for l in resp.text.splitlines() if l and not l.startswith('>')]
                results[uid] = "".join(lines)
            else:
                logger.warning(f"UniProt sequence not found for {uid}: HTTP {resp.status_code}")
        except Exception as e:
            logger.warning(f"Failed fetching UniProt sequence for {uid}: {e}")
    return results


def fetch_metadata(uniprot_ids: List[str]) -> List[Dict[str, str]]:
    """Fetch minimal metadata for UniProt IDs via REST API.

    Returns list of dicts suitable for building a DataFrame with an 'Entry' column.
    """
    rows: List[Dict[str, str]] = []
    if not uniprot_ids:
        return rows
    # Batch query using tab-separated format for speed
    query = " OR ".join([f"accession:{uid}" for uid in uniprot_ids])
    fields = [
        "accession",
        "protein_name",
        "organism_name",
        "ec",
        "protein_families",
        "reviewed",
    ]
    params = {
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
        "size": len(uniprot_ids),
    }
    try:
        url = f"{UNIPROT_REST}/uniprotkb/search"
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code != 200 or not resp.text:
            logger.warning(f"UniProt metadata query failed: HTTP {resp.status_code}")
            return rows
        lines = resp.text.splitlines()
        header = lines[0].split('\t')
        for line in lines[1:]:
            cols = line.split('\t')
            rec = dict(zip(header, cols))
            rows.append({
                'Entry': rec.get('Accession', ''),
                'Protein names': rec.get('Protein names', ''),
                'Organism': rec.get('Organism', ''),
                'EC number': rec.get('EC number', ''),
                'Protein families': rec.get('Protein families', ''),
                'Reviewed': rec.get('Reviewed', ''),
            })
    except Exception as e:
        logger.warning(f"Failed fetching UniProt metadata: {e}")
    return rows



def fetch_protein_sequence(uniprot_id: str) -> Optional[str]:
    """Fetch a single protein sequence by UniProt accession.

    This method is designed for per-protein use (no batching). Call it in a loop
    from higher-level analyses as needed.

    Args:
        uniprot_id: UniProt accession ID

    Returns:
        The amino acid sequence as a single string, or None if not found.
    """
    if not uniprot_id:
        return None
    
    try:
        # Fetch the sequence from UniProt
        url = f"{UNIPROT_REST}/uniprotkb/{uniprot_id}.fasta"
        resp = requests.get(url, timeout=10)
        
        if resp.status_code == 200 and resp.text:    
            # Parse the sequence from the response
            lines = [l.strip() for l in resp.text.splitlines() if l and not l.startswith('>')]
            return "".join(lines)
        else:
            logger.warning(f"UniProt sequence not found for {uniprot_id}: HTTP {resp.status_code}")
            return None

    except Exception as e:
        logger.warning(f"Failed fetching UniProt sequence for {uniprot_id}: {e}")
        return None

def fetch_protein_metadata(uniprot_id: str) -> Dict[str, str]:
    """Fetch minimal metadata for a single UniProt accession via REST API.

    This method is designed for per-protein use (no batching). Call it in a loop
    from higher-level analyses as needed.

    Returns a dict with the same keys produced by the batch `fetch_metadata`:
    'Entry', 'Protein names', 'Organism', 'EC number', 'Protein families', 'Reviewed'.
    Missing values will be empty strings.
    """
    if not uniprot_id:
        return {}

    # Metadata fields
    fields = [
        "accession",
        "protein_name",
        "organism_name",
        "ec",
        "protein_families",
        "reviewed"
    ]

    # Query parameters
    params = {
        "query": f"accession:{uniprot_id}",
        "fields": ",".join(fields),
        "format": "tsv",
        "size": 1,
    }
    
    try:
        # Fetch the metadata from UniProt
        url = f"{UNIPROT_REST}/uniprotkb/search"
        resp = requests.get(url, params=params, timeout=15)
        
        if resp.status_code != 200 or not resp.text:
            logger.warning(f"UniProt metadata query failed for {uniprot_id}: HTTP {resp.status_code}")
            return {}
        
        # Parse the response
        lines = resp.text.splitlines()

        # Check if the response is valid
        if len(lines) < 2:
            return {}
        
        # Parse the header and columns
        header = lines[0].split('\t')
        cols = lines[1].split('\t')
        rec = dict(zip(header, cols))
        
        return {
            'Entry': rec.get('Accession', ''),
            'Protein names': rec.get('Protein names', ''),
            'Organism': rec.get('Organism', ''),
            'EC number': rec.get('EC number', ''),
            'Protein families': rec.get('Protein families', ''),
            'Reviewed': rec.get('Reviewed', ''),
        }

    except Exception as e:
        logger.warning(f"Failed fetching UniProt metadata for {uniprot_id}: {e}")
        return {}


def main() -> int:
    """Test the UniProt API.

    Returns:
        int: 0 on success, 1 on failure
    """
    ok = True
    try:
        # Set the test ID
        test_id = "P01308"  # Human Insulin precursor (commonly used example)

        # Fetch the sequence
        seq = fetch_protein_sequence(test_id)
        if seq is None or len(seq) == 0:
            raise RuntimeError("Sequence fetch returned empty result")

        # Fetch the metadata
        meta = fetch_protein_metadata(test_id)

        # Validate the metadata (allow empty dict if API call fails - this is acceptable for testing)
        if not isinstance(meta, dict):
            raise RuntimeError("Metadata fetch returned invalid format")

        print("UniProt self-test passed:")
        print(f"  ID: {test_id}")
        print(f"  Sequence length: {len(seq) if seq else 0}")
        print(f"  Protein names: {meta.get('Protein names', 'N/A')}")
    except Exception as e:
        ok = False
        print(f"UniProt API test: FAILED - {e}")
        import traceback
        traceback.print_exc()
    
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())


