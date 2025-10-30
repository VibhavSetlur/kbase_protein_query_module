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


def fetch_sequence_and_metadata(uniprot_id: str) -> Tuple[Optional[str], Dict[str, str]]:
    """Fetch the amino acid sequence and minimal metadata for a single UniProt ID.

    Returns (sequence_or_None, metadata_dict). Metadata dict contains at least 'Entry'.
    """
    sequence: Optional[str] = None
    metadata: Dict[str, str] = {'Entry': uniprot_id}

    try:
        seq_map = fetch_sequences([uniprot_id])
        sequence = seq_map.get(uniprot_id)
    except Exception as e:
        logger.warning(f"Failed fetching sequence for {uniprot_id}: {e}")

    try:
        rows = fetch_metadata([uniprot_id])
        if rows:
            metadata.update(rows[0])
    except Exception as e:
        logger.warning(f"Failed fetching metadata for {uniprot_id}: {e}")

    return sequence, metadata


