"""
Protein Existence Checker Module

This module provides a fast method to check if a protein (by UniProt ID)
exists in the storage system, and returns its family and metadata if found.
Integrates with the storage and workflow modules for seamless pipeline use.
"""

import logging
from typing import Optional, Dict, Any
import pandas as pd
from .protein_storage import ProteinStorage, ProteinIDsIndex

logger = logging.getLogger(__name__)

class ProteinExistenceChecker:
    """
    Checks if a protein exists in the storage and returns its family and metadata.
    Uses efficient protein IDs index for fast searching (exact UniProt ID match).
    """
    def __init__(self, storage: Optional[ProteinStorage] = None, base_dir: str = "data"):
        if storage is not None:
            self.storage = storage
        else:
            self.storage = ProteinStorage(base_dir=base_dir)
        # Initialize protein IDs index for efficient searching
        self.protein_ids_index = ProteinIDsIndex(base_dir=base_dir)
        self.family_list = self.storage.get_family_list()
        self.metadata_storage = None
        try:
            from .protein_storage import CompressedMetadataStorage
            self.metadata_storage = CompressedMetadataStorage(metadata_dir=str(self.storage.metadata_dir))
        except Exception as e:
            logger.warning(f"Could not initialize CompressedMetadataStorage: {e}")

    def check_protein_existence(self, uniprot_id: str) -> Dict[str, Any]:
        """
        Check if a protein exists by UniProt ID (exact match only).
        Args:
            uniprot_id: UniProt ID (e.g., P00001)
        Returns:
            Dict with keys: exists (bool), family_id (str or None), metadata (dict or None)
        """
        if not uniprot_id:
            raise ValueError("Must provide a UniProt ID.")
        uniprot_id = uniprot_id.strip()
        if not uniprot_id:
            raise ValueError("UniProt ID cannot be empty or whitespace only.")
        # Use the efficient protein IDs index (exact match)
        index_result = self.protein_ids_index.search_protein(uniprot_id)
        if index_result:
            family_id = index_result['family_id']
            metadata = index_result.get('metadata')
            return {
                "exists": True,
                "family_id": family_id,
                "metadata": metadata
            }
        # Not found
        return {"exists": False, "family_id": None, "metadata": None}
