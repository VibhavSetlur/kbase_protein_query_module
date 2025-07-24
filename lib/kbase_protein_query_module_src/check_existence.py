"""
Protein Existence Checker Module

This module provides a fast method to check if a protein (by UniProt ID)
exists in the storage system, and returns its family and metadata if found.
Integrates with the storage and workflow modules for seamless pipeline use.
"""

import logging
from typing import Optional, Dict, Any
import pandas as pd
from kbase_protein_query_module_src.storage import ProteinStorage

logger = logging.getLogger(__name__)

class ProteinExistenceChecker:
    """
    Checks if a protein exists in the storage and returns its family and metadata.
    """
    def __init__(self, storage: Optional[ProteinStorage] = None, base_dir: str = "data"):
        if storage is not None:
            self.storage = storage
        else:
            self.storage = ProteinStorage(base_dir=base_dir)
        self.family_list = self.storage.get_family_list()
        self.metadata_storage = None
        try:
            from kbase_protein_query_module_src.storage import CompressedMetadataStorage
            self.metadata_storage = CompressedMetadataStorage(metadata_dir=str(self.storage.metadata_dir))
        except Exception as e:
            logger.warning(f"Could not initialize CompressedMetadataStorage: {e}")

    def check_protein_existence(self, protein_id: str) -> Dict[str, Any]:
        """
        Check if a protein exists by UniProt ID.
        Args:
            protein_id: UniProt or other protein ID
        Returns:
            Dict with keys: exists (bool), family_id (str or None), metadata (dict or None)
        """
        if not protein_id:
            raise ValueError("Must provide a protein_id (UniProt ID).")
        for family_id in self.family_list:
            try:
                if self.metadata_storage:
                    try:
                        metadata = self.metadata_storage.load_metadata(family_id=family_id, protein_ids=[protein_id])
                        if not metadata.empty:
                            return {
                                "exists": True,
                                "family_id": family_id,
                                "metadata": metadata.iloc[0].to_dict()
                            }
                    except Exception:
                        pass
                _, protein_ids = self.storage.load_family_embeddings(family_id)
                if protein_id in protein_ids:
                    meta = None
                    if self.metadata_storage:
                        try:
                            meta_df = self.metadata_storage.load_metadata(family_id=family_id, protein_ids=[protein_id])
                            if not meta_df.empty:
                                meta = meta_df.iloc[0].to_dict()
                        except Exception:
                            pass
                    return {
                        "exists": True,
                        "family_id": family_id,
                        "metadata": meta
                    }
            except Exception as e:
                logger.debug(f"Error checking family {family_id}: {e}")
                continue
        return {"exists": False, "family_id": None, "metadata": None} 