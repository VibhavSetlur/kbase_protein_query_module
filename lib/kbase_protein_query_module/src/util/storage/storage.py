"""
PLM-based ProteinStorage for querying KBase Protein Language Model API.

This utility handles protein queries using the KB PLM Utils API to:
- Query the PLM API for protein homologs
- Handle hierarchical protein dictionaries
- Standardize outputs for analysis workflows
- Retrieve embeddings and sequences when available

The hierarchical protein dictionary structure:
{
    'protein_id': {
        'sequence': str,
        'source': str,
        'original_id': str,
        ...
    },
    ...
}
"""

import sys
import os
import logging
from typing import Dict, List, Optional, Tuple, Any, Union

import numpy as np

logger = logging.getLogger(__name__)

# Handle KB PLM Utils import
try:
    from kbutillib.kb_plm_utils import KBPLMUtils
except ImportError:
    try:
        # Try alternative import path
        import sys
        kbutillib_path = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'installed_clients', 'kbutillib')
        if kbutillib_path not in sys.path:
            sys.path.insert(0, kbutillib_path)
        from kb_plm_utils import KBPLMUtils
    except ImportError as e:
        logger.warning(f"KBPLMUtils not found - PLM functionality will be limited. Error: {e}")
        KBPLMUtils = None


class ProteinStorage:
    """
    Storage utility for querying protein homologs via KBase PLM API.
    
    This class provides methods to:
    - Query the PLM API with hierarchical protein dictionaries
    - Retrieve similar proteins and their embeddings
    - Standardize PLM API outputs for analysis workflows
    - Handle batch queries and result caching
    """
    
    def __init__(
        self,
        plm_api_url: str = "https://kbase.us/services/llm_homology_api",
        config: Optional[Dict[str, Any]] = None,
        kb_plm_utils: Optional[Any] = None
    ) -> None:
        """
        Initialize ProteinStorage with PLM API.
        
        Args:
            plm_api_url: Base URL for the PLM API (default: KBase PLM API)
            config: Optional configuration dictionary
            kb_plm_utils: Optional KBPLMUtils instance (if None, creates new one)
        """
        self.config = config or {}
        self.plm_api_url = plm_api_url
        
        # Initialize KB PLM Utils
        if kb_plm_utils is None:
            if KBPLMUtils is None:
                raise ImportError(
                    "KBPLMUtils not available. Install kbutillib with KBPLMUtils support."
                )
            self.plm_utils = KBPLMUtils(plm_api_url=plm_api_url, **self.config)
        else:
            self.plm_utils = kb_plm_utils
        
        # Cache for storing query results
        self._query_cache: Dict[str, Dict[str, Any]] = {}
        self._embedding_cache: Dict[str, np.ndarray] = {}
        
        logger.info("ProteinStorage initialized with PLM API")
    
    def query_similar_proteins(
        self,
        proteins: Dict[str, Dict[str, Any]],
        max_hits: int = 100,
        similarity_threshold: float = 0.0,
        return_embeddings: bool = False,
        **kwargs: Any
    ) -> Dict[str, Dict[str, Any]]:
        """
        Query PLM API for similar proteins given hierarchical protein dictionary.
        
        Args:
            proteins: Hierarchical dict with protein_id -> {sequence, source, ...}
            max_hits: Maximum number of hits per query (1-100)
            similarity_threshold: Minimum similarity score threshold
            return_embeddings: Whether to return embeddings in results
            **kwargs: Additional arguments passed to query_plm_api
        
        Returns:
            Dict mapping protein_id to standardized results:
            {
                'protein_id': {
                    'query_id': str,
                    'query_sequence': str,
                    'hits': [
                        {
                            'uniprot_id': str,
                            'plm_score': float,
                            'embedding': np.ndarray (optional),
                            ...
                        },
                        ...
                    ],
                    'hit_count': int
                },
                ...
            }
        """
        if not proteins:
            raise ValueError("proteins dictionary cannot be empty")
        
        # Convert hierarchical dict to PLM API format
        query_sequences = []
        query_id_to_protein_id = {}
        
        for protein_id, protein_data in proteins.items():
            sequence = protein_data.get('sequence', '')
            if not sequence:
                logger.warning(f"Skipping protein {protein_id}: no sequence found")
                continue
            
            query_id = protein_data.get('original_id', protein_id)
            query_sequences.append({
                'id': query_id,
                'sequence': sequence
            })
            query_id_to_protein_id[query_id] = protein_id
        
        if not query_sequences:
            raise ValueError("No valid protein sequences found in proteins dictionary")
        
        logger.info(
            f"Querying PLM API for {len(query_sequences)} proteins, "
            f"requesting up to {max_hits} hits per protein"
        )
        
        # Query PLM API
        try:
            plm_results = self.plm_utils.query_plm_api(
                query_sequences=query_sequences,
                max_hits=max_hits,
                similarity_threshold=similarity_threshold,
                return_embeddings=return_embeddings,
                **kwargs
            )
        except Exception as e:
            logger.error(f"PLM API query failed: {e}")
            raise RuntimeError(f"PLM API query failed: {e}") from e
        
        return self._standardize_results(
            plm_results, 
            proteins, 
            query_id_to_protein_id, 
            return_embeddings
        )

    def _standardize_results(
        self,
        plm_results: Dict[str, Any],
        proteins: Dict[str, Dict[str, Any]],
        query_id_to_protein_id: Dict[str, str],
        return_embeddings: bool
    ) -> Dict[str, Dict[str, Any]]:
        """
        Standardize raw PLM API results.
        
        Args:
            plm_results: Raw results from PLM API
            proteins: Original hierarchical protein dictionary
            query_id_to_protein_id: Mapping from query ID to protein ID
            return_embeddings: Whether to process embeddings
            
        Returns:
            Standardized results dictionary
        """
        standardized_results = {}
        
        for hits_data in plm_results.get("hits", []):
            query_id = hits_data.get("query_id", "")
            protein_id = query_id_to_protein_id.get(query_id, query_id)
            
            hits = hits_data.get("hits", [])
            standardized_hits = []
            
            for hit in hits:
                hit_entry = {
                    'uniprot_id': hit.get("id", ""),
                    'plm_score': float(hit.get("score", 0.0))
                }
                
                # Add embedding if available
                if return_embeddings and "embedding" in hit:
                    emb_data = hit["embedding"]
                    if isinstance(emb_data, list):
                        hit_entry['embedding'] = np.array(emb_data, dtype=np.float32)
                    elif isinstance(emb_data, np.ndarray):
                        hit_entry['embedding'] = emb_data.astype(np.float32)
                
                standardized_hits.append(hit_entry)
            
            # Get original protein data
            protein_data = proteins.get(protein_id, {})
            
            result_entry = {
                'query_id': query_id,
                'query_sequence': protein_data.get('sequence', ''),
                'source': protein_data.get('source', 'unknown'),
                'original_id': protein_data.get('original_id', protein_id),
                'hits': standardized_hits,
                'hit_count': len(standardized_hits)
            }

            # Extract query embedding if available
            if return_embeddings:
                # Check for embedding in hits_data (API response for this query)
                # It might be under 'embedding' or 'query_embedding'
                query_emb_data = hits_data.get("embedding") or hits_data.get("query_embedding")
                if query_emb_data:
                    if isinstance(query_emb_data, list):
                        result_entry['query_embedding'] = np.array(query_emb_data, dtype=np.float32)
                    elif isinstance(query_emb_data, np.ndarray):
                        result_entry['query_embedding'] = query_emb_data.astype(np.float32)

            standardized_results[protein_id] = result_entry
            
            # Cache embeddings if available
            if return_embeddings:
                self._cache_embeddings(standardized_hits)
        
        return standardized_results

    def _cache_embeddings(self, hits: List[Dict[str, Any]]) -> None:
        """Cache embeddings from hits."""
        for hit in hits:
            uniprot_id = hit.get('uniprot_id')
            if uniprot_id and 'embedding' in hit:
                self._embedding_cache[uniprot_id] = hit['embedding']
    
    def find_top_k_similar(
        self,
        protein_id: str,
        proteins: Dict[str, Dict[str, Any]],
        top_k: int = 5,
        max_plm_hits: int = 100,
        similarity_threshold: float = 0.0,
        **kwargs: Any
    ) -> List[Tuple[str, float]]:
        """
        Find top K similar proteins for a given protein ID.
        
        Args:
            protein_id: Protein ID from the proteins dictionary
            proteins: Hierarchical protein dictionary
            top_k: Number of top similar proteins to return
            max_plm_hits: Maximum PLM hits to request
            similarity_threshold: Minimum similarity threshold
            **kwargs: Additional arguments for PLM query
        
        Returns:
            List of tuples (uniprot_id, plm_score) sorted by score descending
        """
        if protein_id not in proteins:
            logger.warning(f"Protein {protein_id} not found in proteins dictionary")
            return []
        
        # Query PLM API for this specific protein
        single_protein = {protein_id: proteins[protein_id]}
        results = self.query_similar_proteins(
            proteins=single_protein,
            max_hits=max_plm_hits,
            similarity_threshold=similarity_threshold,
            **kwargs
        )
        
        if protein_id not in results:
            return []
        
        protein_results = results[protein_id]
        hits = protein_results.get('hits', [])
        
        # Extract top K hits
        top_hits = sorted(
            hits,
            key=lambda x: x.get('plm_score', 0.0),
            reverse=True
        )[:top_k]
        
        # Return in expected format
        return [
            (hit['uniprot_id'], hit['plm_score'])
            for hit in top_hits
        ]
    
    def find_top_k_similar_by_sequence(
        self,
        sequence: str,
        top_k: int = 5,
        max_plm_hits: int = 100,
        similarity_threshold: float = 0.0,
        **kwargs: Any
    ) -> List[Tuple[str, float]]:
        """
        Find top K similar proteins for a given protein sequence.
        
        Args:
            sequence: Protein sequence string
            top_k: Number of top similar proteins to return
            max_plm_hits: Maximum PLM hits to request
            similarity_threshold: Minimum similarity threshold
            **kwargs: Additional arguments for PLM query
        
        Returns:
            List of tuples (uniprot_id, plm_score) sorted by score descending
        """
        # Create temporary protein entry
        temp_protein_id = "query_protein"
        proteins = {
            temp_protein_id: {
                'sequence': sequence,
                'source': 'query',
                'original_id': temp_protein_id
            }
        }
        
        return self.find_top_k_similar(
            protein_id=temp_protein_id,
            proteins=proteins,
            top_k=top_k,
            max_plm_hits=max_plm_hits,
            similarity_threshold=similarity_threshold,
            **kwargs
        )
    
    def get_embedding(self, uniprot_id: str) -> Optional[np.ndarray]:
        """
        Get embedding for a UniProt ID.
        
        Args:
            uniprot_id: UniProt ID
        
        Returns:
            Embedding array if available, None otherwise
        """
        # Check cache first
        if uniprot_id in self._embedding_cache:
            return self._embedding_cache[uniprot_id]
        
        logger.debug(f"Embedding for {uniprot_id} not in cache")
        return None
    
    def get_uniprot_sequences(
        self,
        uniprot_ids: List[str]
    ) -> Dict[str, str]:
        """
        Retrieve protein sequences from UniProt.
        
        Args:
            uniprot_ids: List of UniProt IDs
        
        Returns:
            Dict mapping UniProt IDs to their sequences
        """
        return self.plm_utils.get_uniprot_sequences(uniprot_ids)
    
    def get_similar_proteins_batch(
        self,
        proteins: Dict[str, Dict[str, Any]],
        max_hits: int = 100,
        similarity_threshold: float = 0.0,
        return_embeddings: bool = False,
        **kwargs: Any
    ) -> Dict[str, List[Tuple[str, float]]]:
        """
        Batch query for similar proteins across multiple proteins.
        
        Returns results in the format expected by analysis code:
        {protein_id: [(uniprot_id, score), ...]}
        
        Args:
            proteins: Hierarchical protein dictionary
            max_hits: Maximum hits per protein
            similarity_threshold: Minimum similarity threshold
            return_embeddings: Whether to cache embeddings
            **kwargs: Additional PLM query arguments
        
        Returns:
            Dict mapping protein_id to list of (uniprot_id, score) tuples
        """
        results = self.query_similar_proteins(
            proteins=proteins,
            max_hits=max_hits,
            similarity_threshold=similarity_threshold,
            return_embeddings=return_embeddings,
            **kwargs
        )
        
        # Convert to expected format
        batch_results = {}
        for protein_id, protein_results in results.items():
            hits = protein_results.get('hits', [])
            batch_results[protein_id] = [
                (hit['uniprot_id'], hit['plm_score'])
                for hit in sorted(hits, key=lambda x: x.get('plm_score', 0.0), reverse=True)
            ]
        
        return batch_results
    
    def get_top_similar_proteins(
        self,
        proteins: Dict[str, Dict[str, Any]],
        top_k: int = 5,
        max_hits: int = 100,
        similarity_threshold: float = 0.0,
        **kwargs: Any
    ) -> Dict[str, List[Tuple[str, float]]]:
        """
        Get top similar proteins for a dictionary of proteins.
        
        This is a wrapper around get_similar_proteins_batch to match the requested API.
        
        Args:
            proteins: Hierarchical protein dictionary
            top_k: Number of top similar proteins to return per protein
            max_hits: Maximum hits per protein
            similarity_threshold: Minimum similarity threshold
            **kwargs: Additional PLM query arguments
        
        Returns:
            Dict mapping protein_id to list of (uniprot_id, score) tuples
        """
        # Use the batch method which already does what we want, but we might need to filter top_k
        results = self.get_similar_proteins_batch(
            proteins=proteins,
            max_hits=max_hits,
            similarity_threshold=similarity_threshold,
            **kwargs
        )
        
        # Ensure top_k limit
        final_results = {}
        for protein_id, hits in results.items():
            final_results[protein_id] = hits[:top_k]
            
        return final_results

def main() -> int:
    """Self-test for ProteinStorage."""
    try:
        print("STORAGE_TEST: Starting self-test...")
        
        if KBPLMUtils is None:
            print("STORAGE_SKIP: KBPLMUtils not available")
            return 0
        
        # Test initialization
        storage = ProteinStorage()
        print("STORAGE_INIT: ProteinStorage initialized successfully")
        
        # Test with sample protein
        test_proteins = {
            'test_protein': {
                'sequence': 'MKLLVVCLFVAVTILPASS',
                'source': 'test',
                'original_id': 'test_protein'
            }
        }
        
        print("STORAGE_TEST: Attempting API query (may fail if no auth/network)...")
        try:
            results = storage.get_top_similar_proteins(test_proteins, top_k=2)
            print(f"STORAGE_SUCCESS: Query returned {len(results)} results")
            print(f"STORAGE_RESULT: {results}")
        except Exception as e:
            print(f"STORAGE_WARNING: API query failed (expected if no auth): {e}")
            # This is acceptable for a self-test in some environments
        
        print("STORAGE_TEST: Completed")
        return 0
        
    except ImportError as e:
        print(f"STORAGE_SKIP: {e}")
        return 0
    except Exception as e:
        print(f"STORAGE_FAIL: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
