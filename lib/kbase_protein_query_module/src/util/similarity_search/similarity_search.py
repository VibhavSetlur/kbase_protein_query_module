"""
Similarity Search Module for KBase Protein Query Module

This module performs similarity search using FAISS float indices for within-family search.
"""

import logging
import time
import numpy as np
from typing import Dict, Any, List, Optional, Union

logger = logging.getLogger(__name__)

class SimilaritySearch:
    """
    Performs similarity search within protein families.
    
    Uses FAISS float indices for efficient similarity search.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.index_dir = config.get('index_dir', 'data/indexes')
        self.top_k = config.get('top_k', 50)
        self.similarity_threshold = config.get('similarity_threshold', 0.75)
        self.hierarchical_index = None
        self._load_hierarchical_index()
    
    def _load_hierarchical_index(self):
        """Load the hierarchical index for similarity search."""
        try:
            from ..similarity_indexing import HierarchicalIndex
            self.hierarchical_index = HierarchicalIndex(base_dir=self.index_dir)
            logger.info("Hierarchical index loaded successfully")
        except Exception as e:
            logger.warning(f"Could not load hierarchical index: {e}")
            self.hierarchical_index = None
    
    def search_similar_proteins(self, 
                               embeddings: Dict[str, np.ndarray],
                               family_assignments: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """
        Search for similar proteins within assigned families.
        
        Args:
            embeddings: Dictionary mapping protein IDs to embeddings
            family_assignments: Dictionary mapping protein IDs to family assignments
            
        Returns:
            Dictionary mapping protein IDs to similarity search results
        """
        if not self.hierarchical_index:
            raise RuntimeError("Hierarchical index not loaded")
        
        try:
            logger.info(f"Searching for similar proteins for {len(embeddings)} proteins")
            
            results = {}
            for protein_id, embedding in embeddings.items():
                try:
                    family_assignment = family_assignments.get(protein_id, {})
                    family_id = family_assignment.get('family_id', 'unknown')
                    
                    if family_id == 'unknown':
                        logger.warning(f"No family assignment for {protein_id}, skipping similarity search")
                        results[protein_id] = {
                            'similar_proteins': [],
                            'similarity_scores': [],
                            'total_found': 0,
                            'family_id': 'unknown',
                            'search_method': 'skipped_no_family'
                        }
                        continue
                    
                    # Search within the assigned family
                    similarities, similar_protein_ids = self.hierarchical_index.search_family_float(
                        family_id, embedding, self.top_k
                    )
                    
                    # Filter by similarity threshold
                    filtered_results = []
                    filtered_scores = []
                    
                    for i, (similarity, similar_id) in enumerate(zip(similarities, similar_protein_ids)):
                        if similarity >= self.similarity_threshold and similar_id != protein_id:
                            filtered_results.append(similar_id)
                            filtered_scores.append(float(similarity))
                    
                    results[protein_id] = {
                        'similar_proteins': filtered_results,
                        'similarity_scores': filtered_scores,
                        'total_found': len(filtered_results),
                        'family_id': family_id,
                        'search_method': 'family_based_float_search',
                        'similarity_threshold': self.similarity_threshold
                    }
                    
                except Exception as e:
                    logger.warning(f"Failed similarity search for {protein_id}: {e}")
                    results[protein_id] = {
                        'similar_proteins': [],
                        'similarity_scores': [],
                        'total_found': 0,
                        'family_id': family_assignments.get(protein_id, {}).get('family_id', 'unknown'),
                        'search_method': 'failed',
                        'error': str(e)
                    }
            
            logger.info(f"Completed similarity search for {len(results)} proteins")
            return results
            
        except Exception as e:
            logger.error(f"Similarity search failed: {e}")
            raise
    
    def search_single_protein(self, 
                            embedding: np.ndarray, 
                            family_id: str) -> Dict[str, Any]:
        """
        Search for similar proteins for a single protein.
        
        Args:
            embedding: Protein embedding vector
            family_id: Family ID to search within
            
        Returns:
            Similarity search results
        """
        if not self.hierarchical_index:
            raise RuntimeError("Hierarchical index not loaded")
        
        try:
            similarities, similar_protein_ids = self.hierarchical_index.search_family_float(
                family_id, embedding, self.top_k
            )
            
            # Filter by similarity threshold
            filtered_results = []
            filtered_scores = []
            
            for similarity, similar_id in zip(similarities, similar_protein_ids):
                if similarity >= self.similarity_threshold:
                    filtered_results.append(similar_id)
                    filtered_scores.append(float(similarity))
            
            return {
                'similar_proteins': filtered_results,
                'similarity_scores': filtered_scores,
                'total_found': len(filtered_results),
                'family_id': family_id,
                'search_method': 'family_based_float_search',
                'similarity_threshold': self.similarity_threshold
            }
            
        except Exception as e:
            logger.error(f"Single protein similarity search failed: {e}")
            return {
                'similar_proteins': [],
                'similarity_scores': [],
                'total_found': 0,
                'family_id': family_id,
                'search_method': 'failed',
                'error': str(e)
            }
    
    def get_index_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the similarity index.
        
        Returns:
            Dictionary with index statistics
        """
        if not self.hierarchical_index:
            return {'error': 'Hierarchical index not loaded'}
        
        try:
            family_stats = self.hierarchical_index.get_family_stats()
            return {
                'index_loaded': True,
                'family_stats': family_stats,
                'index_dir': self.index_dir,
                'top_k': self.top_k,
                'similarity_threshold': self.similarity_threshold
            }
        except Exception as e:
            return {
                'error': str(e),
                'index_loaded': False
            }
    
    def is_available(self) -> bool:
        """Check if similarity search is available."""
        return self.hierarchical_index is not None