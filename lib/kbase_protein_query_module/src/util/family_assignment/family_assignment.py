"""
Family Assignment Module for KBase Protein Query Module
"""

import logging
import os
import time
import numpy as np
from typing import Dict, Any, List, Optional, Union

logger = logging.getLogger(__name__)

class FamilyAssignment:
    """Assigns proteins to families based on embeddings."""
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.centroids_path = config.get('centroids_path', 'data/family_centroids/files/family_centroids_binary.npz')
        self.family_assigner = None
        self._load_family_assigner()
    
    def _load_family_assigner(self):
        """Load the family assigner with centroids."""
        try:
            from ..storage import ProteinFamilyAssigner
            self.family_assigner = ProteinFamilyAssigner()
            self.family_assigner.load_family_centroids(self.centroids_path)
            logger.info("Family assigner loaded successfully")
        except Exception as e:
            logger.warning(f"Could not load family assigner: {e}")
            self.family_assigner = None
    
    def assign_families(self, embeddings: Dict[str, np.ndarray]) -> Dict[str, Dict[str, Any]]:
        """
        Assign proteins to families based on embeddings.
        
        Args:
            embeddings: Dictionary mapping protein IDs to embeddings
            
        Returns:
            Dictionary mapping protein IDs to family assignment results
        """
        if not self.family_assigner:
            raise RuntimeError("Family assigner not loaded")
        
        try:
            logger.info(f"Assigning families for {len(embeddings)} proteins")
            
            results = {}
            for protein_id, embedding in embeddings.items():
                try:
                    assignment = self.family_assigner.assign_family(embedding)
                    results[protein_id] = {
                        'family_id': assignment.get('family_id', 'unknown'),
                        'confidence': assignment.get('confidence', 0.0),
                        'eigenprotein_id': assignment.get('eigenprotein_id', ''),
                        'assignment_method': 'binary_centroid_search'
                    }
                except Exception as e:
                    logger.warning(f"Failed to assign family for {protein_id}: {e}")
                    results[protein_id] = {
                        'family_id': 'unknown',
                        'confidence': 0.0,
                        'eigenprotein_id': '',
                        'assignment_method': 'failed',
                        'error': str(e)
                    }
            
            logger.info(f"Completed family assignment for {len(results)} proteins")
            return results
            
        except Exception as e:
            logger.error(f"Family assignment failed: {e}")
            raise
    
    def assign_family_single(self, embedding: np.ndarray) -> Dict[str, Any]:
        """
        Assign a single protein to a family.
        
        Args:
            embedding: Protein embedding vector
            
        Returns:
            Family assignment result
        """
        if not self.family_assigner:
            raise RuntimeError("Family assigner not loaded")
        
        try:
            assignment = self.family_assigner.assign_family(embedding)
            return {
                'family_id': assignment.get('family_id', 'unknown'),
                'confidence': assignment.get('confidence', 0.0),
                'eigenprotein_id': assignment.get('eigenprotein_id', ''),
                'assignment_method': 'binary_centroid_search'
            }
        except Exception as e:
            logger.error(f"Single family assignment failed: {e}")
            return {
                'family_id': 'unknown',
                'confidence': 0.0,
                'eigenprotein_id': '',
                'assignment_method': 'failed',
                'error': str(e)
            }
    
    def get_family_stats(self) -> Dict[str, Any]:
        """
        Get statistics about available families.
        
        Returns:
            Dictionary with family statistics
        """
        if not self.family_assigner:
            return {'error': 'Family assigner not loaded'}
        
        try:
            families = self.family_assigner.list_families()
            return {
                'total_families': len(families),
                'families': families,
                'centroids_loaded': True
            }
        except Exception as e:
            return {
                'error': str(e),
                'centroids_loaded': False
            }
    
    def is_available(self) -> bool:
        """Check if family assignment is available."""
        return self.family_assigner is not None