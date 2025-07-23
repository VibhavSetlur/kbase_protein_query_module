"""
Protein Family Assignment Module

This module provides fast assignment of protein embeddings to precomputed protein families
using centroid similarity. It is designed for integration with the KBase protein network
analysis pipeline and is compatible with the workflow orchestrator and test suite.
"""

import logging
from typing import Optional, Union, Any, Dict
import numpy as np

logger = logging.getLogger(__name__)

class AssignProteinFamily:
    """
    Assigns protein embeddings to precomputed protein families using centroid similarity.

    This class loads family centroids from a .npz file and provides methods to assign
    a query embedding to the closest family centroid using cosine similarity.
    Compatible with the workflow orchestrator and test suite.
    """
    def __init__(self):
        self.family_ids: Optional[np.ndarray] = None
        self.centroids: Optional[np.ndarray] = None
        self.eigenprotein_ids: Optional[np.ndarray] = None

    def load_family_centroids(self, npz_path: str) -> None:
        """
        Load family centroids, family IDs, and eigenprotein IDs from a .npz file.

        Args:
            npz_path: Path to the .npz file containing 'family_ids', 'centroids', and 'eigenprotein_ids'.
        """
        try:
            data = np.load(npz_path, allow_pickle=True)
            self.family_ids = data['family_ids']
            self.centroids = data['centroids']
            self.eigenprotein_ids = data['eigenprotein_ids']
            logger.info(f"Loaded {len(self.family_ids)} family centroids from {npz_path}")
        except Exception as e:
            logger.error(f"Failed to load family centroids from {npz_path}: {e}")
            raise

    def assign_family(self, embedding: Union[np.ndarray, list, Any]) -> Dict[str, Any]:
        """
        Assign a query embedding to the closest family centroid using Hamming distance (binary).
        Args:
            embedding: Query protein embedding (list or np.ndarray, np.uint8 expected)
        Returns:
            Dict with keys: 'family_id', 'confidence', 'eigenprotein_id'
        """
        if self.centroids is None or self.family_ids is None or self.eigenprotein_ids is None:
            raise ValueError("Family centroids and IDs must be loaded before assignment.")
        emb = np.array(embedding, dtype=np.uint8)
        emb_bits = np.unpackbits(emb)
        centroids_bits = np.array([np.unpackbits(c) for c in self.centroids])
        dists = np.sum(centroids_bits != emb_bits, axis=1)
        idx = int(np.argmin(dists))
        result = {
            'family_id': str(self.family_ids[idx]),
            'confidence': float(-dists[idx]),  # negative Hamming distance as confidence
            'eigenprotein_id': str(self.eigenprotein_ids[idx])
        }
        logger.debug(f"Assigned embedding to family {result['family_id']} (Hamming distance={-result['confidence']:.3f})")
        return result

    def predict_family(self, embedding: Union[np.ndarray, list, Any]) -> (str, float):
        """
        Predict the closest family ID and confidence for a query embedding.
        This is a compatibility alias for assign_family, returning only family_id and confidence.

        Args:
            embedding: Query protein embedding (list or np.ndarray)

        Returns:
            Tuple of (family_id, confidence)
        """
        result = self.assign_family(embedding)
        return result['family_id'], result['confidence'] 