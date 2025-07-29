"""
Protein Family Assignment Module

This module provides fast assignment of protein embeddings to precomputed protein families
using centroid similarity. It is designed for integration with the KBase protein network
analysis pipeline and is compatible with the workflow orchestrator and test suite.
"""

import logging
from typing import Optional, Union, Any, Dict
import numpy as np
import faiss

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
            self.centroids = np.ascontiguousarray(data['centroids'])  # keep as float, do not cast to uint8
            self.eigenprotein_ids = data['eigenprotein_ids']
            self._faiss_index = None
            logger.info(f"Loaded {len(self.family_ids)} family centroids from {npz_path}")
        except Exception as e:
            logger.error(f"Failed to load family centroids from {npz_path}: {e}")
            raise

    def assign_family(self, embedding: Union[np.ndarray, list, Any]) -> Dict[str, Any]:
        """
        Assign a query embedding to the closest family centroid using FAISS binary IVF (Hamming distance).
        Args:
            embedding: Query protein embedding (list or np.ndarray, must be np.float32)
        Returns:
            Dict with keys: 'family_id', 'confidence' (normalized similarity), 'eigenprotein_id'
        """
        if self.centroids is None or self.family_ids is None or self.eigenprotein_ids is None:
            raise ValueError("Family centroids and IDs must be loaded before assignment.")
        emb = np.array(embedding, dtype=np.float32)
        # Binarize: sign(x) > 0 → 1, else 0
        emb_bin = (emb > 0).astype(np.uint8)
        d = self.centroids.shape[1]  # number of features (bits)
        needed_bits = ((d + 7) // 8) * 8  # pad to next byte boundary
        if emb_bin.size < needed_bits:
            emb_bin = np.pad(emb_bin, (0, needed_bits - emb_bin.size), 'constant')
        elif emb_bin.size > needed_bits:
            emb_bin = emb_bin[:needed_bits]
        emb_bin_packed = np.packbits(emb_bin)
        emb_bin_packed = np.ascontiguousarray(emb_bin_packed.reshape(1, -1))
        # Build FAISS binary flat index if not already
        if not hasattr(self, '_faiss_index') or self._faiss_index is None:
            import faiss
            centroids_bin = (self.centroids > 0).astype(np.uint8)
            if centroids_bin.shape[1] < needed_bits:
                centroids_bin = np.pad(centroids_bin, ((0,0),(0, needed_bits - centroids_bin.shape[1])), 'constant')
            elif centroids_bin.shape[1] > needed_bits:
                centroids_bin = centroids_bin[:, :needed_bits]
            centroids_bin_packed = np.packbits(centroids_bin, axis=1)
            centroids_bin_packed = np.ascontiguousarray(centroids_bin_packed)
            self._packed_centroids = centroids_bin_packed
            self._faiss_index = faiss.IndexBinaryFlat(d)
            self._faiss_index.add(self._packed_centroids)
        D, I = self._faiss_index.search(emb_bin_packed, 1)
        idx = int(I[0][0])
        hamming_distance = float(D[0][0])
        # Confidence: normalized similarity (1 - hamming_distance / embedding_length)
        confidence = 1.0 - (hamming_distance / d)
        result = {
            'family_id': str(self.family_ids[idx]),
            'confidence': confidence,
            'eigenprotein_id': str(self.eigenprotein_ids[idx])
        }
        logger.debug(f"Assigned embedding to family {result['family_id']} (confidence={result['confidence']:.3f})")
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