"""
General-purpose similarity search utility using cosine similarity.

This module provides utilities for computing cosine similarity between
embeddings and finding top-K similar items. Designed for general use
cases including protein similarity searches, query comparisons, and
other embedding-based similarity tasks.

The storage module retrieves embeddings (query and nearby embeddings),
and this module performs the similarity computations.
"""

import sys
import os
import logging
from typing import Dict, List, Optional, Tuple, Any, Union

import numpy as np

logger = logging.getLogger(__name__)

# Try importing FaissIndexManager
try:
    # Adjust import based on relative path
    from ..faiss.faiss import FaissIndexManager
except ImportError:
    try:
        # Fallback for different execution contexts
        from lib.kbase_protein_query_module.src.util.faiss.faiss import FaissIndexManager
    except ImportError:
        logger.warning("Could not import FaissIndexManager. FAISS functionality will be disabled.")
        FaissIndexManager = None


class SimilaritySearch:
    """
    General-purpose similarity search using cosine similarity on embeddings.
    
    This class provides methods to:
    - Compute cosine similarity between embeddings
    - Find top-K similar items
    - Perform batch similarity searches
    - Handle query-to-database and query-to-query similarity
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize SimilaritySearch.
        
        Args:
            config: Optional configuration dictionary
        """
        self.config = config or {}
        logger.info("SimilaritySearch initialized")
    
    def cosine_similarity(
        self,
        embedding1: np.ndarray,
        embedding2: np.ndarray
    ) -> float:
        """
        Compute cosine similarity between two embeddings.
        
        Args:
            embedding1: First embedding vector (1D array)
            embedding2: Second embedding vector (1D array)
        
        Returns:
            Cosine similarity score (float, range: -1.0 to 1.0)
        """
        if embedding1.shape != embedding2.shape:
            raise ValueError(
                f"Embedding shapes must match: {embedding1.shape} vs {embedding2.shape}"
            )
        
        # Normalize to unit vectors
        norm1 = np.linalg.norm(embedding1)
        norm2 = np.linalg.norm(embedding2)
        
        if norm1 == 0 or norm2 == 0:
            logger.warning("Zero-norm embedding detected, returning 0.0 similarity")
            return 0.0
        
        # Compute cosine similarity
        similarity = np.dot(embedding1, embedding2) / (norm1 * norm2)
        
        # Ensure result is in valid range (handle floating point errors)
        return float(np.clip(similarity, -1.0, 1.0))
    
    def compute_similarity_matrix(
        self,
        query_embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]],
        target_embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]],
        normalize: bool = True
    ) -> np.ndarray:
        """
        Compute similarity matrix between query and target embeddings.
        
        Args:
            query_embeddings: Query embeddings - can be:
                - 2D array [N_query, D]
                - List of 1D arrays
                - Dict mapping IDs to 1D arrays
            target_embeddings: Target embeddings - can be:
                - 2D array [N_target, D]
                - List of 1D arrays
                - Dict mapping IDs to 1D arrays
            normalize: Whether to normalize embeddings before computation
        
        Returns:
            Similarity matrix [N_query, N_target] with cosine similarities
        """
        # Convert to arrays if needed
        query_array, query_ids = self._to_array_with_ids(query_embeddings)
        target_array, target_ids = self._to_array_with_ids(target_embeddings)
        
        if query_array.ndim != 2 or target_array.ndim != 2:
            raise ValueError(
                f"Embeddings must be 2D arrays, got query: {query_array.ndim}D, "
                f"target: {target_array.ndim}D"
            )
        
        if query_array.shape[1] != target_array.shape[1]:
            raise ValueError(
                f"Embedding dimensions must match: {query_array.shape[1]} vs "
                f"{target_array.shape[1]}"
            )
        
        # Normalize if requested
        if normalize:
            query_norms = np.linalg.norm(query_array, axis=1, keepdims=True)
            target_norms = np.linalg.norm(target_array, axis=1, keepdims=True)
            
            # Avoid division by zero
            query_array = query_array / np.where(query_norms > 0, query_norms, 1.0)
            target_array = target_array / np.where(target_norms > 0, target_norms, 1.0)
        
        # Compute cosine similarity matrix: query @ target.T
        similarity_matrix = np.dot(query_array, target_array.T)
        
        # Clip to valid range
        similarity_matrix = np.clip(similarity_matrix, -1.0, 1.0)
        
        return similarity_matrix.astype(np.float32)
    
    def find_top_k_similar(
        self,
        query_embedding: np.ndarray,
        target_embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]],
        top_k: int = 5,
        similarity_threshold: float = -1.0
    ) -> List[Tuple[str, float]]:
        """
        Find top-K most similar items to a query embedding.
        
        Args:
            query_embedding: Query embedding (1D array)
            target_embeddings: Target embeddings - can be:
                - 2D array [N, D]
                - List of 1D arrays
                - Dict mapping IDs to 1D arrays
            top_k: Number of top similar items to return
            similarity_threshold: Minimum similarity score threshold
        
        Returns:
            List of tuples (target_id, similarity_score) sorted by score descending
        """
        # Convert query to 2D array [1, D]
        if query_embedding.ndim == 1:
            query_array = query_embedding.reshape(1, -1)
        elif query_embedding.ndim == 2:
            if query_embedding.shape[0] != 1:
                raise ValueError("Query embedding must be 1D or 2D with shape [1, D]")
            query_array = query_embedding
        else:
            raise ValueError(f"Invalid query embedding shape: {query_embedding.shape}")
        
        # Convert targets to array with IDs
        target_array, target_ids = self._to_array_with_ids(target_embeddings)
        
        if target_array.ndim != 2:
            raise ValueError(f"Target embeddings must be 2D array, got {target_array.ndim}D")
        
        if query_array.shape[1] != target_array.shape[1]:
            raise ValueError(
                f"Embedding dimensions must match: query {query_array.shape[1]} vs "
                f"target {target_array.shape[1]}"
            )
        
        # Compute similarities
        similarities = self.compute_similarity_matrix(
            query_array,
            target_array,
            normalize=True
        )
        
        # Extract similarities (query is single row)
        similarity_scores = similarities[0]
        
        # Filter by threshold and get top-K
        valid_indices = np.where(similarity_scores >= similarity_threshold)[0]
        valid_scores = similarity_scores[valid_indices]
        
        # Get top-K indices
        top_indices = np.argsort(valid_scores)[::-1][:top_k]
        
        # Build results
        results = []
        for idx in top_indices:
            original_idx = valid_indices[idx]
            target_id = target_ids[original_idx]
            score = float(valid_scores[idx])
            results.append((target_id, score))
        
        return results
    
    def batch_find_top_k_similar(
        self,
        query_embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]],
        target_embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]],
        top_k: int = 5,
        similarity_threshold: float = -1.0
    ) -> Dict[str, List[Tuple[str, float]]]:
        """
        Find top-K similar items for multiple queries.
        
        Args:
            query_embeddings: Query embeddings - can be:
                - 2D array [N_query, D]
                - List of 1D arrays
                - Dict mapping IDs to 1D arrays
            target_embeddings: Target embeddings - can be:
                - 2D array [N_target, D]
                - List of 1D arrays
                - Dict mapping IDs to 1D arrays
            top_k: Number of top similar items to return per query
            similarity_threshold: Minimum similarity score threshold
        
        Returns:
            Dict mapping query_id to list of (target_id, similarity_score) tuples
        """
        # Convert to arrays with IDs
        query_array, query_ids = self._to_array_with_ids(query_embeddings)
        target_array, target_ids = self._to_array_with_ids(target_embeddings)
        
        # Compute similarity matrix
        similarity_matrix = self.compute_similarity_matrix(
            query_array,
            target_array,
            normalize=True
        )
        
        # Process each query
        results = {}
        for i, query_id in enumerate(query_ids):
            similarities = similarity_matrix[i]
            
            # Filter by threshold
            valid_indices = np.where(similarities >= similarity_threshold)[0]
            valid_scores = similarities[valid_indices]
            
            # Get top-K
            top_indices = np.argsort(valid_scores)[::-1][:top_k]
            
            # Build results for this query
            query_results = []
            for idx in top_indices:
                original_idx = valid_indices[idx]
                target_id = target_ids[original_idx]
                score = float(valid_scores[idx])
                query_results.append((target_id, score))
            
            results[query_id] = query_results
        
        return results
    
    def find_similar_within_queries(
        self,
        query_embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]],
        top_k: int = 5,
        similarity_threshold: float = -1.0,
        exclude_self: bool = True
    ) -> Dict[str, List[Tuple[str, float]]]:
        """
        Find similar items within a set of query embeddings.
        
        This is useful for finding similar queries or clustering queries.
        
        Args:
            query_embeddings: Query embeddings - can be:
                - 2D array [N, D]
                - List of 1D arrays
                - Dict mapping IDs to 1D arrays
            top_k: Number of top similar items to return per query
            similarity_threshold: Minimum similarity score threshold
            exclude_self: Whether to exclude each query from its own results
        
        Returns:
            Dict mapping query_id to list of (other_query_id, similarity_score) tuples
        """
        # Convert to arrays with IDs
        query_array, query_ids = self._to_array_with_ids(query_embeddings)
        
        # Compute similarity matrix (query to query)
        similarity_matrix = self.compute_similarity_matrix(
            query_array,
            query_array,
            normalize=True
        )
        
        # Process each query
        results = {}
        for i, query_id in enumerate(query_ids):
            similarities = similarity_matrix[i].copy()
            
            # Exclude self if requested
            if exclude_self:
                similarities[i] = -2.0  # Below any reasonable threshold
            
            # Filter by threshold
            valid_indices = np.where(similarities >= similarity_threshold)[0]
            valid_scores = similarities[valid_indices]
            
            # Get top-K
            top_indices = np.argsort(valid_scores)[::-1][:top_k]
            
            # Build results
            query_results = []
            for idx in top_indices:
                original_idx = valid_indices[idx]
                other_query_id = query_ids[original_idx]
                score = float(valid_scores[idx])
                query_results.append((other_query_id, score))
            
            results[query_id] = query_results
        
        return results
    
    def find_top_k_with_faiss(
        self,
        query_embedding: np.ndarray,
        index_path: Optional[str] = None,
        top_k: int = 5,
        faiss_manager: Optional[Any] = None
    ) -> List[Tuple[str, float]]:
        """
        Find top-K similar items using FAISS index.
        
        Args:
            query_embedding: Query embedding (1D array)
            index_path: Path to the FAISS index (without extension). Required if faiss_manager is None.
            top_k: Number of top similar items to return
            faiss_manager: Optional pre-loaded FaissIndexManager instance
            
        Returns:
            List of tuples (target_id, distance/score)
        """
        if FaissIndexManager is None:
            raise ImportError("FaissIndexManager is not available")
            
        if faiss_manager is None:
            if not index_path:
                raise ValueError("Must provide either index_path or faiss_manager")
            # Load index
            faiss_manager = FaissIndexManager.load_index(index_path)
            
        # Ensure query is 2D [1, D]
        if query_embedding.ndim == 1:
            query_array = query_embedding.reshape(1, -1)
        else:
            query_array = query_embedding
            
        # Search
        distances, indices, external_ids = faiss_manager.search(query_array, k=top_k)
        
        # Format results
        results = []
        # search returns list of lists, we take the first one
        row_ids = external_ids[0]
        row_dists = distances[0]
        
        for ext_id, dist in zip(row_ids, row_dists):
            if ext_id is not None:
                results.append((ext_id, float(dist)))
                
        return results

    def batch_find_top_k_with_faiss(
        self,
        query_embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]],
        index_path: Optional[str] = None,
        top_k: int = 5,
        faiss_manager: Optional[Any] = None
    ) -> Dict[str, List[Tuple[str, float]]]:
        """
        Batch find top-K similar items using FAISS index.
        
        Args:
            query_embeddings: Query embeddings (array, list, or dict)
            index_path: Path to the FAISS index (without extension). Required if faiss_manager is None.
            top_k: Number of top similar items to return
            faiss_manager: Optional pre-loaded FaissIndexManager instance
            
        Returns:
            Dict mapping query_id to list of (target_id, distance) tuples
        """
        if FaissIndexManager is None:
            raise ImportError("FaissIndexManager is not available")
            
        if faiss_manager is None:
            if not index_path:
                raise ValueError("Must provide either index_path or faiss_manager")
            faiss_manager = FaissIndexManager.load_index(index_path)
            
        # Convert queries to array
        query_array, query_ids = self._to_array_with_ids(query_embeddings)
        
        # Search
        distances, indices, external_ids = faiss_manager.search(query_array, k=top_k)
        
        # Format results
        results = {}
        for i, query_id in enumerate(query_ids):
            row_ids = external_ids[i]
            row_dists = distances[i]
            
            query_results = []
            for ext_id, dist in zip(row_ids, row_dists):
                if ext_id is not None:
                    query_results.append((ext_id, float(dist)))
            
            results[query_id] = query_results
            
        return results
    
    def _to_array_with_ids(
        self,
        embeddings: Union[np.ndarray, List[np.ndarray], Dict[str, np.ndarray]]
    ) -> Tuple[np.ndarray, List[str]]:
        """
        Convert various embedding formats to array and ID list.
        
        Args:
            embeddings: Embeddings in various formats
        
        Returns:
            Tuple of (embedding_array [N, D], id_list)
        """
        if isinstance(embeddings, np.ndarray):
            if embeddings.ndim == 1:
                # Single embedding - reshape to [1, D]
                return embeddings.reshape(1, -1), ["item_0"]
            elif embeddings.ndim == 2:
                # Multiple embeddings
                ids = [f"item_{i}" for i in range(embeddings.shape[0])]
                return embeddings, ids
            else:
                raise ValueError(f"Invalid embedding array shape: {embeddings.shape}")
        
        elif isinstance(embeddings, list):
            if not embeddings:
                raise ValueError("Empty embeddings list")
            
            # Convert list to array
            arrays = [np.asarray(emb, dtype=np.float32) for emb in embeddings]
            dim = arrays[0].shape[0]
            
            # Check all have same dimension
            for i, arr in enumerate(arrays):
                if arr.ndim != 1:
                    raise ValueError(f"List item {i} is not 1D array: {arr.shape}")
                if arr.shape[0] != dim:
                    raise ValueError(
                        f"Inconsistent embedding dimensions: item {i} has {arr.shape[0]}, "
                        f"expected {dim}"
                    )
            
            embedding_array = np.vstack(arrays)
            ids = [f"item_{i}" for i in range(len(arrays))]
            return embedding_array, ids
        
        elif isinstance(embeddings, dict):
            if not embeddings:
                raise ValueError("Empty embeddings dictionary")
            
            # Extract IDs and arrays
            ids = list(embeddings.keys())
            
            # Check if values are arrays or dicts with 'embedding' key
            first_val = next(iter(embeddings.values()))
            if isinstance(first_val, dict) and 'embedding' in first_val:
                # Handle nested dict structure: {'id': {'embedding': array, ...}}
                arrays = []
                for key in ids:
                    val = embeddings[key]
                    if 'embedding' not in val:
                        raise ValueError(f"Item {key} missing 'embedding' key")
                    arrays.append(np.asarray(val['embedding'], dtype=np.float32))
            else:
                # Handle direct mapping: {'id': array}
                arrays = [np.asarray(emb, dtype=np.float32) for emb in embeddings.values()]
            
            dim = arrays[0].shape[0]
            
            # Check all have same dimension
            for i, arr in enumerate(arrays):
                if arr.ndim != 1:
                    raise ValueError(f"Dict item {ids[i]} is not 1D array: {arr.shape}")
                if arr.shape[0] != dim:
                    raise ValueError(
                        f"Inconsistent embedding dimensions: {ids[i]} has {arr.shape[0]}, "
                        f"expected {dim}"
                    )
            
            embedding_array = np.vstack(arrays)
            return embedding_array, ids
        
        else:
            raise TypeError(
                f"Unsupported embeddings type: {type(embeddings)}. "
                f"Expected np.ndarray, list, or dict"
            )


def main() -> int:
    """Simple self-test for SimilaritySearch."""
    try:
        # Create test embeddings
        search = SimilaritySearch()
        
        # Test single embedding similarity
        emb1 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        emb2 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        sim = search.cosine_similarity(emb1, emb2)
        if abs(sim - 1.0) > 1e-5:
            raise RuntimeError(f"Expected similarity 1.0 for identical vectors, got {sim}")
        
        # Test top-K search
        query = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        targets = {
            'target_1': np.array([1.0, 0.0, 0.0], dtype=np.float32),
            'target_2': np.array([0.0, 1.0, 0.0], dtype=np.float32),
            'target_3': np.array([0.0, 0.0, 1.0], dtype=np.float32),
        }
        
        results = search.find_top_k_similar(query, targets, top_k=2)
        if len(results) != 2:
            raise RuntimeError(f"Expected 2 results, got {len(results)}")
        if results[0][0] != 'target_1':
            raise RuntimeError(f"Expected 'target_1' as top result, got {results[0][0]}")
        
        # Test query-to-query similarity
        queries = {
            'query_1': np.array([1.0, 0.0, 0.0], dtype=np.float32),
            'query_2': np.array([0.0, 1.0, 0.0], dtype=np.float32),
            'query_3': np.array([0.0, 0.0, 1.0], dtype=np.float32),
        }
        
        query_results = search.find_similar_within_queries(queries, top_k=2, exclude_self=True)
        if 'query_1' not in query_results:
            raise RuntimeError("Missing query_1 in results")
        if len(query_results['query_1']) != 2:
            raise RuntimeError(f"Expected 2 similar queries, got {len(query_results['query_1'])}")
            
        # Test nested dictionary support (Storage compatibility)
        nested_targets = {
            'p1': {'embedding': np.array([1.0, 0.0, 0.0], dtype=np.float32), 'other_data': 'foo'},
            'p2': {'embedding': np.array([0.0, 1.0, 0.0], dtype=np.float32), 'other_data': 'bar'},
        }
        nested_results = search.find_top_k_similar(query, nested_targets, top_k=1)
        if nested_results[0][0] != 'p1':
            raise RuntimeError(f"Nested dict test failed, expected 'p1', got {nested_results[0][0]}")
        
        print("SIMILARITY_SEARCH_OK: All tests passed")
        return 0
        
    except Exception as e:
        print(f"SIMILARITY_SEARCH_FAIL: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

