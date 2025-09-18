"""
Flexible Indexing Strategy Framework

This module provides an extensible framework for implementing different indexing
strategies for protein similarity search. Developers can easily add new indexing
methods by implementing the IndexingStrategy interface.

DEVELOPER GUIDE FOR NEW INDEXING STRATEGIES:
============================================

To implement a new indexing strategy:

1. Create a class inheriting from IndexingStrategy
2. Implement required methods: build_index(), search(), get_index_info()
3. Register your strategy using @register_indexing_strategy decorator
4. Update configuration to use your new strategy

Example Implementation:
----------------------
@register_indexing_strategy("my_custom_index")
class MyCustomIndexingStrategy(IndexingStrategy):
    def build_index(self, embeddings: np.ndarray, metadata: Dict[str, Any]) -> Any:
        # Build your custom index
        return custom_index
    
    def search(self, query: np.ndarray, k: int = 10) -> List[Tuple[int, float]]:
        # Implement search logic
        return [(idx, score), ...]

Classes and Their Locations:
---------------------------
- IndexingStrategy: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:65
- IndexingRegistry: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:180
- FAISSIndexingStrategy: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:220
- AnnoyIndexingStrategy: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:320

Extension Points for New Systems:
---------------------------------
1. Custom distance metrics: Override _calculate_distance() in your strategy
2. Custom quantization: Override _apply_quantization() for memory optimization
3. Custom serialization: Override save_index() and load_index() for persistence
4. Custom metadata handling: Override _process_metadata() for additional features
"""

import os
import logging
import pickle
import json
import numpy as np
from typing import Dict, Any, List, Optional, Tuple, Type, Union
from abc import ABC, abstractmethod
from dataclasses import dataclass
import time

logger = logging.getLogger(__name__)

@dataclass
class IndexingConfig:
    """
    Configuration for indexing strategies.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:50
    USED BY: All indexing strategy implementations
    """
    strategy_name: str
    index_type: str = "faiss"  # "faiss", "annoy", "custom"
    distance_metric: str = "cosine"  # "cosine", "euclidean", "dot_product"
    quantization: str = "none"  # "none", "pq", "sq", "ivf"
    memory_map: bool = False
    use_gpu: bool = False
    batch_size: int = 1000
    cache_size_mb: int = 512
    custom_params: Dict[str, Any] = None

class IndexingStrategy(ABC):
    """
    Abstract base class for all indexing strategies.
    
    This class defines the interface that all indexing implementations must follow.
    New indexing systems should inherit from this class and implement all abstract methods.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:65
    USED BY: ProteinStorage, HierarchicalIndex
    EXTENDS: ABC (Python Abstract Base Class)
    
    IMPLEMENTATION REQUIREMENTS:
    ---------------------------
    1. build_index(): Create the index from embeddings
    2. search(): Perform similarity search
    3. get_index_info(): Return index metadata
    4. save_index(): Persist index to disk
    5. load_index(): Load index from disk
    
    EXTENSION POINTS:
    ----------------
    - Override _preprocess_embeddings() for custom embedding preprocessing
    - Override _calculate_distance() for custom distance metrics
    - Override _apply_quantization() for custom compression
    - Override _optimize_index() for post-build optimization
    """
    
    def __init__(self, config: IndexingConfig):
        self.config = config
        self.logger = logging.getLogger(f"{self.__class__.__module__}.{self.__class__.__name__}")
        self.index = None
        self.metadata = {}
        self.is_built = False
        
        # Performance tracking
        self.build_time = 0.0
        self.search_times: List[float] = []
        
    @abstractmethod
    def build_index(self, embeddings: np.ndarray, metadata: Dict[str, Any] = None) -> Any:
        """
        Build the index from protein embeddings.
        
        Args:
            embeddings: Protein embeddings array (N x D)
            metadata: Optional metadata for embeddings
            
        Returns:
            Built index object
            
        IMPLEMENTATION NOTES:
        - Call self._preprocess_embeddings() before building
        - Store build time in self.build_time
        - Set self.is_built = True when complete
        - Store any relevant metadata in self.metadata
        """
        pass
    
    @abstractmethod
    def search(self, query_embedding: np.ndarray, k: int = 10, **kwargs) -> List[Tuple[int, float]]:
        """
        Search for similar proteins using the index.
        
        Args:
            query_embedding: Query protein embedding (1 x D)
            k: Number of nearest neighbors to return
            **kwargs: Additional search parameters
            
        Returns:
            List of (index, similarity_score) tuples
            
        IMPLEMENTATION NOTES:
        - Track search time and add to self.search_times
        - Handle edge cases (empty index, invalid k, etc.)
        - Apply any post-processing filters
        """
        pass
    
    @abstractmethod
    def get_index_info(self) -> Dict[str, Any]:
        """
        Get information about the current index.
        
        Returns:
            Dictionary with index metadata and statistics
            
        SHOULD INCLUDE:
        - index_size: Number of indexed embeddings
        - memory_usage_mb: Approximate memory usage
        - build_time: Time taken to build index
        - avg_search_time: Average search time
        - strategy_name: Name of indexing strategy
        - distance_metric: Distance metric used
        """
        pass
    
    @abstractmethod
    def save_index(self, filepath: str) -> bool:
        """
        Save index to disk.
        
        Args:
            filepath: Path to save the index
            
        Returns:
            True if successful, False otherwise
        """
        pass
    
    @abstractmethod
    def load_index(self, filepath: str) -> bool:
        """
        Load index from disk.
        
        Args:
            filepath: Path to load the index from
            
        Returns:
            True if successful, False otherwise
        """
        pass
    
    # Extension points - override these for custom behavior
    
    def _preprocess_embeddings(self, embeddings: np.ndarray) -> np.ndarray:
        """
        Preprocess embeddings before indexing.
        
        EXTENSION POINT: Override this method to implement custom preprocessing
        such as normalization, dimensionality reduction, or filtering.
        
        Args:
            embeddings: Raw embeddings
            
        Returns:
            Processed embeddings
        """
        # Default: L2 normalize for cosine similarity
        if self.config.distance_metric == "cosine":
            norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
            norms[norms == 0] = 1  # Avoid division by zero
            embeddings = embeddings / norms
        
        return embeddings.astype(np.float32)  # Use float32 for memory efficiency
    
    def _calculate_distance(self, query: np.ndarray, candidates: np.ndarray) -> np.ndarray:
        """
        Calculate distance between query and candidate embeddings.
        
        EXTENSION POINT: Override this method to implement custom distance metrics.
        
        Args:
            query: Query embedding (1 x D)
            candidates: Candidate embeddings (N x D)
            
        Returns:
            Distance array (N,)
        """
        if self.config.distance_metric == "cosine":
            return 1 - np.dot(candidates, query.T).flatten()
        elif self.config.distance_metric == "euclidean":
            return np.linalg.norm(candidates - query, axis=1)
        elif self.config.distance_metric == "dot_product":
            return -np.dot(candidates, query.T).flatten()  # Negative for descending order
        else:
            raise ValueError(f"Unknown distance metric: {self.config.distance_metric}")
    
    def _apply_quantization(self, embeddings: np.ndarray) -> np.ndarray:
        """
        Apply quantization for memory optimization.
        
        EXTENSION POINT: Override this method to implement custom quantization schemes.
        
        Args:
            embeddings: Input embeddings
            
        Returns:
            Quantized embeddings
        """
        if self.config.quantization == "none":
            return embeddings
        elif self.config.quantization == "int8":
            # Simple int8 quantization
            min_val, max_val = embeddings.min(), embeddings.max()
            scale = (max_val - min_val) / 255
            quantized = ((embeddings - min_val) / scale).astype(np.uint8)
            self.metadata['quantization_scale'] = scale
            self.metadata['quantization_offset'] = min_val
            return quantized
        else:
            logger.warning(f"Unknown quantization method: {self.config.quantization}")
            return embeddings
    
    def _optimize_index(self) -> None:
        """
        Optimize index after building.
        
        EXTENSION POINT: Override this method to implement custom post-build optimizations
        such as index compression, memory mapping, or performance tuning.
        """
        pass

class IndexingRegistry:
    """
    Registry for managing available indexing strategies.
    
    This registry enables dynamic discovery and instantiation of indexing strategies,
    making it easy to add new indexing methods without modifying existing code.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:180
    USED BY: ProteinStorage, HierarchicalIndex
    SINGLETON PATTERN: Use get_indexing_registry() to access global instance
    """
    
    def __init__(self):
        self._strategies: Dict[str, Type[IndexingStrategy]] = {}
        self.logger = logging.getLogger(__name__)
        
        # Register built-in strategies
        self._register_builtin_strategies()
    
    def register_strategy(self, name: str, strategy_class: Type[IndexingStrategy]) -> None:
        """Register a new indexing strategy."""
        if not issubclass(strategy_class, IndexingStrategy):
            raise ValueError(f"Strategy class {strategy_class.__name__} must inherit from IndexingStrategy")
        
        self._strategies[name] = strategy_class
        self.logger.info(f"Registered indexing strategy: {name} ({strategy_class.__name__})")
    
    def get_strategy(self, name: str, config: IndexingConfig) -> IndexingStrategy:
        """Get an instance of a registered indexing strategy."""
        if name not in self._strategies:
            available = list(self._strategies.keys())
            raise ValueError(f"Indexing strategy '{name}' not found. Available: {available}")
        
        strategy_class = self._strategies[name]
        return strategy_class(config)
    
    def list_strategies(self) -> List[str]:
        """Get list of available indexing strategies."""
        return list(self._strategies.keys())
    
    def _register_builtin_strategies(self):
        """Register built-in indexing strategies."""
        # These will be registered when the classes are defined below
        pass

class FAISSIndexingStrategy(IndexingStrategy):
    """
    FAISS-based indexing strategy with advanced features.
    
    Supports various FAISS index types including IVF, PQ, and HNSW for different
    use cases and performance requirements.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:220
    EXTENDS: IndexingStrategy
    USED BY: HierarchicalIndex when strategy="faiss"
    
    CUSTOMIZATION POINTS:
    - Override _select_index_type() to choose different FAISS index types
    - Override _configure_index() to customize index parameters
    - Override _add_embeddings_batch() for custom batch processing
    """
    
    def __init__(self, config: IndexingConfig):
        super().__init__(config)
        
        # Try to import FAISS
        try:
            import faiss
            self.faiss = faiss
        except ImportError:
            raise ImportError("FAISS library is required for FAISSIndexingStrategy")
        
        # FAISS-specific configuration
        self.index_type = config.custom_params.get('index_type', 'IVFFlat') if config.custom_params else 'IVFFlat'
        self.nlist = config.custom_params.get('nlist', 100) if config.custom_params else 100
        self.m = config.custom_params.get('pq_m', 8) if config.custom_params else 8
        self.bits = config.custom_params.get('pq_bits', 8) if config.custom_params else 8
    
    def build_index(self, embeddings: np.ndarray, metadata: Dict[str, Any] = None) -> Any:
        """Build FAISS index from embeddings."""
        start_time = time.time()
        
        # Preprocess embeddings
        processed_embeddings = self._preprocess_embeddings(embeddings)
        
        # Select and configure index type
        dimension = processed_embeddings.shape[1]
        index = self._create_faiss_index(dimension, len(processed_embeddings))
        
        # Add embeddings to index
        if hasattr(index, 'train') and not index.is_trained:
            self.logger.info("Training FAISS index...")
            index.train(processed_embeddings)
        
        self.logger.info(f"Adding {len(processed_embeddings)} embeddings to FAISS index...")
        index.add(processed_embeddings)
        
        # Store index and metadata
        self.index = index
        self.metadata = metadata or {}
        self.metadata.update({
            'index_type': self.index_type,
            'dimension': dimension,
            'num_embeddings': len(processed_embeddings),
            'distance_metric': self.config.distance_metric
        })
        
        self.build_time = time.time() - start_time
        self.is_built = True
        
        self.logger.info(f"FAISS index built in {self.build_time:.2f}s")
        return self.index
    
    def _create_faiss_index(self, dimension: int, num_embeddings: int) -> Any:
        """
        Create appropriate FAISS index based on configuration.
        
        EXTENSION POINT: Override this method to implement custom FAISS index selection
        based on your specific requirements (dataset size, accuracy needs, etc.)
        """
        if self.index_type == "Flat":
            if self.config.distance_metric == "cosine":
                index = self.faiss.IndexFlatIP(dimension)  # Inner product for cosine
            else:
                index = self.faiss.IndexFlatL2(dimension)
                
        elif self.index_type == "IVFFlat":
            # Choose nlist based on dataset size
            nlist = min(self.nlist, max(1, int(np.sqrt(num_embeddings))))
            
            if self.config.distance_metric == "cosine":
                quantizer = self.faiss.IndexFlatIP(dimension)
                index = self.faiss.IndexIVFFlat(quantizer, dimension, nlist)
            else:
                quantizer = self.faiss.IndexFlatL2(dimension)
                index = self.faiss.IndexIVFFlat(quantizer, dimension, nlist)
                
        elif self.index_type == "IVFPQ":
            nlist = min(self.nlist, max(1, int(np.sqrt(num_embeddings))))
            
            if self.config.distance_metric == "cosine":
                quantizer = self.faiss.IndexFlatIP(dimension)
                index = self.faiss.IndexIVFPQ(quantizer, dimension, nlist, self.m, self.bits)
            else:
                quantizer = self.faiss.IndexFlatL2(dimension)
                index = self.faiss.IndexIVFPQ(quantizer, dimension, nlist, self.m, self.bits)
                
        elif self.index_type == "HNSW":
            if self.config.distance_metric == "cosine":
                index = self.faiss.IndexHNSWFlat(dimension, 32)
                index.metric_type = self.faiss.METRIC_INNER_PRODUCT
            else:
                index = self.faiss.IndexHNSWFlat(dimension, 32)
                
        else:
            raise ValueError(f"Unknown FAISS index type: {self.index_type}")
        
        return index
    
    def search(self, query_embedding: np.ndarray, k: int = 10, **kwargs) -> List[Tuple[int, float]]:
        """Search using FAISS index."""
        if not self.is_built or self.index is None:
            raise ValueError("Index not built. Call build_index() first.")
        
        start_time = time.time()
        
        # Preprocess query
        query = self._preprocess_embeddings(query_embedding.reshape(1, -1))
        
        # Perform search
        scores, indices = self.index.search(query, k)
        
        # Convert to list of tuples
        results = []
        for i, (idx, score) in enumerate(zip(indices[0], scores[0])):
            if idx != -1:  # Valid result
                # Convert similarity score based on distance metric
                if self.config.distance_metric == "cosine":
                    similarity = score  # FAISS returns similarity for IP
                else:
                    similarity = 1.0 / (1.0 + score)  # Convert distance to similarity
                
                results.append((int(idx), float(similarity)))
        
        # Track search time
        search_time = time.time() - start_time
        self.search_times.append(search_time)
        
        return results
    
    def get_index_info(self) -> Dict[str, Any]:
        """Get FAISS index information."""
        if not self.is_built:
            return {"status": "not_built"}
        
        return {
            "strategy": "faiss",
            "index_type": self.index_type,
            "dimension": self.metadata.get('dimension', 0),
            "num_embeddings": self.metadata.get('num_embeddings', 0),
            "distance_metric": self.config.distance_metric,
            "build_time": self.build_time,
            "avg_search_time": np.mean(self.search_times) if self.search_times else 0,
            "memory_usage_mb": self._estimate_memory_usage(),
            "is_trained": getattr(self.index, 'is_trained', True) if self.index else False
        }
    
    def save_index(self, filepath: str) -> bool:
        """Save FAISS index to disk."""
        if not self.is_built or self.index is None:
            return False
        
        try:
            # Save FAISS index
            self.faiss.write_index(self.index, f"{filepath}.faiss")
            
            # Save metadata
            with open(f"{filepath}.meta", 'w') as f:
                json.dump({
                    'config': self.config.__dict__,
                    'metadata': self.metadata,
                    'build_time': self.build_time
                }, f, indent=2)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to save FAISS index: {e}")
            return False
    
    def load_index(self, filepath: str) -> bool:
        """Load FAISS index from disk."""
        try:
            # Load FAISS index
            self.index = self.faiss.read_index(f"{filepath}.faiss")
            
            # Load metadata
            with open(f"{filepath}.meta", 'r') as f:
                data = json.load(f)
                self.metadata = data.get('metadata', {})
                self.build_time = data.get('build_time', 0)
            
            self.is_built = True
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to load FAISS index: {e}")
            return False
    
    def _estimate_memory_usage(self) -> float:
        """Estimate memory usage in MB."""
        if not self.index:
            return 0.0
        
        # Rough estimate based on index type and size
        num_embeddings = self.metadata.get('num_embeddings', 0)
        dimension = self.metadata.get('dimension', 0)
        
        if self.index_type == "Flat":
            return (num_embeddings * dimension * 4) / (1024 * 1024)  # float32
        elif self.index_type == "IVFPQ":
            return (num_embeddings * self.m * self.bits / 8) / (1024 * 1024)
        else:
            return (num_embeddings * dimension * 2) / (1024 * 1024)  # Rough estimate

# Global registry instance
_indexing_registry: Optional['IndexingRegistry'] = None

def get_indexing_registry():
    """Get the global indexing registry instance."""
    global _indexing_registry
    if _indexing_registry is None:
        _indexing_registry = IndexingRegistry()
    return _indexing_registry

def register_indexing_strategy(name: str):
    """
    Decorator for registering indexing strategies.
    
    Usage:
    @register_indexing_strategy("my_strategy")
    class MyIndexingStrategy(IndexingStrategy):
        ...
    """
    def decorator(cls: Type[IndexingStrategy]) -> Type[IndexingStrategy]:
        registry = get_indexing_registry()
        registry.register_strategy(name, cls)
        return cls
    return decorator

class DeprecatedIndexingRegistry:
    """Deprecated duplicate registry (kept for backward compatibility). Use get_indexing_registry()."""
    
    def __init__(self):
        self._strategies: Dict[str, Type[IndexingStrategy]] = {}
        self.logger = logging.getLogger(__name__)
    
    def register_strategy(self, name: str, strategy_class: Type[IndexingStrategy]) -> None:
        registry = get_indexing_registry()
        registry.register_strategy(name, strategy_class)
        self.logger.info(f"(Deprecated) Registered indexing strategy: {name}")
    
    def get_strategy(self, name: str, config: IndexingConfig) -> IndexingStrategy:
        registry = get_indexing_registry()
        return registry.get_strategy(name, config)
    
    def list_strategies(self) -> List[str]:
        registry = get_indexing_registry()
        return registry.list_strategies()

# Register built-in FAISS strategy
@register_indexing_strategy("faiss")
class RegisteredFAISSStrategy(FAISSIndexingStrategy):
    """Registered FAISS strategy for the registry."""
    pass
