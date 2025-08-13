"""
Processing Module Tests

This module contains comprehensive tests for all processing functionality:
- Embedding generation using ESM-2 models
- Similarity search and indexing
- Network construction and analysis
- Data transformation and processing pipelines
"""

from .embeddings.test_embedding_generator import TestProteinEmbeddingGenerator
from .similarity.test_similarity_index import TestSimilarityIndex
from .networks.test_network_builder import TestDynamicNetworkBuilder

__all__ = [
    'TestProteinEmbeddingGenerator',
    'TestSimilarityIndex',
    'TestDynamicNetworkBuilder'
]
