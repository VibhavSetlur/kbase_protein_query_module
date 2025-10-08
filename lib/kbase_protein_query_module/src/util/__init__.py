"""
Utility module for KBase Protein Query Module.

This module provides various utility functions organized by category including
embeddings, similarity calculations, storage handling, family assignment,
similarity search, and other common utilities.
"""

from .embeddings import ProteinEmbeddingGenerator
from .similarity_indexing import HierarchicalIndex, StreamingIndex
from .storage import (
    ProteinStorage, MemoryEfficientLoader, 
    ProteinFamilyAssigner, ProteinExistenceChecker, IndexingStrategy
)
from .family_assignment import FamilyAssignment
from .similarity_search import SimilaritySearch

__all__ = [
    # Embeddings
    'ProteinEmbeddingGenerator',
    
    # Similarity
    'HierarchicalIndex',
    'StreamingIndex',
    
    # Storage
    'ProteinStorage',
    'MemoryEfficientLoader',
    'ProteinFamilyAssigner',
    'ProteinExistenceChecker',
    'IndexingStrategy',
    
    # Family Assignment
    'FamilyAssignment',
    
    # Similarity Search
    'SimilaritySearch'
]
