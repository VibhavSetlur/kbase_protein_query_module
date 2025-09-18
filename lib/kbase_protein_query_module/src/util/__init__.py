"""
Utility module for KBase Protein Query Module.

This module provides various utility functions organized by category including
embeddings, similarity calculations, storage handling, and other common utilities.
"""

from .embeddings import *
from .similarity import *
from .storage import *
from .visualization import *
from .documentation import *
from .input_parser import *
from .safe_delete import *

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
    
    # Visualization
    'VisualizationConverter',
    
    # Documentation
    'DocumentationGenerator',
    
    # Input/Output
    'InputParser',
    'SafeDelete'
]
