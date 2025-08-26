"""
Similarity Module for KBase Protein Query Module

This module contains similarity search and indexing functionality.
"""

from .hierarchical_index import HierarchicalIndex, StreamingIndex

# Compatibility alias expected by some tests
class SimilarityIndex(HierarchicalIndex):
    pass

__all__ = ['HierarchicalIndex', 'StreamingIndex', 'SimilarityIndex']
