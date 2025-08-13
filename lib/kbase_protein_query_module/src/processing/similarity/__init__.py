"""
Similarity Module for KBase Protein Query Module

This module contains similarity search and indexing functionality.
"""

from .hierarchical_index import HierarchicalIndex, StreamingIndex

__all__ = ['HierarchicalIndex', 'StreamingIndex']
