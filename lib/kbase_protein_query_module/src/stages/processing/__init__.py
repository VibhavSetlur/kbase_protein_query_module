"""
Processing Stages for KBase Protein Query Module

This module contains processing stages that handle core analysis and computation.
"""

from .embedding_generation import EmbeddingGenerationStage
from .family_assignment import FamilyAssignmentStage
from .similarity_search import SimilaritySearchStage

__all__ = [
    'EmbeddingGenerationStage',
    'FamilyAssignmentStage',
    'SimilaritySearchStage'
]
