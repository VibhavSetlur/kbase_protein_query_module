"""
Analysis components for KBase Protein Query Module

Exposes only the new scalable analyzers.
"""

from .sequence_analysis import SequenceAnalysisStage

__all__ = [
    'SequenceAnalysisStage'
]
