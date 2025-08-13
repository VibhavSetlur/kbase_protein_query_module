"""
Analysis Stages for KBase Protein Query Module

This module contains analysis stages that handle specialized analysis types.
"""

from .sequence_analysis import SequenceAnalysisStage
from .network_analysis import NetworkAnalysisStage
from .bioinformatics_analysis import BioinformaticsAnalysisStage
from .multi_protein_analysis import MultiProteinAnalysisStage

__all__ = [
    'SequenceAnalysisStage',
    'NetworkAnalysisStage',
    'BioinformaticsAnalysisStage',
    'MultiProteinAnalysisStage'
]
