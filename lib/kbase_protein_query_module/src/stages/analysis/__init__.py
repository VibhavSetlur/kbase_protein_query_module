"""
Analysis components for KBase Protein Query Module

This module imports all available analysis implementations to ensure
they are registered with the analysis registry for extensibility.
"""

# Import all analysis implementations to trigger registration
from .sequence_analysis import SequenceAnalysisStage
from .network_analysis import NetworkAnalysisStage  
from .bioinformatics_analysis import BioinformaticsAnalysisStage
from .multi_protein_analysis import MultiProteinAnalysisStage

# Import demo analysis to demonstrate extensibility
# from .example_new_analysis import HydrophobicityAnalysis

__all__ = [
    'SequenceAnalysisStage',
    'NetworkAnalysisStage',
    'BioinformaticsAnalysisStage', 
    'MultiProteinAnalysisStage'
    # 'HydrophobicityAnalysis'
]
