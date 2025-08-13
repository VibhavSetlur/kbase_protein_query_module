"""
Analysis Module for KBase Protein Query Module

This module contains comprehensive protein analysis capabilities.
"""

from .sequence_analyzer import ProteinSequenceAnalyzer
from .network_analyzer import NetworkAnalyzer, create_network_analyzer

# Aliases for backward compatibility
SequenceAnalyzer = ProteinSequenceAnalyzer

__all__ = [
    'ProteinSequenceAnalyzer',
    'SequenceAnalyzer',
    'NetworkAnalyzer',
    'create_network_analyzer'
]
