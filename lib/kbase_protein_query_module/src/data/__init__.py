"""
Data module for KBase Protein Query Module.

This module provides data loading and management functionality for reference data
used in protein analysis, including amino acid properties, sequence motifs,
bioinformatics databases, and physicochemical constants.
"""

from .reference_loader import ReferenceDataLoader

__all__ = ['ReferenceDataLoader']
