"""
Input handling module for KBase Protein Query Module.

This module provides comprehensive input handling including data extraction,
validation, and workspace object management for KBase narrative integration.
"""

from .input_manager import InputManager
from .workspace_object import WorkspaceObjectProcessor
from .protein_sequence import ProteinSequenceProcessor
from .uniprot_ids import UniProtIdsProcessor

__all__ = [
    'InputManager',
    'WorkspaceObjectProcessor',
    'ProteinSequenceProcessor',
    'UniProtIdsProcessor'
]
