"""
Input handling module for KBase Protein Query Module.

This module provides comprehensive input handling including data extraction,
validation, and workspace object management for KBase narrative integration.
"""

from .data_extraction import DataExtractionStage
from .input_validation import InputValidationStage
from .workspace_object import WorkspaceObjectStage

__all__ = [
    'DataExtractionStage',
    'InputValidationStage', 
    'WorkspaceObjectStage'
]
