"""
Outputs module for KBase Protein Query Module.

This module provides comprehensive output management including analysis-specific
output handlers and organized output directory structures.
"""

from .output_manager import OutputManager, ArtifactRecord
from .config import (
    get_output_config,
    get_analysis_output_config,
    get_enabled_output_analyses,
    is_output_enabled_for_analysis,
    get_output_formats_for_analysis,
    get_file_naming_config,
    get_directory_config,
    validate_output_config
)

__all__ = [
    'OutputManager',
    'ArtifactRecord',
    'get_output_config',
    'get_analysis_output_config',
    'get_enabled_output_analyses',
    'is_output_enabled_for_analysis',
    'get_output_formats_for_analysis',
    'get_file_naming_config',
    'get_directory_config',
    'validate_output_config'
]