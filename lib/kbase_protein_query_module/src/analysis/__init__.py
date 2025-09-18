"""
Analysis module for KBase Protein Query Module.

This module provides a comprehensive analysis framework with support for multiple
analysis types including network analysis, sequence analysis, and bioinformatics analysis.
"""

from .analysis_manager import AnalysisManager
from .config import (
    get_enabled_analyses,
    get_analysis_by_category,
    is_analysis_enabled,
    get_analysis_dependencies,
    validate_analysis_config
)

__all__ = [
    'AnalysisManager',
    'get_enabled_analyses',
    'get_analysis_by_category', 
    'is_analysis_enabled',
    'get_analysis_dependencies',
    'validate_analysis_config'
]
