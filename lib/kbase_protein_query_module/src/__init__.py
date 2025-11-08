"""
KBase Protein Query Module - Main Source Package

This package provides a comprehensive protein query and analysis system for KBase,
organized into modular components for input handling, analysis execution, output
generation, and utility functions.
"""

# Import main components
from .input import InputManager, ProteinSequenceProcessor, UniProtIdProcessor
from .analysis import AnalysisManager, get_enabled_analyses
from .output import OutputManager
from .core import WorkflowOrchestrator

# Utilities are imported directly from submodules where needed
# No need to import them here - keeps things simple

__version__ = "1.0.0"
__author__ = "KBase Team"

__all__ = [
    # Input handling
    'InputManager',
    'ProteinSequenceProcessor',
    'UniProtIdProcessor',
    
    # Analysis management
    'AnalysisManager',
    'get_enabled_analyses',
    
    # Output management
    'OutputManager',
    
    # Core workflow
    'WorkflowOrchestrator',
]