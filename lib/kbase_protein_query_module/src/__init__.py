"""
KBase Protein Query Module - Main Source Package

This package provides a comprehensive protein query and analysis system for KBase,
organized into modular components for input handling, analysis execution, output
generation, and utility functions.
"""

# Import main components
from .input import DataExtractionStage, InputValidationStage, WorkspaceObjectStage
from .analysis import AnalysisManager, get_enabled_analyses
from .outputs import (
    OutputManager, ArtifactRecord,
    get_output_config, get_analysis_output_config,
    get_enabled_output_analyses, is_output_enabled_for_analysis
)
from .core import WorkflowOrchestrator, WorkflowResult, PipelineConfig
from .util import (
    ProteinEmbeddingGenerator,
    HierarchicalIndex, StreamingIndex,
    ProteinStorage, MemoryEfficientLoader,
    ProteinFamilyAssigner, ProteinExistenceChecker,
    IndexingStrategy
)

__version__ = "1.0.0"
__author__ = "KBase Team"

__all__ = [
    # Input handling
    'DataExtractionStage',
    'InputValidationStage', 
    'WorkspaceObjectStage',
    
    # Analysis management
    'AnalysisManager',
    'get_enabled_analyses',
    
    # Output management
    'OutputManager',
    'ArtifactRecord',
    'get_output_config',
    'get_analysis_output_config',
    'get_enabled_output_analyses',
    'is_output_enabled_for_analysis',
    
    # Core workflow
    'WorkflowOrchestrator',
    'WorkflowResult',
    'PipelineConfig',
    
    # Utilities
    'ProteinEmbeddingGenerator',
    'HierarchicalIndex',
    'StreamingIndex',
    'ProteinStorage',
    'MemoryEfficientLoader',
    'ProteinFamilyAssigner',
    'ProteinExistenceChecker',
    'IndexingStrategy'
]