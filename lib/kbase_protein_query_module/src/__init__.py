"""
KBase Protein Query Module - Source Package

This package provides comprehensive protein query analysis capabilities for KBase.
The module is organized into logical components for maintainability and scalability.

Module Organization:
- core: Core abstractions, configurations, and base classes
- data: Data models and reference data
- processing: All processing logic (embeddings, similarity, networks)
- analysis: Analysis components and algorithms
- storage: Data storage and persistence layer
- reports: Report generation and visualization
- stages: Pipeline stages (input, processing, output)
- utils: Utility functions and helpers
- workflows: Workflow orchestration and management
"""

# Core module exports
from .core import BaseStage, StageResult, PipelineConfig

# Data module exports
from .data import ReferenceDataLoader

# Processing module exports
from .processing import (
    ProteinEmbeddingGenerator,
    HierarchicalIndex,
    DynamicNetworkBuilder,
    generate_embeddings_from_fasta,
    create_optimized_indexes,
    visualize_interactive_protein_network
)

# Analysis module exports
from .analysis import ProteinSequenceAnalyzer

# Storage module exports
from .storage import (
    ProteinStorage,
    ProteinFamilyAssigner,
    ProteinExistenceChecker,
    MemoryEfficientLoader,
    CompressedMetadataStorage,
    ProteinIDsIndex,
    create_storage_structure
)

# Reports module exports
from .reports import HTMLReportGenerator

# Stages module exports
from .stages import (
    InputValidationStage,
    DataExtractionStage,
    WorkspaceObjectStage,
    EmbeddingGenerationStage,
    FamilyAssignmentStage,
    SimilaritySearchStage,
    SequenceAnalysisStage,
    NetworkAnalysisStage,
    ReportGenerationStage,
    VisualizationStage,
    DataExportStage
)

# Utils module exports
from .utils import InputParser, ProteinRecord

# Workflows module exports
from .workflows import ProteinQueryWorkflow, WorkflowResult

__all__ = [
    # Core
    'BaseStage',
    'StageResult', 
    'PipelineConfig',
    
    # Data
    'ReferenceDataLoader',
    
    # Processing
    'ProteinEmbeddingGenerator',
    'HierarchicalIndex',
    'DynamicNetworkBuilder',
    'generate_embeddings_from_fasta',
    'create_optimized_indexes',
    'visualize_interactive_protein_network',
    
    # Analysis
    'ProteinSequenceAnalyzer',
    
    # Storage
    'ProteinStorage',
    'ProteinFamilyAssigner',
    'ProteinExistenceChecker',
    'MemoryEfficientLoader',
    'CompressedMetadataStorage',
    'ProteinIDsIndex',
    'create_storage_structure',
    
    # Reports
    'HTMLReportGenerator',
    
    # Stages
    'InputValidationStage',
    'DataExtractionStage',
    'WorkspaceObjectStage',
    'EmbeddingGenerationStage',
    'FamilyAssignmentStage',
    'SimilaritySearchStage',
    'SequenceAnalysisStage',
    'NetworkAnalysisStage',
    'ReportGenerationStage',
    'VisualizationStage',
    'DataExportStage',
    
    # Utils
    'InputParser',
    'ProteinRecord',
    
    # Workflows
    'ProteinQueryWorkflow',
    'WorkflowResult'
] 