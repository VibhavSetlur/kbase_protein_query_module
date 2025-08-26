"""
KBase Protein Query Module

A comprehensive protein analysis pipeline that performs sequence analysis,
embedding generation, similarity search, family assignment, network analysis,
and bioinformatics analysis using state-of-the-art machine learning models.
"""

# Core modules
from .core import BaseStage, StageResult, PipelineConfig

# Stage modules
from .stages import (
    STAGE_REGISTRY, STAGE_DEPENDENCIES, get_stage_class, get_stage_dependencies,
    InputValidationStage, DataExtractionStage, WorkspaceObjectStage,
    EmbeddingGenerationStage, FamilyAssignmentStage, SimilaritySearchStage,
    SequenceAnalysisStage, NetworkAnalysisStage, BioinformaticsAnalysisStage,
    ReportGenerationStage, VisualizationStage, DataExportStage
)

# Workflow modules
from .workflows import ProteinQueryWorkflow, WorkflowResult

# Processing modules
from .processing.embeddings.generator import ProteinEmbeddingGenerator
from .processing.similarity.hierarchical_index import HierarchicalIndex
from .processing.networks.builder import DynamicNetworkBuilder

# Storage modules
from .storage import ProteinStorage, ProteinFamilyAssigner, ProteinExistenceChecker

# Utility modules
from .utils import input_parser

__version__ = "2.0.0"
__author__ = "KBase Team"
__email__ = "support@kbase.us"

__all__ = [
    # Core
    'BaseStage',
    'StageResult', 
    'PipelineConfig',
    
    # Stages
    'STAGE_REGISTRY',
    'STAGE_DEPENDENCIES',
    'get_stage_class',
    'get_stage_dependencies',
    'InputValidationStage',
    'DataExtractionStage',
    'WorkspaceObjectStage',
    'EmbeddingGenerationStage',
    'FamilyAssignmentStage',
    'SimilaritySearchStage',
    'SequenceAnalysisStage',
    'NetworkAnalysisStage',
    'BioinformaticsAnalysisStage',
    'ReportGenerationStage',
    'VisualizationStage',
    'DataExportStage',
    
    # Workflows
    'ProteinQueryWorkflow',
    'WorkflowResult',
    
    # Processing
    'ProteinEmbeddingGenerator',
    'HierarchicalIndex',
    'DynamicNetworkBuilder',
    
    # Storage
    'ProteinStorage',
    'ProteinFamilyAssigner',
    'ProteinExistenceChecker',
    
    # Utils
    'input_parser'
] 