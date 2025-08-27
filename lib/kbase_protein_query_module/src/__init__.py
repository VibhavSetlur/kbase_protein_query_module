"""
KBase Protein Query Module

A comprehensive protein analysis pipeline that performs sequence analysis,
embedding generation, similarity search, family assignment, network analysis,
and bioinformatics analysis using state-of-the-art machine learning models.
"""

# Core modules
from .core import BaseStage, StageResult, PipelineConfig
# Scalability and extensibility modules
from .core.analysis_registry import BaseAnalysis, AnalysisRegistry, get_registry, register_analysis
from .core.resource_manager import ResourceManager, ResourceLimits, get_resource_manager
from .core.parallel_processor import ParallelProcessor
from .core.performance_monitor import PerformanceProfiler, get_performance_profiler

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
# Storage extensibility modules
from .storage.indexing_strategy import IndexingStrategy, IndexingConfig, get_indexing_registry, register_indexing_strategy

# Utility modules
from .utils import input_parser
from .utils.documentation_generator import generate_full_documentation

__version__ = "2.0.0"
__author__ = "KBase Team"
__email__ = "support@kbase.us"

__all__ = [
    # Core
    'BaseStage',
    'StageResult', 
    'PipelineConfig',
    # Extensibility Framework
    'BaseAnalysis',
    'AnalysisRegistry', 
    'get_registry',
    'register_analysis',
    'ResourceManager',
    'ResourceLimits',
    'get_resource_manager',
    'ParallelProcessor',
    'PerformanceProfiler',
    'get_performance_profiler',
    # Storage Extensibility
    'IndexingStrategy',
    'IndexingConfig',
    'get_indexing_registry', 
    'register_indexing_strategy',
    # Documentation
    'generate_full_documentation',
    
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