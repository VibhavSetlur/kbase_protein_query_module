"""
Stages module for KBase Protein Query Module

This module contains all pipeline stages organized by functionality.
"""

# Import from reorganized submodules
from .input import InputValidationStage, DataExtractionStage, WorkspaceObjectStage
from .processing import EmbeddingGenerationStage, FamilyAssignmentStage, SimilaritySearchStage
from .analysis import SequenceAnalysisStage, NetworkAnalysisStage, BioinformaticsAnalysisStage, MultiProteinAnalysisStage
from .output import ReportGenerationStage, VisualizationStage, DataExportStage

# Import base stage
from .base_stage import BaseStage, StageResult

# Stage registry and dependencies
STAGE_REGISTRY = {
    'input_validation': InputValidationStage,
    'data_extraction': DataExtractionStage,
    'workspace_object': WorkspaceObjectStage,
    'embedding_generation': EmbeddingGenerationStage,
    'family_assignment': FamilyAssignmentStage,
    'similarity_search': SimilaritySearchStage,
    'sequence_analysis': SequenceAnalysisStage,
    'network_analysis': NetworkAnalysisStage,
    'bioinformatics_analysis': BioinformaticsAnalysisStage,
    'multi_protein_analysis': MultiProteinAnalysisStage,
    'report_generation': ReportGenerationStage,
    'visualization': VisualizationStage,
    'data_export': DataExportStage
}

STAGE_DEPENDENCIES = {
    'input_validation': [],
    'data_extraction': ['input_validation'],
    'workspace_object': ['input_validation'],
    'embedding_generation': ['data_extraction', 'workspace_object'],
    'family_assignment': ['embedding_generation'],
    'similarity_search': ['embedding_generation'],
    'sequence_analysis': ['data_extraction', 'workspace_object'],
    'network_analysis': ['similarity_search', 'family_assignment'],
    'bioinformatics_analysis': ['sequence_analysis'],
    'multi_protein_analysis': ['sequence_analysis', 'data_extraction'],
    'report_generation': ['network_analysis', 'bioinformatics_analysis', 'multi_protein_analysis'],
    'visualization': ['network_analysis', 'multi_protein_analysis'],
    'data_export': ['report_generation', 'visualization']
}

def get_stage_class(stage_name: str):
    """Get stage class by name."""
    if stage_name not in STAGE_REGISTRY:
        raise ValueError(f"Unknown stage: {stage_name}")
    return STAGE_REGISTRY[stage_name]

def get_stage_dependencies(stage_name: str) -> list:
    """Get dependencies for a stage."""
    return STAGE_DEPENDENCIES.get(stage_name, [])

def get_all_stages() -> list:
    """Get all available stage names."""
    return list(STAGE_REGISTRY.keys())

__all__ = [
    # Base classes
    'BaseStage',
    'StageResult',
    
    # Input stages
    'InputValidationStage',
    'DataExtractionStage', 
    'WorkspaceObjectStage',
    
    # Processing stages
    'EmbeddingGenerationStage',
    'FamilyAssignmentStage',
    'SimilaritySearchStage',
    
    # Analysis stages
    'SequenceAnalysisStage',
    'NetworkAnalysisStage',
    'BioinformaticsAnalysisStage',
    'MultiProteinAnalysisStage',
    
    # Output stages
    'ReportGenerationStage',
    'VisualizationStage',
    'DataExportStage',
    
    # Registry and utilities
    'STAGE_REGISTRY',
    'STAGE_DEPENDENCIES',
    'get_stage_class',
    'get_stage_dependencies',
    'get_all_stages'
]
