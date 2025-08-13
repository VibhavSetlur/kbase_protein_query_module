"""
Pipeline Configuration for Protein Query Analysis

This module defines the configuration dataclass for the protein query analysis pipeline.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional


@dataclass
class PipelineConfig:
    """
    Configuration for the protein query analysis pipeline.
    
    This dataclass contains all configuration parameters needed to run
    the protein query analysis workflow.
    """
    
    # Input configuration
    input_proteins: List[str] = field(default_factory=list)
    input_file_path: Optional[str] = None
    workspace_object_ref: Optional[str] = None
    
    # Analysis configuration
    perform_embedding_generation: bool = True
    perform_family_assignment: bool = True
    perform_similarity_search: bool = True
    perform_network_analysis: bool = True
    perform_sequence_analysis: bool = True
    
    # Stage configuration
    enabled_stages: List[str] = field(default_factory=lambda: [
        'input_validation', 'data_extraction', 'embedding_generation', 
        'family_assignment', 'similarity_search', 'sequence_analysis', 
        'network_analysis', 'bioinformatics_analysis', 'report_generation', 
        'visualization', 'data_export'
    ])
    stage_configs: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Embedding configuration
    embedding_model: str = "esm2_t6_8M_UR50D"
    embedding_batch_size: int = 32
    embedding_device: str = "cpu"
    
    # Similarity search configuration
    similarity_threshold: float = 0.8
    max_similar_proteins: int = 100
    
    # Network analysis configuration
    network_min_edges: int = 5
    network_layout: str = "force_directed"
    
    # Output configuration
    generate_html_report: bool = True
    generate_network_visualization: bool = True
    output_format: str = "html"
    
    # Storage configuration
    cache_embeddings: bool = True
    cache_directory: str = "data/cache"
    storage_config: Dict[str, Any] = field(default_factory=dict)
    similarity_config: Dict[str, Any] = field(default_factory=dict)
    
    # Performance configuration
    max_workers: int = 4
    timeout_seconds: int = 300
    
    # KBase specific configuration
    workspace_url: Optional[str] = None
    auth_token: Optional[str] = None
    
    # Custom configuration
    custom_config: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        if not self.input_proteins and not self.input_file_path and not self.workspace_object_ref:
            raise ValueError("Must provide either input_proteins, input_file_path, or workspace_object_ref")
        
        if self.similarity_threshold < 0 or self.similarity_threshold > 1:
            raise ValueError("similarity_threshold must be between 0 and 1")
        
        if self.max_similar_proteins <= 0:
            raise ValueError("max_similar_proteins must be positive")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'input_proteins': self.input_proteins,
            'input_file_path': self.input_file_path,
            'workspace_object_ref': self.workspace_object_ref,
            'perform_embedding_generation': self.perform_embedding_generation,
            'perform_family_assignment': self.perform_family_assignment,
            'perform_similarity_search': self.perform_similarity_search,
            'perform_network_analysis': self.perform_network_analysis,
            'perform_sequence_analysis': self.perform_sequence_analysis,
            'enabled_stages': self.enabled_stages,
            'stage_configs': self.stage_configs,
            'embedding_model': self.embedding_model,
            'embedding_batch_size': self.embedding_batch_size,
            'embedding_device': self.embedding_device,
            'similarity_threshold': self.similarity_threshold,
            'max_similar_proteins': self.max_similar_proteins,
            'network_min_edges': self.network_min_edges,
            'network_layout': self.network_layout,
            'generate_html_report': self.generate_html_report,
            'generate_network_visualization': self.generate_network_visualization,
            'output_format': self.output_format,
            'cache_embeddings': self.cache_embeddings,
            'cache_directory': self.cache_directory,
            'storage_config': self.storage_config,
            'similarity_config': self.similarity_config,
            'max_workers': self.max_workers,
            'timeout_seconds': self.timeout_seconds,
            'workspace_url': self.workspace_url,
            'auth_token': self.auth_token,
            'custom_config': self.custom_config
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PipelineConfig':
        """Create configuration from dictionary."""
        return cls(**config_dict)
