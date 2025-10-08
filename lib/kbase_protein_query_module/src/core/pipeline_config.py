"""
Pipeline Configuration for KBase Protein Query Module

This module provides configuration management for the protein query pipeline,
including input processing, analysis execution, and output generation settings.
"""

import os
import logging
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

@dataclass
class PipelineConfig:
    """
    Configuration for the protein query pipeline.
    
    This class contains all configuration parameters needed to run
    the complete protein analysis workflow.
    """
    
    # Input configuration
    input_proteins: List[str] = field(default_factory=list)
    input_type: str = "protein_input"  # protein_input, uniprot_ids, workspace_object
    workspace_object_ref: Optional[str] = None
    
    # Pipeline stages
    enabled_stages: List[str] = field(default_factory=list)
    stage_configs: Dict[str, Any] = field(default_factory=dict)
    
    # Storage configuration
    storage_config: Dict[str, Any] = field(default_factory=dict)
    similarity_config: Dict[str, Any] = field(default_factory=dict)
    
    # Output configuration
    output_dir: str = "/tmp"
    workspace_name: Optional[str] = None
    workspace_client: Optional[Any] = None
    workspace_url: Optional[str] = None
    auth_token: Optional[str] = None
    
    # Report generation
    generate_html_report: bool = True
    generate_network_visualization: bool = True
    
    # Analysis configuration
    selected_analyses: Optional[List[str]] = None
    analysis_config: Dict[str, Any] = field(default_factory=dict)
    
    # Performance settings
    batch_size: int = 100
    max_memory_gb: float = 4.0
    max_concurrent_tasks: int = 2
    
    def __post_init__(self):
        """Validate and set default values after initialization."""
        # Set default enabled stages if none specified
        if not self.enabled_stages:
            self.enabled_stages = ["input_processing", "analysis", "output_generation"]
        
        # Set default stage configs
        if not self.stage_configs:
            self.stage_configs = {
                "input_processing": {
                    "validate_input": True,
                    "extract_data": True
                },
                "analysis": {
                    "run_enabled_analyses": True
                },
                "output_generation": {
                    "create_reports": True,
                    "zip_outputs": True
                }
            }
        
        # Set default storage config
        if not self.storage_config:
            self.storage_config = {
                "use_compression": True,
                "chunk_size": 1000,
                "max_family_size": 100000
            }
        
        # Set default similarity config
        if not self.similarity_config:
            self.similarity_config = {
                "similarity_threshold": 0.7,
                "top_k_matches": 50,
                "use_family_assignment": True
            }
    
    def get_stage_config(self, stage_name: str) -> Dict[str, Any]:
        """Get configuration for a specific stage."""
        return self.stage_configs.get(stage_name, {})
    
    def is_stage_enabled(self, stage_name: str) -> bool:
        """Check if a stage is enabled."""
        return stage_name in self.enabled_stages
    
    def enable_stage(self, stage_name: str):
        """Enable a stage."""
        if stage_name not in self.enabled_stages:
            self.enabled_stages.append(stage_name)
    
    def disable_stage(self, stage_name: str):
        """Disable a stage."""
        if stage_name in self.enabled_stages:
            self.enabled_stages.remove(stage_name)
    
    def set_stage_config(self, stage_name: str, config: Dict[str, Any]):
        """Set configuration for a specific stage."""
        self.stage_configs[stage_name] = config
    
    def validate(self) -> bool:
        """Validate the configuration."""
        try:
            # Check required fields
            if not self.input_proteins and not self.workspace_object_ref:
                logger.error("Either input_proteins or workspace_object_ref must be specified")
                return False
            
            # Check input type
            valid_input_types = ["protein_input", "uniprot_ids", "workspace_object"]
            if self.input_type not in valid_input_types:
                logger.error(f"Invalid input_type: {self.input_type}. Must be one of: {valid_input_types}")
                return False
            
            # Check workspace object ref for workspace_object input type
            if self.input_type == "workspace_object" and not self.workspace_object_ref:
                logger.error("workspace_object_ref is required for workspace_object input type")
                return False
            
            # Check output directory
            if not os.path.exists(os.path.dirname(self.output_dir)):
                try:
                    os.makedirs(os.path.dirname(self.output_dir), exist_ok=True)
                except Exception as e:
                    logger.error(f"Cannot create output directory: {e}")
                    return False
            
            return True
            
        except Exception as e:
            logger.error(f"Configuration validation failed: {e}")
            return False
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'input_proteins': self.input_proteins,
            'input_type': self.input_type,
            'workspace_object_ref': self.workspace_object_ref,
            'enabled_stages': self.enabled_stages,
            'stage_configs': self.stage_configs,
            'storage_config': self.storage_config,
            'similarity_config': self.similarity_config,
            'output_dir': self.output_dir,
            'workspace_name': self.workspace_name,
            'workspace_url': self.workspace_url,
            'auth_token': self.auth_token,
            'generate_html_report': self.generate_html_report,
            'generate_network_visualization': self.generate_network_visualization,
            'selected_analyses': self.selected_analyses,
            'analysis_config': self.analysis_config,
            'batch_size': self.batch_size,
            'max_memory_gb': self.max_memory_gb,
            'max_concurrent_tasks': self.max_concurrent_tasks
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PipelineConfig':
        """Create PipelineConfig from dictionary."""
        return cls(**config_dict)
    
    def merge(self, other_config: 'PipelineConfig') -> 'PipelineConfig':
        """Merge another configuration into this one."""
        merged_dict = self.to_dict()
        other_dict = other_config.to_dict()
        
        # Merge dictionaries, with other_config taking precedence
        for key, value in other_dict.items():
            if value is not None:
                merged_dict[key] = value
        
        return PipelineConfig.from_dict(merged_dict)