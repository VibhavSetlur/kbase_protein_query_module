"""
Pipeline Configuration for KBase Protein Query Module

Simplified configuration for protein analysis workflows.
"""

import os
import logging
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

@dataclass
class PipelineConfig:
    """Configuration for the protein query pipeline."""
    
    # Input configuration
    input_proteins: List[str] = field(default_factory=list)
    input_type: str = "protein_input"  # "protein_input" or "uniprot_ids"
    
    # Output configuration  
    output_dir: str = "/tmp"
    workspace_name: Optional[str] = None
    workspace_client: Optional[Any] = None
    
    # Analysis configuration
    selected_analyses: Optional[List[str]] = None
    
    def __post_init__(self):
        """Validate and set default values after initialization."""
        # Set default analyses if none specified
        if not self.selected_analyses:
            self.selected_analyses = ["network_analysis"]
    
    
    def validate(self) -> bool:
        """Validate the configuration."""
        try:
            # Check required fields
            if not self.input_proteins:
                logger.error("input_proteins must be specified")
                return False
            
            # Check input type is valid
            valid_input_types = ["protein_input", "uniprot_ids"]
            if self.input_type not in valid_input_types:
                logger.error(f"Invalid input_type: {self.input_type}. Must be one of: {valid_input_types}")
                return False
            
            # Ensure output directory exists and is writable
            try:
                os.makedirs(self.output_dir, exist_ok=True)
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
            'output_dir': self.output_dir,
            'workspace_name': self.workspace_name,
            'workspace_client': self.workspace_client,
            'selected_analyses': self.selected_analyses
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PipelineConfig':
        """Create PipelineConfig from dictionary."""
        return cls(**config_dict)
    
    def merge(self, other_config: 'PipelineConfig') -> 'PipelineConfig':
        """Merge another configuration into this one."""
        merged_dict = self.to_dict()
        other_dict = other_config.to_dict()
        
        for key, value in other_dict.items():
            if value is not None:
                merged_dict[key] = value
        
        return PipelineConfig.from_dict(merged_dict)