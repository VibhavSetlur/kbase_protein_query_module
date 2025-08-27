"""
Pipeline Configuration for Protein Query Analysis

This module defines comprehensive configuration classes for scalable protein query 
analysis pipelines with resource management and performance optimization.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Union
import os
import psutil
import logging

logger = logging.getLogger(__name__)


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
    
    # Performance and scalability configuration
    max_workers: int = 4
    timeout_seconds: int = 300
    max_memory_gb: float = 8.0
    enable_parallel_processing: bool = True
    batch_size_proteins: int = 1000
    chunk_size: int = 1000
    use_streaming: bool = True
    cache_size_mb: int = 512
    enable_resource_monitoring: bool = True
    auto_scale_batch_size: bool = True
    max_concurrent_tasks: int = 4
    gc_threshold_mb: float = 512.0
    
    # KBase specific configuration
    workspace_url: Optional[str] = None
    auth_token: Optional[str] = None
    
    # Custom configuration
    custom_config: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate and auto-configure settings after initialization."""
        # Validate input requirements
        if not self.input_proteins and not self.input_file_path and not self.workspace_object_ref:
            raise ValueError("Must provide either input_proteins, input_file_path, or workspace_object_ref")
        
        # Validate thresholds
        if self.similarity_threshold < 0 or self.similarity_threshold > 1:
            raise ValueError("similarity_threshold must be between 0 and 1")
        
        if self.max_similar_proteins <= 0:
            raise ValueError("max_similar_proteins must be positive")
        
        # Auto-configure scalability parameters based on system resources
        self._auto_configure_scalability()
    
    def _auto_configure_scalability(self):
        """
        Auto-configure scalability parameters for KBase DOE server environments.
        
        Uses percentage-based limits to ensure respectful resource usage on shared servers.
        """
        try:
            # Get system resources
            memory = psutil.virtual_memory()
            cpu_count = psutil.cpu_count()
            
            # Server-aware configuration using percentage-based limits
            total_memory_gb = memory.total / (1024**3)
            
            # Use conservative percentages for shared DOE servers
            max_usable_memory = total_memory_gb * 0.6  # 60% of total memory
            max_usable_cores = max(1, int(cpu_count * 0.7))  # 70% of CPU cores
            
            # Apply server-safe limits
            self.max_memory_gb = min(self.max_memory_gb, max_usable_memory)
            self.max_workers = min(self.max_workers, max_usable_cores)
            self.max_concurrent_tasks = min(self.max_concurrent_tasks, max_usable_cores)
            
            # Conservative batch sizing for server environments
            memory_factor = max_usable_memory / 8.0  # Base factor for 8GB
            server_batch_factor = 0.5  # Additional 50% reduction for server safety
            self.batch_size_proteins = int(self.batch_size_proteins * memory_factor * server_batch_factor)
            self.chunk_size = int(self.chunk_size * memory_factor * server_batch_factor)
            
            # Conservative cache sizing
            cache_factor = min(max_usable_memory / 8.0, 1.5) * 0.8  # 80% of calculated for safety
            self.cache_size_mb = int(self.cache_size_mb * cache_factor)
            
            # Apply absolute minimums for functionality
            self.batch_size_proteins = max(50, self.batch_size_proteins)
            self.chunk_size = max(100, self.chunk_size)
            self.cache_size_mb = max(128, self.cache_size_mb)
            
            # Detect KBase server environment and apply extra restrictions
            if os.environ.get('KBASE_ENDPOINT') or os.environ.get('KB_DEPLOYMENT_CONFIG'):
                self.batch_size_proteins = min(self.batch_size_proteins, 300)
                self.max_workers = min(self.max_workers, 2)
                self.max_concurrent_tasks = min(self.max_concurrent_tasks, 2)
                logger.info("Detected KBase server environment - applying extra conservative limits")
            
            logger.info(f"Server-aware auto-configuration: memory={self.max_memory_gb:.1f}GB "
                       f"({(self.max_memory_gb/total_memory_gb)*100:.1f}% of total), "
                       f"workers={self.max_workers} ({(self.max_workers/cpu_count)*100:.1f}% of cores), "
                       f"batch_size={self.batch_size_proteins}")
                       
        except Exception as e:
            logger.warning(f"Could not auto-configure scalability parameters: {e}")
            # Fallback to very conservative defaults for server safety
            self.max_memory_gb = min(self.max_memory_gb, 4.0)
            self.max_workers = min(self.max_workers, 2)
            self.batch_size_proteins = min(self.batch_size_proteins, 200)
    
    def get_resource_limits(self) -> Dict[str, Any]:
        """Get resource limits for the resource manager."""
        return {
            'max_memory_gb': self.max_memory_gb,
            'max_cpu_percent': 80.0,
            'max_disk_usage_gb': 100.0,
            'batch_size_proteins': self.batch_size_proteins,
            'max_concurrent_tasks': self.max_concurrent_tasks,
            'gc_threshold_mb': self.gc_threshold_mb
        }
    
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
