"""
Core module for KBase Protein Query Module

This module contains the fundamental building blocks for the protein query analysis pipeline:
- Base classes for pipeline stages and analyses
- Data structures for stage results and configuration
- Common interfaces and abstractions
- Extensibility framework for adding new analyses and indexing strategies
- Resource management and performance monitoring
"""

from ..stages.base_stage import BaseStage, StageResult
from .pipeline_config import PipelineConfig
from .analysis_registry import BaseAnalysis, AnalysisRegistry, AnalysisMetadata, get_registry, register_analysis
from .resource_manager import ResourceManager, ResourceLimits, get_resource_manager
from .parallel_processor import ParallelProcessor
from .performance_monitor import PerformanceProfiler, get_performance_profiler

__all__ = [
    # Core pipeline components
    'BaseStage',
    'StageResult', 
    'PipelineConfig',
    
    # Extensibility framework
    'BaseAnalysis',
    'AnalysisRegistry',
    'AnalysisMetadata', 
    'get_registry',
    'register_analysis',
    
    # Resource management and performance
    'ResourceManager',
    'ResourceLimits',
    'get_resource_manager',
    'ParallelProcessor',
    'PerformanceProfiler',
    'get_performance_profiler'
]
