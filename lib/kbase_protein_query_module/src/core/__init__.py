"""
Core module for KBase Protein Query Module.

This module provides core functionality including workflow orchestration,
performance monitoring, and resource management.
"""

from .workflow_orchestrator import WorkflowOrchestrator, WorkflowResult
from .analysis_registry import AnalysisRegistry
from .parallel_processor import ParallelProcessor
from .performance_monitor import PerformanceProfiler, get_performance_profiler, profile_function
from .pipeline_config import PipelineConfig
from .resource_manager import ResourceManager

__all__ = [
    'WorkflowOrchestrator',
    'WorkflowResult',
    'AnalysisRegistry',
    'ParallelProcessor',
    'PerformanceProfiler',
    'get_performance_profiler',
    'profile_function',
    'PipelineConfig',
    'ResourceManager'
]