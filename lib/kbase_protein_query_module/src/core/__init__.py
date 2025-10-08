"""
Core module for KBase Protein Query Module.

This module provides core functionality including workflow orchestration
and pipeline configuration.
"""

from .workflow_orchestrator import WorkflowOrchestrator, WorkflowResult
from .pipeline_config import PipelineConfig

__all__ = [
    'WorkflowOrchestrator',
    'WorkflowResult',
    'PipelineConfig'
]