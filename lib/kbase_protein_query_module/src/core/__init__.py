"""
Core module for KBase Protein Query Module.

This module provides core functionality including workflow orchestration.
"""

from .workflow_orchestrator import WorkflowOrchestrator

__all__ = [
    'WorkflowOrchestrator'
]