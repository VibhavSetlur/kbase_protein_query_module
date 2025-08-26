"""
Workflows module for KBase Protein Query Module

This module provides workflow orchestration for protein query analysis.
"""

from .workflow_orchestrator import ProteinQueryWorkflow, WorkflowResult

__all__ = [
    'ProteinQueryWorkflow',
    'WorkflowResult'
] 
