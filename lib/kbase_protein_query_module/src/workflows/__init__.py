"""
Workflows module for KBase Protein Query Module

This module contains workflow orchestrators that coordinate the execution
of protein query analysis pipelines.
"""

from .workflow_orchestrator import ProteinQueryWorkflow, WorkflowResult

__all__ = [
    'ProteinQueryWorkflow',
    'WorkflowResult'
] 
