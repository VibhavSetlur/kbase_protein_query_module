# -*- coding: utf-8 -*-

"""
KBase Protein Query Module

A comprehensive protein analysis pipeline that performs sequence analysis,
embedding generation, similarity search, family assignment, network analysis,
and bioinformatics analysis using state-of-the-art machine learning models.
"""

# Re-export main workflow and configuration classes (lazy imports)
# from .src.workflows import ProteinQueryWorkflow, WorkflowResult
# from .src.core import PipelineConfig

# Lazy imports to avoid heavy dependencies at startup
def __getattr__(name):
    if name == 'ProteinQueryWorkflow':
        from .src.workflows import ProteinQueryWorkflow
        return ProteinQueryWorkflow
    elif name == 'WorkflowResult':
        from .src.workflows import WorkflowResult
        return WorkflowResult
    elif name == 'PipelineConfig':
        from .src.core import PipelineConfig
        return PipelineConfig
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

__version__ = "2.0.0"
__author__ = "KBase Team"
__email__ = "support@kbase.us"

__all__ = [
    'ProteinQueryWorkflow',
    'WorkflowResult',
    'PipelineConfig'
]


