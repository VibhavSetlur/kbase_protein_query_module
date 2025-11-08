# -*- coding: utf-8 -*-

"""
KBase Protein Query Module

A comprehensive protein analysis pipeline that performs sequence analysis,
embedding generation, similarity search, family assignment, network analysis,
and bioinformatics analysis using state-of-the-art machine learning models.
"""

# Re-export main workflow classes (lazy imports)
# Lazy imports to avoid heavy dependencies at startup
def __getattr__(name):
    if name == 'WorkflowOrchestrator':
        from .src.core import WorkflowOrchestrator
        return WorkflowOrchestrator
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

__version__ = "2.0.0"
__author__ = "KBase Team"
__email__ = "support@kbase.us"

__all__ = [
    'WorkflowOrchestrator'
]


