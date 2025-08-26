# -*- coding: utf-8 -*-

"""
KBase Protein Query Module

A comprehensive protein analysis pipeline that performs sequence analysis,
embedding generation, similarity search, family assignment, network analysis,
and bioinformatics analysis using state-of-the-art machine learning models.
"""

# Re-export main workflow and configuration classes
from .src.workflows import ProteinQueryWorkflow, WorkflowResult
from .src.core import PipelineConfig

__version__ = "2.0.0"
__author__ = "KBase Team"
__email__ = "support@kbase.us"

__all__ = [
    'ProteinQueryWorkflow',
    'WorkflowResult',
    'PipelineConfig'
]


