"""
Core module for KBase Protein Query Module

This module contains the fundamental building blocks for the protein query analysis pipeline:
- Base classes for pipeline stages
- Data structures for stage results and configuration
- Common interfaces and abstractions
"""

from ..stages.base_stage import BaseStage, StageResult
from .pipeline_config import PipelineConfig

__all__ = [
    'BaseStage',
    'StageResult', 
    'PipelineConfig'
]
