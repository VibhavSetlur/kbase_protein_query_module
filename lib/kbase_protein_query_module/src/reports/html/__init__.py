"""
HTML Report Generation Package

This package provides HTML report generation capabilities for the protein analysis pipeline.
"""

from .generator import HTMLReportGenerator
from .visualization_generator import VisualizationGenerator

__all__ = ['HTMLReportGenerator', 'VisualizationGenerator']
