"""
HTML Report Components Module

This module contains components for generating different parts of HTML reports.
"""

from .dashboard import DashboardGenerator
from .charts import ChartGenerator
from .network_viz import NetworkVisualizationGenerator
from .protein_details import ProteinDetailsGenerator

__all__ = [
    'DashboardGenerator',
    'ChartGenerator',
    'NetworkVisualizationGenerator',
    'ProteinDetailsGenerator'
]
