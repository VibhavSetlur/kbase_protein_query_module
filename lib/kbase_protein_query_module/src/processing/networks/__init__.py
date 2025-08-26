"""
Networks Module for KBase Protein Query Module

This module contains network construction and analysis functionality.
"""

from .builder import DynamicNetworkBuilder

# Compatibility alias expected by some tests
class NetworkBuilder(DynamicNetworkBuilder):
    pass

__all__ = ['DynamicNetworkBuilder', 'NetworkBuilder']
