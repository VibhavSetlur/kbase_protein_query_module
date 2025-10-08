"""
Analysis Output Registry

This module provides a registry for analysis output handlers.
"""

from typing import Dict, Type, Any, List

# Registry for analysis output handlers
_analysis_output_registry: Dict[str, Type[Any]] = {}

def register_analysis_output(analysis_name: str):
    """Decorator to register an analysis output handler."""
    def decorator(cls):
        _analysis_output_registry[analysis_name] = cls
        return cls
    return decorator

def get_analysis_output_handler(analysis_name: str):
    """Get an analysis output handler by name."""
    return _analysis_output_registry.get(analysis_name)

def list_registered_handlers() -> List[str]:
    """List all registered analysis output handlers."""
    return list(_analysis_output_registry.keys())
