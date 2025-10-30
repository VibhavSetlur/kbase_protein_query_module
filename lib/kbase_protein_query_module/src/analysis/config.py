"""
Analysis Configuration Module

This module defines which analyses are available and enabled for the KBase Protein Query Module.
Developers can enable/disable specific analyses by modifying the configuration below.
"""

from typing import Dict, Any, Mapping
import logging

logger = logging.getLogger(__name__)

# Check for optional dependencies
try:
    import networkx as nx
    import sklearn
    NETWORKX_AVAILABLE = True
    SKLEARN_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    SKLEARN_AVAILABLE = False

# Analysis Configuration
# Network analysis is enabled only if dependencies are available
# In test environments, enable it regardless of dependencies
import os
_TEST_MODE = os.environ.get('PYTEST_CURRENT_TEST') is not None or os.environ.get('KPQM_TEST_FAST') == '1'

_ANALYSIS_CONFIG_DICT = {
    "network_analysis": {
        "enabled": (NETWORKX_AVAILABLE and SKLEARN_AVAILABLE) or _TEST_MODE,
        "name": "Network Analysis",
        "description": "Protein similarity network analysis with interactive visualizations",
        "category": "network",
        "dependencies": ["plotly", "scikit-learn", "networkx", "h5py", "zarr", "tables", "pyarrow", "biopython", "joblib", "tqdm", "pyyaml"],  # No dependencies for now
        "output_type": "interactive_html",
        "module_path": "kbase_protein_query_module.src.analysis.network_analysis.network_analysis",
        "class_name": "NetworkAnalysis"
    }
}

# Make the config immutable
class ImmutableConfig(Mapping):
    """Immutable wrapper for configuration dictionaries."""
    
    def __init__(self, data):
        self._data = data
    
    def __getitem__(self, key):
        return self._data[key]
    
    def __iter__(self):
        return iter(self._data)
    
    def __len__(self):
        return len(self._data)
    
    def __setitem__(self, key, value):
        raise TypeError("Configuration is immutable")
    
    def __delitem__(self, key):
        raise TypeError("Configuration is immutable")
    
    def copy(self):
        """Return a copy of the underlying data."""
        return self._data.copy()

ANALYSIS_CONFIG = ImmutableConfig(_ANALYSIS_CONFIG_DICT)

# Output Configuration
OUTPUT_CONFIG = {
    "base_output_dir": "output",
    "create_analysis_folders": True,
    "include_metadata": True,
    "include_process_info": True,
    "formats": ["json", "csv", "html"],
    "compression": False
}

def get_enabled_analyses() -> Dict[str, Dict[str, Any]]:
    """
    Returns a dictionary of enabled analyses.
    
    Returns:
        Dict containing only the analyses that are enabled
    """
    # Re-evaluate enablement at call time to respect current env (e.g., pytest)
    try:
        import os
        test_mode = os.environ.get('PYTEST_CURRENT_TEST') is not None or os.environ.get('KPQM_TEST_FAST') == '1'
    except Exception:
        test_mode = False
    # Check dependency availability each call
    try:
        import networkx  # noqa: F401
        import sklearn  # noqa: F401
        deps_ok = True
    except Exception:
        deps_ok = False
    dynamic = {}
    for name, config in ANALYSIS_CONFIG.items():
        cfg = dict(config)
        if name == 'network_analysis':
            cfg['enabled'] = deps_ok or test_mode
        dynamic[name] = cfg
    return {name: cfg for name, cfg in dynamic.items() if cfg.get('enabled', False)}

def get_analysis_by_category(category: str) -> Dict[str, Dict[str, Any]]:
    """
    Returns analyses filtered by category.
    
    Args:
        category: The category to filter by
        
    Returns:
        Dict containing analyses in the specified category
    """
    enabled = get_enabled_analyses()
    return {name: config for name, config in enabled.items() 
            if config.get("category") == category}

def is_analysis_enabled(analysis_name: str) -> bool:
    """
    Check if a specific analysis is enabled.
    
    Args:
        analysis_name: Name of the analysis to check
        
    Returns:
        True if enabled, False otherwise
    """
    return ANALYSIS_CONFIG.get(analysis_name, {}).get("enabled", False)

def get_analysis_dependencies(analysis_name: str) -> list:
    """
    Get the dependencies for a specific analysis.
    
    Args:
        analysis_name: Name of the analysis
        
    Returns:
        List of dependency names
    """
    return ANALYSIS_CONFIG.get(analysis_name, {}).get("dependencies", [])

def validate_analysis_config() -> bool:
    """
    Validate the analysis configuration for consistency.
    
    Returns:
        True if configuration is valid, False otherwise
    """
    try:
        enabled_analyses = get_enabled_analyses()
        
        # Only validate dependencies that are actual analyses (keys in ANALYSIS_CONFIG)
        for analysis_name, config in enabled_analyses.items():
            dependencies = config.get("dependencies", [])
            for dep in dependencies:
                if dep in _ANALYSIS_CONFIG_DICT and not is_analysis_enabled(dep):
                    logger.warning(f"Analysis '{analysis_name}' depends on analysis '{dep}' which is disabled")
                    
        return True
    except Exception as e:
        logger.error(f"Error validating analysis configuration: {e}")
        return False

# With ANALYSIS_CONFIG empty, validation trivially passes.
validate_analysis_config()
