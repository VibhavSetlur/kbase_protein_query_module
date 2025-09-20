"""
Analysis Configuration Module

This module defines which analyses are available and enabled for the KBase Protein Query Module.
Developers can enable/disable specific analyses by modifying the configuration below.
"""

from typing import Dict, Any
import logging

logger = logging.getLogger(__name__)

# Analysis Configuration
# Set to True to enable an analysis, False to disable
ANALYSIS_CONFIG = {
    "network_analysis": {
        "enabled": True,
        "name": "Network Analysis",
        "description": "Constructs protein interaction networks using similarity-based connections",
        "category": "network",
        "dependencies": ["similarity_search", "embeddings"],
        "output_type": "network_data"
    },
    "sequence_analysis": {
        "enabled": True,
        "name": "Sequence Analysis", 
        "description": "Performs comprehensive sequence-based analysis including motifs and patterns",
        "category": "sequence",
        "dependencies": ["embeddings"],
        "output_type": "sequence_data"
    },
    "bioinformatics_analysis": {
        "enabled": True,
        "name": "Bioinformatics Analysis",
        "description": "Advanced bioinformatics analysis including domain prediction and functional annotation",
        "category": "bioinformatics", 
        "dependencies": ["embeddings", "similarity_search"],
        "output_type": "bioinformatics_data"
    },
    "multi_protein_analysis": {
        "enabled": True,
        "name": "Multi-Protein Analysis",
        "description": "Comparative analysis across multiple proteins with statistical insights",
        "category": "comparative",
        "dependencies": ["embeddings", "similarity_search"],
        "output_type": "comparative_data"
    },
    "similarity_search": {
        "enabled": True,
        "name": "Similarity Search",
        "description": "Finds similar proteins using embedding-based similarity metrics",
        "category": "search",
        "dependencies": ["embeddings"],
        "output_type": "similarity_data"
    },
    "family_assignment": {
        "enabled": True,
        "name": "Family Assignment",
        "description": "Assigns proteins to functional families based on sequence similarity",
        "category": "classification",
        "dependencies": ["embeddings", "similarity_search"],
        "output_type": "family_data"
    },
    "embeddings": {
        "enabled": True,
        "name": "Protein Embeddings",
        "description": "Generates protein sequence embeddings using deep learning models",
        "category": "preprocessing",
        "dependencies": [],
        "output_type": "embedding_data"
    }
}

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
    return {name: config for name, config in ANALYSIS_CONFIG.items() 
            if config.get("enabled", False)}

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
        
        # Check that all dependencies of enabled analyses are also enabled
        for analysis_name, config in enabled_analyses.items():
            dependencies = config.get("dependencies", [])
            for dep in dependencies:
                if not is_analysis_enabled(dep):
                    logger.warning(f"Analysis '{analysis_name}' depends on '{dep}' which is disabled")
                    
        return True
    except Exception as e:
        logger.error(f"Error validating analysis configuration: {e}")
        return False

# Validate configuration on import
if not validate_analysis_config():
    logger.warning("Analysis configuration validation failed. Some analyses may not work correctly.")
