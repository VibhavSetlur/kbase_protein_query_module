"""
Output Configuration Module

This module defines output settings and configurations for the KBase Protein Query Module.
It manages output formats, directory structures, and output-specific settings.
"""

from typing import Dict, Any, List, Optional
import logging

logger = logging.getLogger(__name__)

# Output Configuration
OUTPUT_CONFIG = {
    "base_output_dir": "output",
    "create_analysis_folders": True,
    "include_metadata": True,
    "include_process_info": True,
    "formats": ["json", "csv", "html"],
    "compression": False,
    "timestamp_format": "%Y%m%d_%H%M%S",
    "max_file_size_mb": 100,
    "cleanup_temp_files": True
}

# Analysis-specific output configurations
ANALYSIS_OUTPUT_CONFIG = {
    "network_analysis": {
        "enabled": True,
        "output_formats": ["json", "html", "csv"],
        "include_visualizations": True,
        "include_network_data": True,
        "include_statistics": True,
        "max_nodes_display": 1000
    },
    "sequence_analysis": {
        "enabled": True,
        "output_formats": ["json", "csv"],
        "include_alignments": True,
        "include_motifs": True,
        "include_statistics": True
    },
    "bioinformatics_analysis": {
        "enabled": True,
        "output_formats": ["json", "csv", "html"],
        "include_domain_predictions": True,
        "include_functional_annotations": True,
        "include_statistics": True
    },
    "multi_protein_analysis": {
        "enabled": True,
        "output_formats": ["json", "csv", "html"],
        "include_comparisons": True,
        "include_statistics": True,
        "include_visualizations": True
    },
    "similarity_search": {
        "enabled": True,
        "output_formats": ["json", "csv"],
        "include_similarity_scores": True,
        "include_metadata": True,
        "max_results": 1000
    },
    "family_assignment": {
        "enabled": True,
        "output_formats": ["json", "csv"],
        "include_confidence_scores": True,
        "include_family_metadata": True,
        "include_statistics": True
    }
}

# File naming conventions
FILE_NAMING_CONFIG = {
    "results_file": "results.json",
    "summary_file": "summary.txt",
    "metadata_file": "metadata.json",
    "manifest_file": "manifest.json",
    "process_info_file": "process_info.json",
    "final_output_file": "final_output.json"
}

# Directory structure configuration
DIRECTORY_CONFIG = {
    "analysis_subdir": "analysis",
    "metadata_subdir": "metadata", 
    "process_info_subdir": "process_info",
    "logs_subdir": "logs",
    "temp_subdir": "temp"
}

def get_output_config() -> Dict[str, Any]:
    """
    Get the main output configuration.
    
    Returns:
        Dictionary containing output configuration
    """
    return OUTPUT_CONFIG.copy()

def get_analysis_output_config(analysis_name: str) -> Dict[str, Any]:
    """
    Get output configuration for a specific analysis.
    
    Args:
        analysis_name: Name of the analysis
        
    Returns:
        Dictionary containing analysis-specific output configuration
    """
    return ANALYSIS_OUTPUT_CONFIG.get(analysis_name, {}).copy()

def get_enabled_output_analyses() -> Dict[str, Dict[str, Any]]:
    """
    Get configurations for analyses that have output enabled.
    
    Returns:
        Dictionary containing enabled analysis output configurations
    """
    return {name: config for name, config in ANALYSIS_OUTPUT_CONFIG.items() 
            if config.get("enabled", False)}

def is_output_enabled_for_analysis(analysis_name: str) -> bool:
    """
    Check if output is enabled for a specific analysis.
    
    Args:
        analysis_name: Name of the analysis
        
    Returns:
        True if output is enabled, False otherwise
    """
    return ANALYSIS_OUTPUT_CONFIG.get(analysis_name, {}).get("enabled", False)

def get_output_formats_for_analysis(analysis_name: str) -> List[str]:
    """
    Get supported output formats for a specific analysis.
    
    Args:
        analysis_name: Name of the analysis
        
    Returns:
        List of supported output formats
    """
    return ANALYSIS_OUTPUT_CONFIG.get(analysis_name, {}).get("output_formats", ["json"])

def get_file_naming_config() -> Dict[str, str]:
    """
    Get file naming configuration.
    
    Returns:
        Dictionary containing file naming conventions
    """
    return FILE_NAMING_CONFIG.copy()

def get_directory_config() -> Dict[str, str]:
    """
    Get directory structure configuration.
    
    Returns:
        Dictionary containing directory structure settings
    """
    return DIRECTORY_CONFIG.copy()

def validate_output_config() -> bool:
    """
    Validate the output configuration for consistency.
    
    Returns:
        True if configuration is valid, False otherwise
    """
    try:
        # Check that all analysis output configs have required fields
        for analysis_name, config in ANALYSIS_OUTPUT_CONFIG.items():
            if not isinstance(config, dict):
                logger.error(f"Invalid config for analysis {analysis_name}: not a dictionary")
                return False
                
            if "enabled" not in config:
                logger.warning(f"Analysis {analysis_name} missing 'enabled' field")
                
            if "output_formats" not in config:
                logger.warning(f"Analysis {analysis_name} missing 'output_formats' field")
        
        # Check that file naming config has required files
        required_files = ["results_file", "summary_file", "metadata_file", "manifest_file"]
        for required_file in required_files:
            if required_file not in FILE_NAMING_CONFIG:
                logger.error(f"Missing required file naming config: {required_file}")
                return False
        
        return True
        
    except Exception as e:
        logger.error(f"Error validating output configuration: {e}")
        return False

# Validate configuration on import
if not validate_output_config():
    logger.warning("Output configuration validation failed. Some outputs may not work correctly.")
