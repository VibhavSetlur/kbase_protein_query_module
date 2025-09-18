"""
Analysis Manager Module

This module manages the execution of all available analyses in the KBase Protein Query Module.
It coordinates between different analysis types and ensures proper dependency management.
"""

import os
import logging
import importlib
from typing import Dict, Any, List, Optional, Union
from pathlib import Path

from .config import (
    get_enabled_analyses, 
    get_analysis_dependencies, 
    is_analysis_enabled,
    validate_analysis_config
)

logger = logging.getLogger(__name__)

class AnalysisManager:
    """
    Manages the execution of protein analysis workflows.
    
    This class coordinates different analysis types, handles dependencies,
    and ensures proper integration with the output system.
    """
    
    def __init__(self, output_manager=None):
        """
        Initialize the Analysis Manager.
        
        Args:
            output_manager: Instance of OutputManager for handling results
        """
        self.output_manager = output_manager
        self.analysis_modules = {}
        self.results = {}
        
        # Validate configuration
        if not validate_analysis_config():
            logger.warning("Analysis configuration validation failed")
            
        # Load available analysis modules
        self._load_analysis_modules()
    
    def _load_analysis_modules(self):
        """Load all available analysis modules dynamically."""
        enabled_analyses = get_enabled_analyses()
        
        for analysis_name, config in enabled_analyses.items():
            try:
                # Import the analysis module
                module_path = f"analysis.{analysis_name}"
                module = importlib.import_module(module_path)
                
                # Look for the main analysis class (convention: AnalysisName + "Analysis")
                class_name = f"{analysis_name.title().replace('_', '')}Analysis"
                if hasattr(module, class_name):
                    self.analysis_modules[analysis_name] = getattr(module, class_name)
                    logger.info(f"Loaded analysis module: {analysis_name}")
                else:
                    logger.warning(f"Could not find analysis class {class_name} in {module_path}")
                    
            except ImportError as e:
                logger.error(f"Failed to import analysis module {analysis_name}: {e}")
            except Exception as e:
                logger.error(f"Error loading analysis module {analysis_name}: {e}")
    
    def get_available_analyses(self) -> Dict[str, Dict[str, Any]]:
        """
        Get list of available analyses.
        
        Returns:
            Dictionary of available analyses with their configurations
        """
        return get_enabled_analyses()
    
    def get_analysis_dependencies(self, analysis_name: str) -> List[str]:
        """
        Get dependencies for a specific analysis.
        
        Args:
            analysis_name: Name of the analysis
            
        Returns:
            List of dependency names
        """
        return get_analysis_dependencies(analysis_name)
    
    def can_run_analysis(self, analysis_name: str, available_data: Dict[str, Any]) -> bool:
        """
        Check if an analysis can be run with the available data.
        
        Args:
            analysis_name: Name of the analysis to check
            available_data: Dictionary of available data types
            
        Returns:
            True if analysis can be run, False otherwise
        """
        if not is_analysis_enabled(analysis_name):
            return False
            
        if analysis_name not in self.analysis_modules:
            return False
            
        # Check if required data is available
        dependencies = self.get_analysis_dependencies(analysis_name)
        for dep in dependencies:
            if dep not in available_data:
                logger.warning(f"Analysis {analysis_name} requires {dep} which is not available")
                return False
                
        return True
    
    def run_analysis(self, analysis_name: str, input_data: Dict[str, Any], 
                    output_dir: str, **kwargs) -> Dict[str, Any]:
        """
        Run a specific analysis.
        
        Args:
            analysis_name: Name of the analysis to run
            input_data: Input data for the analysis
            output_dir: Directory to save outputs
            **kwargs: Additional parameters for the analysis
            
        Returns:
            Dictionary containing analysis results
        """
        if not self.can_run_analysis(analysis_name, input_data):
            raise ValueError(f"Cannot run analysis {analysis_name} with available data")
        
        try:
            logger.info(f"Starting analysis: {analysis_name}")
            
            # Get the analysis class
            analysis_class = self.analysis_modules[analysis_name]
            
            # Initialize the analysis
            analysis_instance = analysis_class()
            
            # Run the analysis
            result = analysis_instance.run(input_data, output_dir, **kwargs)
            
            # Store results
            self.results[analysis_name] = result
            
            # If output manager is available, handle outputs
            if self.output_manager:
                self.output_manager.save_analysis_output(analysis_name, result, output_dir)
            
            logger.info(f"Completed analysis: {analysis_name}")
            return result
            
        except Exception as e:
            logger.error(f"Error running analysis {analysis_name}: {e}")
            raise
    
    def run_multiple_analyses(self, analysis_names: List[str], input_data: Dict[str, Any],
                            output_dir: str, **kwargs) -> Dict[str, Any]:
        """
        Run multiple analyses in dependency order.
        
        Args:
            analysis_names: List of analysis names to run
            input_data: Input data for the analyses
            output_dir: Directory to save outputs
            **kwargs: Additional parameters
            
        Returns:
            Dictionary containing results from all analyses
        """
        results = {}
        current_data = input_data.copy()
        
        # Sort analyses by dependencies
        ordered_analyses = self._order_analyses_by_dependencies(analysis_names)
        
        for analysis_name in ordered_analyses:
            try:
                # Run the analysis
                result = self.run_analysis(analysis_name, current_data, output_dir, **kwargs)
                results[analysis_name] = result
                
                # Add results to current data for dependent analyses
                current_data[analysis_name] = result
                
            except Exception as e:
                logger.error(f"Failed to run analysis {analysis_name}: {e}")
                results[analysis_name] = {"error": str(e)}
        
        return results
    
    def _order_analyses_by_dependencies(self, analysis_names: List[str]) -> List[str]:
        """
        Order analyses by their dependencies.
        
        Args:
            analysis_names: List of analysis names
            
        Returns:
            Ordered list of analysis names
        """
        ordered = []
        remaining = set(analysis_names)
        
        while remaining:
            # Find analyses with no unmet dependencies
            ready = []
            for analysis in remaining:
                dependencies = set(self.get_analysis_dependencies(analysis))
                if dependencies.issubset(set(ordered)):
                    ready.append(analysis)
            
            if not ready:
                # Circular dependency or missing dependency
                logger.warning(f"Could not resolve dependencies for: {remaining}")
                ordered.extend(remaining)
                break
                
            # Add ready analyses to ordered list
            ordered.extend(ready)
            remaining -= set(ready)
        
        return ordered
    
    def get_analysis_results(self, analysis_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Get results from completed analyses.
        
        Args:
            analysis_name: Specific analysis name, or None for all results
            
        Returns:
            Dictionary containing analysis results
        """
        if analysis_name:
            return self.results.get(analysis_name, {})
        return self.results.copy()
    
    def clear_results(self):
        """Clear all stored analysis results."""
        self.results.clear()
        logger.info("Cleared all analysis results")
