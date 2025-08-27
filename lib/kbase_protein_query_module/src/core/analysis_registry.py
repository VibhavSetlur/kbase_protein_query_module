"""
Analysis Registry for Modular Extension Framework

This module provides a registry system for dynamically adding new analysis types
to the protein query pipeline. Developers can easily extend the module by
registering new analysis classes that follow the standardized interface.

DEVELOPER GUIDE:
================

To add a new analysis type:

1. Create your analysis class inheriting from BaseAnalysis
2. Implement required methods: analyze(), get_output_files(), get_metadata()
3. Register your analysis using @register_analysis decorator
4. Add HTML template for your analysis results
5. Update pipeline configuration to include your analysis

Example:
--------
@register_analysis("my_new_analysis")
class MyNewAnalysis(BaseAnalysis):
    def analyze(self, proteins: List[ProteinRecord]) -> Dict[str, Any]:
        # Your analysis logic here
        return results
    
    def get_output_files(self) -> List[str]:
        return ["my_analysis.html", "my_data.csv"]

Classes Used:
- BaseAnalysis: lib/kbase_protein_query_module/src/core/analysis_registry.py:45
- AnalysisRegistry: lib/kbase_protein_query_module/src/core/analysis_registry.py:120
- AnalysisMetadata: lib/kbase_protein_query_module/src/core/analysis_registry.py:25
"""

import os
import logging
import inspect
from typing import Dict, Any, List, Optional, Type, Callable
from dataclasses import dataclass
from abc import ABC, abstractmethod
import importlib
import pkgutil

logger = logging.getLogger(__name__)

@dataclass
class AnalysisMetadata:
    """Metadata for registered analyses."""
    name: str
    description: str
    version: str
    author: str
    output_files: List[str]
    dependencies: List[str]
    category: str
    computational_complexity: str  # "low", "medium", "high"
    memory_requirements: str  # "low", "medium", "high" 
    supports_batch: bool = True
    supports_streaming: bool = False

class BaseAnalysis(ABC):
    """
    Base class for all protein analyses.
    
    All new analysis types MUST inherit from this class and implement
    the required abstract methods.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/core/analysis_registry.py:45
    USED BY: All analysis implementations in src/stages/analysis/
    EXTENDS: ABC (Python Abstract Base Class)
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.logger = logging.getLogger(f"{self.__class__.__module__}.{self.__class__.__name__}")
        self.metadata = self._get_metadata()
    
    @abstractmethod
    def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
        """
        Perform the analysis on protein data.
        
        Args:
            proteins: List of protein records or sequences
            **kwargs: Additional analysis parameters
            
        Returns:
            Dictionary containing analysis results
        """
        pass
    
    @abstractmethod
    def get_output_files(self, output_dir: str) -> List[str]:
        """
        Get list of output files this analysis will create.
        
        Args:
            output_dir: Base output directory
            
        Returns:
            List of file paths that will be created
        """
        pass
    
    @abstractmethod
    def get_metadata(self) -> AnalysisMetadata:
        """Get metadata for this analysis type."""
        pass
    
    def _get_metadata(self) -> AnalysisMetadata:
        """Internal method to get metadata - override get_metadata() instead."""
        return self.get_metadata()
    
    def validate_input(self, proteins: List[Any]) -> bool:
        """
        Validate input data for this analysis.
        
        Args:
            proteins: Input protein data
            
        Returns:
            True if input is valid, False otherwise
        """
        return len(proteins) > 0
    
    def estimate_resources(self, proteins: List[Any]) -> Dict[str, Any]:
        """
        Estimate computational resources needed for this analysis.
        
        Args:
            proteins: Input protein data
            
        Returns:
            Dictionary with resource estimates
        """
        protein_count = len(proteins)
        
        # Base estimates - override in subclasses for more accurate estimates
        complexity_multipliers = {"low": 1, "medium": 5, "high": 20}
        memory_multipliers = {"low": 1, "medium": 10, "high": 50}
        
        complexity = complexity_multipliers.get(self.metadata.computational_complexity, 5)
        memory = memory_multipliers.get(self.metadata.memory_requirements, 10)
        
        return {
            "estimated_time_seconds": protein_count * complexity * 0.1,
            "estimated_memory_mb": protein_count * memory,
            "cpu_cores_recommended": min(4, max(1, protein_count // 100)),
            "supports_parallel": self.metadata.supports_batch
        }

class AnalysisRegistry:
    """
    Registry for managing available analysis types.
    
    This registry allows dynamic discovery and registration of analysis types,
    making it easy to extend the module with new analyses.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/core/analysis_registry.py:120
    USED BY: WorkflowOrchestrator, ReportGenerationStage
    SINGLETON PATTERN: Use get_registry() to access global instance
    """
    
    def __init__(self):
        self._analyses: Dict[str, Type[BaseAnalysis]] = {}
        self._metadata: Dict[str, AnalysisMetadata] = {}
        self.logger = logging.getLogger(__name__)
        
        # Auto-discover existing analyses
        self._discover_analyses()
    
    def register_analysis(self, name: str, analysis_class: Type[BaseAnalysis]) -> None:
        """
        Register a new analysis type.
        
        Args:
            name: Unique name for the analysis
            analysis_class: Class implementing the analysis
        """
        if not issubclass(analysis_class, BaseAnalysis):
            raise ValueError(f"Analysis class {analysis_class.__name__} must inherit from BaseAnalysis")
        
        # Get metadata from the class
        try:
            instance = analysis_class()
            metadata = instance.get_metadata()
            
            self._analyses[name] = analysis_class
            self._metadata[name] = metadata
            
            self.logger.info(f"Registered analysis: {name} ({analysis_class.__name__})")
            
        except Exception as e:
            self.logger.error(f"Failed to register analysis {name}: {e}")
            raise
    
    def get_analysis(self, name: str, config: Dict[str, Any] = None) -> BaseAnalysis:
        """
        Get an instance of a registered analysis.
        
        Args:
            name: Name of the analysis
            config: Configuration for the analysis
            
        Returns:
            Instance of the analysis class
        """
        if name not in self._analyses:
            available = list(self._analyses.keys())
            raise ValueError(f"Analysis '{name}' not found. Available: {available}")
        
        analysis_class = self._analyses[name]
        return analysis_class(config)
    
    def list_analyses(self) -> Dict[str, AnalysisMetadata]:
        """Get all registered analyses and their metadata."""
        return self._metadata.copy()
    
    def get_analyses_by_category(self, category: str) -> Dict[str, AnalysisMetadata]:
        """Get analyses filtered by category."""
        return {
            name: metadata 
            for name, metadata in self._metadata.items()
            if metadata.category == category
        }
    
    def _discover_analyses(self):
        """Auto-discover analysis classes in the stages/analysis directory."""
        try:
            # Import analysis modules to trigger registration
            analysis_modules = [
                'lib.kbase_protein_query_module.src.stages.analysis.sequence_analysis',
                'lib.kbase_protein_query_module.src.stages.analysis.network_analysis',
                'lib.kbase_protein_query_module.src.stages.analysis.multi_protein_analysis',
                'lib.kbase_protein_query_module.src.stages.analysis.bioinformatics_analysis'
            ]
            
            for module_name in analysis_modules:
                try:
                    importlib.import_module(module_name)
                    self.logger.debug(f"Loaded analysis module: {module_name}")
                except ImportError as e:
                    self.logger.debug(f"Could not load analysis module {module_name}: {e}")
                    
        except Exception as e:
            self.logger.warning(f"Error during analysis discovery: {e}")
    
    def generate_analysis_documentation(self) -> str:
        """Generate comprehensive documentation for all analyses."""
        doc = """
# Available Analysis Types

This document lists all available analysis types and their characteristics.

"""
        
        for name, metadata in self._metadata.items():
            analysis_class = self._analyses[name]
            doc += f"""
## {metadata.name}

**Class:** `{analysis_class.__name__}`  
**Location:** `{inspect.getfile(analysis_class)}`  
**Category:** {metadata.category}  
**Version:** {metadata.version}  
**Author:** {metadata.author}  

**Description:** {metadata.description}

**Computational Requirements:**
- Complexity: {metadata.computational_complexity}
- Memory: {metadata.memory_requirements}
- Batch Processing: {'✅' if metadata.supports_batch else '❌'}
- Streaming: {'✅' if metadata.supports_streaming else '❌'}

**Output Files:**
{chr(10).join(f'- {file}' for file in metadata.output_files)}

**Dependencies:**
{chr(10).join(f'- {dep}' for dep in metadata.dependencies)}

---
"""
        
        return doc

# Global registry instance
_registry: Optional[AnalysisRegistry] = None

def get_registry() -> AnalysisRegistry:
    """Get the global analysis registry instance."""
    global _registry
    if _registry is None:
        _registry = AnalysisRegistry()
    return _registry

def register_analysis(name: str):
    """
    Decorator for registering analysis classes.
    
    Usage:
    @register_analysis("my_analysis")
    class MyAnalysis(BaseAnalysis):
        ...
    """
    def decorator(cls: Type[BaseAnalysis]) -> Type[BaseAnalysis]:
        registry = get_registry()
        registry.register_analysis(name, cls)
        return cls
    return decorator

def list_available_analyses() -> Dict[str, AnalysisMetadata]:
    """Convenience function to list all available analyses."""
    registry = get_registry()
    return registry.list_analyses()

def create_analysis(name: str, config: Dict[str, Any] = None) -> BaseAnalysis:
    """Convenience function to create an analysis instance."""
    registry = get_registry()
    return registry.get_analysis(name, config)
