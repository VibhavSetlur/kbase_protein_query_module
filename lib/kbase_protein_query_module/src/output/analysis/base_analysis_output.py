"""
Base Analysis Output Handler

This module provides the base class for all analysis output handlers.
"""

from typing import Dict, Any, List
from dataclasses import dataclass

@dataclass
class AnalysisOutputResult:
    """Result of analysis output generation."""
    success: bool
    output_files: List[str]
    metadata: Dict[str, Any]
    summary: str = ""
    error_message: str = ""

class BaseAnalysisOutput:
    """Base class for all analysis output handlers."""
    
    def __init__(self, output_manager=None):
        """Initialize the base analysis output handler."""
        self.output_manager = output_manager
    
    def generate_outputs(self, analysis_data: Dict[str, Any], stage_name: str) -> AnalysisOutputResult:
        """Generate outputs for the analysis data."""
        raise NotImplementedError("Subclasses must implement generate_outputs")
    
    def get_output_description(self) -> str:
        """Get description of the outputs."""
        return "Analysis output"
    
    def get_supported_formats(self) -> List[str]:
        """Get supported output formats."""
        return ['json']
    
    def validate_analysis_data(self, analysis_data: Dict[str, Any]) -> bool:
        """Validate that analysis data contains required fields."""
        return True

