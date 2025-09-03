"""
Data Formatter Utility

This module provides utilities for formatting data for HTML reports.
"""

import logging
from typing import Dict, Any, List, Union
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

class DataFormatter:
    """Formats data for HTML report generation."""
    
    def __init__(self):
        """Initialize the data formatter."""
        pass
    
    def format_protein_data(self, protein_data: Dict[str, Any]) -> Dict[str, Any]:
        """Format protein data for display."""
        formatted = {}
        for key, value in protein_data.items():
            if isinstance(value, np.ndarray):
                formatted[key] = value.tolist()
            elif isinstance(value, float):
                formatted[key] = round(value, 3)
            else:
                formatted[key] = value
        return formatted
    
    def format_similarity_scores(self, similarities: List[float]) -> List[float]:
        """Format similarity scores for display."""
        return [round(score, 3) for score in similarities]
    
    def format_family_data(self, family_data: Dict[str, Any]) -> Dict[str, Any]:
        """Format family data for display."""
        formatted = {}
        for family_id, data in family_data.items():
            if isinstance(data, dict):
                formatted[family_id] = self.format_protein_data(data)
            else:
                formatted[family_id] = data
        return formatted
    
    def create_summary_table(self, pipeline_results: Dict[str, Any]) -> pd.DataFrame:
        """Create a summary table from pipeline results."""
        summary_data = []
        
        # Extract protein information
        proteins = pipeline_results.get('proteins', [])
        for protein in proteins:
            summary_data.append({
                'Protein ID': protein.get('protein_id', 'Unknown'),
                'Family': protein.get('family', 'Unknown'),
                'Similarity Score': round(protein.get('similarity_score', 0), 3),
                'Sequence Length': protein.get('sequence_length', 0)
            })
        
        return pd.DataFrame(summary_data)
