"""
Example of Adding New Analysis

This demonstrates how easy it is to add a new analysis type to the module.
This analysis calculates protein hydrophobicity patterns using the Kyte-Doolittle scale.

CLASS LOCATION: lib/kbase_protein_query_module/src/stages/analysis/demo_hydrophobicity_analysis.py:12
EXTENDS: BaseAnalysis
USED BY: WorkflowOrchestrator when "hydrophobicity_analysis" is requested
REGISTERED AS: "hydrophobicity_analysis"
"""

import os
import time
import logging
from typing import Dict, Any, List, Tuple
import numpy as np

from ...core.analysis_registry import BaseAnalysis, AnalysisMetadata, register_analysis

logger = logging.getLogger(__name__)

@register_analysis("hydrophobicity_analysis")  # This single decorator makes it available!
class HydrophobicityAnalysis(BaseAnalysis):
    """
    Protein hydrophobicity analysis using the Kyte-Doolittle scale.
    
    This analysis demonstrates the extensibility framework - it was added with
    just this single file and is automatically integrated into the entire pipeline.
    
    Features:
    - Calculates average hydrophobicity per protein
    - Identifies hydrophobic and hydrophilic regions
    - Generates hydrophobicity profiles
    - Creates CSV data exports and HTML visualizations
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        
        # Kyte-Doolittle hydrophobicity scale
        self.hydro_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        # Configuration parameters
        self.window_size = config.get('hydrophobic_window_size', 7) if config else 7
        self.hydrophobic_threshold = config.get('hydrophobic_threshold', 1.0) if config else 1.0
        self.hydrophilic_threshold = config.get('hydrophilic_threshold', -1.0) if config else -1.0
    
    def get_metadata(self) -> AnalysisMetadata:
        """Define analysis metadata for the registry."""
        return AnalysisMetadata(
            name="Hydrophobicity Analysis",
            description="Analyzes protein hydrophobicity patterns using the Kyte-Doolittle scale to identify hydrophobic and hydrophilic regions",
            version="1.0.0",
            author="KBase Development Team",
            output_files=["hydrophobicity_analysis.html", "hydrophobicity_data.csv", "hydrophobicity_summary.json"],
            dependencies=["numpy"],
            category="sequence_properties",
            computational_complexity="low",
            memory_requirements="low",
            supports_batch=True,
            supports_streaming=True
        )
    
    def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
        """
        Perform hydrophobicity analysis on protein sequences.
        
        Args:
            proteins: List of protein records with 'id' and 'sequence' fields
            **kwargs: Additional parameters (pipeline_results, analysis_config, etc.)
            
        Returns:
            Dictionary containing hydrophobicity analysis results
        """
        self.logger.info(f"Starting hydrophobicity analysis for {len(proteins)} proteins")
        start_time = time.time()
        
        results = {}
        summary_stats = {
            'total_proteins': len(proteins),
            'successful_analyses': 0,
            'avg_hydrophobicity': 0.0,
            'hydrophobic_proteins': 0,
            'hydrophilic_proteins': 0
        }
        
        all_hydrophobicity_values = []
        
        for protein in proteins:
            protein_id = protein.get('id', 'unknown')
            sequence = protein.get('sequence', '')
            
            if sequence:
                try:
                    # Calculate hydrophobicity profile
                    hydro_profile = self._calculate_hydrophobicity_profile(sequence)
                    
                    # Find regions of interest
                    hydrophobic_regions = self._find_hydrophobic_regions(hydro_profile)
                    hydrophilic_regions = self._find_hydrophilic_regions(hydro_profile)
                    
                    # Calculate statistics
                    avg_hydro = np.mean(hydro_profile)
                    all_hydrophobicity_values.append(avg_hydro)
                    
                    # Classify protein
                    if avg_hydro > 0.5:
                        summary_stats['hydrophobic_proteins'] += 1
                        classification = 'hydrophobic'
                    elif avg_hydro < -0.5:
                        summary_stats['hydrophilic_proteins'] += 1
                        classification = 'hydrophilic'
                    else:
                        classification = 'mixed'
                    
                    results[protein_id] = {
                        'avg_hydrophobicity': float(avg_hydro),
                        'max_hydrophobicity': float(np.max(hydro_profile)),
                        'min_hydrophobicity': float(np.min(hydro_profile)),
                        'std_hydrophobicity': float(np.std(hydro_profile)),
                        'classification': classification,
                        'hydrophobic_regions': hydrophobic_regions,
                        'hydrophilic_regions': hydrophilic_regions,
                        'sequence_length': len(sequence),
                        'hydrophobicity_profile': hydro_profile.tolist()  # For detailed analysis
                    }
                    
                    summary_stats['successful_analyses'] += 1
                    
                except Exception as e:
                    self.logger.warning(f"Failed to analyze protein {protein_id}: {e}")
                    results[protein_id] = {
                        'error': str(e),
                        'status': 'failed'
                    }
        
        # Calculate overall statistics
        if all_hydrophobicity_values:
            summary_stats['avg_hydrophobicity'] = float(np.mean(all_hydrophobicity_values))
            summary_stats['std_hydrophobicity'] = float(np.std(all_hydrophobicity_values))
        
        execution_time = time.time() - start_time
        
        self.logger.info(f"Hydrophobicity analysis completed in {execution_time:.2f}s")
        
        return {
            'analysis_type': 'hydrophobicity_analysis',
            'results': results,
            'summary': f"Analyzed hydrophobicity patterns for {summary_stats['successful_analyses']}/{len(proteins)} proteins",
            'statistics': summary_stats,
            'execution_time': execution_time,
            'metadata': {
                'window_size': self.window_size,
                'hydrophobic_threshold': self.hydrophobic_threshold,
                'hydrophilic_threshold': self.hydrophilic_threshold,
                'scale_used': 'Kyte-Doolittle'
            }
        }
    
    def get_output_files(self, output_dir: str) -> List[str]:
        """Define output files this analysis will create."""
        return [
            os.path.join(output_dir, "hydrophobicity_analysis.html"),
            os.path.join(output_dir, "hydrophobicity_data.csv"),
            os.path.join(output_dir, "hydrophobicity_summary.json")
        ]
    
    def validate_input(self, proteins: List[Any]) -> bool:
        """Validate input data for hydrophobicity analysis."""
        if not proteins:
            return False
        
        # Check that we have sequences to analyze
        valid_proteins = 0
        for protein in proteins:
            sequence = protein.get('sequence', '')
            if sequence and len(sequence) > 5:  # Minimum sequence length
                valid_proteins += 1
        
        return valid_proteins > 0
    
    def estimate_resources(self, proteins: List[Any]) -> Dict[str, Any]:
        """Estimate computational resources for this analysis."""
        protein_count = len(proteins)
        avg_sequence_length = 300  # Typical protein length
        
        return {
            "estimated_time_seconds": protein_count * 0.01,  # Very fast analysis
            "estimated_memory_mb": protein_count * 0.1,      # Minimal memory usage
            "cpu_cores_recommended": 1,                      # Single-threaded is fine
            "supports_parallel": True                        # Can be parallelized
        }
    
    def _calculate_hydrophobicity_profile(self, sequence: str) -> np.ndarray:
        """Calculate hydrophobicity values for each amino acid in the sequence."""
        hydro_values = []
        
        for aa in sequence.upper():
            if aa in self.hydro_scale:
                hydro_values.append(self.hydro_scale[aa])
            else:
                # Unknown amino acid - use neutral value
                hydro_values.append(0.0)
        
        return np.array(hydro_values)
    
    def _find_hydrophobic_regions(self, hydro_profile: np.ndarray) -> List[Tuple[int, int, float]]:
        """
        Find hydrophobic regions using sliding window analysis.
        
        Returns:
            List of (start, end, avg_hydrophobicity) tuples
        """
        regions = []
        
        if len(hydro_profile) < self.window_size:
            return regions
        
        for i in range(len(hydro_profile) - self.window_size + 1):
            window = hydro_profile[i:i + self.window_size]
            window_avg = np.mean(window)
            
            if window_avg > self.hydrophobic_threshold:
                # Check if this extends an existing region
                if regions and regions[-1][1] >= i - 1:
                    # Extend the last region
                    start = regions[-1][0]
                    regions[-1] = (start, i + self.window_size, 
                                 np.mean(hydro_profile[start:i + self.window_size]))
                else:
                    # New region
                    regions.append((i, i + self.window_size, window_avg))
        
        return regions
    
    def _find_hydrophilic_regions(self, hydro_profile: np.ndarray) -> List[Tuple[int, int, float]]:
        """
        Find hydrophilic regions using sliding window analysis.
        
        Returns:
            List of (start, end, avg_hydrophobicity) tuples
        """
        regions = []
        
        if len(hydro_profile) < self.window_size:
            return regions
        
        for i in range(len(hydro_profile) - self.window_size + 1):
            window = hydro_profile[i:i + self.window_size]
            window_avg = np.mean(window)
            
            if window_avg < self.hydrophilic_threshold:
                # Check if this extends an existing region
                if regions and regions[-1][1] >= i - 1:
                    # Extend the last region
                    start = regions[-1][0]
                    regions[-1] = (start, i + self.window_size,
                                 np.mean(hydro_profile[start:i + self.window_size]))
                else:
                    # New region
                    regions.append((i, i + self.window_size, window_avg))
        
        return regions
