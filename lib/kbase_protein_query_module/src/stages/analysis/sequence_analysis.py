"""
Sequence Analysis Stage for KBase Protein Query Module

This stage performs comprehensive sequence analysis on proteins using the
ProteinSequenceAnalyzer for detailed characterization.
"""

import time
import logging
from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from ...analysis.sequence_analyzer import ProteinSequenceAnalyzer

logger = logging.getLogger(__name__)

class SequenceAnalysisStage(BaseStage):
    """
    Comprehensive sequence analysis stage that performs detailed protein characterization.
    
    This stage analyzes protein sequences using:
    - Amino acid composition analysis
    - Physicochemical properties calculation
    - Secondary structure prediction
    - Sequence motif identification
    - Bioinformatics database links generation
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """Initialize the sequence analysis stage."""
        super().__init__(config)
        self.analyzer = ProteinSequenceAnalyzer()
    
    def get_stage_name(self) -> str:
        return "sequence_analysis"
    
    def get_required_inputs(self) -> List[str]:
        return ['protein_records']
    
    def get_optional_inputs(self) -> List[str]:
        return ['analysis_config']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data for sequence analysis."""
        if 'protein_records' not in input_data:
            return False
        
        protein_records = input_data['protein_records']
        if not protein_records:
            logger.warning("No protein records provided for sequence analysis")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'sequence_analysis': {
                'type': 'object',
                'properties': {
                    'analysis_results': {'type': 'object'},
                    'metadata': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Execute comprehensive sequence analysis."""
        start_time = time.time()
        
        try:
            protein_records = input_data['protein_records']
            analysis_results = {}
            
            logger.info(f"Starting sequence analysis for {len(protein_records)} proteins")
            
            # Analyze each protein sequence
            for record in protein_records:
                try:
                    protein_id = record.protein_id
                    sequence = record.sequence
                    
                    logger.debug(f"Analyzing sequence for protein: {protein_id}")
                    
                    # Perform comprehensive analysis
                    analysis = self.analyzer.analyze_sequence(sequence, protein_id)
                    analysis_results[protein_id] = analysis
                    
                except Exception as e:
                    logger.error(f"Failed to analyze protein {record.protein_id}: {str(e)}")
                    analysis_results[record.protein_id] = {
                        'error': str(e),
                        'protein_id': record.protein_id,
                        'sequence': record.sequence
                    }
            
            # Calculate summary statistics
            summary_stats = self._calculate_summary_statistics(analysis_results)
            
            execution_time = time.time() - start_time
            
            output_data = {
                'sequence_analysis': {
                    'analysis_results': analysis_results,
                    'summary_statistics': summary_stats,
                    'metadata': {
                        'num_proteins_analyzed': len(analysis_results),
                        'successful_analyses': len([r for r in analysis_results.values() if 'error' not in r]),
                        'failed_analyses': len([r for r in analysis_results.values() if 'error' in r]),
                        'execution_time': execution_time
                    }
                }
            }
            
            logger.info(f"Sequence analysis completed in {execution_time:.2f}s")
            
            return StageResult(
                success=True,
                output_data=output_data,
                metadata={'stage': 'sequence_analysis', 'execution_time': execution_time},
                execution_time=execution_time
            )
            
        except Exception as e:
            logger.error(f"Sequence analysis stage failed: {str(e)}")
            return StageResult(
                success=False,
                output_data={},
                metadata={'error': str(e)},
                execution_time=time.time() - start_time
            )
    
    def _calculate_summary_statistics(self, analysis_results: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate summary statistics across all analyzed proteins."""
        successful_analyses = [r for r in analysis_results.values() if 'error' not in r]
        
        if not successful_analyses:
            return {'error': 'No successful analyses to summarize'}
        
        # Collect statistics
        lengths = [r['length'] for r in successful_analyses]
        molecular_weights = [r['physicochemical_properties']['molecular_weight'] for r in successful_analyses]
        isoelectric_points = [r['physicochemical_properties']['isoelectric_point'] for r in successful_analyses]
        hydrophobicity_scores = [r['physicochemical_properties']['average_hydrophobicity'] for r in successful_analyses]
        
        # Calculate summary statistics
        summary = {
            'sequence_length': {
                'mean': sum(lengths) / len(lengths),
                'min': min(lengths),
                'max': max(lengths),
                'median': sorted(lengths)[len(lengths)//2]
            },
            'molecular_weight': {
                'mean': sum(molecular_weights) / len(molecular_weights),
                'min': min(molecular_weights),
                'max': max(molecular_weights),
                'median': sorted(molecular_weights)[len(molecular_weights)//2]
            },
            'isoelectric_point': {
                'mean': sum(isoelectric_points) / len(isoelectric_points),
                'min': min(isoelectric_points),
                'max': max(isoelectric_points),
                'median': sorted(isoelectric_points)[len(isoelectric_points)//2]
            },
            'hydrophobicity': {
                'mean': sum(hydrophobicity_scores) / len(hydrophobicity_scores),
                'min': min(hydrophobicity_scores),
                'max': max(hydrophobicity_scores),
                'median': sorted(hydrophobicity_scores)[len(hydrophobicity_scores)//2]
            },
            'total_proteins': len(analysis_results),
            'successful_analyses': len(successful_analyses),
            'failed_analyses': len(analysis_results) - len(successful_analyses)
        }
        
        return summary
