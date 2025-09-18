"""
Bioinformatics Analysis Stage for KBase Protein Query Module

This stage performs bioinformatics analysis on proteins.
"""

from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from .sequence_analysis import SequenceAnalysisStage
import time

class BioinformaticsAnalysisStage(BaseStage):
    """Bioinformatics analysis stage leveraging SequenceAnalysisStage links and motifs."""
    
    def get_stage_name(self) -> str:
        return "bioinformatics_analysis"
    
    def get_required_inputs(self) -> List[str]:
        return ['sequence_analysis']
    
    def get_optional_inputs(self) -> List[str]:
        return ['bioinformatics_config']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'bioinformatics_results': {
                'type': 'object',
                'properties': {
                    'analysis_results': {'type': 'object'},
                    'metadata': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data, workspace_client=None):
        start = time.time()
        try:
            seq_analysis = input_data['sequence_analysis']
            analysis_results = seq_analysis.get('analysis_results', {})
            # Aggregate bioinformatics links and motif summaries
            bio_links = {}
            motif_summary = {}
            for pid, res in analysis_results.items():
                bio_links[pid] = res.get('bioinformatics_links', {})
                motifs = res.get('sequence_motifs', {})
                motif_summary[pid] = {k: len(v or []) for k, v in motifs.items()}
            output = {
                'bioinformatics_results': {
                    'links': bio_links,
                    'motifs': motif_summary
                }
            }
            return StageResult(success=True, output_data=output, metadata={'count': len(bio_links)}, execution_time=time.time()-start)
        except Exception as e:
            return StageResult(success=False, output_data={}, metadata={}, execution_time=time.time()-start, error_message=str(e))
