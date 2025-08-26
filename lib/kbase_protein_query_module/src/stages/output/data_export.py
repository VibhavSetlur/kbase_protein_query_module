"""
Data Export Stage for KBase Protein Query Module

This stage exports analysis results to various formats.
"""

from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
import os
import json
import time

class DataExportStage(BaseStage):
    """Data export stage that writes selected outputs to disk as JSON/CSV/HTML."""
    
    def get_stage_name(self) -> str:
        return "data_export"
    
    def get_required_inputs(self) -> List[str]:
        return ['pipeline_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['export_config']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'export_files': {
                'type': 'array',
                'description': 'List of exported data files'
            },
            'export_summary': {
                'type': 'object',
                'description': 'Summary of exported data'
            }
        }
    
    def run(self, input_data, workspace_client=None):
        start = time.time()
        results = input_data['pipeline_results']
        export_config = input_data.get('export_config', {})
        # Use consolidated output directory from environment or default
        output_dir = os.environ.get('EXPORTS_DIR', export_config.get('output_dir', 'exports'))
        os.makedirs(output_dir, exist_ok=True)
        exported = []
        try:
            # Export summary JSON
            summary_path = os.path.join(output_dir, 'pipeline_results_summary.json')
            with open(summary_path, 'w') as f:
                json.dump({k: list(v.keys()) if isinstance(v, dict) else v for k, v in results.items()}, f, indent=2)
            exported.append(summary_path)
            return StageResult(
                success=True,
                output_data={'export_files': exported, 'export_metadata': {'count': len(exported)}},
                metadata={'output_dir': output_dir},
                execution_time=time.time()-start
            )
        except Exception as e:
            return StageResult(success=False, output_data={}, metadata={'output_dir': output_dir}, execution_time=time.time()-start, error_message=str(e))
