"""
Report Generation Stage for KBase Protein Query Module

This stage generates analysis reports.
"""

from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from ...reports.html.generator import HTMLReportGenerator
import time

class ReportGenerationStage(BaseStage):
    """Report generation stage that produces comprehensive HTML reports."""
    
    def get_stage_name(self) -> str:
        return "report_generation"
    
    def get_required_inputs(self) -> List[str]:
        return ['pipeline_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['report_config']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'report_file': {
                'type': 'string',
                'description': 'Path to generated report file'
            },
            'report_summary': {
                'type': 'object',
                'description': 'Summary of report contents'
            }
        }
    
    def run(self, input_data, workspace_client=None):
        start = time.time()
        try:
            results = input_data['pipeline_results']
            protein_id = input_data.get('protein_id')
            generator = HTMLReportGenerator()
            report_info = generator.generate_comprehensive_report(results, protein_id=protein_id)
            exec_time = time.time() - start
            return StageResult(
                success=True,
                output_data={'report': report_info},
                metadata={'report_path': report_info.get('html_path')},
                execution_time=exec_time
            )
        except Exception as e:
            return StageResult(success=False, output_data={}, metadata={}, execution_time=time.time()-start, error_message=str(e))
