"""
Multi-Protein Analysis Output Handler

Handles output generation for multi-protein analysis results.
"""

from typing import Dict, Any, List
from ...base_analysis_output import BaseAnalysisOutput, AnalysisOutputResult
from ...analysis_output_registry import register_analysis_output


@register_analysis_output("multi_protein_analysis")
class MultiProteinAnalysisOutput(BaseAnalysisOutput):
    """Output handler for multi-protein analysis."""

    def generate_outputs(self, analysis_data: Dict[str, Any], stage_name: str) -> AnalysisOutputResult:
        output_files: List[str] = []
        try:
            cohort = analysis_data.get('cohort', [])
            aggregates = analysis_data.get('aggregates', {})
            comparison = analysis_data.get('comparison', {})

            summary_path = self.output_manager.write_json(
                'multi_protein_analysis',
                'multi_protein_summary.json',
                {
                    'cohort_size': len(cohort),
                    'aggregates': aggregates,
                    'comparison': comparison,
                    'stage_name': stage_name
                },
                stage=stage_name,
                description='Multi-protein analysis summary and aggregates'
            )
            output_files.append(summary_path)

            return AnalysisOutputResult(
                success=True,
                output_files=output_files,
                metadata={'analysis_type': 'multi_protein_analysis', 'output_files': len(output_files)},
                summary='Multi-protein analysis outputs generated'
            )
        except Exception as e:
            return AnalysisOutputResult(
                success=False,
                output_files=[],
                metadata={},
                summary=f'Error generating multi-protein analysis outputs: {e}',
                error_message=str(e)
            )

    def get_output_description(self) -> str:
        return 'Multi-protein aggregates and comparative analysis summary'

    def get_supported_formats(self) -> List[str]:
        return ['json']


