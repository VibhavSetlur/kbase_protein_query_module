"""
Bioinformatics Analysis Output Handler

Handles output generation for bioinformatics analysis results.
"""

from typing import Dict, Any, List
from ...base_analysis_output import BaseAnalysisOutput, AnalysisOutputResult
from ...analysis_output_registry import register_analysis_output


@register_analysis_output("bioinformatics_analysis")
class BioinformaticsAnalysisOutput(BaseAnalysisOutput):
    """Output handler for bioinformatics analysis."""

    def generate_outputs(self, analysis_data: Dict[str, Any], stage_name: str) -> AnalysisOutputResult:
        output_files: List[str] = []
        try:
            summary_path = self.output_manager.write_json(
                'bioinformatics_analysis',
                'bioinformatics_summary.json',
                {
                    'metrics': analysis_data.get('metrics', {}),
                    'enrichment': analysis_data.get('enrichment', {}),
                    'annotations': analysis_data.get('annotations', {}),
                    'stage_name': stage_name
                },
                stage=stage_name,
                description='Bioinformatics analysis summary'
            )
            output_files.append(summary_path)

            return AnalysisOutputResult(
                success=True,
                output_files=output_files,
                metadata={'analysis_type': 'bioinformatics_analysis', 'output_files': len(output_files)},
                summary='Bioinformatics analysis outputs generated'
            )
        except Exception as e:
            return AnalysisOutputResult(
                success=False,
                output_files=[],
                metadata={},
                summary=f'Error generating bioinformatics analysis outputs: {e}',
                error_message=str(e)
            )

    def get_output_description(self) -> str:
        return 'Bioinformatics metrics, enrichment results, and annotations summary'

    def get_supported_formats(self) -> List[str]:
        return ['json']


