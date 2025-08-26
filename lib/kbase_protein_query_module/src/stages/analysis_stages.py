"""
Integrated Analysis Stages for KBase Protein Query Module

This module provides an integrated wrapper that orchestrates all analysis stages
without performing analysis itself. It delegates to modular analysis implementations
under stages/analysis/.
"""

import logging
import time
from typing import Dict, Any, List, Optional, Union

from .base_stage import BaseStage, StageResult
from .analysis.sequence_analysis import SequenceAnalysisStage
from .analysis.network_analysis import NetworkAnalysisStage
from .analysis.bioinformatics_analysis import BioinformaticsAnalysisStage

logger = logging.getLogger(__name__)


class IntegratedAnalysisStage(BaseStage):
    """
    Integrated wrapper stage that orchestrates all analysis sub-stages without performing
    any analysis itself. It delegates to:
      - `SequenceAnalysisStage` (from analysis package)
      - `NetworkAnalysisStage` (from analysis package) 
      - `BioinformaticsAnalysisStage` (from analysis package)

    This provides a single entrypoint for running analysis in sequence.
    """

    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self._sequence_stage_cls = SequenceAnalysisStage
        self._network_stage_cls = NetworkAnalysisStage
        self._bio_stage_cls = BioinformaticsAnalysisStage

    def get_stage_name(self) -> str:
        return "analysis"

    def get_required_inputs(self) -> List[str]:
        # Minimum needed to start sequence analysis; downstream stages handle
        # their own validation based on what is available
        return ['protein_records']

    def get_optional_inputs(self) -> List[str]:
        return [
            'analysis_config', 'workspace_client',
            # For network analysis
            'embeddings', 'similarity_results', 'network_config',
            # For bioinformatics analysis
            'bioinformatics_config'
        ]

    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        return 'protein_records' in input_data and isinstance(input_data['protein_records'], list) and len(input_data['protein_records']) > 0

    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'analysis': {
                'type': 'object',
                'properties': {
                    'sequence_analyses': {'type': 'object'},
                    'network_analysis': {'type': 'object'},
                    'bioinformatics_analyses': {'type': 'object'}
                }
            }
        }

    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        start_time = time.time()
        aggregate_output: Dict[str, Any] = {}
        metadata: Dict[str, Any] = {}
        warnings: List[str] = []

        try:
            # 1) Sequence analysis
            seq_stage = self._sequence_stage_cls(self.config.get('sequence_analysis', {}))
            seq_result = seq_stage.execute(input_data, workspace_client)
            aggregate_output.update(seq_result.output_data)
            metadata['sequence_analysis'] = seq_result.metadata
            if not seq_result.success:
                warnings.append(f"Sequence analysis reported failure: {seq_result.error_message}")

            # Prepare inputs for downstream analysis
            current_data = input_data.copy()
            current_data.update(seq_result.output_data)

            # 2) Network analysis (only if similarity results and embeddings are provided)
            if 'similarity_results' in current_data and 'embeddings' in current_data:
                net_stage = self._network_stage_cls(self.config.get('network_analysis', {}))
                net_result = net_stage.execute(current_data, workspace_client)
                aggregate_output.update(net_result.output_data)
                metadata['network_analysis'] = net_result.metadata
                if not net_result.success:
                    warnings.append(f"Network analysis reported failure: {net_result.error_message}")

            # 3) Bioinformatics analysis (only if sequence analyses present)
            if 'sequence_analyses' in current_data:
                bio_stage = self._bio_stage_cls(self.config.get('bioinformatics_analysis', {}))
                bio_result = bio_stage.execute(current_data, workspace_client)
                aggregate_output.update(bio_result.output_data)
                metadata['bioinformatics_analysis'] = bio_result.metadata
                if not bio_result.success:
                    warnings.append(f"Bioinformatics analysis reported failure: {bio_result.error_message}")

            return StageResult(
                success=('sequence_analyses' in aggregate_output),
                output_data=aggregate_output,
                metadata=metadata,
                execution_time=time.time() - start_time,
                warnings=warnings
            )
        except Exception as e:
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=time.time() - start_time,
                error_message=str(e)
            )
