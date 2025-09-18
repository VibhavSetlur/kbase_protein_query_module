"""
Sequence Analysis Output Handler

Handles output generation for sequence analysis results.
"""

import json
from typing import Dict, Any, List
from ...base_analysis_output import BaseAnalysisOutput, AnalysisOutputResult
from ...analysis_output_registry import register_analysis_output


@register_analysis_output("sequence_analysis")
class SequenceAnalysisOutput(BaseAnalysisOutput):
    """Output handler for sequence analysis."""
    
    def generate_outputs(self, analysis_data: Dict[str, Any], stage_name: str) -> AnalysisOutputResult:
        """Generate sequence analysis outputs."""
        output_files = []
        metadata = {}
        
        try:
            # Extract sequence analysis results
            sequence_features = analysis_data.get('sequence_features', {})
            conserved_regions = analysis_data.get('conserved_regions', [])
            motifs = analysis_data.get('motifs', [])
            domain_annotations = analysis_data.get('domain_annotations', {})
            
            # Generate sequence features summary
            features_data = {
                'sequence_features': sequence_features,
                'conserved_regions': conserved_regions,
                'motifs': motifs,
                'domain_annotations': domain_annotations,
                'analysis_timestamp': analysis_data.get('timestamp'),
                'stage_name': stage_name
            }
            
            features_path = self.output_manager.write_json(
                'sequence_analysis',
                'sequence_features.json',
                features_data,
                stage=stage_name,
                description='Sequence features, motifs, and domain annotations'
            )
            output_files.append(features_path)
            
            # Generate conserved regions analysis
            if conserved_regions:
                conserved_data = {
                    'conserved_regions': conserved_regions,
                    'region_count': len(conserved_regions),
                    'analysis_method': analysis_data.get('analysis_method', 'unknown'),
                    'stage_name': stage_name
                }
                conserved_path = self.output_manager.write_json(
                    'sequence_analysis',
                    'conserved_regions.json',
                    conserved_data,
                    stage=stage_name,
                    description='Conserved sequence regions analysis'
                )
                output_files.append(conserved_path)
            
            # Generate motif analysis
            if motifs:
                motif_data = {
                    'motifs': motifs,
                    'motif_count': len(motifs),
                    'motif_types': list(set(m.get('type', 'unknown') for m in motifs)),
                    'stage_name': stage_name
                }
                motif_path = self.output_manager.write_json(
                    'sequence_analysis',
                    'motif_analysis.json',
                    motif_data,
                    stage=stage_name,
                    description='Protein motif analysis results'
                )
                output_files.append(motif_path)
            
            # Generate domain annotations
            if domain_annotations:
                domain_path = self.output_manager.write_json(
                    'sequence_analysis',
                    'domain_annotations.json',
                    domain_annotations,
                    stage=stage_name,
                    description='Protein domain annotations'
                )
                output_files.append(domain_path)
            
            metadata = {
                'features_analyzed': len(sequence_features),
                'conserved_regions': len(conserved_regions),
                'motifs_found': len(motifs),
                'domains_annotated': len(domain_annotations),
                'output_files': len(output_files),
                'analysis_type': 'sequence_analysis'
            }
            
            summary = f"Sequence analysis completed: {len(sequence_features)} features, {len(conserved_regions)} conserved regions, {len(motifs)} motifs"
            
            return AnalysisOutputResult(
                success=True,
                output_files=output_files,
                metadata=metadata,
                summary=summary
            )
            
        except Exception as e:
            return AnalysisOutputResult(
                success=False,
                output_files=[],
                metadata={},
                summary=f"Error generating sequence analysis outputs: {str(e)}",
                error_message=str(e)
            )
    
    def get_output_description(self) -> str:
        """Get description of sequence analysis outputs."""
        return "Sequence features, conserved regions, motifs, and domain annotations"
    
    def get_supported_formats(self) -> List[str]:
        """Get supported output formats."""
        return ['json']
    
    def validate_analysis_data(self, analysis_data: Dict[str, Any]) -> bool:
        """Validate sequence analysis data."""
        # At least one type of analysis should be present
        analysis_types = ['sequence_features', 'conserved_regions', 'motifs', 'domain_annotations']
        return any(field in analysis_data for field in analysis_types)
