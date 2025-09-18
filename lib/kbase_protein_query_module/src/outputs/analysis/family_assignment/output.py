"""
Family Assignment Analysis Output Handler

Handles output generation for protein family assignment analysis results.
"""

import json
from typing import Dict, Any, List
from ...base_analysis_output import BaseAnalysisOutput, AnalysisOutputResult
from ...analysis_output_registry import register_analysis_output


@register_analysis_output("family_assignment")
class FamilyAssignmentOutput(BaseAnalysisOutput):
    """Output handler for family assignment analysis."""
    
    def generate_outputs(self, analysis_data: Dict[str, Any], stage_name: str) -> AnalysisOutputResult:
        """Generate family assignment outputs."""
        output_files = []
        metadata = {}
        
        try:
            # Extract family assignment results
            family_assignments = analysis_data.get('family_assignments', {})
            protein_families = analysis_data.get('protein_families', [])
            confidence_scores = analysis_data.get('confidence_scores', {})
            
            # Generate family assignment summary
            summary_data = {
                'total_proteins': len(family_assignments),
                'families_found': len(set(family_assignments.values())),
                'family_assignments': family_assignments,
                'confidence_scores': confidence_scores,
                'analysis_timestamp': analysis_data.get('timestamp'),
                'stage_name': stage_name
            }
            
            json_path = self.output_manager.write_json(
                'family_assignment',
                'family_assignments.json',
                summary_data,
                stage=stage_name,
                description='Protein family assignments with confidence scores'
            )
            output_files.append(json_path)
            
            # Generate family statistics
            family_stats = self._generate_family_statistics(family_assignments, confidence_scores)
            stats_path = self.output_manager.write_json(
                'family_assignment',
                'family_statistics.json',
                family_stats,
                stage=stage_name,
                description='Family assignment statistics and distributions'
            )
            output_files.append(stats_path)
            
            # Generate protein-to-family mapping
            if protein_families:
                mapping_data = {
                    'protein_families': protein_families,
                    'assignment_method': analysis_data.get('assignment_method', 'unknown'),
                    'stage_name': stage_name
                }
                mapping_path = self.output_manager.write_json(
                    'family_assignment',
                    'protein_family_mapping.json',
                    mapping_data,
                    stage=stage_name,
                    description='Detailed protein to family mapping'
                )
                output_files.append(mapping_path)
            
            metadata = {
                'total_proteins': len(family_assignments),
                'families_found': len(set(family_assignments.values())),
                'output_files': len(output_files),
                'analysis_type': 'family_assignment'
            }
            
            summary = f"Family assignment completed: {len(family_assignments)} proteins assigned to {len(set(family_assignments.values()))} families"
            
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
                summary=f"Error generating family assignment outputs: {str(e)}",
                error_message=str(e)
            )
    
    def _generate_family_statistics(self, family_assignments: Dict[str, str], confidence_scores: Dict[str, float]) -> Dict[str, Any]:
        """Generate family assignment statistics."""
        from collections import Counter
        
        # Count proteins per family
        family_counts = Counter(family_assignments.values())
        
        # Calculate confidence statistics
        confidences = list(confidence_scores.values())
        avg_confidence = sum(confidences) / len(confidences) if confidences else 0.0
        
        return {
            'family_distribution': dict(family_counts),
            'total_families': len(family_counts),
            'average_confidence': avg_confidence,
            'confidence_range': {
                'min': min(confidences) if confidences else 0.0,
                'max': max(confidences) if confidences else 0.0
            },
            'high_confidence_assignments': len([c for c in confidences if c > 0.8]),
            'low_confidence_assignments': len([c for c in confidences if c < 0.5])
        }
    
    def get_output_description(self) -> str:
        """Get description of family assignment outputs."""
        return "Protein family assignments with confidence scores, statistics, and detailed mappings"
    
    def get_supported_formats(self) -> List[str]:
        """Get supported output formats."""
        return ['json']
    
    def validate_analysis_data(self, analysis_data: Dict[str, Any]) -> bool:
        """Validate family assignment data."""
        required_fields = ['family_assignments']
        return all(field in analysis_data for field in required_fields)
