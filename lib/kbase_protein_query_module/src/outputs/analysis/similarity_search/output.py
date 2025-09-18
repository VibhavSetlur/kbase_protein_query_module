"""
Similarity Search Analysis Output Handler

Handles output generation for similarity search analysis results.
"""

import json
import csv
from typing import Dict, Any, List
from ...base_analysis_output import BaseAnalysisOutput, AnalysisOutputResult
from ...analysis_output_registry import register_analysis_output


@register_analysis_output("similarity_search")
class SimilaritySearchOutput(BaseAnalysisOutput):
    """Output handler for similarity search analysis."""
    
    def generate_outputs(self, analysis_data: Dict[str, Any], stage_name: str) -> AnalysisOutputResult:
        """Generate similarity search outputs."""
        output_files = []
        metadata = {}
        
        try:
            # Extract similarity results
            top_matches = analysis_data.get('top_matches', [])
            query_proteins = analysis_data.get('query_proteins', [])
            similarity_scores = analysis_data.get('similarity_scores', {})
            
            # Generate JSON summary
            summary_data = {
                'query_count': len(query_proteins),
                'matches_found': len(top_matches),
                'top_matches': top_matches[:10],  # Top 10 for summary
                'similarity_threshold': analysis_data.get('similarity_threshold', 0.7),
                'analysis_timestamp': analysis_data.get('timestamp'),
                'stage_name': stage_name
            }
            
            json_path = self.output_manager.write_json(
                'similarity_search',
                'similarity_summary.json',
                summary_data,
                stage=stage_name,
                description='Similarity search summary with top matches'
            )
            output_files.append(json_path)
            
            # Generate detailed CSV results
            if top_matches:
                csv_data = self._generate_similarity_csv(top_matches, similarity_scores)
                csv_path = self.output_manager.write_csv(
                    'similarity_search',
                    'similarity_results.csv',
                    csv_data,
                    stage=stage_name,
                    description='Detailed similarity search results'
                )
                output_files.append(csv_path)
            
            # Generate protein list for downstream analysis
            if query_proteins:
                protein_list = {
                    'query_proteins': query_proteins,
                    'analysis_type': 'similarity_search',
                    'stage_name': stage_name
                }
                protein_path = self.output_manager.write_json(
                    'similarity_search',
                    'query_proteins.json',
                    protein_list,
                    stage=stage_name,
                    description='Query proteins used in similarity search'
                )
                output_files.append(protein_path)
            
            metadata = {
                'query_count': len(query_proteins),
                'matches_found': len(top_matches),
                'output_files': len(output_files),
                'analysis_type': 'similarity_search'
            }
            
            summary = f"Similarity search completed: {len(query_proteins)} queries, {len(top_matches)} matches found"
            
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
                summary=f"Error generating similarity search outputs: {str(e)}",
                error_message=str(e)
            )
    
    def _generate_similarity_csv(self, top_matches: List[Dict], similarity_scores: Dict[str, float]) -> str:
        """Generate CSV data for similarity results."""
        import io
        output = io.StringIO()
        writer = csv.writer(output)
        
        # Header
        writer.writerow(['query_protein', 'match_protein', 'similarity_score', 'rank', 'family', 'description'])
        
        # Data rows
        for i, match in enumerate(top_matches, 1):
            query_protein = match.get('query_protein', '')
            match_protein = match.get('protein_id', '')
            score = match.get('similarity_score', 0.0)
            family = match.get('family', '')
            description = match.get('description', '')
            
            writer.writerow([query_protein, match_protein, score, i, family, description])
        
        return output.getvalue()
    
    def get_output_description(self) -> str:
        """Get description of similarity search outputs."""
        return "Similarity search results including top matches, CSV data, and query protein lists"
    
    def get_supported_formats(self) -> List[str]:
        """Get supported output formats."""
        return ['json', 'csv']
    
    def validate_analysis_data(self, analysis_data: Dict[str, Any]) -> bool:
        """Validate similarity search data."""
        required_fields = ['top_matches', 'query_proteins']
        return all(field in analysis_data for field in required_fields)
