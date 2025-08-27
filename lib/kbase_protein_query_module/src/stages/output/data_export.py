"""
Data Export Stage for KBase Protein Query Module

This stage exports analysis results to CSV, JSON, and other user-friendly formats.
"""

from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
import os
import json
import time
import pandas as pd
import numpy as np

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
        
        # Use the same output directory as report generation for consistency
        output_dir = input_data.get('output_directory')
        if not output_dir:
            # Use scratch space for KBase Narrative integration, test/outputs for local testing
            base_dir = os.environ.get('EXPORTS_DIR')
            if not base_dir:
                if os.path.exists('test/outputs'):
                    base_dir = 'test/outputs'
                else:
                    # For KBase Narrative, use scratch space
                    base_dir = os.environ.get('SCRATCH_DIR', '/kb/module/work/tmp')
            output_dir = base_dir
        
        os.makedirs(output_dir, exist_ok=True)
        
        exported_files = []
        
        try:
            # 1. Export comprehensive protein metadata CSV
            metadata_csv = self._create_protein_metadata_csv(results, input_data, output_dir)
            if metadata_csv:
                exported_files.append(metadata_csv)
            
            # 2. Export top similarity matches CSV with input protein correlation
            matches_csv = self._create_top_matches_csv(results, output_dir, input_data)
            if matches_csv:
                exported_files.append(matches_csv)
            
            # 3. Export detailed pipeline results JSON
            pipeline_json = self._create_detailed_pipeline_json(results, input_data, output_dir)
            if pipeline_json:
                exported_files.append(pipeline_json)
            
            # 4. Export analysis summary CSV for research use
            summary_csv = self._create_analysis_summary_csv(results, input_data, output_dir)
            if summary_csv:
                exported_files.append(summary_csv)
            
            return StageResult(
                success=True,
                output_data={
                    'export_files': exported_files, 
                    'export_metadata': {
                        'count': len(exported_files),
                        'output_directory': output_dir,
                        'file_types': ['CSV', 'JSON']
                    }
                },
                metadata={'output_dir': output_dir, 'files_exported': len(exported_files)},
                execution_time=time.time()-start
            )
        except Exception as e:
            return StageResult(success=False, output_data={}, metadata={'output_dir': output_dir}, execution_time=time.time()-start, error_message=str(e))
    
    def _create_protein_metadata_csv(self, results: Dict[str, Any], input_data: Dict[str, Any], output_dir: str) -> str:
        """Create comprehensive protein metadata CSV file."""
        try:
            # Extract input proteins
            input_proteins = input_data.get('input_proteins', [])
            if isinstance(input_proteins, list) and input_proteins:
                if isinstance(input_proteins[0], str):
                    proteins_data = [{'protein_id': f'protein_{i+1}', 'sequence': seq, 'source': 'user_input'} 
                                   for i, seq in enumerate(input_proteins)]
                else:
                    proteins_data = input_proteins
            else:
                proteins_data = []
            
            metadata_rows = []
            
            for protein in proteins_data:
                protein_id = protein.get('protein_id', 'unknown')
                sequence = protein.get('sequence', '')
                
                row = {
                    'protein_id': protein_id,
                    'sequence_length': len(sequence),
                    'source': protein.get('source', 'user_input'),
                    'input_type': 'user_input'
                }
                
                # Add sequence analysis data if available
                if 'sequence_analysis' in results:
                    seq_data = results['sequence_analysis'].get('output_data', {})
                    row.update({
                        'molecular_weight': seq_data.get('molecular_weight'),
                        'isoelectric_point': seq_data.get('isoelectric_point'),
                        'hydropathy': seq_data.get('hydropathy')
                    })
                
                # Add family assignment data if available
                if 'family_assignment' in results:
                    family_data = results['family_assignment'].get('output_data', {})
                    row.update({
                        'assigned_family': family_data.get('family_id'),
                        'family_confidence': family_data.get('confidence_score'),
                        'family_description': family_data.get('family_description')
                    })
                
                metadata_rows.append(row)
            
            # Create DataFrame and save as CSV
            if metadata_rows:
                df = pd.DataFrame(metadata_rows)
                csv_path = os.path.join(output_dir, 'protein_metadata.csv')
                df.to_csv(csv_path, index=False)
                return csv_path
            
        except Exception as e:
            print(f"Error creating protein metadata CSV: {e}")
        
        return None
    
    def _create_top_matches_csv(self, results: Dict[str, Any], output_dir: str, input_data: Dict[str, Any] = None) -> str:
        """Create comprehensive CSV file with top N proteins including metadata and their corresponding input proteins."""
        try:
            # Extract input proteins for reference
            input_proteins = input_data.get('input_proteins', []) if input_data else []
            
            # Collect top matches from similarity search
            matches_data = []
            
            if 'similarity_search' in results:
                sim_data = results['similarity_search'].get('output_data', {})
                top_matches = sim_data.get('top_matches', [])
                
                for i, match in enumerate(top_matches[:50]):  # Top 50 matches
                    match_entry = {
                        'rank': i + 1,
                        'protein_id': match.get('protein_id', f'match_{i+1}'),
                        'similarity_score': round(match.get('similarity_score', 0), 4),
                        'family': match.get('family', 'Unknown'),
                        'description': match.get('description', 'No description available'),
                        'database': match.get('database', 'Internal'),
                        'source_type': 'database_match',
                        'match_type': 'similarity_search',
                        'sequence_length': match.get('sequence_length', 'N/A'),
                        'molecular_weight': match.get('molecular_weight', 'N/A'),
                        'isoelectric_point': match.get('isoelectric_point', 'N/A'),
                        'input_protein_index': match.get('query_index', 0),  # Which input protein this matches
                        'input_protein_sequence': input_proteins[match.get('query_index', 0)] if match.get('query_index', 0) < len(input_proteins) else 'N/A'
                    }
                    matches_data.append(match_entry)
            
            # Add input proteins as reference entries
            for i, protein_seq in enumerate(input_proteins):
                input_entry = {
                    'rank': 'INPUT',
                    'protein_id': f'input_protein_{i+1}',
                    'similarity_score': 1.0,  # Perfect match to itself
                    'family': results.get('family_assignment', {}).get('output_data', {}).get('family_id', 'Unknown'),
                    'description': 'User input protein',
                    'database': 'User Input',
                    'source_type': 'user_input',
                    'match_type': 'query_protein',
                    'sequence_length': len(protein_seq) if isinstance(protein_seq, str) else 'N/A',
                    'molecular_weight': 'N/A',  # Could be calculated if needed
                    'isoelectric_point': 'N/A',  # Could be calculated if needed
                    'input_protein_index': i,
                    'input_protein_sequence': protein_seq if isinstance(protein_seq, str) else str(protein_seq)
                }
                matches_data.append(input_entry)
            
            if not matches_data:
                # Create minimal entry if no data available
                matches_data.append({
                    'rank': 1,
                    'protein_id': 'no_results',
                    'similarity_score': 0.0,
                    'family': 'Unknown',
                    'description': 'No matches found',
                    'database': 'N/A',
                    'source_type': 'no_data',
                    'match_type': 'no_results',
                    'sequence_length': 'N/A',
                    'molecular_weight': 'N/A',
                    'isoelectric_point': 'N/A',
                    'input_protein_index': 0,
                    'input_protein_sequence': 'N/A'
                })
            
            df = pd.DataFrame(matches_data)
            
            # Sort by rank (keeping INPUT entries at top, then by similarity score)
            df['sort_key'] = df.apply(lambda row: (0 if row['rank'] == 'INPUT' else 1, -row['similarity_score'] if isinstance(row['similarity_score'], (int, float)) else 999), axis=1)
            df = df.sort_values('sort_key').drop('sort_key', axis=1)
            
            csv_path = os.path.join(output_dir, 'top_proteins_with_metadata.csv')
            df.to_csv(csv_path, index=False)
            
            # Also create a simplified version for quick reference
            simple_df = df[['rank', 'protein_id', 'similarity_score', 'family', 'source_type', 'input_protein_index']].copy()
            simple_csv_path = os.path.join(output_dir, 'top_proteins_summary.csv')
            simple_df.to_csv(simple_csv_path, index=False)
            
            return csv_path
            
        except Exception as e:
            print(f"Error creating top matches CSV: {e}")
            # Create minimal fallback file
            try:
                fallback_data = [{
                    'rank': 1,
                    'protein_id': 'error_occurred',
                    'similarity_score': 0.0,
                    'family': 'Error',
                    'description': f'Error generating results: {str(e)}',
                    'database': 'N/A',
                    'source_type': 'error',
                    'match_type': 'error',
                    'sequence_length': 'N/A',
                    'molecular_weight': 'N/A',
                    'isoelectric_point': 'N/A',
                    'input_protein_index': 0,
                    'input_protein_sequence': 'N/A'
                }]
                df = pd.DataFrame(fallback_data)
                csv_path = os.path.join(output_dir, 'top_proteins_with_metadata.csv')
                df.to_csv(csv_path, index=False)
                return csv_path
            except:
                pass
        
        return None
    
    def _create_detailed_pipeline_json(self, results: Dict[str, Any], input_data: Dict[str, Any], output_dir: str) -> str:
        """Create detailed JSON file with all pipeline execution details."""
        try:
            detailed_results = {
                'pipeline_metadata': {
                    'execution_timestamp': time.strftime("%Y-%m-%d %H:%M:%S"),
                    'input_parameters': input_data,
                    'stages_executed': list(results.keys()),
                    'pipeline_version': '1.0.0'
                },
                'stage_details': {},
                'raw_results': results
            }
            
            # Process each stage with detailed information
            for stage_name, stage_result in results.items():
                if isinstance(stage_result, dict):
                    detailed_results['stage_details'][stage_name] = {
                        'success': stage_result.get('success', False),
                        'execution_time': stage_result.get('execution_time', 0),
                        'metadata': stage_result.get('metadata', {}),
                        'output_data': stage_result.get('output_data', {}),
                        'error_message': stage_result.get('error_message')
                    }
            
            json_path = os.path.join(output_dir, 'detailed_pipeline_results.json')
            with open(json_path, 'w') as f:
                json.dump(detailed_results, f, indent=2, default=str)
            
            return json_path
            
        except Exception as e:
            print(f"Error creating detailed pipeline JSON: {e}")
        
        return None
    
    def _create_analysis_summary_csv(self, results: Dict[str, Any], input_data: Dict[str, Any], output_dir: str) -> str:
        """Create analysis summary CSV for research use."""
        try:
            summary_data = []
            
            # Extract basic info
            input_proteins = input_data.get('input_proteins', [])
            protein_count = len(input_proteins) if input_proteins else 0
            
            summary_row = {
                'analysis_id': f"analysis_{int(time.time())}",
                'protein_count': protein_count,
                'analysis_type': 'multi_protein' if protein_count > 1 else 'single_protein',
                'timestamp': time.strftime("%Y-%m-%d %H:%M:%S"),
                'stages_completed': len([s for s, r in results.items() if isinstance(r, dict) and r.get('success', False)]),
                'total_stages': len(results)
            }
            
            # Add stage-specific metrics
            if 'embedding_generation' in results:
                emb_data = results['embedding_generation'].get('output_data', {})
                summary_row.update({
                    'embeddings_generated': len(emb_data.get('embeddings', [])),
                    'embedding_dimensions': emb_data.get('embedding_dimensions')
                })
            
            if 'family_assignment' in results:
                family_data = results['family_assignment'].get('output_data', {})
                summary_row.update({
                    'family_assigned': family_data.get('family_id'),
                    'family_confidence': family_data.get('confidence_score')
                })
            
            if 'similarity_search' in results:
                sim_data = results['similarity_search'].get('output_data', {})
                top_matches = sim_data.get('top_matches', [])
                summary_row.update({
                    'similarity_matches_found': len(top_matches),
                    'best_similarity_score': max([m.get('similarity_score', 0) for m in top_matches]) if top_matches else 0,
                    'avg_similarity_score': np.mean([m.get('similarity_score', 0) for m in top_matches]) if top_matches else 0
                })
            
            if 'multi_protein_analysis' in results:
                multi_data = results['multi_protein_analysis'].get('output_data', {})
                summary_row.update({
                    'alignment_length': multi_data.get('alignment_length'),
                    'conservation_score': multi_data.get('conservation_score'),
                    'identity_score': multi_data.get('identity_score')
                })
            
            if 'network_analysis' in results:
                net_data = results['network_analysis'].get('output_data', {})
                summary_row.update({
                    'network_nodes': net_data.get('nodes_count', 0),
                    'network_edges': net_data.get('edges_count', 0),
                    'network_density': net_data.get('network_density')
                })
            
            summary_data.append(summary_row)
            
            df = pd.DataFrame(summary_data)
            csv_path = os.path.join(output_dir, 'analysis_summary.csv')
            df.to_csv(csv_path, index=False)
            return csv_path
            
        except Exception as e:
            print(f"Error creating analysis summary CSV: {e}")
        
        return None
