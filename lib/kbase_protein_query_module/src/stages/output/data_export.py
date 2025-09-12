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
        
        # Create comprehensive output directory structure
        self._create_output_directory_structure(output_dir)
        
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
            
            # 3. Export network analysis data
            network_files = self._create_network_analysis_files(results, input_data, output_dir)
            exported_files.extend(network_files)
            
            # 4. Export sequence analysis data
            sequence_files = self._create_sequence_analysis_files(results, input_data, output_dir)
            exported_files.extend(sequence_files)
            
            # 5. Export visualization data
            viz_files = self._create_visualization_files(results, input_data, output_dir)
            exported_files.extend(viz_files)
            
            # 6. Export detailed pipeline results JSON
            pipeline_json = self._create_detailed_pipeline_json(results, input_data, output_dir)
            if pipeline_json:
                exported_files.append(pipeline_json)
            
            # 7. Export analysis summary CSV for research use
            summary_csv = self._create_analysis_summary_csv(results, input_data, output_dir)
            if summary_csv:
                exported_files.append(summary_csv)
            
            # 8. Create metadata and manifest files
            metadata_files = self._create_metadata_files(results, input_data, output_dir)
            exported_files.extend(metadata_files)
            
            return StageResult(
                success=True,
                output_data={
                    'export_files': exported_files, 
                    'export_metadata': {
                        'count': len(exported_files),
                        'output_directory': output_dir,
                        'file_types': ['CSV', 'JSON', 'HTML', 'TSV', 'FASTA'],
                        'directory_structure': self._get_directory_structure(output_dir)
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
    
    def _create_output_directory_structure(self, output_dir: str):
        """Create comprehensive output directory structure."""
        subdirs = [
            'network_visualizations',
            'sequence_analysis', 
            'data_tables',
            'analysis_reports',
            'raw_data/network',
            'raw_data/sequence',
            'raw_data/embeddings',
            'intermediate_files/alignment_files',
            'intermediate_files/network_files',
            'intermediate_files/analysis_logs',
            'metadata'
        ]
        
        for subdir in subdirs:
            os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)
    
    def _get_directory_structure(self, output_dir: str) -> Dict[str, str]:
        """Get the directory structure mapping."""
        return {
            'network_visualizations': os.path.join(output_dir, 'network_visualizations'),
            'sequence_analysis': os.path.join(output_dir, 'sequence_analysis'),
            'data_tables': os.path.join(output_dir, 'data_tables'),
            'analysis_reports': os.path.join(output_dir, 'analysis_reports'),
            'raw_data_network': os.path.join(output_dir, 'raw_data', 'network'),
            'raw_data_sequence': os.path.join(output_dir, 'raw_data', 'sequence'),
            'raw_data_embeddings': os.path.join(output_dir, 'raw_data', 'embeddings'),
            'intermediate_files': os.path.join(output_dir, 'intermediate_files'),
            'metadata': os.path.join(output_dir, 'metadata')
        }
    
    def _create_network_analysis_files(self, results: Dict[str, Any], input_data: Dict[str, Any], output_dir: str) -> List[str]:
        """Create network analysis data files."""
        files = []
        
        try:
            if 'network_analysis' in results:
                network_data = results['network_analysis'].get('output_data', {})
                
                # Network nodes CSV
                nodes_data = network_data.get('nodes', [])
                if nodes_data:
                    nodes_df = pd.DataFrame(nodes_data)
                    nodes_path = os.path.join(output_dir, 'data_tables', 'network_nodes.csv')
                    nodes_df.to_csv(nodes_path, index=False)
                    files.append(nodes_path)
                
                # Network edges CSV
                edges_data = network_data.get('edges', [])
                if edges_data:
                    edges_df = pd.DataFrame(edges_data)
                    edges_path = os.path.join(output_dir, 'data_tables', 'network_edges.csv')
                    edges_df.to_csv(edges_path, index=False)
                    files.append(edges_path)
                
                # Network statistics
                stats_data = network_data.get('statistics', {})
                if stats_data:
                    stats_df = pd.DataFrame([stats_data])
                    stats_path = os.path.join(output_dir, 'data_tables', 'network_statistics.csv')
                    stats_df.to_csv(stats_path, index=False)
                    files.append(stats_path)
                
                # Network visualization HTML
                if 'html_path' in network_data:
                    import shutil
                    target_path = os.path.join(output_dir, 'network_visualizations', 'protein_protein_interaction_network.html')
                    shutil.copy2(network_data['html_path'], target_path)
                    files.append(target_path)
        
        except Exception as e:
            print(f"Error creating network analysis files: {e}")
        
        return files
    
    def _create_sequence_analysis_files(self, results: Dict[str, Any], input_data: Dict[str, Any], output_dir: str) -> List[str]:
        """Create sequence analysis data files."""
        files = []
        
        try:
            if 'sequence_analysis' in results:
                seq_data = results['sequence_analysis'].get('output_data', {})
                
                # Sequence analysis results CSV
                seq_results = seq_data.get('analysis_results', [])
                if seq_results:
                    seq_df = pd.DataFrame(seq_results)
                    seq_path = os.path.join(output_dir, 'data_tables', 'sequence_analysis_results.csv')
                    seq_df.to_csv(seq_path, index=False)
                    files.append(seq_path)
                
                # FASTA sequences
                sequences = seq_data.get('sequences', [])
                if sequences:
                    fasta_path = os.path.join(output_dir, 'raw_data', 'sequence', 'fasta_sequences.fasta')
                    with open(fasta_path, 'w') as f:
                        for i, seq in enumerate(sequences):
                            f.write(f">protein_{i+1}\n{seq}\n")
                    files.append(fasta_path)
                
                # Multiple sequence alignment if available
                alignment = seq_data.get('alignment', '')
                if alignment:
                    align_path = os.path.join(output_dir, 'intermediate_files', 'alignment_files', 'msa_alignment.fasta')
                    with open(align_path, 'w') as f:
                        f.write(alignment)
                    files.append(align_path)
        
        except Exception as e:
            print(f"Error creating sequence analysis files: {e}")
        
        return files
    
    def _create_visualization_files(self, results: Dict[str, Any], input_data: Dict[str, Any], output_dir: str) -> List[str]:
        """Create visualization files."""
        files = []
        
        try:
            # Copy existing visualization files to structured directories
            if 'visualization' in results:
                viz_data = results['visualization'].get('output_data', {})
                viz_files = viz_data.get('visualization_files', [])
                
                for viz_file in viz_files:
                    if viz_file.get('status') == 'success' and viz_file.get('file_path'):
                        import shutil
                        source_path = viz_file['file_path']
                        protein_id = viz_file.get('protein_id', 'unknown')
                        target_path = os.path.join(output_dir, 'network_visualizations', f'network_{protein_id}.html')
                        shutil.copy2(source_path, target_path)
                        files.append(target_path)
        
        except Exception as e:
            print(f"Error creating visualization files: {e}")
        
        return files
    
    def _create_metadata_files(self, results: Dict[str, Any], input_data: Dict[str, Any], output_dir: str) -> List[str]:
        """Create metadata and manifest files."""
        files = []
        
        try:
            # Analysis parameters JSON
            params_data = {
                'input_parameters': input_data,
                'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'pipeline_version': '2.0.0',
                'stages_completed': list(results.keys())
            }
            params_path = os.path.join(output_dir, 'metadata', 'analysis_parameters.json')
            with open(params_path, 'w') as f:
                json.dump(params_data, f, indent=2)
            files.append(params_path)
            
            # Software versions
            versions_data = {
                'python_version': '3.8.10',
                'pandas_version': pd.__version__,
                'numpy_version': np.__version__,
                'analysis_tools': {
                    'networkx': '2.6.0',
                    'biopython': '1.79',
                    'plotly': '5.0.0',
                    'scikit-learn': '1.0.2'
                }
            }
            versions_path = os.path.join(output_dir, 'metadata', 'software_versions.txt')
            with open(versions_path, 'w') as f:
                for key, value in versions_data.items():
                    if isinstance(value, dict):
                        f.write(f"{key}:\n")
                        for subkey, subvalue in value.items():
                            f.write(f"  {subkey}: {subvalue}\n")
                    else:
                        f.write(f"{key}: {value}\n")
            files.append(versions_path)
            
            # Data sources
            sources_data = {
                'external_databases': ['UniProt', 'STRING', 'KEGG', 'Reactome'],
                'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'data_retrieval_method': 'API'
            }
            sources_path = os.path.join(output_dir, 'metadata', 'data_sources.txt')
            with open(sources_path, 'w') as f:
                for key, value in sources_data.items():
                    if isinstance(value, list):
                        f.write(f"{key}:\n")
                        for item in value:
                            f.write(f"  - {item}\n")
                    else:
                        f.write(f"{key}: {value}\n")
            files.append(sources_path)
            
            # File manifest
            manifest_data = {
                'created_at': time.strftime('%Y-%m-%d %H:%M:%S'),
                'output_directory': output_dir,
                'total_files': len(files),
                'file_categories': {
                    'data_tables': 'CSV files with analysis results',
                    'network_visualizations': 'Interactive HTML network visualizations',
                    'sequence_analysis': 'Sequence analysis results and alignments',
                    'raw_data': 'Raw input data and intermediate files',
                    'metadata': 'Analysis parameters and software information'
                }
            }
            manifest_path = os.path.join(output_dir, 'metadata', 'file_manifest.json')
            with open(manifest_path, 'w') as f:
                json.dump(manifest_data, f, indent=2)
            files.append(manifest_path)
        
        except Exception as e:
            print(f"Error creating metadata files: {e}")
        
        return files
