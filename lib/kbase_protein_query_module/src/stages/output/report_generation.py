"""
Report Generation Stage for KBase Protein Query Module

This stage organizes raw data files and visualizations for user output.
No fancy HTML reports - just file organization and manifests.
"""

from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from .file_organizer import FileOrganizer, create_simple_file_list
import time
import os
import logging

logger = logging.getLogger(__name__)

class ReportGenerationStage(BaseStage):
    """Simple report generation stage that organizes raw data files for output."""
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config or {})
        self.output_dir = self.config.get('output_dir', 'test/outputs')
    
    def get_stage_name(self) -> str:
        return "report_generation"
    
    def get_required_inputs(self) -> List[str]:
        return ['pipeline_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['pipeline_metadata', 'protein_id']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        return 'pipeline_results' in input_data and isinstance(input_data['pipeline_results'], dict)
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'organized_files': {
                'type': 'object',
                'description': 'Dictionary of organized output files by type'
            },
            'files_manifest': {
                'type': 'string',
                'description': 'Path to files manifest JSON'
            },
            'file_list': {
                'type': 'string', 
                'description': 'Path to simple text file listing all outputs'
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        start_time = time.time()
        
        try:
            pipeline_results = input_data['pipeline_results']
            pipeline_metadata = input_data.get('pipeline_metadata', {})
            protein_id = input_data.get('protein_id', 'analysis')
            
            # Initialize file organizer
            organizer = FileOrganizer(self.output_dir)
            
            organized_files = {
                'visualizations': {},
                'data_files': {},
                'metadata_files': {}
            }
            
            # Process network analysis results
            if 'network_analysis' in pipeline_results:
                network_data = pipeline_results['network_analysis']
                
                # Register network visualization HTML
                if 'html_path' in network_data:
                    html_path = network_data['html_path']
                    new_path = organizer.register_visualization('network_visualization', html_path, 'html')
                    organized_files['visualizations']['network'] = new_path
                
                # Register network CSV files
                if 'csv_files' in network_data:
                    for csv_type, csv_path in network_data['csv_files'].items():
                        new_path = organizer.register_data_file(f'network_{csv_type}', csv_path, 
                                                   f'Network analysis {csv_type} data', 'network')
                        organized_files['data_files'][f'network_{csv_type}'] = new_path
            
            # Process multi-protein network analysis if available
            if 'multi_protein_network_analysis' in pipeline_results:
                multi_data = pipeline_results['multi_protein_network_analysis']
                
                # Register individual network visualizations
                if 'individual_networks' in multi_data:
                    for protein_id, network_result in multi_data['individual_networks'].items():
                        if 'html_path' in network_result:
                            organizer.register_visualization(f'network_{protein_id}', 
                                                           network_result['html_path'], 'html')
                            organized_files['visualizations'][f'network_{protein_id}'] = network_result['html_path']
                
                # Register multi-protein CSV files
                if 'csv_data' in multi_data:
                    for csv_type, csv_path in multi_data['csv_data'].items():
                        organizer.register_data_file(f'multi_protein_{csv_type}', csv_path,
                                                   f'Multi-protein {csv_type} data')
                        organized_files['data_files'][f'multi_protein_{csv_type}'] = csv_path
            
            # Process sequence analysis results
            if 'sequence_analysis' in pipeline_results:
                seq_data = pipeline_results['sequence_analysis']
                
                if 'csv_files' in seq_data:
                    for csv_type, csv_path in seq_data['csv_files'].items():
                        new_path = organizer.register_data_file(f'sequence_{csv_type}', csv_path,
                                                   f'Sequence analysis {csv_type} data', 'sequence')
                        organized_files['data_files'][f'sequence_{csv_type}'] = new_path
            
            # Process embedding results
            if 'embedding_generation' in pipeline_results:
                emb_data = pipeline_results['embedding_generation']
                
                if 'h5_path' in emb_data:
                    new_path = organizer.register_data_file('protein_embeddings', emb_data['h5_path'],
                                               'Protein embeddings and sequences (HDF5)', 'embeddings')
                    organized_files['data_files']['embeddings'] = new_path
            
            # Create comprehensive pipeline summary (includes all metadata)
            pipeline_summary_path = organizer.create_comprehensive_pipeline_summary(
                pipeline_results=pipeline_results,
                pipeline_metadata={
                    'protein_id': protein_id,
                    'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'pipeline_version': '2.0.0',
                    **pipeline_metadata
                }
            )
            organized_files['metadata_files']['pipeline_summary'] = pipeline_summary_path
            
            # Save manifest and file list
            manifest_path = organizer.save_manifest()
            file_list_path = create_simple_file_list(self.output_dir, {
                **organized_files['visualizations'],
                **organized_files['data_files'], 
                **organized_files['metadata_files']
            })
            
            execution_time = time.time() - start_time
            
            output_data = {
                'organized_files': organized_files,
                'files_manifest': manifest_path,
                'file_list': file_list_path,
                'file_summary': organizer.get_file_summary()
            }
            
            logger.info(f"Organized {sum(len(files) for files in organized_files.values())} files in {execution_time:.2f}s")
            
            return StageResult(
                success=True,
                output_data=output_data,
                metadata={
                    'total_files': sum(len(files) for files in organized_files.values()),
                    'output_directory': self.output_dir,
                    'manifest_path': manifest_path
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            logger.error(f"Report generation failed: {str(e)}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=time.time() - start_time,
                error_message=str(e)
            )