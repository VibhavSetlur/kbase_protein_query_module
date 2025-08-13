"""
Visualization Stage for KBase Protein Query Module

This stage generates interactive visualizations including network visualizations.
"""

import logging
import time
import os
import numpy as np
import pandas as pd
from typing import Dict, Any, List, Optional
from ..base_stage import BaseStage, StageResult
from ...processing.networks.builder import visualize_interactive_protein_network, DynamicNetworkBuilder

logger = logging.getLogger(__name__)

class VisualizationStage(BaseStage):
    """Visualization stage for generating interactive network visualizations."""
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.output_dir = config.get('output_dir', 'html_reports') if config else 'html_reports'
        self.k_neighbors = config.get('k_neighbors', 8) if config else 8
        self.similarity_threshold = config.get('similarity_threshold', 0.5) if config else 0.5
        self.include_metadata = config.get('include_metadata', True) if config else True
        
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
    
    def get_stage_name(self) -> str:
        return "visualization"
    
    def get_required_inputs(self) -> List[str]:
        return ['embeddings', 'network_results', 'protein_ids']
    
    def get_optional_inputs(self) -> List[str]:
        return ['visualization_config', 'metadata_df', 'output_dir']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'visualization_files': {
                'type': 'object',
                'description': 'Generated visualization files',
                'properties': {
                    'network_visualizations': {'type': 'array'},
                    'interactive_charts': {'type': 'array'},
                    'summary_plots': {'type': 'array'}
                }
            }
        }
    
    def run(self, input_data, workspace_client=None):
        """Generate interactive visualizations."""
        start_time = time.time()
        
        try:
            embeddings = input_data['embeddings']
            network_results = input_data['network_results']
            protein_ids = input_data['protein_ids']
            metadata_df = input_data.get('metadata_df', pd.DataFrame())
            
            # Generate visualizations for each protein
            visualization_files = []
            
            for i, embedding in enumerate(embeddings):
                protein_id = protein_ids[i] if i < len(protein_ids) else f"protein_{i}"
                
                try:
                    # Get network results for this protein
                    network_result = network_results.get(protein_id, {})
                    
                    if network_result.get('status') != 'success' or network_result.get('network') is None:
                        logger.warning(f"No network results for {protein_id}, skipping visualization")
                        continue
                    
                    # Get network and similar proteins
                    network = network_result['network']
                    similar_proteins = network_result.get('similar_proteins', [])
                    
                    if not similar_proteins:
                        logger.warning(f"No similar proteins for {protein_id}, skipping visualization")
                        continue
                    
                    # Create embeddings array for visualization
                    # We need embeddings for all proteins in the network
                    network_proteins = list(network.nodes())
                    
                    # For now, we'll create a simplified visualization
                    # In a full implementation, you'd load the actual embeddings for all network proteins
                    
                    # Generate interactive network visualization
                    try:
                        # Create a simple metadata DataFrame for network proteins
                        if metadata_df.empty:
                            # Create minimal metadata
                            viz_metadata = pd.DataFrame(index=network_proteins)
                            viz_metadata['Protein names'] = [f'Protein {pid}' for pid in network_proteins]
                            viz_metadata['Organism'] = ['Unknown'] * len(network_proteins)
                            viz_metadata['Function [CC]'] = ['No metadata available'] * len(network_proteins)
                        else:
                            # Use provided metadata, filtering for network proteins
                            viz_metadata = metadata_df.loc[metadata_df.index.intersection(network_proteins)]
                        
                        # Create embeddings array (simplified - in reality you'd load actual embeddings)
                        # For visualization purposes, we'll create dummy embeddings
                        network_embeddings = np.random.rand(len(network_proteins), 320).astype(np.float32)
                        
                        # Generate interactive visualization
                        output_file = os.path.join(self.output_dir, f"network_{protein_id}_{int(time.time())}.html")
                        
                        fig, G = visualize_interactive_protein_network(
                            embeddings=network_embeddings,
                            protein_ids=network_proteins,
                            metadata_df=viz_metadata,
                            k_neighbors=self.k_neighbors,
                            similarity_threshold=self.similarity_threshold,
                            query_embedding=embedding,
                            query_protein_id=protein_id,
                            output_file=output_file
                        )
                        
                        if output_file and os.path.exists(output_file):
                            visualization_files.append({
                                'protein_id': protein_id,
                                'file_path': output_file,
                                'file_type': 'interactive_network',
                                'status': 'success'
                            })
                            
                            logger.info(f"Generated network visualization for {protein_id}: {output_file}")
                        
                    except Exception as e:
                        logger.warning(f"Failed to generate visualization for {protein_id}: {e}")
                        visualization_files.append({
                            'protein_id': protein_id,
                            'file_path': None,
                            'file_type': 'interactive_network',
                            'status': 'error',
                            'error': str(e)
                        })
                    
                except Exception as e:
                    logger.warning(f"Failed to process visualization for {protein_id}: {e}")
                    visualization_files.append({
                        'protein_id': protein_id,
                        'file_path': None,
                        'file_type': 'interactive_network',
                        'status': 'error',
                        'error': str(e)
                    })
            
            execution_time = time.time() - start_time
            
            logger.info(f"Visualization generation completed for {len(embeddings)} proteins in {execution_time:.2f}s")
            
            return StageResult(
                success=True,
                output_data={'visualization_files': visualization_files},
                metadata={
                    'num_proteins': len(embeddings),
                    'successful_visualizations': len([v for v in visualization_files if v['status'] == 'success']),
                    'output_directory': self.output_dir,
                    'execution_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Visualization stage failed: {e}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={'execution_time': execution_time},
                execution_time=execution_time,
                error_message=str(e)
            )
