"""
Network Analysis Stage for KBase Protein Query Module

This stage performs comprehensive network analysis using the new NetworkAnalyzer.
"""

import logging
import time
import numpy as np
import pandas as pd
from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from ...analysis.network_analyzer import NetworkAnalyzer

logger = logging.getLogger(__name__)

class NetworkAnalysisStage(BaseStage):
    """Network analysis stage using the new NetworkAnalyzer for comprehensive analysis."""
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.k_neighbors = config.get('k_neighbors', 8) if config else 8
        self.similarity_threshold = config.get('similarity_threshold', 0.1) if config else 0.1
        self.mutual_knn = config.get('mutual_knn', True) if config else True
        self.min_network_size = config.get('min_network_size', 5) if config else 5
        self.max_network_size = config.get('max_network_size', 100) if config else 100
        
        # Initialize the network analyzer
        self.network_analyzer = NetworkAnalyzer(
            k_neighbors=self.k_neighbors,
            similarity_threshold=self.similarity_threshold,
            mutual_knn=self.mutual_knn,
            min_network_size=self.min_network_size,
            max_network_size=self.max_network_size
        )
    
    def get_stage_name(self) -> str:
        return "network_analysis"
    
    def get_required_inputs(self) -> List[str]:
        return ['embeddings', 'similarity_results', 'family_assignments']
    
    def get_optional_inputs(self) -> List[str]:
        return ['network_config', 'k_neighbors', 'similarity_threshold', 'mutual_knn']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'network_data': {
                'type': 'object',
                'description': 'Network analysis results for each protein',
                'properties': {
                    'network': {'type': 'object'},
                    'network_properties': {'type': 'object'},
                    'node_centrality': {'type': 'object'},
                    'community_detection': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data, workspace_client=None):
        """Execute comprehensive network analysis using the new NetworkAnalyzer."""
        start_time = time.time()
        
        try:
            embeddings = input_data['embeddings']
            similarity_results = input_data['similarity_results']
            protein_ids = input_data.get('protein_ids', [])
            metadata_df = input_data.get('metadata_df', pd.DataFrame())
            
            # Perform network analysis for each protein
            network_results = {}
            
            for i, embedding in enumerate(embeddings):
                protein_id = protein_ids[i] if i < len(protein_ids) else f"protein_{i}"
                
                try:
                    # Get similarity results for this protein
                    similarity_result = similarity_results.get(protein_id, {})
                    
                    if similarity_result.get('status') != 'success' or not similarity_result.get('similar_proteins'):
                        logger.warning(f"No similarity results for {protein_id}, skipping network analysis")
                        network_results[protein_id] = self._create_empty_result(protein_id, 'No similarity results')
                        continue
                    
                    # Get similar proteins
                    similar_proteins = similarity_result['similar_proteins']
                    
                    if not similar_proteins:
                        network_results[protein_id] = self._create_empty_result(protein_id, 'No similar proteins found')
                        continue
                    
                    # Perform comprehensive network analysis using the new NetworkAnalyzer
                    try:
                        analysis_results = self.network_analyzer.analyze_similarity_search_results(
                            similar_proteins=similar_proteins,
                            embeddings=embeddings,
                            protein_ids=protein_ids,
                            metadata_df=metadata_df,
                            query_embedding=embedding,
                            query_protein_id=protein_id
                        )
                        
                        network_results[protein_id] = {
                            'network_analysis_results': analysis_results,
                            'similarity_table': analysis_results.get('similarity_table'),
                            'top_similar_proteins': analysis_results.get('top_similar_proteins'),
                            'network_visualization': analysis_results.get('network_visualization'),
                            'network_properties': analysis_results.get('network_properties'),
                            'network_statistics': analysis_results.get('network_statistics'),
                            'clustering_analysis': analysis_results.get('clustering_analysis'),
                            'status': 'success'
                        }
                        
                    except Exception as e:
                        logger.warning(f"Network analysis failed for {protein_id}: {e}")
                        network_results[protein_id] = self._create_empty_result(protein_id, str(e))
                    
                except Exception as e:
                    logger.warning(f"Failed to perform network analysis for {protein_id}: {e}")
                    network_results[protein_id] = self._create_empty_result(protein_id, str(e))
            
            execution_time = time.time() - start_time
            
            logger.info(f"Network analysis completed for {len(embeddings)} proteins in {execution_time:.2f}s")
            
            return StageResult(
                success=True,
                output_data={'network_results': network_results},
                metadata={
                    'num_proteins': len(embeddings),
                    'successful_networks': len([r for r in network_results.values() if r['status'] == 'success']),
                    'execution_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Network analysis stage failed: {e}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={'execution_time': execution_time},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _create_empty_result(self, protein_id: str, error_message: str) -> Dict:
        """Create an empty result for failed network analysis."""
        return {
            'network_analysis_results': {
                'query_protein_id': protein_id,
                'similarity_table': pd.DataFrame(),
                'top_similar_proteins': [],
                'network_visualization': None,
                'network_properties': {},
                'network_statistics': pd.DataFrame(),
                'clustering_analysis': {}
            },
            'similarity_table': pd.DataFrame(),
            'top_similar_proteins': [],
            'network_visualization': None,
            'network_properties': {},
            'network_statistics': pd.DataFrame(),
            'clustering_analysis': {},
            'status': 'error',
            'error': error_message
        }
    

