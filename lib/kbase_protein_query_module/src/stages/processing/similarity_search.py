"""
Similarity Search Stage for KBase Protein Query Module

This stage performs similarity search using FAISS float indices for within-family search.
"""

import logging
import time
import numpy as np
from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from ...processing.similarity.hierarchical_index import HierarchicalIndex

logger = logging.getLogger(__name__)

class SimilaritySearchStage(BaseStage):
    """Similarity search stage using FAISS float indices."""
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.hierarchical_index = None
        self.index_dir = config.get('index_dir', 'data/indexes') if config else 'data/indexes'
        self.top_k = config.get('top_k', 50) if config else 50
        self.similarity_threshold = config.get('similarity_threshold', 0.1) if config else 0.1
    
    def get_stage_name(self) -> str:
        return "similarity_search"
    
    def get_required_inputs(self) -> List[str]:
        return ['embeddings', 'family_assignments']
    
    def get_optional_inputs(self) -> List[str]:
        return ['similarity_config', 'index_dir', 'top_k', 'similarity_threshold']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'similarity_results': {
                'type': 'object',
                'description': 'Similarity search results for each protein',
                'properties': {
                    'similar_proteins': {'type': 'array'},
                    'similarity_scores': {'type': 'array'},
                    'total_found': {'type': 'number'}
                }
            }
        }
    
    def run(self, input_data, workspace_client=None):
        """Execute similarity search using FAISS float indices."""
        start_time = time.time()
        
        try:
            embeddings = input_data['embeddings']
            family_assignments = input_data['family_assignments']
            protein_ids = input_data.get('protein_ids', [])
            
            # Initialize hierarchical index if not already done
            if self.hierarchical_index is None:
                self.hierarchical_index = HierarchicalIndex(base_dir=self.index_dir)
            
            # Perform similarity search for each protein
            similarity_results = {}
            
            for i, embedding in enumerate(embeddings):
                protein_id = protein_ids[i] if i < len(protein_ids) else f"protein_{i}"
                
                try:
                    # Get family assignment for this protein
                    family_assignment = family_assignments.get(protein_id, {})
                    family_id = family_assignment.get('family_id', 'unknown')
                    
                    if family_id == 'unknown':
                        logger.warning(f"No family assignment for {protein_id}, skipping similarity search")
                        similarity_results[protein_id] = {
                            'similar_proteins': [],
                            'similarity_scores': [],
                            'total_found': 0,
                            'status': 'error',
                            'error': 'No family assignment'
                        }
                        continue
                    
                    # Ensure embedding is float32
                    if embedding.dtype != np.float32:
                        embedding = embedding.astype(np.float32)
                    
                    # Search within family using FAISS float index
                    try:
                        distances, similar_protein_ids = self.hierarchical_index.search_family_float(
                            family_id, embedding, top_k=self.top_k
                        )
                        
                        # Convert distances to similarities and filter by threshold
                        similar_proteins = []
                        similarity_scores = []
                        
                        for j, (dist, similar_id) in enumerate(zip(distances, similar_protein_ids)):
                            # Convert L2 distance to similarity (1 - normalized_distance)
                            similarity = max(0.0, 1.0 - abs(dist))
                            
                            if similarity >= self.similarity_threshold and similar_id != protein_id:
                                similar_proteins.append({
                                    'protein_id': similar_id,
                                    'similarity': similarity,
                                    'rank': j + 1
                                })
                                similarity_scores.append(similarity)
                        
                        similarity_results[protein_id] = {
                            'similar_proteins': similar_proteins,
                            'similarity_scores': similarity_scores,
                            'total_found': len(similar_proteins),
                            'family_id': family_id,
                            'status': 'success'
                        }
                        
                    except Exception as e:
                        logger.warning(f"FAISS search failed for {protein_id} in family {family_id}: {e}")
                        # Fallback to empty results
                        similarity_results[protein_id] = {
                            'similar_proteins': [],
                            'similarity_scores': [],
                            'total_found': 0,
                            'family_id': family_id,
                            'status': 'error',
                            'error': str(e)
                        }
                    
                except Exception as e:
                    logger.warning(f"Failed to perform similarity search for {protein_id}: {e}")
                    similarity_results[protein_id] = {
                        'similar_proteins': [],
                        'similarity_scores': [],
                        'total_found': 0,
                        'status': 'error',
                        'error': str(e)
                    }
            
            execution_time = time.time() - start_time
            
            logger.info(f"Similarity search completed for {len(embeddings)} proteins in {execution_time:.2f}s")
            
            return StageResult(
                success=True,
                output_data={'similarity_results': similarity_results},
                metadata={
                    'num_proteins': len(embeddings),
                    'successful_searches': len([r for r in similarity_results.values() if r['status'] == 'success']),
                    'total_similar_proteins': sum(r['total_found'] for r in similarity_results.values()),
                    'execution_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Similarity search stage failed: {e}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={'execution_time': execution_time},
                execution_time=execution_time,
                error_message=str(e)
            )
