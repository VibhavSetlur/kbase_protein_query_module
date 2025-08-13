"""
Family Assignment Stage for KBase Protein Query Module

This stage assigns proteins to families based on embeddings using FAISS binary centroid search.
"""

import logging
import time
from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from ...storage import ProteinFamilyAssigner

logger = logging.getLogger(__name__)

class FamilyAssignmentStage(BaseStage):
    """Family assignment stage using FAISS binary centroid search."""
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.family_assigner = None
        self.centroids_path = config.get('centroids_path', 'data/family_centroids/files/family_centroids_binary.npz') if config else 'data/family_centroids/files/family_centroids_binary.npz'
    
    def get_stage_name(self) -> str:
        return "family_assignment"
    
    def get_required_inputs(self) -> List[str]:
        return ['embeddings']
    
    def get_optional_inputs(self) -> List[str]:
        return ['family_config', 'centroids_path']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'family_assignments': {
                'type': 'object',
                'description': 'Protein family assignments with confidence scores',
                'properties': {
                    'family_id': {'type': 'string'},
                    'confidence': {'type': 'number'},
                    'eigenprotein_id': {'type': 'string'}
                }
            }
        }
    
    def run(self, input_data, workspace_client=None):
        """Execute family assignment using FAISS binary centroid search."""
        start_time = time.time()
        
        try:
            embeddings = input_data['embeddings']
            protein_ids = input_data.get('protein_ids', [])
            
            # Initialize family assigner if not already done
            if self.family_assigner is None:
                self.family_assigner = ProteinFamilyAssigner()
                self.family_assigner.load_family_centroids(self.centroids_path)
            
            # Perform family assignments
            family_assignments = {}
            
            for i, embedding in enumerate(embeddings):
                protein_id = protein_ids[i] if i < len(protein_ids) else f"protein_{i}"
                
                try:
                    # Assign family using FAISS binary centroid search
                    assignment = self.family_assigner.assign_family(embedding)
                    
                    family_assignments[protein_id] = {
                        'family_id': assignment['family_id'],
                        'confidence': assignment['confidence'],
                        'eigenprotein_id': assignment['eigenprotein_id'],
                        'status': 'success'
                    }
                    
                except Exception as e:
                    logger.warning(f"Failed to assign family for {protein_id}: {e}")
                    family_assignments[protein_id] = {
                        'family_id': 'unknown',
                        'confidence': 0.0,
                        'eigenprotein_id': None,
                        'status': 'error',
                        'error': str(e)
                    }
            
            execution_time = time.time() - start_time
            
            logger.info(f"Family assignment completed for {len(embeddings)} proteins in {execution_time:.2f}s")
            
            return StageResult(
                success=True,
                output_data={'family_assignments': family_assignments},
                metadata={
                    'num_proteins': len(embeddings),
                    'successful_assignments': len([a for a in family_assignments.values() if a['status'] == 'success']),
                    'execution_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Family assignment stage failed: {e}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={'execution_time': execution_time},
                execution_time=execution_time,
                error_message=str(e)
            )
