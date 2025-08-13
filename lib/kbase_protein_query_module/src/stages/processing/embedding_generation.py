"""
Embedding Generation Stage for KBase Protein Query Module

This stage generates protein embeddings using deep learning models.
"""

import logging
import time
import numpy as np
from typing import Dict, Any, List, Optional, Union

from ..base_stage import BaseStage, StageResult
from ...processing.embeddings.generator import ProteinEmbeddingGenerator

logger = logging.getLogger(__name__)

class EmbeddingGenerationStage(BaseStage):
    """
    Generates protein embeddings using deep learning models.
    
    Handles:
    - Protein sequence embedding generation
    - Model loading and caching
    - Batch processing
    - Embedding storage
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.model_name = config.get('model_name', 'esm2_t6_8M_UR50D') if config else 'esm2_t6_8M_UR50D'
        self.device = config.get('device', 'cpu') if config else 'cpu'
        self.batch_size = config.get('batch_size', 32) if config else 32
        self.cache_embeddings = config.get('cache_embeddings', True) if config else True
        
        # Initialize embedding generator
        self.embedding_generator = ProteinEmbeddingGenerator(
            model_name=self.model_name,
            device=self.device
        )
    
    def get_stage_name(self) -> str:
        return "embedding_generation"
    
    def get_required_inputs(self) -> List[str]:
        return ['extracted_data']
    
    def get_optional_inputs(self) -> List[str]:
        return ['embedding_config']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input parameters."""
        if 'extracted_data' not in input_data:
            logger.error("Missing extracted_data")
            return False
        
        extracted_data = input_data['extracted_data']
        if 'protein_records' not in extracted_data:
            logger.error("Missing protein_records in extracted_data")
            return False
        
        protein_records = extracted_data['protein_records']
        if not protein_records:
            logger.error("No protein records to process")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'embeddings': {
                'type': 'object',
                'properties': {
                    'protein_embeddings': {'type': 'array'},
                    'embedding_stats': {'type': 'object'},
                    'metadata': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Generate embeddings for protein sequences."""
        start_time = time.time()
        
        try:
            extracted_data = input_data['extracted_data']
            embedding_config = input_data.get('embedding_config', {})
            
            protein_records = extracted_data['protein_records']
            
            # Generate embeddings
            embeddings_data = self._generate_embeddings(protein_records, embedding_config)
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=True,
                output_data={
                    'embeddings': embeddings_data
                },
                metadata={
                    'model_name': self.model_name,
                    'device': self.device,
                    'generation_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Embedding generation failed: {str(e)}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _generate_embeddings(self, protein_records: List, embedding_config: Dict[str, Any]) -> Dict[str, Any]:
        """Generate embeddings for protein records."""
        protein_embeddings = []
        embedding_stats = {
            'total_proteins': len(protein_records),
            'successfully_embedded': 0,
            'failed_embeddings': 0,
            'errors': []
        }
        
        # Process in batches
        for i in range(0, len(protein_records), self.batch_size):
            batch = protein_records[i:i + self.batch_size]
            
            for record in batch:
                try:
                    # Generate embedding
                    embedding = self.embedding_generator.generate_embedding(record.sequence)
                    
                    protein_embeddings.append({
                        'protein_id': record.id,
                        'embedding': embedding,
                        'sequence_length': len(record.sequence),
                        'metadata': record.metadata or {}
                    })
                    
                    embedding_stats['successfully_embedded'] += 1
                    
                except Exception as e:
                    embedding_stats['failed_embeddings'] += 1
                    embedding_stats['errors'].append(f"Failed to embed {record.id}: {str(e)}")
        
        return {
            'protein_embeddings': protein_embeddings,
            'embedding_stats': embedding_stats,
            'metadata': {
                'model_name': self.model_name,
                'embedding_dimension': protein_embeddings[0]['embedding'].shape[0] if protein_embeddings else 0,
                'device': self.device
            }
        }
