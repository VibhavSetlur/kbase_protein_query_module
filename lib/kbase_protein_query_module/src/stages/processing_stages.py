"""
Processing Stages for KBase Protein Query Module

This module contains all processing stages that transform and analyze protein data.
"""

import logging
import time
import numpy as np
import pandas as pd
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

from .base_stage import BaseStage, StageResult
from ..processing.embeddings.generator import ProteinEmbeddingGenerator
from ..storage import ProteinFamilyAssigner
from ..processing.similarity.hierarchical_index import HierarchicalIndex

logger = logging.getLogger(__name__)

class EmbeddingGenerationStage(BaseStage):
    """
    Generates protein embeddings using ESM-2 models.
    
    Handles:
    - ESM-2 model loading and initialization
    - Batch embedding generation
    - Embedding validation and quality control
    - Memory-efficient processing
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.model_name = config.get('model_name', 'esm2_t6_8M_UR50D') if config else 'esm2_t6_8M_UR50D'
        self.device = config.get('device', 'auto') if config else 'auto'
        self.batch_size = config.get('batch_size', 8) if config else 8
        self.embedding_generator = None
    
    def get_stage_name(self) -> str:
        return "embedding_generation"
    
    def get_required_inputs(self) -> List[str]:
        return ['protein_records']
    
    def get_optional_inputs(self) -> List[str]:
        return ['embedding_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['data_extraction']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'protein_records' not in input_data:
            logger.error("Missing protein_records")
            return False
        
        protein_records = input_data['protein_records']
        if not isinstance(protein_records, list) or len(protein_records) == 0:
            logger.error("protein_records must be a non-empty list")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'embeddings': {
                'type': 'object',
                'description': 'Dictionary mapping protein IDs to embeddings'
            },
            'embedding_stats': {
                'type': 'object',
                'properties': {
                    'total_proteins': {'type': 'integer'},
                    'successful_embeddings': {'type': 'integer'},
                    'failed_embeddings': {'type': 'integer'},
                    'embedding_dimension': {'type': 'integer'},
                    'model_info': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Generate embeddings for protein sequences."""
        start_time = time.time()
        
        try:
            protein_records = input_data['protein_records']
            embedding_config = input_data.get('embedding_config', {})
            
            # Initialize embedding generator if not already done
            if self.embedding_generator is None:
                self.embedding_generator = ProteinEmbeddingGenerator(
                    model_name=self.model_name,
                    device=self.device
                )
            
            # Extract sequences and IDs
            sequences = [record.sequence for record in protein_records]
            protein_ids = [record.protein_id for record in protein_records]
            
            # Generate embeddings
            embeddings_dict = self.embedding_generator.generate_embeddings_batch(
                sequences=sequences,
                protein_ids=protein_ids,
                batch_size=self.batch_size
            )
            
            # Calculate statistics
            total_proteins = len(protein_records)
            successful_embeddings = len(embeddings_dict)
            failed_embeddings = total_proteins - successful_embeddings
            
            # Get model information
            model_info = {
                'model_name': self.model_name,
                'device': str(self.embedding_generator.device),
                'embedding_dimension': self.embedding_generator.embedding_dim,
                'batch_size': self.batch_size
            }
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_embeddings > 0,
                output_data={
                    'embeddings': embeddings_dict,
                    'embedding_stats': {
                        'total_proteins': total_proteins,
                        'successful_embeddings': successful_embeddings,
                        'failed_embeddings': failed_embeddings,
                        'embedding_dimension': self.embedding_generator.embedding_dim,
                        'model_info': model_info
                    }
                },
                metadata={
                    'model_name': self.model_name,
                    'embedding_dimension': self.embedding_generator.embedding_dim,
                    'success_rate': successful_embeddings / total_proteins if total_proteins > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Embedding generation failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class FamilyAssignmentStage(BaseStage):
    """
    Assigns proteins to families using centroid similarity.
    
    Handles:
    - Family centroid loading
    - Binary FAISS indexing
    - Confidence scoring
    - Family assignment validation
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.centroids_path = config.get('centroids_path', 'data/family_centroids/family_centroids_binary.npz') if config else 'data/family_centroids/family_centroids_binary.npz'
        self.confidence_threshold = config.get('confidence_threshold', 0.1) if config else 0.1
        self.family_assigner = None
    
    def get_stage_name(self) -> str:
        return "family_assignment"
    
    def get_required_inputs(self) -> List[str]:
        return ['embeddings']
    
    def get_optional_inputs(self) -> List[str]:
        return ['assignment_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['embedding_generation']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'embeddings' not in input_data:
            logger.error("Missing embeddings")
            return False
        
        embeddings = input_data['embeddings']
        if not isinstance(embeddings, dict) or len(embeddings) == 0:
            logger.error("embeddings must be a non-empty dictionary")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'family_assignments': {
                'type': 'object',
                'description': 'Dictionary mapping protein IDs to family assignment results'
            },
            'assignment_stats': {
                'type': 'object',
                'properties': {
                    'total_proteins': {'type': 'integer'},
                    'assigned_proteins': {'type': 'integer'},
                    'unassigned_proteins': {'type': 'integer'},
                    'family_distribution': {'type': 'object'},
                    'confidence_stats': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Assign proteins to families."""
        start_time = time.time()
        
        try:
            embeddings = input_data['embeddings']
            assignment_config = input_data.get('assignment_config', {})
            
            # Initialize family assigner if not already done
            if self.family_assigner is None:
                self.family_assigner = ProteinFamilyAssigner()
                
                # Try multiple possible paths for centroids
                possible_paths = [
                    self.centroids_path,
                    "data/family_centroids_binary.npz",
                    "/kb/module/data/family_centroids/family_centroids_binary.npz",
                    "/kb/module/data/family_centroids_binary.npz"
                ]
                
                centroids_loaded = False
                for path in possible_paths:
                    if os.path.exists(path):
                        try:
                            self.family_assigner.load_family_centroids(path)
                            centroids_loaded = True
                            logger.info(f"Loaded family centroids from {path}")
                            break
                        except Exception as e:
                            logger.warning(f"Failed to load centroids from {path}: {e}")
                            continue
                
                if not centroids_loaded:
                    raise RuntimeError("Could not load family centroids from any expected location")
            
            # Assign families
            family_assignments = {}
            protein_ids = list(embeddings.keys())
            
            for protein_id in protein_ids:
                try:
                    embedding = embeddings[protein_id]
                    assignment_result = self.family_assigner.assign_family_with_monitoring(embedding)
                    family_assignments[protein_id] = assignment_result
                except Exception as e:
                    logger.warning(f"Failed to assign family for {protein_id}: {e}")
                    family_assignments[protein_id] = {
                        'family_id': None,
                        'confidence': 0.0,
                        'eigenprotein_id': None,
                        'status': 'error',
                        'error': str(e)
                    }
            
            # Calculate statistics
            total_proteins = len(protein_ids)
            assigned_proteins = sum(1 for result in family_assignments.values() 
                                  if result.get('status') == 'success' and result.get('family_id'))
            unassigned_proteins = total_proteins - assigned_proteins
            
            # Family distribution
            family_distribution = {}
            confidence_scores = []
            
            for result in family_assignments.values():
                if result.get('status') == 'success' and result.get('family_id'):
                    family_id = result['family_id']
                    family_distribution[family_id] = family_distribution.get(family_id, 0) + 1
                    confidence_scores.append(result.get('confidence', 0.0))
            
            # Confidence statistics
            confidence_stats = {}
            if confidence_scores:
                confidence_stats = {
                    'mean': np.mean(confidence_scores),
                    'median': np.median(confidence_scores),
                    'min': np.min(confidence_scores),
                    'max': np.max(confidence_scores),
                    'std': np.std(confidence_scores)
                }
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=assigned_proteins > 0,
                output_data={
                    'family_assignments': family_assignments,
                    'assignment_stats': {
                        'total_proteins': total_proteins,
                        'assigned_proteins': assigned_proteins,
                        'unassigned_proteins': unassigned_proteins,
                        'family_distribution': family_distribution,
                        'confidence_stats': confidence_stats
                    }
                },
                metadata={
                    'centroids_path': self.centroids_path,
                    'confidence_threshold': self.confidence_threshold,
                    'assignment_rate': assigned_proteins / total_proteins if total_proteins > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Family assignment failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class SimilaritySearchStage(BaseStage):
    """
    Performs similarity search within protein families.
    
    Handles:
    - FAISS index loading and management
    - Similarity search within families
    - Top-k retrieval
    - Similarity score calculation
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.top_k = config.get('top_k', 50) if config else 50
        self.similarity_threshold = config.get('similarity_threshold', 0.1) if config else 0.1
        self.max_families = config.get('max_families', 10) if config else 10
        self.hierarchical_index = None
        self.storage = None
    
    def get_stage_name(self) -> str:
        return "similarity_search"
    
    def get_required_inputs(self) -> List[str]:
        return ['embeddings', 'family_assignments']
    
    def get_optional_inputs(self) -> List[str]:
        return ['search_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['family_assignment']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        required = ['embeddings', 'family_assignments']
        for field in required:
            if field not in input_data:
                logger.error(f"Missing required input: {field}")
                return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'similarity_results': {
                'type': 'object',
                'description': 'Dictionary mapping protein IDs to similarity search results'
            },
            'search_stats': {
                'type': 'object',
                'properties': {
                    'total_searches': {'type': 'integer'},
                    'successful_searches': {'type': 'integer'},
                    'failed_searches': {'type': 'integer'},
                    'average_similarity': {'type': 'number'},
                    'family_coverage': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Perform similarity search for proteins."""
        start_time = time.time()
        
        try:
            embeddings = input_data['embeddings']
            family_assignments = input_data['family_assignments']
            search_config = input_data.get('search_config', {})
            
            # Initialize components if not already done
            if self.hierarchical_index is None:
                self.hierarchical_index = HierarchicalIndex()
                self.storage = ProteinStorage()
            
            # Perform similarity search
            similarity_results = {}
            protein_ids = list(embeddings.keys())
            
            for protein_id in protein_ids:
                try:
                    embedding = embeddings[protein_id]
                    family_assignment = family_assignments.get(protein_id, {})
                    
                    if family_assignment.get('status') == 'success' and family_assignment.get('family_id'):
                        family_id = family_assignment['family_id']
                        
                        # Search within family
                        try:
                            distances, similar_protein_ids = self.hierarchical_index.search_family_float(
                                family_id=family_id,
                                query_embedding=embedding,
                                top_k=self.top_k
                            )
                            
                            # Convert distances to similarities and filter by threshold
                            similarities = 1.0 - distances
                            similar_proteins = []
                            
                            for i, (similarity, similar_id) in enumerate(zip(similarities, similar_protein_ids)):
                                if similarity >= self.similarity_threshold and similar_id != protein_id:
                                    similar_proteins.append({
                                        'protein_id': similar_id,
                                        'similarity': float(similarity),
                                        'family_id': family_id,
                                        'rank': i + 1
                                    })
                            
                            similarity_results[protein_id] = {
                                'status': 'success',
                                'family_id': family_id,
                                'similar_proteins': similar_proteins,
                                'total_found': len(similar_proteins),
                                'search_stats': {
                                    'max_similarity': max([p['similarity'] for p in similar_proteins]) if similar_proteins else 0.0,
                                    'min_similarity': min([p['similarity'] for p in similar_proteins]) if similar_proteins else 0.0,
                                    'avg_similarity': np.mean([p['similarity'] for p in similar_proteins]) if similar_proteins else 0.0
                                }
                            }
                            
                        except Exception as e:
                            logger.warning(f"Failed to search family {family_id} for {protein_id}: {e}")
                            similarity_results[protein_id] = {
                                'status': 'error',
                                'family_id': family_id,
                                'error': str(e),
                                'similar_proteins': []
                            }
                    else:
                        similarity_results[protein_id] = {
                            'status': 'no_family',
                            'family_id': None,
                            'similar_proteins': [],
                            'error': 'No family assigned'
                        }
                        
                except Exception as e:
                    logger.warning(f"Failed to perform similarity search for {protein_id}: {e}")
                    similarity_results[protein_id] = {
                        'status': 'error',
                        'family_id': None,
                        'error': str(e),
                        'similar_proteins': []
                    }
            
            # Calculate statistics
            total_searches = len(protein_ids)
            successful_searches = sum(1 for result in similarity_results.values() 
                                    if result.get('status') == 'success')
            failed_searches = total_searches - successful_searches
            
            # Average similarity across all successful searches
            all_similarities = []
            family_coverage = {}
            
            for result in similarity_results.values():
                if result.get('status') == 'success':
                    family_id = result['family_id']
                    family_coverage[family_id] = family_coverage.get(family_id, 0) + 1
                    
                    for protein in result['similar_proteins']:
                        all_similarities.append(protein['similarity'])
            
            average_similarity = np.mean(all_similarities) if all_similarities else 0.0
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_searches > 0,
                output_data={
                    'similarity_results': similarity_results,
                    'search_stats': {
                        'total_searches': total_searches,
                        'successful_searches': successful_searches,
                        'failed_searches': failed_searches,
                        'average_similarity': average_similarity,
                        'family_coverage': family_coverage
                    }
                },
                metadata={
                    'top_k': self.top_k,
                    'similarity_threshold': self.similarity_threshold,
                    'search_success_rate': successful_searches / total_searches if total_searches > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Similarity search failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
