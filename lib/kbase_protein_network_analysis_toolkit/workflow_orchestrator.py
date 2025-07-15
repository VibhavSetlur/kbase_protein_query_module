"""
Protein Network Analysis Workflow Orchestrator

This module provides a workflow orchestrator that integrates with
the hierarchical storage system for massive datasets (250M+ proteins) while
maintaining backward compatibility with existing workflows.
"""

import numpy as np
import pandas as pd
import networkx as nx
import logging
import os
import yaml
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
import time
import gc

from kbase_protein_network_analysis_toolkit.embedding_generator import ProteinEmbeddingGenerator
from kbase_protein_network_analysis_toolkit.family_classifier import ProteinFamilyClassifier
from kbase_protein_network_analysis_toolkit.similarity_index import FamilySpecificIndex, SimilarityIndex, HierarchicalIndex, StreamingIndex
from kbase_protein_network_analysis_toolkit.network_builder import DynamicNetworkBuilder
from kbase_protein_network_analysis_toolkit.storage import ProteinStorage, MemoryEfficientLoader

logger = logging.getLogger(__name__)


class ProteinNetworkWorkflow:
    """
    Workflow orchestrator for massive protein datasets.
    
    This class provides the same interface as the original workflow but uses
    storage and indexing for massive datasets (250M+ proteins).
    """
    
    def __init__(self, config_file: str = "config_optimized.yaml"):
        """
        Initialize the workflow with configuration.
        
        Args:
            config_file: Path to configuration file
        """
        self.config = self._load_config(config_file)
        self._setup_logging()
        
        # Initialize components
        self.embedding_generator = None
        self.family_classifier = None
        
        # Initialize storage components
        self.storage = None
        self.hierarchical_index = None
        self.streaming_index = None
        self.memory_loader = None
        
        # Load pre-computed data if available
        self._load_precomputed_data()
        
        # Performance monitoring
        self.performance_metrics = {}
    
    def _load_config(self, config_file: str) -> Dict:
        """Load configuration from YAML file."""
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config
    
    def _setup_logging(self):
        """Setup logging configuration."""
        log_config = self.config.get('logging', {})
        log_level = getattr(logging, log_config.get('level', 'INFO'))
        
        # Ensure log directory exists
        log_file = log_config.get('log_file', 'logs/optimized_protein_network.log')
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
        
        # Configure logging
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler() if log_config.get('console_output', True) else logging.NullHandler()
            ]
        )
    
    def _load_precomputed_data(self):
        """Load pre-computed data using storage."""
        storage_config = self.config.get('storage', {})
        
        # Initialize storage
        optimized_storage_dir = storage_config.get('optimized_storage_dir', 'data')
        if os.path.exists(optimized_storage_dir):
            logger.info(f"Loading storage from {optimized_storage_dir}")
            self.storage = ProteinStorage(base_dir=optimized_storage_dir)
            self.memory_loader = MemoryEfficientLoader(self.storage)
            
            # Get available families
            self.available_families = self.storage.get_family_list()
            logger.info(f"Found {len(self.available_families)} families in storage")
            
            # Load family statistics
            self.family_stats = self.storage.get_family_stats()
            
        else:
            logger.warning(f"Storage not found at {optimized_storage_dir}")
            logger.info("Falling back to legacy storage system")
            self._load_legacy_data()
            return
        
        # Initialize hierarchical index
        index_storage_dir = storage_config.get('index_storage_dir', 'data/indexes')
        if os.path.exists(index_storage_dir):
            logger.info(f"Loading hierarchical index from {index_storage_dir}")
            index_config = self.config.get('similarity_search', {})
            self.hierarchical_index = HierarchicalIndex(
                base_dir=index_storage_dir,
                index_type=index_config.get('index_type', 'faiss'),
                quantization=index_config.get('faiss', {}).get('quantization', 'pq'),
                cache_size=index_config.get('cache_size', 10)
            )
        else:
            logger.warning(f"Hierarchical index not found at {index_storage_dir}")
            self.hierarchical_index = None
        
        # Initialize streaming index for memory-efficient processing
        streaming_dir = storage_config.get('streaming_dir', 'data/streaming')
        self.streaming_index = StreamingIndex(
            storage_dir=streaming_dir,
            batch_size=storage_config.get('streaming_batch_size', 10000),
            max_memory_gb=storage_config.get('max_memory_gb', 8.0)
        )
        
        # Set configuration attributes for compatibility
        embedding_config = self.config.get('embedding', {})
        self.model_name = embedding_config.get('model_name', 'esm2_t48_15B_UR50D')
        self.device = embedding_config.get('device', 'auto')
        
        # Get embedding dimension from first family
        if self.family_stats:
            first_family = list(self.family_stats.keys())[0]
            self.embedding_dim = self.family_stats[first_family]['embedding_dim']
        else:
            self.embedding_dim = 5120  # Default for ESM-2 15B
        
        logger.info(f"Configured model: {self.model_name}")
        logger.info(f"Expected embedding dimension: {self.embedding_dim}")
    
    def _load_legacy_data(self):
        """Load data using legacy storage system for backward compatibility."""
        logger.info("Loading legacy data for backward compatibility")
        
        # This would load data using the original workflow system
        # For now, we'll just set up the basic structure
        self.available_families = []
        self.family_stats = {}
        self.embedding_dim = 5120
        
        logger.warning("Legacy data loading not fully implemented - please migrate to optimized storage")
    
    def generate_query_embedding(self, query_sequence: str, 
                               query_protein_id: str = "QUERY_PROTEIN") -> np.ndarray:
        """
        Generate embedding for a query protein sequence.
        
        Args:
            query_sequence: Amino acid sequence string
            query_protein_id: Query protein ID
            
        Returns:
            Query protein embedding
        """
        logger.info(f"Generating embedding for query protein: {query_protein_id}")
        
        # Initialize embedding generator if not already done
        if self.embedding_generator is None:
            embedding_config = self.config.get('embedding', {})
            self.embedding_generator = ProteinEmbeddingGenerator(
                model_name=embedding_config.get('model_name', 'esm2_t48_15B_UR50D'),
                device=embedding_config.get('device', 'auto')
            )
        
        # Generate embedding
        query_embedding = self.embedding_generator.generate_embedding(
            query_sequence, 
            pooling_method="mean"
        )
        
        logger.info(f"Generated embedding with shape: {query_embedding.shape}")
        return query_embedding
    
    def classify_query_family(self, query_embedding: np.ndarray) -> Tuple[str, float]:
        """
        Classify query protein into a family using storage.
        
        Args:
            query_embedding: Query protein embedding
            
        Returns:
            Tuple of (family_id, confidence_score)
        """
        if not self.available_families:
            raise ValueError("No families available in storage")
        
        logger.info("Classifying query protein into family...")
        
        # Normalize query embedding
        query_norm = query_embedding / (np.linalg.norm(query_embedding) + 1e-8)
        
        # Search across all families to find the best match
        best_family = None
        best_similarity = -1.0
        
        # Use hierarchical index if available
        if self.hierarchical_index:
            logger.info("Using hierarchical index for family classification")
            
            # Search all families
            family_results = self.hierarchical_index.search_all_families(
                query_norm,
                top_k=1,
                max_families=100  # Limit for performance
            )
            
            if family_results:
                # Get the best result
                best_family, similarities, protein_ids = family_results[0]
                best_similarity = similarities[0] if len(similarities) > 0 else 0.0
                
        else:
            logger.info("Using streaming search for family classification")
            
            # Use streaming search as fallback
            # This is more memory-intensive but works without hierarchical index
            for family_id in self.available_families[:50]:  # Limit for performance
                try:
                    family_embeddings, family_protein_ids = self.storage.load_family_embeddings(family_id)
                    
                    # Compute similarity to family centroid
                    family_centroid = np.mean(family_embeddings, axis=0)
                    family_centroid_norm = family_centroid / (np.linalg.norm(family_centroid) + 1e-8)
                    
                    similarity = np.dot(query_norm, family_centroid_norm)
                    
                    if similarity > best_similarity:
                        best_similarity = similarity
                        best_family = family_id
                        
                except Exception as e:
                    logger.warning(f"Error processing family {family_id}: {e}")
                    continue
        
        if best_family is None:
            raise ValueError("Could not classify query into any family")
        
        logger.info(f"Query classified into family {best_family} with confidence {best_similarity:.3f}")
        return best_family, best_similarity
    
    def load_family_subset(self, family_id: str) -> Tuple[np.ndarray, List[str], pd.DataFrame]:
        """
        Load family-specific subset using storage.
        
        Args:
            family_id: Family ID to load
            
        Returns:
            Tuple of (embeddings, protein_ids, metadata)
        """
        if not self.storage:
            raise ValueError("Storage not available")
        
        logger.info(f"Loading family {family_id} subset from storage...")
        
        # Load family embeddings and protein IDs
        family_embeddings, family_protein_ids = self.storage.load_family_embeddings(family_id)
        
        # Load family metadata
        try:
            family_metadata = self.storage.load_metadata(family_id=family_id)
        except FileNotFoundError:
            logger.warning(f"No metadata found for family {family_id}, creating empty metadata")
            family_metadata = pd.DataFrame(index=family_protein_ids)
        
        logger.info(f"Loaded {len(family_protein_ids)} proteins from family {family_id}")
        return family_embeddings, family_protein_ids, family_metadata
    
    def perform_optimized_similarity_search(self, query_embedding: np.ndarray,
                                          family_id: str,
                                          k: int = 50) -> List[Dict]:
        """
        Perform optimized similarity search within family.
        
        Args:
            query_embedding: Query protein embedding
            family_id: Family ID to search within
            k: Number of similar proteins to retrieve
            
        Returns:
            List of similar proteins with metadata
        """
        logger.info(f"Performing optimized similarity search for top {k} proteins in family {family_id}...")
        
        # Use hierarchical index if available
        if self.hierarchical_index:
            try:
                similarities, protein_ids = self.hierarchical_index.search_family(
                    family_id, query_embedding, top_k=k
                )
                
                # Convert to result format
                similar_proteins = []
                for i, (similarity, protein_id) in enumerate(zip(similarities, protein_ids)):
                    protein_result = {
                        'protein_id': protein_id,
                        'similarity_score': float(similarity),
                        'rank': i + 1
                    }
                    similar_proteins.append(protein_result)
                
                logger.info(f"Found {len(similar_proteins)} similar proteins using hierarchical index")
                return similar_proteins
                
            except Exception as e:
                logger.warning(f"Hierarchical index search failed: {e}, falling back to streaming search")
        
        # Fallback to streaming search
        logger.info("Using streaming search as fallback")
        
        # Load family embeddings
        family_embeddings, family_protein_ids = self.storage.load_family_embeddings(family_id)
        
        # Normalize query embedding
        query_norm = query_embedding / (np.linalg.norm(query_embedding) + 1e-8)
        
        # Compute similarities
        similarities = np.dot(family_embeddings, query_norm)
        
        # Get top k results
        top_indices = np.argsort(similarities)[::-1][:k]
        
        similar_proteins = []
        for i, idx in enumerate(top_indices):
            protein_result = {
                'protein_id': family_protein_ids[idx],
                'similarity_score': float(similarities[idx]),
                'rank': i + 1
            }
            similar_proteins.append(protein_result)
        
        logger.info(f"Found {len(similar_proteins)} similar proteins using streaming search")
        return similar_proteins
    
    def build_optimized_network(self, query_embedding: np.ndarray,
                               query_protein_id: str,
                               similar_proteins: List[Dict],
                               family_id: str,
                               network_method: str = "mutual_knn") -> Tuple[nx.Graph, Dict]:
        """
        Build optimized network using family-specific data.
        
        Args:
            query_embedding: Query protein embedding
            query_protein_id: Query protein ID
            similar_proteins: List of similar proteins
            family_id: Family ID
            network_method: Network construction method
            
        Returns:
            Tuple of (NetworkX graph, network properties)
        """
        logger.info(f"Building optimized network for family {family_id}...")
        
        # Load family embeddings for network construction
        family_embeddings, family_protein_ids = self.storage.load_family_embeddings(family_id)
        
        # Get embeddings for similar proteins
        similar_protein_ids = [p['protein_id'] for p in similar_proteins]
        similar_indices = [family_protein_ids.index(pid) for pid in similar_protein_ids if pid in family_protein_ids]
        similar_embeddings = family_embeddings[similar_indices]
        
        # Create network builder
        network_config = self.config.get('network', {})
        builder = DynamicNetworkBuilder(**network_config)
        
        # Build network
        G, properties = builder.build_network_from_similar_proteins(
            similar_proteins=similar_proteins,
            embeddings=similar_embeddings,
            protein_ids=similar_protein_ids,
            query_embedding=query_embedding,
            query_protein_id=query_protein_id,
            method=network_method
        )
        
        logger.info(f"Built network with {len(G.nodes())} nodes and {len(G.edges())} edges")
        return G, properties
    
    def run_optimized_workflow(self, query_sequence: str,
                              query_protein_id: str = "QUERY_PROTEIN",
                              k_similar: int = 50,
                              network_method: str = "mutual_knn",
                              save_results: bool = True) -> Dict:
        """
        Run the complete optimized workflow.
        
        Args:
            query_sequence: Query protein sequence
            query_protein_id: Query protein ID
            k_similar: Number of similar proteins to retrieve
            network_method: Network construction method
            save_results: Whether to save results
            
        Returns:
            Dictionary with workflow results
        """
        start_time = time.time()
        
        # Initialize results dictionary
        results = {
            'status': 'running',
            'query_protein_id': query_protein_id,
            'query_sequence': query_sequence,
            'workflow_steps': {},
            'timing': {},
            'performance_metrics': {}
        }
        
        try:
            # Step 1: Generate query embedding
            step_start = time.time()
            query_embedding = self.generate_query_embedding(query_sequence, query_protein_id)
            results['workflow_steps']['embedding_generation'] = {
                'status': 'success',
                'embedding_shape': query_embedding.shape
            }
            results['timing']['embedding_generation'] = time.time() - step_start
            
            # Step 2: Classify into family
            step_start = time.time()
            family_id, confidence = self.classify_query_family(query_embedding)
            results['workflow_steps']['family_classification'] = {
                'status': 'success',
                'family_id': family_id,
                'confidence': confidence
            }
            results['timing']['family_classification'] = time.time() - step_start
            
            # Step 3: Load family subset
            step_start = time.time()
            family_embeddings, family_protein_ids, family_metadata = self.load_family_subset(family_id)
            results['workflow_steps']['family_subset_loading'] = {
                'status': 'success',
                'family_size': len(family_protein_ids)
            }
            results['timing']['family_subset_loading'] = time.time() - step_start
            
            # Step 4: Optimized similarity search
            step_start = time.time()
            similar_proteins = self.perform_optimized_similarity_search(
                query_embedding, family_id, k_similar
            )
            results['workflow_steps']['similarity_search'] = {
                'status': 'success',
                'num_similar_proteins': len(similar_proteins)
            }
            results['timing']['similarity_search'] = time.time() - step_start
            
            # Step 5: Build optimized network
            step_start = time.time()
            G, network_properties = self.build_optimized_network(
                query_embedding, query_protein_id, similar_proteins, family_id, network_method
            )
            results['workflow_steps']['network_construction'] = {
                'status': 'success',
                'network_properties': network_properties
            }
            results['timing']['network_construction'] = time.time() - step_start
            
            # Store results
            results['query_embedding'] = query_embedding
            results['family_id'] = family_id
            results['similar_proteins'] = similar_proteins
            results['network'] = G
            results['network_properties'] = network_properties
            
            # Add performance metrics
            results['performance_metrics'] = {
                'memory_usage_mb': self._get_memory_usage(),
                'total_proteins_processed': len(family_protein_ids),
                'storage_efficiency': self._calculate_storage_efficiency()
            }
            
            # Save results if requested
            if save_results:
                self._save_workflow_results(results, query_protein_id)
            
            total_time = time.time() - start_time
            results['timing']['total'] = total_time
            results['status'] = 'success'
            
            logger.info(f"Optimized workflow completed successfully in {total_time:.2f} seconds")
            
        except Exception as e:
            logger.error(f"Optimized workflow failed: {e}")
            results['status'] = 'error'
            results['error'] = str(e)
            
        finally:
            # Clean up memory
            gc.collect()
        
        return results
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB."""
        try:
            import psutil
            process = psutil.Process()
            return process.memory_info().rss / 1024 / 1024
        except ImportError:
            return 0.0
    
    def _calculate_storage_efficiency(self) -> float:
        """Calculate storage efficiency ratio."""
        if not self.family_stats:
            return 1.0
        
        # Calculate compression ratio based on family statistics
        total_size = sum(stats['file_size_mb'] for stats in self.family_stats.values())
        total_proteins = sum(stats['num_proteins'] for stats in self.family_stats.values())
        
        # Estimate original size (assuming 4 bytes per dimension per protein)
        estimated_original_size = total_proteins * self.embedding_dim * 4 / 1024 / 1024
        
        if estimated_original_size > 0:
            return total_size / estimated_original_size
        else:
            return 1.0
    
    def _save_workflow_results(self, results: Dict, query_protein_id: str):
        """Save workflow results to disk."""
        output_dir = self.config.get('storage', {}).get('output_dir', 'results')
        os.makedirs(output_dir, exist_ok=True)
        
        # Save results as JSON
        import json
        results_file = os.path.join(output_dir, f"{query_protein_id}_optimized_results.json")
        
        # Convert numpy arrays to lists for JSON serialization
        serializable_results = self._make_json_serializable(results)
        
        with open(results_file, 'w') as f:
            json.dump(serializable_results, f, indent=2)
        
        logger.info(f"Saved optimized workflow results to {results_file}")
    
    def _make_json_serializable(self, obj):
        """Convert object to JSON serializable format."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: self._make_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_json_serializable(item) for item in obj]
        else:
            return obj
    
    def get_system_info(self) -> Dict:
        """Get system information and capabilities."""
        return {
            'storage_available': self.storage is not None,
            'hierarchical_index_available': self.hierarchical_index is not None,
            'available_families': len(self.available_families),
            'total_proteins': sum(stats['num_proteins'] for stats in self.family_stats.values()),
            'storage_efficiency': self._calculate_storage_efficiency(),
            'memory_usage_mb': self._get_memory_usage(),
            'embedding_dimension': self.embedding_dim,
            'model_name': self.model_name
        }


# Backward compatibility wrapper
class LegacyProteinNetworkWorkflow(ProteinNetworkWorkflow):
    """
    Backward compatibility wrapper for the original workflow interface.
    
    This class provides the same interface as the original workflow but uses
    storage and indexing internally.
    """
    
    def __init__(self, config_file: str = "config.yaml"):
        """Initialize with backward compatibility."""
        # Try optimized config first, fall back to original
        optimized_config = "config_optimized.yaml"
        if os.path.exists(optimized_config):
            super().__init__(optimized_config)
        else:
            super().__init__(config_file)
    
    def run_complete_workflow(self, query_sequence: str,
                            query_protein_id: str = "QUERY_PROTEIN",
                            k_similar: int = 50,
                            network_method: str = "mutual_knn",
                            save_results: bool = True) -> Dict:
        """Backward compatibility method."""
        return self.run_optimized_workflow(
            query_sequence, query_protein_id, k_similar, network_method, save_results
        ) 
