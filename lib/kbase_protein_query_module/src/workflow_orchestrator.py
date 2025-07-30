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

from .embedding_generator import ProteinEmbeddingGenerator
from .assign_protein_family import AssignProteinFamily
from .similarity_index import HierarchicalIndex, StreamingIndex
from .network_builder import DynamicNetworkBuilder
from .storage import ProteinStorage, MemoryEfficientLoader

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
        self.assign_protein_family = AssignProteinFamily()
        
        # Initialize storage components
        self.storage = None
        self.hierarchical_index = None
        self.streaming_index = None
        self.memory_loader = None
        
        # Load pre-computed data if available
        self._load_precomputed_data()
        
        # Load centroids for assignment
        centroids_path = self.config.get('storage', {}).get('centroids_path', None)
        if centroids_path is None:
            # Try multiple possible paths for centroids
            possible_paths = [
                str(self.storage.base_dir / "family_centroids" / "family_centroids_binary.npz") if self.storage else None,
                str(self.storage.base_dir / "family_centroids_binary.npz") if self.storage else None,
                "data/family_centroids/family_centroids_binary.npz",
                "data/family_centroids_binary.npz",
                "/kb/module/data/family_centroids/family_centroids_binary.npz",
                "/kb/module/data/family_centroids_binary.npz"
            ]
            centroids_path = None
            for path in possible_paths:
                if path and os.path.exists(path):
                    centroids_path = path
                    break
        
        if centroids_path and os.path.exists(centroids_path):
            logger.info(f"Loading family centroids from: {centroids_path}")
            self.assign_protein_family.load_family_centroids(centroids_path)
        else:
            logger.warning("Family centroids file not found. Family assignment will not be available.")
            # Try to find centroids by searching
            import glob
            search_paths = ["data", "/kb/module/data", "."]
            for search_path in search_paths:
                if os.path.exists(search_path):
                    centroids_files = glob.glob(os.path.join(search_path, "**/*family_centroids*.npz"), recursive=True)
                    if centroids_files:
                        centroids_path = centroids_files[0]
                        logger.info(f"Found centroids file: {centroids_path}")
                        self.assign_protein_family.load_family_centroids(centroids_path)
                        break
        
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
        # Also check the families directory where indexes are actually stored
        families_index_dir = os.path.join(optimized_storage_dir, 'families')
        
        if os.path.exists(index_storage_dir):
            logger.info(f"Loading hierarchical index from {index_storage_dir}")
            index_config = self.config.get('similarity_search', {})
            self.hierarchical_index = HierarchicalIndex(
                base_dir=index_storage_dir,
                index_type=index_config.get('index_type', 'faiss'),
                quantization=index_config.get('faiss', {}).get('quantization', 'pq'),
                cache_size=index_config.get('cache_size', 10)
            )
        elif os.path.exists(families_index_dir):
            logger.info(f"Loading hierarchical index from {families_index_dir}")
            index_config = self.config.get('similarity_search', {})
            self.hierarchical_index = HierarchicalIndex(
                base_dir=families_index_dir,
                index_type=index_config.get('index_type', 'faiss'),
                quantization=index_config.get('faiss', {}).get('quantization', 'pq'),
                cache_size=index_config.get('cache_size', 10)
            )
        else:
            logger.warning(f"Hierarchical index not found at {index_storage_dir} or {families_index_dir}")
            self.hierarchical_index = None
        
        # Create FAISS indexes if they don't exist
        if self.hierarchical_index is None and self.storage is not None:
            logger.info("Creating FAISS indexes for available families...")
            self._create_faiss_indexes()
        
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
    
    def _create_faiss_indexes(self):
        """Create FAISS indexes for all available families."""
        if not self.storage or not self.available_families:
            logger.warning("No storage or families available for index creation")
            return
        
        logger.info("Creating FAISS indexes for all families...")
        created_count = 0
        
        for family_id in self.available_families:
            try:
                # Check if index already exists
                index_file = self.storage.base_dir / "indexes" / "families" / f"{family_id}.faiss"
                if index_file.exists():
                    logger.debug(f"FAISS index already exists for family {family_id}")
                    continue
                
                # Create the index
                self.storage.create_family_index_float(family_id)
                created_count += 1
                logger.info(f"Created FAISS index for family {family_id}")
                
            except Exception as e:
                logger.error(f"Failed to create FAISS index for family {family_id}: {e}")
                continue
        
        logger.info(f"Created {created_count} FAISS indexes")
        
        # Reinitialize hierarchical index after creating indexes
        index_storage_dir = self.storage.base_dir / "indexes"
        if index_storage_dir.exists():
            index_config = self.config.get('similarity_search', {})
            self.hierarchical_index = HierarchicalIndex(
                base_dir=str(index_storage_dir),
                index_type=index_config.get('index_type', 'faiss'),
                quantization=index_config.get('faiss', {}).get('quantization', 'pq'),
                cache_size=index_config.get('cache_size', 10)
            )
            logger.info("Reinitialized hierarchical index with created FAISS indexes")
    
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
        query_embedding = self.embedding_generator.generate_embedding(query_sequence)
        
        logger.info(f"Generated embedding with shape: {query_embedding.shape}")
        return query_embedding
    
    def classify_query_family(self, query_embedding: np.ndarray) -> Tuple[str, float]:
        """
        Classify query protein into a family using binary FAISS IVF (Hamming search on binary centroids).
        Args:
            query_embedding: Query protein embedding (must be np.float32)
        Returns:
            Tuple of (family_id, confidence_score)
        """
        logger.info("Classifying query protein into family using binary FAISS IVF centroids...")
        result = self.assign_protein_family.assign_family(query_embedding)
        return result['family_id'], result['confidence']
    
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
        Perform optimized similarity search within family using FAISS IVF float32 only.
        Args:
            query_embedding: Query protein embedding (np.float32, shape [D])
            family_id: Family ID to search within
            k: Number of similar proteins to retrieve
        Returns:
            List of similar proteins with metadata
        """
        logger.info(f"Performing optimized similarity search for top {k} proteins in family {family_id} using FAISS IVF float32...")
        if self.hierarchical_index is None:
            raise RuntimeError("FAISS hierarchical index not available for similarity search.")
        
        if query_embedding.dtype != np.float32:
            raise ValueError("Query embedding must be np.float32 for float FAISS search.")
        try:
            D, protein_ids = self.hierarchical_index.search_family_float(
                family_id, query_embedding, top_k=k
            )
            similar_proteins = []
            for i, (dist, protein_id) in enumerate(zip(D, protein_ids)):
                protein_result = {
                    'protein_id': protein_id,
                    'similarity_score': float(-dist),  # negative L2 distance for compatibility
                    'rank': i + 1
                }
                similar_proteins.append(protein_result)
            logger.info(f"Found {len(similar_proteins)} similar proteins using FAISS IVF float32 index")
            return similar_proteins
        except Exception as e:
            logger.error(f"FAISS IVF float32 index search failed: {e}")
            raise RuntimeError(f"FAISS IVF float32 index search failed: {e}")
    
    def build_optimized_network(self, query_embedding: np.ndarray,
                               query_protein_id: str,
                               similar_proteins: List[Dict],
                               family_id: str,
                               network_method: str = "mutual_knn") -> Tuple[nx.Graph, Dict]:
        """
        Build optimized network using family-specific data and FAISS IVF float32.
        Args:
            query_embedding: Query protein embedding
            query_protein_id: Query protein ID
            similar_proteins: List of similar proteins
            family_id: Family ID
            network_method: Network construction method
        Returns:
            Tuple of (NetworkX graph, network properties)
        """
        logger.info(f"Building optimized network for family {family_id} using FAISS IVF float32...")
        
        if self.hierarchical_index is None:
            raise RuntimeError("FAISS hierarchical index is required for network building.")
        
        family_embeddings, family_protein_ids = self.storage.load_family_embeddings(family_id)
        from .network_builder import DynamicNetworkBuilder
        builder = DynamicNetworkBuilder(**self.config.get('network', {}))
        G = None
        if network_method == "mutual_knn":
            G = builder.build_mutual_knn_network(family_embeddings, family_protein_ids, family_id=family_id, index=self.hierarchical_index)
        elif network_method == "threshold":
            G = builder.build_threshold_network(family_embeddings, family_protein_ids, family_id=family_id, index=self.hierarchical_index)
        elif network_method == "hybrid":
            G = builder.build_hybrid_network(family_embeddings, family_protein_ids, family_id=family_id, index=self.hierarchical_index)
        else:
            raise ValueError(f"Unknown network method: {network_method}")
        properties = builder.analyze_network_properties(G)
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

            # --- Generate interactive network HTML file ---
            # Use a dedicated directory for HTML reports
            html_dir = self.config.get('storage', {}).get('output_dir', 'results')
            html_filename = f"network_{query_protein_id}_{int(time.time())}.html"
            # Generate the HTML file using DynamicNetworkBuilder and storage
            builder = DynamicNetworkBuilder(**self.config.get('network', {}))
            fig, _ = builder.create_interactive_visualization(
                embeddings=family_embeddings,
                protein_ids=family_protein_ids,
                metadata_df=family_metadata,
                query_embedding=query_embedding,
                query_protein_id=query_protein_id,
                output_file=os.path.join(html_dir, html_filename),
                storage=self.storage
            )
            # Get the full path from storage (guaranteed to exist)
            html_file_path = self.storage.base_dir / html_dir / html_filename if not os.path.isabs(html_dir) else Path(html_dir) / html_filename
            results['network_html_file'] = str(html_file_path)
            results['family_metadata'] = family_metadata
            
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
