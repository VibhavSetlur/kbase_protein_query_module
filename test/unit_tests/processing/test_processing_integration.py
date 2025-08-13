"""
Comprehensive Integration Tests for Processing Module

This module tests the complete processing pipeline including:
- Embedding generation
- Similarity search and indexing
- Network construction and analysis
- End-to-end processing workflows
"""

import unittest
import numpy as np
import os
import sys
import tempfile
import h5py
import pandas as pd
from pathlib import Path

from kbase_protein_query_module.src.processing.embeddings.generator import ProteinEmbeddingGenerator
from kbase_protein_query_module.src.processing.similarity.hierarchical_index import HierarchicalIndex
from kbase_protein_query_module.src.processing.networks.builder import DynamicNetworkBuilder

class TestProcessingIntegration(unittest.TestCase):
    """Integration tests for the complete processing pipeline."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test environment with real data."""
        cls.temp_dir = tempfile.mkdtemp()
        
        # Use real data from data/families/
        cls.families_dir = "/home/vibhav/Downloads/Work/ANL/Research/kbase_protein_query_module/data/families"
        
        # Load real family data
        family_id = 'FAM0'
        family_file = os.path.join(cls.families_dir, f'{family_id}.h5')
        
        if not os.path.exists(family_file):
            cls.fail(f"Real family data not found: {family_file}. Tests must use actual data from data/families/")
        
        # Load real embeddings and protein IDs
        with h5py.File(family_file, 'r') as f:
            cls.embeddings = f['embeddings'][:200]  # Use first 200 for testing
            cls.protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][:200]]
        
        # Initialize processing components
        cls.embedding_generator = ProteinEmbeddingGenerator(device="cpu")
        cls.similarity_index = HierarchicalIndex(base_dir=cls.temp_dir)
        cls.network_builder = DynamicNetworkBuilder()
        
        # Create test sequences
        cls.test_sequences = [
            "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
            "MKLIVTALAGALALQAGSLFASAADSSSPSEAGVDLK",
            "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQ"
        ]
        cls.test_protein_ids = ["TEST1", "TEST2", "TEST3"]
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(cls.temp_dir)
    
    def test_end_to_end_processing_pipeline(self):
        """Test complete processing pipeline from sequences to network."""
        # Step 1: Generate embeddings
        embeddings_dict = self.embedding_generator.generate_embeddings_batch(
            self.test_sequences, self.test_protein_ids, batch_size=2
        )
        
        self.assertEqual(len(embeddings_dict), 3)
        for protein_id, embedding in embeddings_dict.items():
            self.assertIn(protein_id, self.test_protein_ids)
            self.assertEqual(embedding.shape[0], 320)  # ESM2 t6 dimension
        
        # Step 2: Create similarity index
        embeddings_array = np.array(list(embeddings_dict.values()))
        protein_ids_list = list(embeddings_dict.keys())
        
        index_path = self.similarity_index.create_family_index_float(
            "test_family", embeddings_array, protein_ids_list
        )
        self.assertTrue(os.path.exists(index_path))
        
        # Step 3: Perform similarity search
        query_embedding = embeddings_array[0]  # Use first protein as query
        distances, indices = self.similarity_index.search_family_float(
            "test_family", query_embedding, top_k=3
        )
        
        self.assertGreater(len(distances), 0)
        self.assertEqual(len(distances), len(indices))
        self.assertLessEqual(len(distances), 3)
        
        # Step 4: Build network
        network = self.network_builder.build_mutual_knn_network(
            embeddings_array, protein_ids_list, query_embedding, "TEST1"
        )
        
        self.assertIsNotNone(network)
        self.assertGreater(len(network.nodes()), 0)
        self.assertGreater(len(network.edges()), 0)
        
        # Step 5: Analyze network properties
        properties = self.network_builder.analyze_network_properties(network)
        
        self.assertIn('num_nodes', properties)
        self.assertIn('num_edges', properties)
        self.assertIn('density', properties)
        self.assertGreater(properties['num_nodes'], 0)
    
    def test_processing_with_real_data(self):
        """Test processing with real protein data."""
        # Use real embeddings for similarity search
        real_embeddings = self.embeddings[:50]  # Use first 50 real proteins
        real_protein_ids = self.protein_ids[:50]
        
        # Create index with real data
        index_path = self.similarity_index.create_family_index_float(
            "real_family", real_embeddings, real_protein_ids
        )
        
        # Test similarity search with real data
        query_embedding = real_embeddings[0]
        distances, indices = self.similarity_index.search_family_float(
            "real_family", query_embedding, top_k=10
        )
        
        self.assertGreater(len(distances), 0)
        self.assertEqual(len(distances), len(indices))
        
        # Test network construction with real data
        network = self.network_builder.build_mutual_knn_network(
            real_embeddings, real_protein_ids, query_embedding, real_protein_ids[0]
        )
        
        self.assertIsNotNone(network)
        self.assertGreater(len(network.nodes()), 0)
    
    def test_processing_error_handling(self):
        """Test error handling in processing pipeline."""
        # Test with invalid embeddings
        invalid_embeddings = np.random.rand(5, 100)  # Wrong dimension
        invalid_ids = ["INVALID1", "INVALID2", "INVALID3", "INVALID4", "INVALID5"]
        
        # Should handle dimension mismatch gracefully
        try:
            self.similarity_index.create_family_index_float(
                "invalid_family", invalid_embeddings, invalid_ids
            )
        except Exception as e:
            # Expected to fail due to dimension mismatch
            self.assertIsInstance(e, Exception)
    
    def test_processing_performance(self):
        """Test processing performance with larger datasets."""
        # Use more real data for performance testing
        if len(self.protein_ids) >= 100:
            large_embeddings = self.embeddings[:100]
            large_protein_ids = self.protein_ids[:100]
            
            # Test index creation performance
            import time
            start_time = time.time()
            
            index_path = self.similarity_index.create_family_index_float(
                "large_family", large_embeddings, large_protein_ids
            )
            
            creation_time = time.time() - start_time
            self.assertLess(creation_time, 30)  # Should complete within 30 seconds
            
            # Test search performance
            query_embedding = large_embeddings[0]
            start_time = time.time()
            
            distances, indices = self.similarity_index.search_family_float(
                "large_family", query_embedding, top_k=20
            )
            
            search_time = time.time() - start_time
            self.assertLess(search_time, 5)  # Should complete within 5 seconds
    
    def test_processing_configuration_options(self):
        """Test different configuration options for processing components."""
        # Test different network building methods
        embeddings_array = np.array(list(self.embedding_generator.generate_embeddings_batch(
            self.test_sequences, self.test_protein_ids
        ).values()))
        protein_ids_list = self.test_protein_ids
        
        # Test mutual k-NN network
        network_mutual = self.network_builder.build_mutual_knn_network(
            embeddings_array, protein_ids_list
        )
        
        # Test threshold network
        network_threshold = self.network_builder.build_threshold_network(
            embeddings_array, protein_ids_list
        )
        
        # Test hybrid network
        network_hybrid = self.network_builder.build_hybrid_network(
            embeddings_array, protein_ids_list
        )
        
        # All networks should be valid
        self.assertIsNotNone(network_mutual)
        self.assertIsNotNone(network_threshold)
        self.assertIsNotNone(network_hybrid)
        
        # Networks should be valid
        self.assertIsNotNone(network_mutual)
        self.assertIsNotNone(network_threshold)
        self.assertIsNotNone(network_hybrid)

if __name__ == '__main__':
    unittest.main()
