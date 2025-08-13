#!/usr/bin/env python3
"""
Comprehensive test suite for all query types

Tests all supported query types including protein existence checks,
similarity searches, family assignments, and network analysis queries
for the KBase Protein Query Module.
"""

import unittest
import tempfile
import os
import json
import numpy as np
from unittest.mock import patch, MagicMock
from pathlib import Path

# Add lib directory to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'lib'))

from kbase_protein_query_module.src.check_existence import ProteinExistenceChecker
from kbase_protein_query_module.src.similarity_index import HierarchicalIndex
from kbase_protein_query_module.src.assign_protein_family import AssignProteinFamily
from kbase_protein_query_module.src.network_builder import DynamicNetworkBuilder
from kbase_protein_query_module.src.embedding_generator import ProteinEmbeddingGenerator

class TestQueryTypes(unittest.TestCase):
    """Test cases for all query types."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Initialize query components
        self.check_existence = CheckExistence()
        self.similarity_index = SimilarityIndex()
        self.family_assigner = AssignProteinFamily()
        self.network_builder = NetworkBuilder()
        self.embedding_generator = EmbeddingGenerator()
        
        # Test data
        self.test_proteins = [
            {"protein_id": "P12345", "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"},
            {"protein_id": "P67890", "sequence": "MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"},
            {"protein_id": "P11111", "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"}
        ]
        
        self.test_embeddings = {
            "P12345": np.random.rand(1280).astype(np.float32),
            "P67890": np.random.rand(1280).astype(np.float32),
            "P11111": np.random.rand(1280).astype(np.float32)
        }
        
        # Create test index files
        self.test_index_file = os.path.join(self.temp_dir, 'test_index.faiss')
        self.test_centroids_file = os.path.join(self.temp_dir, 'test_centroids.npz')
        
        # Create test centroids file
        family_ids = np.array(['FAM1', 'FAM2', 'FAM3'])
        centroids = np.random.rand(3, 1280).astype(np.float32)
        eigenprotein_ids = np.array(['P12345', 'P67890', 'P11111'])
        np.savez(self.test_centroids_file, family_ids=family_ids, centroids=centroids, eigenprotein_ids=eigenprotein_ids)
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)

class TestExistenceQueries(unittest.TestCase):
    """Test cases for protein existence queries."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.check_existence = ProteinExistenceChecker()
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test database files
        self.test_db_file = os.path.join(self.temp_dir, 'test_db.fasta')
        with open(self.test_db_file, 'w') as f:
            f.write(">P12345|Test protein 1\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
            f.write(">P67890|Test protein 2\n")
            f.write("MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_single_protein_existence(self):
        """Test checking existence of a single protein."""
        protein_id = "P12345"
        
        with patch.object(self.check_existence, 'load_database') as mock_load:
            mock_load.return_value = {"P12345": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"}
            
            result = self.check_existence.check_protein_exists(protein_id)
            
            self.assertIsInstance(result, dict)
            self.assertIn('exists', result)
            self.assertIn('protein_id', result)
            self.assertTrue(result['exists'])
            self.assertEqual(result['protein_id'], protein_id)
    
    def test_multiple_protein_existence(self):
        """Test checking existence of multiple proteins."""
        protein_ids = ["P12345", "P67890", "P99999"]
        
        with patch.object(self.check_existence, 'load_database') as mock_load:
            mock_load.return_value = {
                "P12345": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
                "P67890": "MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
            }
            
            results = self.check_existence.check_multiple_proteins(protein_ids)
            
            self.assertIsInstance(results, list)
            self.assertEqual(len(results), 3)
            
            # Check existing proteins
            existing_proteins = [r for r in results if r['exists']]
            self.assertEqual(len(existing_proteins), 2)
            
            # Check non-existing protein
            non_existing = [r for r in results if not r['exists']]
            self.assertEqual(len(non_existing), 1)
            self.assertEqual(non_existing[0]['protein_id'], "P99999")
    
    def test_sequence_based_existence(self):
        """Test checking existence based on sequence similarity."""
        query_sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        
        with patch.object(self.check_existence, 'load_database') as mock_load:
            mock_load.return_value = {
                "P12345": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
                "P67890": "MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
            }
            
            result = self.check_existence.check_sequence_exists(query_sequence, similarity_threshold=0.9)
            
            self.assertIsInstance(result, dict)
            self.assertIn('exists', result)
            self.assertIn('matches', result)
            self.assertIsInstance(result['matches'], list)
    
    def test_existence_with_metadata(self):
        """Test existence check with metadata retrieval."""
        protein_id = "P12345"
        
        with patch.object(self.check_existence, 'load_database_with_metadata') as mock_load:
            mock_load.return_value = {
                "P12345": {
                    "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
                    "description": "Test protein 1",
                    "organism": "Escherichia coli",
                    "length": 40
                }
            }
            
            result = self.check_existence.check_protein_with_metadata(protein_id)
            
            self.assertIsInstance(result, dict)
            self.assertIn('exists', result)
            self.assertIn('metadata', result)
            self.assertIn('sequence', result['metadata'])
            self.assertIn('description', result['metadata'])
            self.assertIn('organism', result['metadata'])
            self.assertIn('length', result['metadata'])

class TestSimilarityQueries(unittest.TestCase):
    """Test cases for similarity search queries."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.similarity_index = SimilarityIndex()
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test embeddings
        self.test_embeddings = {
            "P12345": np.random.rand(1280).astype(np.float32),
            "P67890": np.random.rand(1280).astype(np.float32),
            "P11111": np.random.rand(1280).astype(np.float32),
            "P22222": np.random.rand(1280).astype(np.float32),
            "P33333": np.random.rand(1280).astype(np.float32)
        }
        
        # Create test index
        self.test_index_file = os.path.join(self.temp_dir, 'test_index.faiss')
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_similarity_search_single_query(self):
        """Test similarity search with single query embedding."""
        query_embedding = np.random.rand(1280).astype(np.float32)
        
        with patch.object(self.similarity_index, 'load_index') as mock_load:
            mock_load.return_value = MagicMock()
            
            with patch.object(self.similarity_index, 'search') as mock_search:
                mock_search.return_value = {
                    'distances': np.array([[0.1, 0.2, 0.3]]),
                    'indices': np.array([[0, 1, 2]])
                }
                
                result = self.similarity_index.search_similar_proteins(
                    query_embedding, top_k=5, similarity_threshold=0.5
                )
                
                self.assertIsInstance(result, dict)
                self.assertIn('matches', result)
                self.assertIn('scores', result)
                self.assertIsInstance(result['matches'], list)
                self.assertIsInstance(result['scores'], list)
    
    def test_similarity_search_multiple_queries(self):
        """Test similarity search with multiple query embeddings."""
        query_embeddings = {
            "query1": np.random.rand(1280).astype(np.float32),
            "query2": np.random.rand(1280).astype(np.float32)
        }
        
        with patch.object(self.similarity_index, 'load_index') as mock_load:
            mock_load.return_value = MagicMock()
            
            with patch.object(self.similarity_index, 'batch_search') as mock_search:
                mock_search.return_value = {
                    'query1': {
                        'matches': ['P12345', 'P67890'],
                        'scores': [0.95, 0.85]
                    },
                    'query2': {
                        'matches': ['P11111', 'P22222'],
                        'scores': [0.92, 0.78]
                    }
                }
                
                result = self.similarity_index.batch_similarity_search(
                    query_embeddings, top_k=5, similarity_threshold=0.5
                )
                
                self.assertIsInstance(result, dict)
                self.assertIn('query1', result)
                self.assertIn('query2', result)
                self.assertIn('matches', result['query1'])
                self.assertIn('scores', result['query1'])
    
    def test_similarity_search_with_filters(self):
        """Test similarity search with additional filters."""
        query_embedding = np.random.rand(1280).astype(np.float32)
        
        filters = {
            'organism': 'Escherichia coli',
            'length_min': 100,
            'length_max': 500,
            'family': 'FAM1'
        }
        
        with patch.object(self.similarity_index, 'load_index') as mock_load:
            mock_load.return_value = MagicMock()
            
            with patch.object(self.similarity_index, 'filtered_search') as mock_search:
                mock_search.return_value = {
                    'matches': ['P12345', 'P67890'],
                    'scores': [0.95, 0.85],
                    'metadata': [
                        {'organism': 'Escherichia coli', 'length': 150},
                        {'organism': 'Escherichia coli', 'length': 200}
                    ]
                }
                
                result = self.similarity_index.search_with_filters(
                    query_embedding, filters, top_k=5
                )
                
                self.assertIsInstance(result, dict)
                self.assertIn('matches', result)
                self.assertIn('scores', result)
                self.assertIn('metadata', result)
    
    def test_similarity_threshold_search(self):
        """Test similarity search with different thresholds."""
        query_embedding = np.random.rand(1280).astype(np.float32)
        
        with patch.object(self.similarity_index, 'load_index') as mock_load:
            mock_load.return_value = MagicMock()
            
            # Test high threshold
            with patch.object(self.similarity_index, 'search') as mock_search:
                mock_search.return_value = {
                    'distances': np.array([[0.1, 0.2, 0.3]]),
                    'indices': np.array([[0, 1, 2]])
                }
                
                result_high = self.similarity_index.search_similar_proteins(
                    query_embedding, top_k=5, similarity_threshold=0.9
                )
                
                # Test low threshold
                result_low = self.similarity_index.search_similar_proteins(
                    query_embedding, top_k=5, similarity_threshold=0.3
                )
                
                # High threshold should return fewer or equal matches
                self.assertLessEqual(len(result_high['matches']), len(result_low['matches']))

class TestFamilyAssignmentQueries(unittest.TestCase):
    """Test cases for protein family assignment queries."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.family_assigner = AssignProteinFamily()
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test centroids file
        self.test_centroids_file = os.path.join(self.temp_dir, 'test_centroids.npz')
        family_ids = np.array(['FAM1', 'FAM2', 'FAM3'])
        centroids = np.random.rand(3, 1280).astype(np.float32)
        eigenprotein_ids = np.array(['P12345', 'P67890', 'P11111'])
        np.savez(self.test_centroids_file, family_ids=family_ids, centroids=centroids, eigenprotein_ids=eigenprotein_ids)
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_single_protein_family_assignment(self):
        """Test family assignment for a single protein."""
        protein_embedding = np.random.rand(1280).astype(np.float32)
        
        self.family_assigner.load_family_centroids(self.test_centroids_file)
        
        result = self.family_assigner.assign_family(protein_embedding)
        
        self.assertIsInstance(result, dict)
        self.assertIn('family_id', result)
        self.assertIn('confidence', result)
        self.assertIn('eigenprotein_id', result)
        self.assertIsInstance(result['family_id'], str)
        self.assertIsInstance(result['confidence'], float)
        self.assertGreaterEqual(result['confidence'], 0.0)
        self.assertLessEqual(result['confidence'], 1.0)
    
    def test_multiple_protein_family_assignment(self):
        """Test family assignment for multiple proteins."""
        protein_embeddings = {
            "P12345": np.random.rand(1280).astype(np.float32),
            "P67890": np.random.rand(1280).astype(np.float32),
            "P11111": np.random.rand(1280).astype(np.float32)
        }
        
        self.family_assigner.load_family_centroids(self.test_centroids_file)
        
        results = self.family_assigner.assign_families_batch(protein_embeddings)
        
        self.assertIsInstance(results, dict)
        self.assertEqual(len(results), 3)
        
        for protein_id, result in results.items():
            self.assertIn('family_id', result)
            self.assertIn('confidence', result)
            self.assertIn('eigenprotein_id', result)
    
    def test_family_assignment_with_confidence_threshold(self):
        """Test family assignment with confidence threshold."""
        protein_embedding = np.random.rand(1280).astype(np.float32)
        
        self.family_assigner.load_family_centroids(self.test_centroids_file)
        
        # Test with high confidence threshold
        result_high = self.family_assigner.assign_family_with_threshold(
            protein_embedding, confidence_threshold=0.9
        )
        
        # Test with low confidence threshold
        result_low = self.family_assigner.assign_family_with_threshold(
            protein_embedding, confidence_threshold=0.3
        )
        
        # High threshold might return no assignment
        if result_high['family_id'] is not None:
            self.assertGreaterEqual(result_high['confidence'], 0.9)
        
        # Low threshold should always return an assignment
        self.assertIsNotNone(result_low['family_id'])
        self.assertGreaterEqual(result_low['confidence'], 0.3)
    
    def test_family_assignment_with_metadata(self):
        """Test family assignment with additional metadata."""
        protein_embedding = np.random.rand(1280).astype(np.float32)
        protein_metadata = {
            'protein_id': 'P12345',
            'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ',
            'organism': 'Escherichia coli'
        }
        
        self.family_assigner.load_family_centroids(self.test_centroids_file)
        
        result = self.family_assigner.assign_family_with_metadata(
            protein_embedding, protein_metadata
        )
        
        self.assertIsInstance(result, dict)
        self.assertIn('family_id', result)
        self.assertIn('confidence', result)
        self.assertIn('eigenprotein_id', result)
        self.assertIn('metadata', result)
        self.assertEqual(result['metadata']['protein_id'], 'P12345')

class TestNetworkAnalysisQueries(unittest.TestCase):
    """Test cases for network analysis queries."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.network_builder = NetworkBuilder()
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test similarity matrix
        self.test_similarity_matrix = np.array([
            [1.0, 0.8, 0.6, 0.4],
            [0.8, 1.0, 0.7, 0.5],
            [0.6, 0.7, 1.0, 0.9],
            [0.4, 0.5, 0.9, 1.0]
        ])
        
        self.test_protein_ids = ['P12345', 'P67890', 'P11111', 'P22222']
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_network_construction(self):
        """Test network construction from similarity matrix."""
        result = self.network_builder.build_network(
            self.test_similarity_matrix,
            self.test_protein_ids,
            k_neighbors=2,
            similarity_threshold=0.5
        )
        
        self.assertIsInstance(result, dict)
        self.assertIn('nodes', result)
        self.assertIn('edges', result)
        self.assertIn('network_metrics', result)
        
        # Check nodes
        self.assertIsInstance(result['nodes'], list)
        self.assertEqual(len(result['nodes']), 4)
        
        # Check edges
        self.assertIsInstance(result['edges'], list)
        self.assertGreater(len(result['edges']), 0)
        
        # Check metrics
        self.assertIsInstance(result['network_metrics'], dict)
        self.assertIn('density', result['network_metrics'])
        self.assertIn('clustering_coefficient', result['network_metrics'])
    
    def test_network_analysis_with_different_methods(self):
        """Test network analysis with different construction methods."""
        # Test KNN method
        result_knn = self.network_builder.build_network(
            self.test_similarity_matrix,
            self.test_protein_ids,
            network_method='knn',
            k_neighbors=2
        )
        
        # Test mutual KNN method
        result_mutual_knn = self.network_builder.build_network(
            self.test_similarity_matrix,
            self.test_protein_ids,
            network_method='mutual_knn',
            k_neighbors=2
        )
        
        # Test threshold method
        result_threshold = self.network_builder.build_network(
            self.test_similarity_matrix,
            self.test_protein_ids,
            network_method='threshold',
            similarity_threshold=0.7
        )
        
        # All should return valid networks
        for result in [result_knn, result_mutual_knn, result_threshold]:
            self.assertIn('nodes', result)
            self.assertIn('edges', result)
            self.assertIn('network_metrics', result)
    
    def test_network_community_detection(self):
        """Test community detection in protein networks."""
        result = self.network_builder.build_network(
            self.test_similarity_matrix,
            self.test_protein_ids,
            k_neighbors=2,
            similarity_threshold=0.5
        )
        
        communities = self.network_builder.detect_communities(result)
        
        self.assertIsInstance(communities, dict)
        self.assertIn('communities', communities)
        self.assertIn('modularity', communities)
        self.assertIsInstance(communities['communities'], list)
        self.assertIsInstance(communities['modularity'], float)
    
    def test_network_centrality_analysis(self):
        """Test centrality analysis of protein networks."""
        result = self.network_builder.build_network(
            self.test_similarity_matrix,
            self.test_protein_ids,
            k_neighbors=2,
            similarity_threshold=0.5
        )
        
        centrality = self.network_builder.calculate_centrality(result)
        
        self.assertIsInstance(centrality, dict)
        self.assertIn('degree_centrality', centrality)
        self.assertIn('betweenness_centrality', centrality)
        self.assertIn('closeness_centrality', centrality)
        self.assertIn('eigenvector_centrality', centrality)
        
        # Check that centrality values are calculated for each node
        for centrality_type, values in centrality.items():
            self.assertIsInstance(values, dict)
            self.assertEqual(len(values), 4)  # 4 proteins
    
    def test_network_export(self):
        """Test network export to different formats."""
        result = self.network_builder.build_network(
            self.test_similarity_matrix,
            self.test_protein_ids,
            k_neighbors=2,
            similarity_threshold=0.5
        )
        
        # Test JSON export
        json_file = os.path.join(self.temp_dir, 'network.json')
        self.network_builder.export_network(result, json_file, format='json')
        self.assertTrue(os.path.exists(json_file))
        
        # Test CSV export
        csv_file = os.path.join(self.temp_dir, 'network.csv')
        self.network_builder.export_network(result, csv_file, format='csv')
        self.assertTrue(os.path.exists(csv_file))
        
        # Test GML export
        gml_file = os.path.join(self.temp_dir, 'network.gml')
        self.network_builder.export_network(result, gml_file, format='gml')
        self.assertTrue(os.path.exists(gml_file))

if __name__ == '__main__':
    unittest.main()
