#!/usr/bin/env python3
"""
Comprehensive integration test suite for the complete workflow pipeline

Tests the entire workflow from input to output, including all stages,
error handling, and edge cases for the KBase Protein Query Module.
"""

import unittest
import tempfile
import os
import json
import numpy as np
import time
from unittest.mock import patch, MagicMock
from pathlib import Path

# Add lib directory to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'lib'))

from kbase_protein_query_module.src.unified_workflow_orchestrator import UnifiedProteinQueryWorkflow, PipelineConfig
from kbase_protein_query_module.src.workflow_orchestrator import ProteinNetworkWorkflow
from kbase_protein_query_module.src.input_parser import InputParser
from kbase_protein_query_module.src.embedding_generator import ProteinEmbeddingGenerator
from kbase_protein_query_module.src.similarity_index import HierarchicalIndex
from kbase_protein_query_module.src.assign_protein_family import AssignProteinFamily
from kbase_protein_query_module.src.network_builder import DynamicNetworkBuilder
from kbase_protein_query_module.src.html_report_generator import HTMLReportGenerator

class TestWorkflowIntegration(unittest.TestCase):
    """Test cases for complete workflow integration."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test data files
        self.test_fasta_file = os.path.join(self.temp_dir, 'test_proteins.fasta')
        self.test_csv_file = os.path.join(self.temp_dir, 'test_proteins.csv')
        self.test_json_file = os.path.join(self.temp_dir, 'test_proteins.json')
        
        # Create test FASTA file
        with open(self.test_fasta_file, 'w') as f:
            f.write(">protein1|P12345|Test protein 1\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
            f.write(">protein2|P67890|Test protein 2\n")
            f.write("MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
            f.write(">protein3|P11111|Test protein 3\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
        
        # Create test CSV file
        with open(self.test_csv_file, 'w') as f:
            f.write("protein_id,sequence,description\n")
            f.write("P12345,MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ,Test protein 1\n")
            f.write("P67890,MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ,Test protein 2\n")
        
        # Create test JSON file
        test_data = {
            "proteins": [
                {
                    "protein_id": "P12345",
                    "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
                    "description": "Test protein 1"
                },
                {
                    "protein_id": "P67890",
                    "sequence": "MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
                    "description": "Test protein 2"
                }
            ]
        }
        with open(self.test_json_file, 'w') as f:
            json.dump(test_data, f)
        
        # Initialize workflow orchestrators
        self.workflow_orchestrator = WorkflowOrchestrator()
        
        # Test configuration
        self.test_config = {
            'model_name': 'esm2_t6_8M_UR50D',
            'embedding_dim': 1280,
            'pooling_method': 'mean',
            'k_neighbors': 8,
            'similarity_threshold': 0.5,
            'network_method': 'mutual_knn',
            'output_dir': self.temp_dir
        }
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_complete_fasta_workflow(self):
        """Test complete workflow with FASTA input."""
        input_params = {
            'input_type': 'fasta',
            'input_path': self.test_fasta_file,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'protein1': np.random.rand(1280).astype(np.float32),
                'protein2': np.random.rand(1280).astype(np.float32),
                'protein3': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    # Check result structure
                    self.assertIn('status', result)
                    self.assertIn('results', result)
                    self.assertIn('output_files', result)
                    
                    # Check status
                    self.assertEqual(result['status'], 'success')
                    
                    # Check results
                    self.assertIn('proteins', result['results'])
                    self.assertIn('embeddings', result['results'])
                    self.assertIn('family_assignments', result['results'])
                    self.assertIn('network_data', result['results'])
                    
                    # Check output files
                    self.assertIsInstance(result['output_files'], list)
    
    def test_complete_csv_workflow(self):
        """Test complete workflow with CSV input."""
        input_params = {
            'input_type': 'csv',
            'input_path': self.test_csv_file,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'P12345': np.random.rand(1280).astype(np.float32),
                'P67890': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    self.assertEqual(result['status'], 'success')
                    self.assertIn('proteins', result['results'])
                    self.assertEqual(len(result['results']['proteins']), 2)
    
    def test_complete_json_workflow(self):
        """Test complete workflow with JSON input."""
        input_params = {
            'input_type': 'json',
            'input_path': self.test_json_file,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'P12345': np.random.rand(1280).astype(np.float32),
                'P67890': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    self.assertEqual(result['status'], 'success')
                    self.assertIn('proteins', result['results'])
    
    def test_sequence_string_workflow(self):
        """Test workflow with sequence string input."""
        input_params = {
            'input_type': 'sequence',
            'sequences': [
                "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
                "MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
            ],
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'sequence_0': np.random.rand(1280).astype(np.float32),
                'sequence_1': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    self.assertEqual(result['status'], 'success')
                    self.assertIn('proteins', result['results'])
                    self.assertEqual(len(result['results']['proteins']), 2)
    
    def test_protein_id_workflow(self):
        """Test workflow with protein ID input."""
        input_params = {
            'input_type': 'protein_id',
            'protein_ids': ['P12345', 'P67890'],
            'workspace_name': 'test_workspace',
            'config': self.test_config
        }
        
        with patch.object(self.unified_orchestrator, 'fetch_proteins_from_workspace') as mock_fetch:
            mock_fetch.return_value = [
                {'protein_id': 'P12345', 'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'},
                {'protein_id': 'P67890', 'sequence': 'MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'}
            ]
            
            with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
                mock_embeddings.return_value = {
                    'P12345': np.random.rand(1280).astype(np.float32),
                    'P67890': np.random.rand(1280).astype(np.float32)
                }
                
                with patch.object(SimilarityIndex, 'build_index') as mock_index:
                    mock_index.return_value = MagicMock()
                    
                    with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                        mock_family.return_value = {
                            'family_id': 'FAM1',
                            'confidence': 0.95,
                            'eigenprotein_id': 'P12345'
                        }
                        
                        result = self.unified_orchestrator.run_workflow(input_params)
                        
                        self.assertEqual(result['status'], 'success')
                        self.assertIn('proteins', result['results'])
    
    def test_workflow_with_different_configurations(self):
        """Test workflow with different configuration options."""
        # Test with different model
        config1 = self.test_config.copy()
        config1['model_name'] = 'esm2_t12_35M_UR50D'
        
        # Test with different pooling method
        config2 = self.test_config.copy()
        config2['pooling_method'] = 'cls'
        
        # Test with different network parameters
        config3 = self.test_config.copy()
        config3['k_neighbors'] = 5
        config3['similarity_threshold'] = 0.7
        config3['network_method'] = 'knn'
        
        input_params = {
            'input_type': 'fasta',
            'input_path': self.test_fasta_file,
            'config': config1
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'protein1': np.random.rand(1280).astype(np.float32),
                'protein2': np.random.rand(1280).astype(np.float32),
                'protein3': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    self.assertEqual(result['status'], 'success')
    
    def test_workflow_error_handling(self):
        """Test workflow error handling."""
        # Test invalid input file
        input_params = {
            'input_type': 'fasta',
            'input_path': '/nonexistent/file.fasta',
            'config': self.test_config
        }
        
        result = self.unified_orchestrator.run_workflow(input_params)
        self.assertEqual(result['status'], 'error')
        self.assertIn('error_message', result)
        
        # Test invalid input type
        input_params = {
            'input_type': 'invalid_type',
            'input_path': self.test_fasta_file,
            'config': self.test_config
        }
        
        result = self.unified_orchestrator.run_workflow(input_params)
        self.assertEqual(result['status'], 'error')
        self.assertIn('error_message', result)
        
        # Test embedding generation error
        input_params = {
            'input_type': 'fasta',
            'input_path': self.test_fasta_file,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings', side_effect=Exception("Embedding error")):
            result = self.unified_orchestrator.run_workflow(input_params)
            self.assertEqual(result['status'], 'error')
            self.assertIn('error_message', result)
    
    def test_workflow_performance(self):
        """Test workflow performance with large datasets."""
        # Create large test file
        large_fasta_file = os.path.join(self.temp_dir, 'large_test.fasta')
        with open(large_fasta_file, 'w') as f:
            for i in range(100):  # 100 proteins
                f.write(f">protein{i}|P{i:05d}|Test protein {i}\n")
                f.write("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
        
        input_params = {
            'input_type': 'fasta',
            'input_path': large_fasta_file,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                f'protein{i}': np.random.rand(1280).astype(np.float32)
                for i in range(100)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    start_time = time.time()
                    result = self.unified_orchestrator.run_workflow(input_params)
                    end_time = time.time()
                    
                    self.assertEqual(result['status'], 'success')
                    self.assertIn('proteins', result['results'])
                    self.assertEqual(len(result['results']['proteins']), 100)
                    
                    # Check performance (should complete within reasonable time)
                    execution_time = end_time - start_time
                    self.assertLess(execution_time, 60)  # Should complete within 60 seconds
    
    def test_workflow_output_validation(self):
        """Test workflow output validation."""
        input_params = {
            'input_type': 'fasta',
            'input_path': self.test_fasta_file,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'protein1': np.random.rand(1280).astype(np.float32),
                'protein2': np.random.rand(1280).astype(np.float32),
                'protein3': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    # Validate output structure
                    self.assertIn('status', result)
                    self.assertIn('results', result)
                    self.assertIn('output_files', result)
                    
                    # Validate results structure
                    results = result['results']
                    self.assertIn('proteins', results)
                    self.assertIn('embeddings', results)
                    self.assertIn('family_assignments', results)
                    self.assertIn('network_data', results)
                    
                    # Validate protein data
                    proteins = results['proteins']
                    self.assertIsInstance(proteins, list)
                    self.assertGreater(len(proteins), 0)
                    
                    for protein in proteins:
                        self.assertIn('protein_id', protein)
                        self.assertIn('sequence', protein)
                        self.assertIsInstance(protein['protein_id'], str)
                        self.assertIsInstance(protein['sequence'], str)
                    
                    # Validate embeddings
                    embeddings = results['embeddings']
                    self.assertIsInstance(embeddings, dict)
                    self.assertGreater(len(embeddings), 0)
                    
                    for protein_id, embedding in embeddings.items():
                        self.assertIsInstance(protein_id, str)
                        self.assertIsInstance(embedding, np.ndarray)
                        self.assertEqual(embedding.shape, (1280,))
                        self.assertEqual(embedding.dtype, np.float32)
                    
                    # Validate family assignments
                    family_assignments = results['family_assignments']
                    self.assertIsInstance(family_assignments, dict)
                    
                    # Validate network data
                    network_data = results['network_data']
                    self.assertIsInstance(network_data, dict)
                    self.assertIn('nodes', network_data)
                    self.assertIn('edges', network_data)
    
    def test_workflow_edge_cases(self):
        """Test workflow edge cases."""
        # Test single protein
        single_protein_fasta = os.path.join(self.temp_dir, 'single.fasta')
        with open(single_protein_fasta, 'w') as f:
            f.write(">protein1|P12345|Single protein\n")
            f.write("M\n")  # Single amino acid
        
        input_params = {
            'input_type': 'fasta',
            'input_path': single_protein_fasta,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'protein1': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    self.assertEqual(result['status'], 'success')
                    self.assertIn('proteins', result['results'])
                    self.assertEqual(len(result['results']['proteins']), 1)
        
        # Test very long sequence
        long_sequence_fasta = os.path.join(self.temp_dir, 'long.fasta')
        with open(long_sequence_fasta, 'w') as f:
            f.write(">protein1|P12345|Long protein\n")
            f.write("M" * 10000)  # Very long sequence
        
        input_params = {
            'input_type': 'fasta',
            'input_path': long_sequence_fasta,
            'config': self.test_config
        }
        
        with patch.object(EmbeddingGenerator, 'generate_embeddings') as mock_embeddings:
            mock_embeddings.return_value = {
                'protein1': np.random.rand(1280).astype(np.float32)
            }
            
            with patch.object(SimilarityIndex, 'build_index') as mock_index:
                mock_index.return_value = MagicMock()
                
                with patch.object(AssignProteinFamily, 'assign_family') as mock_family:
                    mock_family.return_value = {
                        'family_id': 'FAM1',
                        'confidence': 0.95,
                        'eigenprotein_id': 'P12345'
                    }
                    
                    result = self.unified_orchestrator.run_workflow(input_params)
                    
                    self.assertEqual(result['status'], 'success')
                    self.assertIn('proteins', result['results'])
                    self.assertEqual(len(result['results']['proteins']), 1)

if __name__ == '__main__':
    unittest.main()
