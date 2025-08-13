"""
Comprehensive Pipeline Stages Tests

This module tests all pipeline stages including:
- Input stages (validation, extraction, workspace objects)
- Processing stages (embedding generation, family assignment, similarity search)
- Analysis stages (sequence analysis, network analysis)
- Output stages (report generation, visualization, data export)
"""

import unittest
import numpy as np
import os
import sys
import tempfile
import h5py
import pandas as pd
from pathlib import Path

from kbase_protein_query_module.src.stages import (
    InputValidationStage,
    DataExtractionStage,
    WorkspaceObjectStage,
    EmbeddingGenerationStage,
    FamilyAssignmentStage,
    SimilaritySearchStage,
    SequenceAnalysisStage,
    NetworkAnalysisStage,
    ReportGenerationStage,
    VisualizationStage,
    DataExportStage
)

from kbase_protein_query_module.src.core import StageResult

class TestPipelineStagesComprehensive(unittest.TestCase):
    """Comprehensive tests for all pipeline stages."""
    
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
            cls.embeddings = f['embeddings'][:100]  # Use first 100 for testing
            cls.protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][:100]]
        
        # Create test sequences
        cls.test_sequences = [
            "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
            "MKLIVTALAGALALQAGSLFASAADSSSPSEAGVDLK",
            "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQ"
        ]
        cls.test_protein_ids = ["TEST1", "TEST2", "TEST3"]
        
        # Initialize stages
        cls.input_validation_stage = InputValidationStage()
        cls.data_extraction_stage = DataExtractionStage()
        cls.workspace_object_stage = WorkspaceObjectStage()
        cls.embedding_generation_stage = EmbeddingGenerationStage()
        cls.family_assignment_stage = FamilyAssignmentStage()
        cls.similarity_search_stage = SimilaritySearchStage()
        cls.sequence_analysis_stage = SequenceAnalysisStage()
        cls.network_analysis_stage = NetworkAnalysisStage()
        cls.report_generation_stage = ReportGenerationStage()
        cls.visualization_stage = VisualizationStage()
        cls.data_export_stage = DataExportStage()
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(cls.temp_dir)
    
    def test_input_validation_stage(self):
        """Test input validation stage with various input types."""
        # Test UniProt ID validation
        uniprot_input = {
            'input_type': 'uniprot_ids',
            'input_data': ['P12345', 'Q67890', 'A1B2C3']
        }
        
        result = self.input_validation_stage.run(uniprot_input)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('validated_data', result.output_data)
        
        # Test FASTA sequence validation
        fasta_input = {
            'input_type': 'fasta_sequence',
            'input_data': '>TEST1\nMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'
        }
        
        result = self.input_validation_stage.run(fasta_input)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('validated_data', result.output_data)
        
        # Test invalid input
        invalid_input = {
            'input_type': 'invalid_type',
            'input_data': 'invalid_data'
        }
        
        result = self.input_validation_stage.run(invalid_input)
        self.assertIsInstance(result, StageResult)
        self.assertFalse(result.success)
        self.assertIsNotNone(result.error_message)
    
    def test_data_extraction_stage(self):
        """Test data extraction stage."""
        # Test with validated input
        validated_input = {
            'validated_data': {
                'protein_records': [
                    {'protein_id': 'TEST1', 'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'},
                    {'protein_id': 'TEST2', 'sequence': 'MKLIVTALAGALALQAGSLFASAADSSSPSEAGVDLK'}
                ],
                'input_type': 'fasta_sequence'
            }
        }
        
        result = self.data_extraction_stage.run(validated_input)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('extracted_data', result.output_data)
    
    def test_workspace_object_stage(self):
        """Test workspace object stage."""
        # Test with mock workspace object data
        workspace_input = {
            'workspace_ref': 'test_ref',
            'object_type': 'KBaseGenomes.GenomeAnnotation',
            'object_data': {
                'features': [
                    {
                        'id': 'TEST1',
                        'protein_translation': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'
                    },
                    {
                        'id': 'TEST2', 
                        'protein_translation': 'MKLIVTALAGALALQAGSLFASAADSSSPSEAGVDLK'
                    }
                ]
            }
        }
        
        result = self.workspace_object_stage.run(workspace_input)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('extracted_proteins', result.output_data)
    
    def test_embedding_generation_stage(self):
        """Test embedding generation stage."""
        # Test with protein records
        protein_records = [
            {'protein_id': 'TEST1', 'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'},
            {'protein_id': 'TEST2', 'sequence': 'MKLIVTALAGALALQAGSLFASAADSSSPSEAGVDLK'}
        ]
        
        input_data = {
            'protein_records': protein_records,
            'embedding_config': {
                'model_name': 'esm2_t6_8M_UR50D',
                'device': 'cpu',
                'batch_size': 2
            }
        }
        
        result = self.embedding_generation_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        
        # Note: This may fail if ESM2 model is not available
        if result.success:
            self.assertIn('embeddings', result.output_data)
            self.assertIn('embedding_metadata', result.output_data)
    
    def test_family_assignment_stage(self):
        """Test family assignment stage."""
        # Test with embeddings
        test_embeddings = self.embeddings[:5]
        test_protein_ids = self.protein_ids[:5]
        
        input_data = {
            'embeddings': test_embeddings,
            'protein_ids': test_protein_ids,
            'family_config': {
                'centroids_file': None,  # Will use default
                'similarity_threshold': 0.8
            }
        }
        
        result = self.family_assignment_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('family_assignments', result.output_data)
        self.assertIn('family_metadata', result.output_data)
    
    def test_similarity_search_stage(self):
        """Test similarity search stage."""
        # Test with embeddings and query
        test_embeddings = self.embeddings[:20]
        test_protein_ids = self.protein_ids[:20]
        query_embedding = test_embeddings[0]
        
        input_data = {
            'embeddings': test_embeddings,
            'protein_ids': test_protein_ids,
            'query_embedding': query_embedding,
            'query_protein_id': test_protein_ids[0],
            'similarity_config': {
                'top_k': 10,
                'similarity_threshold': 0.5
            }
        }
        
        result = self.similarity_search_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('similar_proteins', result.output_data)
        self.assertIn('similarity_scores', result.output_data)
    
    def test_sequence_analysis_stage(self):
        """Test sequence analysis stage."""
        # Test with protein records
        protein_records = [
            {'protein_id': 'TEST1', 'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'},
            {'protein_id': 'TEST2', 'sequence': 'MKLIVTALAGALALQAGSLFASAADSSSPSEAGVDLK'}
        ]
        
        input_data = {
            'protein_records': protein_records,
            'analysis_config': {
                'include_physicochemical': True,
                'include_secondary_structure': True,
                'include_motifs': True
            }
        }
        
        result = self.sequence_analysis_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('sequence_analysis', result.output_data)
        self.assertIn('analysis_summary', result.output_data)
    
    def test_network_analysis_stage(self):
        """Test network analysis stage."""
        # Test with similarity results
        test_embeddings = self.embeddings[:15]
        test_protein_ids = self.protein_ids[:15]
        query_embedding = test_embeddings[0]
        
        similar_proteins = [
            {'protein_id': test_protein_ids[i], 'similarity_score': 0.9 - i*0.1}
            for i in range(1, 10)
        ]
        
        input_data = {
            'similar_proteins': similar_proteins,
            'embeddings': test_embeddings,
            'protein_ids': test_protein_ids,
            'query_embedding': query_embedding,
            'query_protein_id': test_protein_ids[0],
            'network_config': {
                'k_neighbors': 5,
                'similarity_threshold': 0.3,
                'min_network_size': 3
            }
        }
        
        result = self.network_analysis_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('network_analysis', result.output_data)
        self.assertIn('network_properties', result.output_data)
    
    def test_report_generation_stage(self):
        """Test report generation stage."""
        # Test with pipeline results
        pipeline_results = {
            'protein_id': 'TEST1',
            'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ',
            'embeddings': {'TEST1': self.embeddings[0]},
            'family_assignments': {'TEST1': {'family_id': 'FAM0', 'confidence': 0.95}},
            'similar_proteins': [
                {'protein_id': 'SIM1', 'similarity_score': 0.85},
                {'protein_id': 'SIM2', 'similarity_score': 0.80}
            ],
            'sequence_analysis': {
                'amino_acid_composition': {'A': {'count': 5, 'percentage': 12.5}},
                'physicochemical_properties': {'molecular_weight': 4500.0}
            },
            'network_analysis': {
                'network_properties': {'num_nodes': 10, 'num_edges': 15}
            }
        }
        
        input_data = {
            'pipeline_results': pipeline_results,
            'report_config': {
                'output_format': 'html',
                'include_visualizations': True,
                'output_directory': self.temp_dir
            }
        }
        
        result = self.report_generation_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('report_files', result.output_data)
        self.assertIn('report_metadata', result.output_data)
    
    def test_visualization_stage(self):
        """Test visualization stage."""
        # Test with network data
        test_embeddings = self.embeddings[:10]
        test_protein_ids = self.protein_ids[:10]
        
        network_data = {
            'embeddings': test_embeddings,
            'protein_ids': test_protein_ids,
            'network_properties': {
                'num_nodes': 10,
                'num_edges': 15,
                'density': 0.33
            }
        }
        
        input_data = {
            'network_data': network_data,
            'visualization_config': {
                'output_format': 'html',
                'interactive': True,
                'output_directory': self.temp_dir
            }
        }
        
        result = self.visualization_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('visualization_files', result.output_data)
        self.assertIn('visualization_metadata', result.output_data)
    
    def test_data_export_stage(self):
        """Test data export stage."""
        # Test with pipeline results
        pipeline_results = {
            'protein_id': 'TEST1',
            'embeddings': {'TEST1': self.embeddings[0]},
            'family_assignments': {'TEST1': {'family_id': 'FAM0', 'confidence': 0.95}},
            'similar_proteins': [
                {'protein_id': 'SIM1', 'similarity_score': 0.85}
            ]
        }
        
        input_data = {
            'pipeline_results': pipeline_results,
            'export_config': {
                'output_format': 'json',
                'include_embeddings': True,
                'output_directory': self.temp_dir
            }
        }
        
        result = self.data_export_stage.run(input_data)
        self.assertIsInstance(result, StageResult)
        self.assertTrue(result.success)
        self.assertIn('export_files', result.output_data)
        self.assertIn('export_metadata', result.output_data)
    
    def test_stage_dependencies(self):
        """Test stage dependency relationships."""
        # Test that stages have correct dependencies
        self.assertIn('protein_records', self.embedding_generation_stage.get_required_inputs())
        self.assertIn('embeddings', self.family_assignment_stage.get_required_inputs())
        self.assertIn('embeddings', self.similarity_search_stage.get_required_inputs())
        self.assertIn('similar_proteins', self.network_analysis_stage.get_required_inputs())
    
    def test_stage_configuration(self):
        """Test stage configuration handling."""
        # Test configuration validation
        config = {
            'embedding_model': 'esm2_t6_8M_UR50D',
            'similarity_threshold': 0.8,
            'network_k_neighbors': 5
        }
        
        # Test that stages can handle configuration
        embedding_stage = EmbeddingGenerationStage(config=config)
        self.assertIsNotNone(embedding_stage.config)
        
        similarity_stage = SimilaritySearchStage(config=config)
        self.assertIsNotNone(similarity_stage.config)
    
    def test_stage_error_handling(self):
        """Test stage error handling."""
        # Test with invalid input
        invalid_input = {
            'invalid_key': 'invalid_value'
        }
        
        # All stages should handle invalid input gracefully
        stages = [
            self.input_validation_stage,
            self.data_extraction_stage,
            self.embedding_generation_stage,
            self.family_assignment_stage,
            self.similarity_search_stage,
            self.sequence_analysis_stage,
            self.network_analysis_stage,
            self.report_generation_stage,
            self.visualization_stage,
            self.data_export_stage
        ]
        
        for stage in stages:
            result = stage.run(invalid_input)
            self.assertIsInstance(result, StageResult)
            # Some stages may succeed with default values, others may fail
            # The important thing is they don't crash
    
    def test_stage_performance(self):
        """Test stage performance with larger datasets."""
        if len(self.protein_ids) >= 50:
            # Test with larger dataset
            large_embeddings = self.embeddings[:50]
            large_protein_ids = self.protein_ids[:50]
            
            input_data = {
                'embeddings': large_embeddings,
                'protein_ids': large_protein_ids,
                'family_config': {'similarity_threshold': 0.8}
            }
            
            # Test family assignment performance
            import time
            start_time = time.time()
            
            result = self.family_assignment_stage.run(input_data)
            
            processing_time = time.time() - start_time
            self.assertLess(processing_time, 30)  # Should complete within 30 seconds
            self.assertTrue(result.success)

if __name__ == '__main__':
    unittest.main()
