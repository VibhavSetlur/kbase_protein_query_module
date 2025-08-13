#!/usr/bin/env python3
"""
Comprehensive test suite for workflow stages

Tests all stage modules including input, processing, analysis, and output stages
for the KBase Protein Query Module workflow pipeline.
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

from kbase_protein_query_module.src.core import BaseStage, StageResult
from kbase_protein_query_module.src.stages import (
    InputValidationStage, DataExtractionStage, WorkspaceObjectStage,
    EmbeddingGenerationStage, SimilaritySearchStage, FamilyAssignmentStage,
    NetworkAnalysisStage, SequenceAnalysisStage, BioinformaticsAnalysisStage,
    ReportGenerationStage, VisualizationStage, DataExportStage
)

class TestBaseStage(unittest.TestCase):
    """Test cases for BaseStage abstract class."""
    
    def test_base_stage_abstract(self):
        """Test that BaseStage cannot be instantiated directly."""
        with self.assertRaises(TypeError):
            BaseStage()

class TestInputStages(unittest.TestCase):
    """Test cases for input stage classes."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.input_stage = InputValidationStage()
        self.data_extraction_stage = DataExtractionStage()
        self.workspace_stage = WorkspaceObjectStage()
        
        # Create test files
        self.test_fasta_file = os.path.join(self.temp_dir, 'test.fasta')
        with open(self.test_fasta_file, 'w') as f:
            f.write(">P12345|Test protein\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
        
        self.test_csv_file = os.path.join(self.temp_dir, 'test.csv')
        with open(self.test_csv_file, 'w') as f:
            f.write("protein_id,sequence\n")
            f.write("P12345,MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_input_validation_stage_initialization(self):
        """Test input validation stage initialization."""
        self.assertIsInstance(self.input_stage, BaseStage)
        self.assertEqual(self.input_stage.get_stage_name(), "input_validation")
    
    def test_data_extraction_stage_initialization(self):
        """Test data extraction stage initialization."""
        self.assertIsInstance(self.data_extraction_stage, BaseStage)
        self.assertEqual(self.data_extraction_stage.get_stage_name(), "data_extraction")
    
    def test_workspace_object_stage_initialization(self):
        """Test workspace object stage initialization."""
        self.assertIsInstance(self.workspace_stage, BaseStage)
        self.assertEqual(self.workspace_stage.get_stage_name(), "workspace_object")
    
    def test_input_validation_stage_validation(self):
        """Test input validation stage validation."""
        # Test valid input
        valid_input = {
            "input_type": "fasta",
            "input_data": {"sequences": ["MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"]}
        }
        self.assertTrue(self.input_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.input_stage.validate_input(invalid_input))
    
    def test_data_extraction_stage_validation(self):
        """Test data extraction stage validation."""
        # Create a mock validated_input object
        mock_validated_input = MagicMock()
        mock_validated_input.input_type = "fasta"
        
        # Test valid input
        valid_input = {"validated_input": mock_validated_input}
        self.assertTrue(self.data_extraction_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.data_extraction_stage.validate_input(invalid_input))
    
    def test_workspace_object_stage_validation(self):
        """Test workspace object stage validation."""
        # Create a mock validated_input object
        mock_validated_input = MagicMock()
        mock_validated_input.input_type = "workspace_object"
        mock_validated_input.input_source = "workspace_ref"
        
        # Test valid input
        valid_input = {"validated_input": mock_validated_input}
        self.assertTrue(self.workspace_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.workspace_stage.validate_input(invalid_input))

class TestProcessingStages(unittest.TestCase):
    """Test cases for processing stage classes."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.embedding_stage = EmbeddingGenerationStage()
        self.similarity_stage = SimilaritySearchStage()
        self.family_stage = FamilyAssignmentStage()
        
        self.test_proteins = [
            {"protein_id": "P12345", "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"},
            {"protein_id": "P67890", "sequence": "MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"}
        ]
    
    def test_embedding_stage_initialization(self):
        """Test embedding generation stage initialization."""
        self.assertIsInstance(self.embedding_stage, BaseStage)
        self.assertEqual(self.embedding_stage.get_stage_name(), "embedding_generation")
    
    def test_similarity_stage_initialization(self):
        """Test similarity search stage initialization."""
        self.assertIsInstance(self.similarity_stage, BaseStage)
        self.assertEqual(self.similarity_stage.get_stage_name(), "similarity_search")
    
    def test_family_stage_initialization(self):
        """Test family assignment stage initialization."""
        self.assertIsInstance(self.family_stage, BaseStage)
        self.assertEqual(self.family_stage.get_stage_name(), "family_assignment")
    
    def test_embedding_stage_validation(self):
        """Test embedding generation stage validation."""
        # Test valid input - the stage expects extracted_data to be a dict with protein_records
        valid_input = {"extracted_data": {"protein_records": self.test_proteins}}
        self.assertTrue(self.embedding_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.embedding_stage.validate_input(invalid_input))
    
    def test_similarity_stage_validation(self):
        """Test similarity search stage validation."""
        # Test valid input
        valid_input = {"embeddings": {"P12345": np.random.rand(1280).astype(np.float32)}}
        self.assertTrue(self.similarity_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.similarity_stage.validate_input(invalid_input))
    
    def test_family_stage_validation(self):
        """Test family assignment stage validation."""
        # Test valid input
        valid_input = {"embeddings": {"P12345": np.random.rand(1280).astype(np.float32)}}
        self.assertTrue(self.family_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.family_stage.validate_input(invalid_input))

class TestAnalysisStages(unittest.TestCase):
    """Test cases for analysis stage classes."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.family_assignment_stage = FamilyAssignmentStage()
        self.network_analysis_stage = NetworkAnalysisStage()
        self.sequence_analysis_stage = SequenceAnalysisStage()
        self.bioinformatics_analysis_stage = BioinformaticsAnalysisStage()
        
        self.test_embeddings = {
            "P12345": np.random.rand(1280).astype(np.float32),
            "P67890": np.random.rand(1280).astype(np.float32)
        }
        
        self.test_similarity_matrix = np.array([
            [1.0, 0.8],
            [0.8, 1.0]
        ])
    
    def test_sequence_analysis_stage_initialization(self):
        """Test sequence analysis stage initialization."""
        self.assertIsInstance(self.sequence_analysis_stage, BaseStage)
        self.assertEqual(self.sequence_analysis_stage.get_stage_name(), "sequence_analysis")
    
    def test_network_analysis_stage_initialization(self):
        """Test network analysis stage initialization."""
        self.assertIsInstance(self.network_analysis_stage, BaseStage)
        self.assertEqual(self.network_analysis_stage.get_stage_name(), "network_analysis")
    
    def test_bioinformatics_analysis_stage_initialization(self):
        """Test bioinformatics analysis stage initialization."""
        self.assertIsInstance(self.bioinformatics_analysis_stage, BaseStage)
        self.assertEqual(self.bioinformatics_analysis_stage.get_stage_name(), "bioinformatics_analysis")
    
    def test_sequence_analysis_stage_validation(self):
        """Test sequence analysis stage validation."""
        # Test valid input
        valid_input = {"extracted_data": {"proteins": [{"protein_id": "P12345", "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"}]}}
        self.assertTrue(self.sequence_analysis_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.sequence_analysis_stage.validate_input(invalid_input))
    
    def test_network_analysis_stage_validation(self):
        """Test network analysis stage validation."""
        # Test valid input
        valid_input = {
            "similarity_matrix": self.test_similarity_matrix,
            "protein_ids": ["P12345", "P67890"]
        }
        self.assertTrue(self.network_analysis_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.network_analysis_stage.validate_input(invalid_input))
    
    def test_bioinformatics_analysis_stage_validation(self):
        """Test bioinformatics analysis stage validation."""
        # Test valid input
        valid_input = {"sequence_analysis": {"proteins": [{"protein_id": "P12345", "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"}]}}
        self.assertTrue(self.bioinformatics_analysis_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.bioinformatics_analysis_stage.validate_input(invalid_input))

class TestOutputStages(unittest.TestCase):
    """Test cases for output stage classes."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.report_generation_stage = ReportGenerationStage()
        self.visualization_stage = VisualizationStage()
        self.data_export_stage = DataExportStage()
        
        self.temp_dir = tempfile.mkdtemp()
        
        self.test_results = {
            "proteins": [
                {"protein_id": "P12345", "sequence": "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"}
            ],
            "embeddings": {
                "P12345": np.random.rand(1280).astype(np.float32)
            },
            "family_assignments": {
                "P12345": {"family_id": "FAM1", "confidence": 0.95}
            },
            "network_data": {
                "nodes": [{"id": "P12345", "label": "P12345"}],
                "edges": []
            }
        }
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_report_generation_stage_initialization(self):
        """Test report generation stage initialization."""
        self.assertIsInstance(self.report_generation_stage, BaseStage)
        self.assertEqual(self.report_generation_stage.get_stage_name(), "report_generation")
    
    def test_visualization_stage_initialization(self):
        """Test visualization stage initialization."""
        self.assertIsInstance(self.visualization_stage, BaseStage)
        self.assertEqual(self.visualization_stage.get_stage_name(), "visualization")
    
    def test_data_export_stage_initialization(self):
        """Test data export stage initialization."""
        self.assertIsInstance(self.data_export_stage, BaseStage)
        self.assertEqual(self.data_export_stage.get_stage_name(), "data_export")
    
    def test_report_generation_stage_validation(self):
        """Test report generation stage validation."""
        # Test valid input
        valid_input = {"results": self.test_results}
        self.assertTrue(self.report_generation_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.report_generation_stage.validate_input(invalid_input))
    
    def test_visualization_stage_validation(self):
        """Test visualization stage validation."""
        # Test valid input
        valid_input = {"results": self.test_results}
        self.assertTrue(self.visualization_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.visualization_stage.validate_input(invalid_input))
    
    def test_data_export_stage_validation(self):
        """Test data export stage validation."""
        # Test valid input
        valid_input = {"results": self.test_results}
        self.assertTrue(self.data_export_stage.validate_input(valid_input))
        
        # Test invalid input
        invalid_input = {"invalid": "data"}
        self.assertFalse(self.data_export_stage.validate_input(invalid_input))

class TestStageIntegration(unittest.TestCase):
    """Test cases for stage integration and workflow."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create stage instances
        self.input_stage = InputValidationStage()
        self.embedding_stage = EmbeddingGenerationStage()
        self.analysis_stage = SequenceAnalysisStage()
        self.output_stage = ReportGenerationStage()
        
        # Create test data
        self.test_fasta_file = os.path.join(self.temp_dir, 'test.fasta')
        with open(self.test_fasta_file, 'w') as f:
            f.write(">P12345|Test protein\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_stage_output_schemas(self):
        """Test that all stages have valid output schemas."""
        stages = [
            self.input_stage,
            self.embedding_stage,
            self.analysis_stage,
            self.output_stage
        ]
        
        for stage in stages:
            schema = stage.get_output_schema()
            self.assertIsInstance(schema, dict)
            # Check that schema has at least one key
            self.assertGreater(len(schema), 0)
    
    def test_stage_required_inputs(self):
        """Test that all stages have valid required inputs lists."""
        stages = [
            self.input_stage,
            self.embedding_stage,
            self.analysis_stage,
            self.output_stage
        ]
        
        for stage in stages:
            required_inputs = stage.get_required_inputs()
            self.assertIsInstance(required_inputs, list)
    
    def test_stage_optional_inputs(self):
        """Test that all stages have valid optional inputs lists."""
        stages = [
            self.input_stage,
            self.embedding_stage,
            self.analysis_stage,
            self.output_stage
        ]
        
        for stage in stages:
            optional_inputs = stage.get_optional_inputs()
            self.assertIsInstance(optional_inputs, list)

if __name__ == '__main__':
    unittest.main()
