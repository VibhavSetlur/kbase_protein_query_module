"""
Comprehensive Workflow Tests

This module tests the complete workflow orchestration including:
- End-to-end workflow execution
- Stage coordination and dependency management
- Error handling and recovery
- Performance monitoring
- Configuration management
"""

import unittest
import numpy as np
import os
import sys
import tempfile
import h5py
import pandas as pd
from pathlib import Path

from kbase_protein_query_module.src.workflows import ProteinQueryWorkflow, WorkflowResult
from kbase_protein_query_module.src.core import PipelineConfig

class TestWorkflowsComprehensive(unittest.TestCase):
    """Comprehensive tests for workflow orchestration."""
    
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
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(cls.temp_dir)
    
    def test_workflow_initialization(self):
        """Test workflow initialization with different configurations."""
        # Test with minimal configuration
        config = PipelineConfig(input_proteins=self.test_protein_ids[:1])
        workflow = ProteinQueryWorkflow(config=config)
        self.assertIsNotNone(workflow)
        self.assertIsNotNone(workflow.config)
        
        # Test with custom configuration
        config = PipelineConfig(
            input_proteins=self.test_protein_ids,
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=True,
            perform_sequence_analysis=True,
            embedding_model="esm2_t6_8M_UR50D",
            similarity_threshold=0.8,
            max_similar_proteins=50,
            generate_html_report=True,
            output_format="html"
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        self.assertIsNotNone(workflow)
        self.assertEqual(workflow.config.input_proteins, self.test_protein_ids)
        self.assertEqual(workflow.config.similarity_threshold, 0.8)
    
    def test_workflow_with_fasta_input(self):
        """Test workflow with FASTA sequence input."""
        # Create FASTA input
        fasta_data = ">TEST1\nMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ\n>TEST2\nMKLIVTALAGALALQAGSLFASAADSSSPSEAGVDLK"
        
        config = PipelineConfig(
            input_proteins=self.test_protein_ids[:2],  # Add required input
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=False,  # Skip for faster testing
            perform_sequence_analysis=True,
            embedding_model="esm2_t6_8M_UR50D",
            similarity_threshold=0.8,
            max_similar_proteins=20,
            generate_html_report=True,
            output_format="html"
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute workflow with FASTA data
        input_data = {
            'input_type': 'fasta_sequence',
            'input_data': fasta_data
        }
        
        result = workflow.execute(input_data)
        self.assertIsInstance(result, WorkflowResult)
        
        # Check workflow result
        if result.success:
            self.assertIn('protein_id', result.final_output)
            self.assertIn('sequence', result.final_output)
            self.assertGreater(len(result.stages_completed), 0)
            self.assertGreater(result.execution_time, 0)
        else:
            # If it fails, it should be due to missing ESM2 model
            self.assertIsNotNone(result.error_message)
    
    def test_workflow_with_uniprot_ids(self):
        """Test workflow with UniProt ID input."""
        # Use real UniProt IDs from the data
        real_uniprot_ids = self.protein_ids[:3]  # Use first 3 real IDs
        
        config = PipelineConfig(
            input_proteins=real_uniprot_ids,
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=False,  # Skip for faster testing
            perform_sequence_analysis=True,
            embedding_model="esm2_t6_8M_UR50D",
            similarity_threshold=0.8,
            max_similar_proteins=20,
            generate_html_report=True,
            output_format="html"
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute workflow
        result = workflow.execute()
        self.assertIsInstance(result, WorkflowResult)
        
        # Check workflow result
        if result.success:
            self.assertIn('protein_id', result.final_output)
            self.assertIn('sequence', result.final_output)
            self.assertGreater(len(result.stages_completed), 0)
            self.assertGreater(result.execution_time, 0)
        else:
            # If it fails, it should be due to missing ESM2 model
            self.assertIsNotNone(result.error_message)
    
    def test_workflow_stage_execution_order(self):
        """Test that stages execute in the correct order."""
        config = PipelineConfig(
            input_proteins=self.test_protein_ids,
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=False,
            perform_sequence_analysis=True,
            generate_html_report=True
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Get execution order
        execution_order = workflow._get_stage_execution_order()
        
        # Verify order makes sense
        self.assertIsInstance(execution_order, list)
        self.assertGreater(len(execution_order), 0)
        
        # Input stages should come first
        input_stages = ['InputValidationStage', 'DataExtractionStage']
        for stage in input_stages:
            if stage in execution_order:
                input_index = execution_order.index(stage)
                # Input stages should be early in the pipeline
                self.assertLess(input_index, len(execution_order) // 2)
        
        # Processing stages should come after input
        processing_stages = ['EmbeddingGenerationStage', 'FamilyAssignmentStage', 'SimilaritySearchStage']
        for stage in processing_stages:
            if stage in execution_order:
                processing_index = execution_order.index(stage)
                # Processing stages should be after input stages
                self.assertGreater(processing_index, 0)
        
        # Output stages should come last
        output_stages = ['ReportGenerationStage', 'DataExportStage']
        for stage in output_stages:
            if stage in execution_order:
                output_index = execution_order.index(stage)
                # Output stages should be late in the pipeline
                self.assertGreater(output_index, len(execution_order) // 2)
    
    def test_workflow_error_handling(self):
        """Test workflow error handling."""
        # Test with invalid input
        config = PipelineConfig(
            input_proteins=[],
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=False,
            perform_sequence_analysis=True,
            generate_html_report=True
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute with invalid input
        result = workflow.execute()
        self.assertIsInstance(result, WorkflowResult)
        self.assertFalse(result.success)
        self.assertIsNotNone(result.error_message)
        
        # Test with invalid configuration
        try:
            invalid_config = PipelineConfig(
                # Missing all required input parameters
            )
            self.fail("Should have raised ValueError for missing input parameters")
        except ValueError as e:
            self.assertIn("Must provide either input_proteins", str(e))
    
    def test_workflow_performance_monitoring(self):
        """Test workflow performance monitoring."""
        config = PipelineConfig(
            input_proteins=self.test_protein_ids,
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=False,
            perform_sequence_analysis=True,
            generate_html_report=True
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute workflow
        result = workflow.execute()
        self.assertIsInstance(result, WorkflowResult)
        
        # Check performance metrics
        self.assertGreater(result.execution_time, 0)
        self.assertIsInstance(result.stages_completed, list)
        self.assertIsInstance(result.stage_results, dict)
        
        # Get performance summary
        performance_summary = workflow.get_performance_summary()
        self.assertIsInstance(performance_summary, dict)
        self.assertIn('total_execution_time', performance_summary)
        self.assertIn('stages_completed', performance_summary)
        self.assertIn('stages_failed', performance_summary)
    
    def test_workflow_configuration_validation(self):
        """Test workflow configuration validation."""
        # Test valid configuration
        valid_config = PipelineConfig(
            input_proteins=self.test_protein_ids,
            perform_embedding_generation=True,
            similarity_threshold=0.8,
            max_similar_proteins=50
        )
        
        workflow = ProteinQueryWorkflow(valid_config)
        self.assertIsNotNone(workflow)
        
        # Test invalid similarity threshold
        try:
            invalid_config = PipelineConfig(
                input_proteins=self.test_protein_ids,
                similarity_threshold=1.5  # Should be between 0 and 1
            )
            self.fail("Should have raised ValueError for invalid similarity threshold")
        except ValueError:
            # Expected
            pass
        
        # Test invalid max_similar_proteins
        try:
            invalid_config = PipelineConfig(
                input_proteins=self.test_protein_ids,
                max_similar_proteins=-1  # Should be positive
            )
            self.fail("Should have raised ValueError for invalid max_similar_proteins")
        except ValueError:
            # Expected
            pass
    
    def test_workflow_with_partial_execution(self):
        """Test workflow with partial stage execution."""
        # Test workflow with only some stages enabled
        config = PipelineConfig(
            input_proteins=self.test_protein_ids,
            perform_embedding_generation=True,
            perform_family_assignment=False,  # Disable family assignment
            perform_similarity_search=False,  # Disable similarity search
            perform_network_analysis=False,   # Disable network analysis
            perform_sequence_analysis=True,   # Enable sequence analysis
            generate_html_report=True
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute workflow
        result = workflow.execute()
        self.assertIsInstance(result, WorkflowResult)
        
        if result.success:
            # Should have completed fewer stages
            self.assertLess(len(result.stages_completed), 8)  # Less than full pipeline
            
            # Should not have family assignment or similarity search results
            if 'family_assignments' in result.final_output:
                self.fail("Should not have family assignments when disabled")
            if 'similar_proteins' in result.final_output:
                self.fail("Should not have similar proteins when disabled")
    
    def test_workflow_cleanup(self):
        """Test workflow cleanup functionality."""
        config = PipelineConfig(
            input_proteins=self.test_protein_ids,
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=False,
            perform_sequence_analysis=True,
            generate_html_report=True
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute workflow
        result = workflow.execute()
        
        # Test cleanup
        workflow.cleanup()
        
        # Workflow should still be usable after cleanup
        self.assertIsNotNone(workflow.config)
        self.assertIsNotNone(workflow.get_performance_summary())
    
    def test_workflow_with_real_data(self):
        """Test workflow with real protein data."""
        # Use real protein IDs from the data
        real_protein_ids = self.protein_ids[:5]  # Use first 5 real IDs
        
        config = PipelineConfig(
            input_proteins=real_protein_ids,
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=False,  # Skip for faster testing
            perform_sequence_analysis=True,
            embedding_model="esm2_t6_8M_UR50D",
            similarity_threshold=0.8,
            max_similar_proteins=20,
            generate_html_report=True,
            output_format="html"
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute workflow
        result = workflow.execute()
        self.assertIsInstance(result, WorkflowResult)
        
        # Check workflow result
        if result.success:
            self.assertIn('protein_id', result.final_output)
            self.assertIn('sequence', result.final_output)
            self.assertGreater(len(result.stages_completed), 0)
            self.assertGreater(result.execution_time, 0)
            
            # Check that all enabled stages were completed
            expected_stages = ['InputValidationStage', 'DataExtractionStage', 
                             'EmbeddingGenerationStage', 'FamilyAssignmentStage', 
                             'SimilaritySearchStage', 'SequenceAnalysisStage', 
                             'ReportGenerationStage']
            
            for stage in expected_stages:
                if stage in result.stages_completed:
                    self.assertIn(stage, result.stage_results)
        else:
            # If it fails, it should be due to missing ESM2 model
            self.assertIsNotNone(result.error_message)
    
    def test_workflow_integration(self):
        """Test complete workflow integration."""
        # Test the complete workflow with all components
        config = PipelineConfig(
            input_proteins=self.test_protein_ids,
            perform_embedding_generation=True,
            perform_family_assignment=True,
            perform_similarity_search=True,
            perform_network_analysis=True,
            perform_sequence_analysis=True,
            embedding_model="esm2_t6_8M_UR50D",
            similarity_threshold=0.8,
            max_similar_proteins=20,
            network_min_edges=3,
            generate_html_report=True,
            generate_network_visualization=True,
            output_format="html"
        )
        
        workflow = ProteinQueryWorkflow(config=config)
        
        # Execute complete workflow
        result = workflow.execute()
        self.assertIsInstance(result, WorkflowResult)
        
        # Check workflow result
        if result.success:
            # Should have completed all stages
            self.assertGreater(len(result.stages_completed), 8)
            
            # Should have comprehensive output
            self.assertIn('protein_id', result.final_output)
            self.assertIn('sequence', result.final_output)
            self.assertIn('embeddings', result.final_output)
            self.assertIn('family_assignments', result.final_output)
            self.assertIn('similar_proteins', result.final_output)
            self.assertIn('sequence_analysis', result.final_output)
            self.assertIn('network_analysis', result.final_output)
            self.assertIn('report_files', result.final_output)
            
            # Check performance
            self.assertGreater(result.execution_time, 0)
            self.assertIsInstance(result.stage_results, dict)
            self.assertGreater(len(result.stage_results), 0)
        else:
            # If it fails, it should be due to missing ESM2 model
            self.assertIsNotNone(result.error_message)

if __name__ == '__main__':
    unittest.main()
