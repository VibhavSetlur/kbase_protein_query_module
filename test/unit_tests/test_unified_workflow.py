"""
Test suite for the Unified Workflow Orchestrator

This module tests the unified workflow orchestrator functionality including:
- Pipeline configuration
- Network analysis capabilities
- Storage and indexing integration
- Error handling
"""

import pytest
import numpy as np
import pandas as pd
import tempfile
import os
from unittest.mock import Mock, patch, MagicMock

from lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator import (
    UnifiedProteinQueryWorkflow,
    PipelineConfig,
    PipelineResult,
    create_basic_pipeline,
    create_full_pipeline,
    create_embedding_only_pipeline
)

class TestPipelineConfig:
    """Test PipelineConfig dataclass."""
    
    def test_pipeline_config_creation(self):
        """Test creating a basic pipeline configuration."""
        config = PipelineConfig(
            input_type="FASTA",
            input_data="test_sequence.fasta"
        )
        
        assert config.input_type == "FASTA"
        assert config.input_data == "test_sequence.fasta"
        assert config.model_name == "esm2_t6_8M_UR50D"
        assert config.top_k == 50
        assert config.network_method == "mutual_knn"
        assert config.similarity_threshold == 0.5
    
    def test_pipeline_config_with_custom_values(self):
        """Test creating a pipeline configuration with custom values."""
        config = PipelineConfig(
            input_type="ProteinSequenceSet",
            input_data="test_proteins",
            model_name="esm2_t30_150M_UR50D",
            top_k=100,
            network_method="threshold",
            similarity_threshold=0.7,
            max_concurrency=8,
            include_interactive_charts=False
        )
        
        assert config.model_name == "esm2_t30_150M_UR50D"
        assert config.top_k == 100
        assert config.network_method == "threshold"
        assert config.similarity_threshold == 0.7
        assert config.max_concurrency == 8
        assert config.include_interactive_charts == False

class TestPipelineResult:
    """Test PipelineResult dataclass."""
    
    def test_pipeline_result_creation(self):
        """Test creating a pipeline result."""
        result = PipelineResult(
            success=True,
            run_id="test_run_123",
            stages_completed=["input_validation", "data_extraction"],
            stage_results={},
            final_output={},
            execution_time=1.5
        )
        
        assert result.success == True
        assert result.run_id == "test_run_123"
        assert result.stages_completed == ["input_validation", "data_extraction"]
        assert result.execution_time == 1.5
        assert result.error_message is None
        assert result.warnings == []

class TestUnifiedWorkflowOrchestrator:
    """Test the Unified Workflow Orchestrator."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.config = PipelineConfig(
            input_type="FASTA",
            input_data="test_sequence.fasta"
        )
        
        # Mock workspace client
        self.mock_workspace = Mock()
        self.config.workspace_client = self.mock_workspace
    
    @patch('lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator.ProteinStorage')
    @patch('lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator.ProteinEmbeddingGenerator')
    @patch('lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator.assign_protein_family.AssignProteinFamily')
    @patch('lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator.DynamicNetworkBuilder')
    def test_workflow_initialization(self, mock_network_builder, mock_family_assigner, 
                                   mock_embedding_generator, mock_storage):
        """Test workflow initialization."""
        # Mock the components
        mock_storage.return_value = Mock()
        mock_embedding_generator.return_value = Mock()
        mock_family_assigner.return_value = Mock()
        mock_network_builder.return_value = Mock()
        
        workflow = UnifiedProteinQueryWorkflow(self.config)
        
        assert workflow.config == self.config
        assert workflow.run_id is not None
        assert isinstance(workflow.run_id, str)
        assert len(workflow.run_id) > 0
    
    def test_workflow_initialization_without_components(self):
        """Test workflow initialization when components fail to load."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        
        # Should still initialize even if components fail
        assert workflow.config == self.config
        assert workflow.run_id is not None
    
    @patch('lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator.ProteinEmbeddingGenerator')
    def test_generate_query_embedding(self, mock_embedding_generator):
        """Test query embedding generation."""
        # Mock embedding generator
        mock_generator = Mock()
        mock_generator.generate_embedding.return_value = np.random.rand(320)
        mock_embedding_generator.return_value = mock_generator
        
        workflow = UnifiedProteinQueryWorkflow(self.config)
        workflow.embedding_generator = mock_generator
        
        query_sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        embedding = workflow.generate_query_embedding(query_sequence, "TEST_PROTEIN")
        
        assert embedding.shape == (320,)
        mock_generator.generate_embedding.assert_called_once_with(query_sequence)
    
    def test_generate_query_embedding_no_generator(self):
        """Test query embedding generation when generator is not initialized."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        workflow.embedding_generator = None
        
        with pytest.raises(ValueError, match="Embedding generator not initialized"):
            workflow.generate_query_embedding("TEST_SEQUENCE", "TEST_PROTEIN")
    
    @patch('lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator.assign_protein_family.AssignProteinFamily')
    def test_classify_query_family(self, mock_family_assigner):
        """Test family classification."""
        # Mock family assigner
        mock_assigner = Mock()
        mock_assigner.assign_family.return_value = {
            'family_id': 'test_family',
            'confidence': 0.85
        }
        mock_family_assigner.return_value = mock_assigner
        
        workflow = UnifiedProteinQueryWorkflow(self.config)
        workflow.assign_protein_family = mock_assigner
        
        query_embedding = np.random.rand(320).astype(np.float32)
        family_id, confidence = workflow.classify_query_family(query_embedding)
        
        assert family_id == 'test_family'
        assert confidence == 0.85
        mock_assigner.assign_family.assert_called_once_with(query_embedding)
    
    def test_classify_query_family_no_assigner(self):
        """Test family classification when assigner is not initialized."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        workflow.assign_protein_family = None
        
        with pytest.raises(ValueError, match="Family assignment not initialized"):
            workflow.classify_query_family(np.random.rand(320))
    
    @patch('lib.kbase_protein_query_module.src.workflows.unified_workflow_orchestrator.ProteinStorage')
    def test_load_family_subset(self, mock_storage):
        """Test loading family subset."""
        # Mock storage
        mock_storage_instance = Mock()
        mock_storage_instance.load_family_data.return_value = (
            np.random.rand(100, 320),  # embeddings
            [f"protein_{i}" for i in range(100)],  # protein_ids
            pd.DataFrame({'family': ['test_family'] * 100})  # metadata
        )
        mock_storage.return_value = mock_storage_instance
        
        workflow = UnifiedProteinQueryWorkflow(self.config)
        workflow.storage = mock_storage_instance
        
        embeddings, protein_ids, metadata = workflow.load_family_subset("test_family")
        
        assert embeddings.shape == (100, 320)
        assert len(protein_ids) == 100
        assert len(metadata) == 100
        mock_storage_instance.load_family_data.assert_called_once_with("test_family")
    
    def test_load_family_subset_test_family_fallback(self):
        """Test loading family subset with test_family fallback."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        
        # Should raise FileNotFoundError for test_family when no real data is available
        with pytest.raises(FileNotFoundError, match="Family data not found for test_family"):
            workflow.load_family_subset("test_family")
    
    def test_perform_optimized_similarity_search(self):
        """Test optimized similarity search."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        
        # Mock the load_family_subset method
        workflow.load_family_subset = Mock(return_value=(
            np.random.rand(100, 320),  # embeddings
            [f"protein_{i}" for i in range(100)],  # protein_ids
            pd.DataFrame({'family': ['test_family'] * 100})  # metadata
        ))
        
        query_embedding = np.random.rand(320)
        similar_proteins = workflow.perform_optimized_similarity_search(
            query_embedding, "test_family", k_similar=10
        )
        
        assert len(similar_proteins) == 10
        assert all('protein_id' in protein for protein in similar_proteins)
        assert all('similarity_score' in protein for protein in similar_proteins)
        assert all('metadata' in protein for protein in similar_proteins)
    
    def test_build_optimized_network(self):
        """Test optimized network construction."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        
        # Mock the load_family_subset method
        workflow.load_family_subset = Mock(return_value=(
            np.random.rand(100, 320),  # embeddings
            [f"protein_{i}" for i in range(100)],  # protein_ids
            pd.DataFrame({'family': ['test_family'] * 100})  # metadata
        ))
        
        # Mock network builder
        workflow.network_builder = Mock()
        workflow.network_builder.create_interactive_visualization.return_value = (
            Mock(),  # Graph object
            {'num_nodes': 101, 'num_edges': 100}  # network properties
        )
        
        query_embedding = np.random.rand(320)
        similar_proteins = [
            {'protein_id': f'protein_{i}', 'similarity_score': 0.8 + i*0.01}
            for i in range(10)
        ]
        
        G, network_properties = workflow.build_optimized_network(
            query_embedding, "TEST_PROTEIN", similar_proteins, "test_family"
        )
        
        assert network_properties['num_nodes'] == 101
        assert network_properties['num_edges'] == 100
        workflow.network_builder.create_interactive_visualization.assert_called_once()
    
    def test_run_network_analysis_workflow(self):
        """Test complete network analysis workflow."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        
        # Mock all the methods
        workflow.generate_query_embedding = Mock(return_value=np.random.rand(320))
        workflow.classify_query_family = Mock(return_value=("test_family", 0.85))
        workflow.load_family_subset = Mock(return_value=(
            np.random.rand(100, 320),
            [f"protein_{i}" for i in range(100)],
            pd.DataFrame({'family': ['test_family'] * 100})
        ))
        workflow.perform_optimized_similarity_search = Mock(return_value=[
            {'protein_id': f'protein_{i}', 'similarity_score': 0.8 + i*0.01}
            for i in range(10)
        ])
        workflow.build_optimized_network = Mock(return_value=(
            Mock(),  # Graph
            {'num_nodes': 101, 'num_edges': 100}  # properties
        ))
        workflow._save_workflow_results = Mock()
        
        # Mock network builder for HTML generation
        workflow.network_builder = Mock()
        workflow.network_builder.create_interactive_visualization.return_value = (Mock(), Mock())
        
        results = workflow.run_network_analysis_workflow(
            "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG",
            "TEST_PROTEIN",
            k_similar=10,
            network_method="mutual_knn"
        )
        
        assert results['status'] == 'success'
        assert results['query_protein_id'] == 'TEST_PROTEIN'
        assert 'workflow_steps' in results
        assert 'timing' in results
        assert 'performance_metrics' in results
        assert len(results['workflow_steps']) == 6  # 6 steps in the workflow
    
    def test_get_system_info(self):
        """Test getting system information."""
        workflow = UnifiedProteinQueryWorkflow(self.config)
        
        # Mock storage and family stats
        workflow.storage = Mock()
        workflow.available_families = ['family1', 'family2']
        workflow.family_stats = {
            'family1': {'num_proteins': 1000},
            'family2': {'num_proteins': 2000}
        }
        
        system_info = workflow.get_system_info()
        
        assert 'storage_available' in system_info
        assert 'available_families' in system_info
        assert 'total_proteins' in system_info
        assert 'memory_usage_mb' in system_info
        assert 'model_name' in system_info
        assert system_info['available_families'] == 2
        assert system_info['total_proteins'] == 3000

class TestConvenienceFunctions:
    """Test convenience functions for creating pipelines."""
    
    def test_create_basic_pipeline(self):
        """Test creating a basic pipeline."""
        pipeline = create_basic_pipeline("FASTA", "test.fasta")
        
        assert isinstance(pipeline, UnifiedProteinQueryWorkflow)
        assert pipeline.config.input_type == "FASTA"
        assert pipeline.config.input_data == "test.fasta"
        assert pipeline.config.enabled_stages == ['input_validation', 'data_extraction', 'embedding_generation']
    
    def test_create_full_pipeline(self):
        """Test creating a full pipeline."""
        pipeline = create_full_pipeline("FASTA", "test.fasta")
        
        assert isinstance(pipeline, UnifiedProteinQueryWorkflow)
        assert pipeline.config.input_type == "FASTA"
        assert pipeline.config.input_data == "test.fasta"
        assert 'network_analysis' in pipeline.config.enabled_stages
        assert 'report_generation' in pipeline.config.enabled_stages
    
    def test_create_embedding_only_pipeline(self):
        """Test creating an embedding-only pipeline."""
        pipeline = create_embedding_only_pipeline("FASTA", "test.fasta")
        
        assert isinstance(pipeline, UnifiedProteinQueryWorkflow)
        assert pipeline.config.input_type == "FASTA"
        assert pipeline.config.input_data == "test.fasta"
        assert pipeline.config.stop_after_stage == 'embedding_generation'
        assert pipeline.config.enabled_stages == ['input_validation', 'data_extraction', 'embedding_generation']

class TestErrorHandling:
    """Test error handling in the workflow orchestrator."""
    
    def test_workflow_initialization_with_invalid_config(self):
        """Test workflow initialization with invalid configuration."""
        # Should not raise an exception, just log warnings
        config = PipelineConfig(
            input_type="",  # Invalid empty input type
            input_data=""   # Invalid empty input data
        )
        
        workflow = UnifiedProteinQueryWorkflow(config)
        assert workflow.config == config
    
    def test_network_analysis_workflow_error_handling(self):
        """Test error handling in network analysis workflow."""
        workflow = UnifiedProteinQueryWorkflow(PipelineConfig("FASTA", "test.fasta"))
        
        # Mock generate_query_embedding to raise an exception
        workflow.generate_query_embedding = Mock(side_effect=Exception("Test error"))
        
        results = workflow.run_network_analysis_workflow(
            "TEST_SEQUENCE",
            "TEST_PROTEIN"
        )
        
        assert results['status'] == 'error'
        assert 'Test error' in results['error']
        assert 'timing' in results
        assert 'total' in results['timing']

if __name__ == "__main__":
    pytest.main([__file__])
