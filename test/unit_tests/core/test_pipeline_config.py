"""
Unit tests for PipelineConfig.

Tests configuration validation, merging, and serialization.
"""

import pytest
import sys
import os
from unittest.mock import patch, Mock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.core.pipeline_config import PipelineConfig


class TestPipelineConfig:
    """Test cases for PipelineConfig."""
    
    def test_default_initialization(self):
        """Test default initialization."""
        config = PipelineConfig()
        
        assert config.input_proteins == []
        assert config.input_type == "protein_input"
        assert config.workspace_object_ref is None
        assert config.enabled_stages == ["input_processing", "analysis", "output_generation"]
        assert config.output_dir == "/tmp"
        assert config.workspace_name is None
        assert config.generate_html_report is True
        assert config.generate_network_visualization is True
        assert config.batch_size == 100
        assert config.max_memory_gb == 4.0
        assert config.max_concurrent_tasks == 2
    
    def test_custom_initialization(self):
        """Test initialization with custom values."""
        config = PipelineConfig(
            input_proteins=['P12345', 'P67890'],
            input_type='uniprot_ids',
            workspace_object_ref='1/2/3',
            output_dir='/custom/output',
            workspace_name='test_workspace',
            batch_size=50
        )
        
        assert config.input_proteins == ['P12345', 'P67890']
        assert config.input_type == 'uniprot_ids'
        assert config.workspace_object_ref == '1/2/3'
        assert config.output_dir == '/custom/output'
        assert config.workspace_name == 'test_workspace'
        assert config.batch_size == 50
    
    def test_post_init_default_stages(self):
        """Test that default stages are set correctly."""
        config = PipelineConfig()
        
        assert "input_processing" in config.enabled_stages
        assert "analysis" in config.enabled_stages
        assert "output_generation" in config.enabled_stages
    
    def test_post_init_default_stage_configs(self):
        """Test that default stage configs are set."""
        config = PipelineConfig()
        
        assert "input_processing" in config.stage_configs
        assert "analysis" in config.stage_configs
        assert "output_generation" in config.stage_configs
        
        # Check specific stage configs
        input_config = config.stage_configs["input_processing"]
        assert input_config["validate_input"] is True
        assert input_config["extract_data"] is True
        
        analysis_config = config.stage_configs["analysis"]
        assert analysis_config["run_enabled_analyses"] is True
        
        output_config = config.stage_configs["output_generation"]
        assert output_config["create_reports"] is True
        assert output_config["zip_outputs"] is True
    
    def test_post_init_default_storage_config(self):
        """Test that default storage config is set."""
        config = PipelineConfig()
        
        storage_config = config.storage_config
        assert storage_config["use_compression"] is True
        assert storage_config["chunk_size"] == 1000
        assert storage_config["max_family_size"] == 100000
    
    def test_post_init_default_similarity_config(self):
        """Test that default similarity config is set."""
        config = PipelineConfig()
        
        similarity_config = config.similarity_config
        assert similarity_config["similarity_threshold"] == 0.7
        assert similarity_config["top_k_matches"] == 50
        assert similarity_config["use_family_assignment"] is True
    
    def test_get_stage_config_existing(self):
        """Test getting existing stage config."""
        config = PipelineConfig()
        
        stage_config = config.get_stage_config("input_processing")
        
        assert isinstance(stage_config, dict)
        assert stage_config["validate_input"] is True
    
    def test_get_stage_config_nonexistent(self):
        """Test getting non-existent stage config."""
        config = PipelineConfig()
        
        stage_config = config.get_stage_config("nonexistent_stage")
        
        assert stage_config == {}
    
    def test_is_stage_enabled_true(self):
        """Test checking enabled stage."""
        config = PipelineConfig()
        
        assert config.is_stage_enabled("input_processing") is True
        assert config.is_stage_enabled("analysis") is True
    
    def test_is_stage_enabled_false(self):
        """Test checking disabled stage."""
        config = PipelineConfig()
        
        assert config.is_stage_enabled("nonexistent_stage") is False
    
    def test_enable_stage_new(self):
        """Test enabling new stage."""
        config = PipelineConfig()
        
        config.enable_stage("new_stage")
        
        assert "new_stage" in config.enabled_stages
    
    def test_enable_stage_existing(self):
        """Test enabling existing stage."""
        config = PipelineConfig()
        
        original_length = len(config.enabled_stages)
        config.enable_stage("input_processing")  # Already enabled
        
        assert len(config.enabled_stages) == original_length
        assert "input_processing" in config.enabled_stages
    
    def test_disable_stage_existing(self):
        """Test disabling existing stage."""
        config = PipelineConfig()
        
        config.disable_stage("input_processing")
        
        assert "input_processing" not in config.enabled_stages
    
    def test_disable_stage_nonexistent(self):
        """Test disabling non-existent stage."""
        config = PipelineConfig()
        
        original_length = len(config.enabled_stages)
        config.disable_stage("nonexistent_stage")
        
        assert len(config.enabled_stages) == original_length
    
    def test_set_stage_config(self):
        """Test setting stage config."""
        config = PipelineConfig()
        
        new_config = {"custom_param": "custom_value"}
        config.set_stage_config("test_stage", new_config)
        
        assert config.get_stage_config("test_stage") == new_config
    
    def test_validate_valid_config(self):
        """Test validating valid configuration."""
        config = PipelineConfig(
            input_proteins=['P12345'],
            input_type='protein_input'
        )
        
        assert config.validate() is True
    
    def test_validate_with_workspace_object(self):
        """Test validating configuration with workspace object."""
        config = PipelineConfig(
            input_type='workspace_object',
            workspace_object_ref='1/2/3'
        )
        
        assert config.validate() is True
    
    def test_validate_missing_input_proteins(self):
        """Test validating configuration without input proteins."""
        config = PipelineConfig(
            input_type='protein_input'
            # Missing input_proteins
        )
        
        assert config.validate() is False
    
    def test_validate_missing_workspace_object_ref(self):
        """Test validating configuration with workspace_object input but no ref."""
        config = PipelineConfig(
            input_type='workspace_object'
            # Missing workspace_object_ref
        )
        
        assert config.validate() is False
    
    def test_validate_invalid_input_type(self):
        """Test validating configuration with invalid input type."""
        config = PipelineConfig(
            input_type='invalid_type',
            input_proteins=['P12345']
        )
        
        assert config.validate() is False
    
    def test_validate_invalid_output_directory(self):
        """Test validating configuration with invalid output directory."""
        config = PipelineConfig(
            input_proteins=['P12345'],
            output_dir='/invalid/path/that/does/not/exist'
        )
        
        # Should return False as it cannot create the directory
        assert config.validate() is False
    
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_validate_output_directory_creation(self, mock_makedirs, mock_exists):
        """Test that output directory is created if it doesn't exist."""
        mock_exists.return_value = False
        
        config = PipelineConfig(
            input_proteins=['P12345'],
            output_dir='/tmp/test_output'
        )
        
        result = config.validate()
        
        assert result is True
        mock_makedirs.assert_called_once_with('/tmp', exist_ok=True)
    
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_validate_output_directory_creation_error(self, mock_makedirs, mock_exists):
        """Test handling of output directory creation error."""
        mock_exists.return_value = False
        mock_makedirs.side_effect = Exception("Permission denied")
        
        config = PipelineConfig(
            input_proteins=['P12345'],
            output_dir='/invalid/path'
        )
        
        assert config.validate() is False
    
    def test_validate_exception_handling(self):
        """Test exception handling during validation."""
        config = PipelineConfig()
        
        # Mock a method to raise exception
        with patch.object(config, 'input_type', side_effect=Exception("Validation error")):
            assert config.validate() is False
    
    def test_to_dict(self):
        """Test converting configuration to dictionary."""
        config = PipelineConfig(
            input_proteins=['P12345'],
            input_type='protein_input',
            output_dir='/tmp/test'
        )
        
        config_dict = config.to_dict()
        
        assert isinstance(config_dict, dict)
        assert config_dict['input_proteins'] == ['P12345']
        assert config_dict['input_type'] == 'protein_input'
        assert config_dict['output_dir'] == '/tmp/test'
        assert 'enabled_stages' in config_dict
        assert 'stage_configs' in config_dict
        assert 'storage_config' in config_dict
        assert 'similarity_config' in config_dict
    
    def test_from_dict(self):
        """Test creating configuration from dictionary."""
        config_dict = {
            'input_proteins': ['P12345', 'P67890'],
            'input_type': 'uniprot_ids',
            'output_dir': '/tmp/test',
            'batch_size': 50
        }
        
        config = PipelineConfig.from_dict(config_dict)
        
        assert config.input_proteins == ['P12345', 'P67890']
        assert config.input_type == 'uniprot_ids'
        assert config.output_dir == '/tmp/test'
        assert config.batch_size == 50
        # Should still have defaults for other fields
        assert config.enabled_stages == ["input_processing", "analysis", "output_generation"]
    
    def test_merge_configs(self):
        """Test merging two configurations."""
        config1 = PipelineConfig(
            input_proteins=['P12345'],
            input_type='protein_input',
            batch_size=100
        )
        
        config2 = PipelineConfig(
            input_proteins=['P67890'],
            output_dir='/tmp/merged',
            batch_size=50
        )
        
        merged = config1.merge(config2)
        
        # config2 should take precedence
        assert merged.input_proteins == ['P67890']
        assert merged.output_dir == '/tmp/merged'
        assert merged.batch_size == 50
        # Should inherit from config1 where config2 has None
        assert merged.input_type == 'protein_input'
    
    def test_merge_configs_with_none_values(self):
        """Test merging configurations where one has None values."""
        config1 = PipelineConfig(
            input_proteins=['P12345'],
            workspace_name='workspace1'
        )
        
        config2 = PipelineConfig(
            input_proteins=None,  # None value
            workspace_name='workspace2'
        )
        
        merged = config1.merge(config2)
        
        # None values should not override
        assert merged.input_proteins == ['P12345']
        assert merged.workspace_name == 'workspace2'
    
    def test_configuration_immutability_after_merge(self):
        """Test that original configurations are not modified after merge."""
        config1 = PipelineConfig(input_proteins=['P12345'])
        config2 = PipelineConfig(output_dir='/tmp/test')
        
        original_input_proteins = config1.input_proteins.copy()
        original_output_dir = config2.output_dir
        
        merged = config1.merge(config2)
        
        # Original configs should be unchanged
        assert config1.input_proteins == original_input_proteins
        assert config2.output_dir == original_output_dir
    
    def test_comprehensive_configuration(self):
        """Test comprehensive configuration with all fields."""
        config = PipelineConfig(
            input_proteins=['P12345', 'P67890'],
            input_type='uniprot_ids',
            workspace_object_ref=None,
            enabled_stages=['custom_stage'],
            stage_configs={'custom_stage': {'param': 'value'}},
            storage_config={'custom_storage': 'value'},
            similarity_config={'custom_similarity': 'value'},
            output_dir='/comprehensive/output',
            workspace_name='comprehensive_workspace',
            workspace_client=Mock(),
            workspace_url='https://test.kbase.us/services/ws',
            auth_token='test_token',
            generate_html_report=False,
            generate_network_visualization=False,
            selected_analyses=['analysis1', 'analysis2'],
            analysis_config={'custom_analysis': 'value'},
            batch_size=200,
            max_memory_gb=8.0,
            max_concurrent_tasks=4
        )
        
        # Test all fields
        assert config.input_proteins == ['P12345', 'P67890']
        assert config.input_type == 'uniprot_ids'
        assert config.enabled_stages == ['custom_stage']
        assert config.output_dir == '/comprehensive/output'
        assert config.workspace_name == 'comprehensive_workspace'
        assert config.generate_html_report is False
        assert config.batch_size == 200
        assert config.max_memory_gb == 8.0
        assert config.max_concurrent_tasks == 4
        
        # Test validation - mock directory creation to avoid permission issues
        with patch('os.path.exists', return_value=False), \
             patch('os.makedirs', return_value=None):
            assert config.validate() is True

