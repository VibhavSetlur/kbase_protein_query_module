"""
Unit tests for InputManager.

Tests input coordination, validation, and routing to appropriate processors.
"""

import pytest
import sys
import os
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.input.input_manager import InputManager, InputResult


class TestInputManager:
    """Test cases for InputManager."""
    
    def test_initialization(self, test_config, mock_kb_util):
        """Test InputManager initialization."""
        manager = InputManager(test_config, mock_kb_util)
        
        assert manager.config == test_config
        assert manager.kb_util == mock_kb_util
        assert manager.workspace_processor is not None
        assert manager.protein_sequence_processor is not None
        assert manager.uniprot_processor is not None
        assert 'protein_input' in manager.enabled_input_types
        assert 'uniprot_ids' in manager.enabled_input_types
        assert 'workspace_object' in manager.enabled_input_types
    
    def test_get_supported_input_types(self, test_config, mock_kb_util):
        """Test getting supported input types."""
        manager = InputManager(test_config, mock_kb_util)
        supported_types = manager.get_supported_input_types()
        
        assert isinstance(supported_types, list)
        assert 'protein_input' in supported_types
        assert 'uniprot_ids' in supported_types
        assert 'workspace_object' in supported_types
    
    def test_is_input_type_enabled(self, test_config, mock_kb_util):
        """Test checking if input type is enabled."""
        manager = InputManager(test_config, mock_kb_util)
        
        assert manager.is_input_type_enabled('protein_input') is True
        assert manager.is_input_type_enabled('uniprot_ids') is True
        assert manager.is_input_type_enabled('workspace_object') is True
        assert manager.is_input_type_enabled('invalid_type') is False
    
    def test_enable_input_type(self, test_config, mock_kb_util):
        """Test enabling input type."""
        # Create a modified config with only protein_input enabled
        modified_config = test_config.copy()
        modified_config['enabled_input_types'] = ['protein_input']
        manager = InputManager(modified_config, mock_kb_util)
        
        # Initially disabled
        assert 'uniprot_ids' not in manager.enabled_input_types
        
        # Enable it
        manager.enable_input_type('uniprot_ids')
        assert 'uniprot_ids' in manager.enabled_input_types
    
    def test_disable_input_type(self, test_config, mock_kb_util):
        """Test disabling input type."""
        manager = InputManager(test_config, mock_kb_util)
        
        # Initially enabled
        assert 'protein_input' in manager.enabled_input_types
        
        # Disable it
        manager.disable_input_type('protein_input')
        assert 'protein_input' not in manager.enabled_input_types
    
    def test_validate_input_data_missing_type(self, test_config, mock_kb_util):
        """Test validation with missing input type."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        }
        
        # Should fail without input_type
        assert manager.validate_input_data(input_data) is False
    
    def test_validate_input_data_invalid_type(self, test_config, mock_kb_util):
        """Test validation with invalid input type."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'invalid_type',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        }
        
        assert manager.validate_input_data(input_data) is False
    
    def test_validate_input_data_protein_input_valid(self, test_config, mock_kb_util):
        """Test validation with valid protein input."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        }
        
        assert manager.validate_input_data(input_data) is True
    
    def test_validate_input_data_protein_input_missing_data(self, test_config, mock_kb_util):
        """Test validation with missing protein input data."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'protein_input'
            # Missing protein_input field
        }
        
        assert manager.validate_input_data(input_data) is False
    
    def test_validate_input_data_uniprot_ids_valid(self, test_config, mock_kb_util):
        """Test validation with valid UniProt IDs."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'uniprot_ids',
            'uniprot_ids': ['P12345', 'P67890']
        }
        
        assert manager.validate_input_data(input_data) is True
    
    def test_validate_input_data_uniprot_ids_missing_data(self, test_config, mock_kb_util):
        """Test validation with missing UniProt IDs data."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'uniprot_ids'
            # Missing uniprot_ids field
        }
        
        assert manager.validate_input_data(input_data) is False
    
    def test_validate_input_data_workspace_object_valid(self, test_config, mock_kb_util):
        """Test validation with valid workspace object."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'workspace_object',
            'workspace_object': '1/2/3'
        }
        
        assert manager.validate_input_data(input_data) is True
    
    def test_validate_input_data_workspace_object_missing_data(self, test_config, mock_kb_util):
        """Test validation with missing workspace object data."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'workspace_object'
            # Missing workspace_object field
        }
        
        assert manager.validate_input_data(input_data) is False
    
    @patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor')
    def test_process_input_protein_input(self, mock_processor_class, test_config, mock_kb_util):
        """Test processing protein input."""
        # Setup mock
        mock_processor = Mock()
        mock_processor.process.return_value = {
            'success': True,
            'proteins': [{'protein_id': 'test_protein', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'}],
            'input_type': 'protein_input',
            'processing_time': 0.1
        }
        mock_processor_class.return_value = mock_processor
        
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        }
        
        result = manager.process_input(input_data)
        
        assert result['success'] is True
        assert result['input_type'] == 'protein_input'
        assert len(result['proteins']) == 1
        mock_processor.process.assert_called_once_with(input_data)
    
    @patch('kbase_protein_query_module.src.input.input_manager.UniProtIdsProcessor')
    def test_process_input_uniprot_ids(self, mock_processor_class, test_config, mock_kb_util):
        """Test processing UniProt IDs input."""
        # Setup mock
        mock_processor = Mock()
        mock_processor.process.return_value = {
            'success': True,
            'proteins': [{'protein_id': 'P12345', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'}],
            'input_type': 'uniprot_ids',
            'processing_time': 0.1
        }
        mock_processor_class.return_value = mock_processor
        
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'uniprot_ids',
            'uniprot_ids': ['P12345']
        }
        
        result = manager.process_input(input_data)
        
        assert result['success'] is True
        assert result['input_type'] == 'uniprot_ids'
        assert len(result['proteins']) == 1
        mock_processor.process.assert_called_once_with(input_data)
    
    @patch('kbase_protein_query_module.src.input.input_manager.WorkspaceObjectProcessor')
    def test_process_input_workspace_object(self, mock_processor_class, test_config, mock_kb_util):
        """Test processing workspace object input."""
        # Setup mock
        mock_processor = Mock()
        mock_processor.process.return_value = {
            'success': True,
            'proteins': [{'protein_id': 'P12345', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'}],
            'input_type': 'workspace_object',
            'processing_time': 0.1
        }
        mock_processor_class.return_value = mock_processor
        
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'workspace_object',
            'workspace_object': '1/2/3'
        }
        
        result = manager.process_input(input_data)
        
        assert result['success'] is True
        assert result['input_type'] == 'workspace_object'
        assert len(result['proteins']) == 1
        mock_processor.process.assert_called_once_with(input_data)
    
    def test_process_input_invalid_type(self, test_config, mock_kb_util):
        """Test processing with invalid input type."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'invalid_type',
            'data': 'some_data'
        }
        
        result = manager.process_input(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
        assert 'not enabled' in result['error_message']
    
    def test_process_input_validation_failure(self, test_config, mock_kb_util):
        """Test processing with validation failure."""
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'protein_input'
            # Missing required protein_input field
        }
        
        result = manager.process_input(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
        assert 'required for protein sequence input type' in result['error_message']
    
    @patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor')
    def test_process_input_processor_error(self, mock_processor_class, test_config, mock_kb_util):
        """Test processing when processor raises an error."""
        # Setup mock to raise exception
        mock_processor = Mock()
        mock_processor.process.side_effect = Exception("Processor error")
        mock_processor_class.return_value = mock_processor
        
        manager = InputManager(test_config, mock_kb_util)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        }
        
        result = manager.process_input(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
        assert 'Processor error' in result['error_message']
    
    def test_configuration_override(self):
        """Test initialization with custom configuration."""
        custom_config = {
            'enabled_input_types': ['protein_input'],
            'custom_setting': 'custom_value'
        }
        
        manager = InputManager(custom_config, None)
        
        assert manager.config == custom_config
        assert manager.enabled_input_types == ['protein_input']
        assert 'uniprot_ids' not in manager.enabled_input_types
        assert 'workspace_object' not in manager.enabled_input_types

