"""
Unit tests for WorkspaceObjectProcessor.

Tests KBase workspace object retrieval, validation, and processing.
"""

import pytest
import sys
import os
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.input.workspace_object.input import WorkspaceObjectProcessor, WorkspaceData


class TestWorkspaceObjectProcessor:
    """Test cases for WorkspaceObjectProcessor."""
    
    def test_initialization(self, test_config, mock_kb_util):
        """Test WorkspaceObjectProcessor initialization."""
        processor = WorkspaceObjectProcessor(test_config, mock_kb_util)
        
        assert processor.config == test_config
        assert processor.kb_util == mock_kb_util
    
    def test_initialization_without_kb_util(self, test_config):
        """Test initialization without KBase utility."""
        processor = WorkspaceObjectProcessor(test_config)
        
        assert processor.config == test_config
        assert processor.kb_util is None
    
    def test_process_valid_workspace_object(self, test_config, mock_workspace_client, sample_workspace_object_ref):
        """Test processing valid workspace object."""
        processor = WorkspaceObjectProcessor(test_config, kb_util=None)
        
        # Mock the workspace client methods
        with patch.object(processor, '_get_workspace_client', return_value=mock_workspace_client):
            input_data = {
                'input_type': 'workspace_object',
                'workspace_object': sample_workspace_object_ref
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert result['input_type'] == 'workspace_object'
            assert len(result['proteins']) > 0
            assert all('protein_id' in protein for protein in result['proteins'])
            assert all('sequence' in protein for protein in result['proteins'])
            assert all(protein['source'] == 'workspace_object' for protein in result['proteins'])
    
    def test_process_invalid_workspace_ref(self, test_config):
        """Test processing with invalid workspace reference."""
        processor = WorkspaceObjectProcessor(test_config)
        
        input_data = {
            'input_type': 'workspace_object',
            'workspace_object': 'invalid_ref'
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
        assert 'Invalid workspace reference' in result['error_message']
    
    def test_process_missing_workspace_object(self, test_config):
        """Test processing with missing workspace object field."""
        processor = WorkspaceObjectProcessor(test_config)
        
        input_data = {
            'input_type': 'workspace_object'
            # Missing workspace_object field
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
        assert 'workspace_object' in result['error_message']
    
    def test_process_workspace_client_error(self, test_config):
        """Test processing when workspace client raises error."""
        processor = WorkspaceObjectProcessor(test_config)
        
        # Mock workspace client to raise exception
        mock_client = Mock()
        mock_client.get_objects2.side_effect = Exception("Workspace error")
        
        with patch.object(processor, '_get_workspace_client', return_value=mock_client):
            input_data = {
                'input_type': 'workspace_object',
                'workspace_object': '1/2/3'
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is False
            assert 'error_message' in result
            assert 'Workspace error' in result['error_message']
    
    def test_process_object_not_found(self, test_config):
        """Test processing when workspace object is not found."""
        processor = WorkspaceObjectProcessor(test_config)
        
        # Mock workspace client to return empty data
        mock_client = Mock()
        mock_client.get_objects2.return_value = {'data': []}
        
        with patch.object(processor, '_get_workspace_client', return_value=mock_client):
            input_data = {
                'input_type': 'workspace_object',
                'workspace_object': '1/2/3'
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is False
            assert 'error_message' in result
            assert 'Object not found' in result['error_message']
    
    def test_is_valid_workspace_ref_valid(self, test_config):
        """Test workspace reference validation with valid references."""
        processor = WorkspaceObjectProcessor(test_config)
        
        valid_refs = ['1/2/3', '123/456/789', 'workspace/object/version']
        
        for ref in valid_refs:
            assert processor._is_valid_workspace_ref(ref) is True
    
    def test_is_valid_workspace_ref_invalid(self, test_config):
        """Test workspace reference validation with invalid references."""
        processor = WorkspaceObjectProcessor(test_config)
        
        invalid_refs = ['invalid', '1/2', '1/2/3/4', '', None]
        
        for ref in invalid_refs:
            assert processor._is_valid_workspace_ref(ref) is False
    
    def test_get_workspace_client_with_kb_util(self, test_config, mock_kb_util):
        """Test getting workspace client with KBase utility."""
        processor = WorkspaceObjectProcessor(test_config, mock_kb_util)
        
        # Mock the KBase utility workspace client creation
        mock_client = Mock()
        mock_kb_util.get_client.return_value = mock_client
        
        client = processor._get_workspace_client()
        
        assert client is not None
        mock_kb_util.get_client.assert_called_once()
    
    def test_get_workspace_client_without_kb_util(self, test_config):
        """Test getting workspace client without KBase utility."""
        processor = WorkspaceObjectProcessor(test_config)
        
        # Should raise an error when no kb_util is available
        with pytest.raises(Exception):
            processor._get_workspace_client()
    
    def test_retrieve_object_data_success(self, test_config, mock_workspace_client):
        """Test successful object data retrieval."""
        processor = WorkspaceObjectProcessor(test_config)
        
        object_data = processor._retrieve_object_data(mock_workspace_client, '1/2/3')
        
        assert object_data is not None
        mock_workspace_client.get_objects2.assert_called_once()
    
    def test_retrieve_object_data_error(self, test_config):
        """Test object data retrieval with error."""
        processor = WorkspaceObjectProcessor(test_config)
        
        mock_client = Mock()
        mock_client.get_objects2.side_effect = Exception("Retrieval error")
        
        with pytest.raises(Exception, match="Retrieval error"):
            processor._retrieve_object_data(mock_client, '1/2/3')
    
    def test_parse_workspace_object_genome_data(self, test_config):
        """Test parsing workspace object with genome data."""
        processor = WorkspaceObjectProcessor(test_config)
        
        object_data = {
            'data': {
                'proteins': [
                    {
                        'id': 'P12345',
                        'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
                    },
                    {
                        'id': 'P67890',
                        'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
                    }
                ]
            },
            'info': [1, 'test_object', 'KBaseGenomes.Genome-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]
        }
        
        proteins = processor._parse_workspace_object(object_data)
        
        assert len(proteins) == 2
        assert proteins[0]['protein_id'] == 'P12345'
        assert proteins[1]['protein_id'] == 'P67890'
        assert all('sequence' in protein for protein in proteins)
        assert all(protein['source'] == 'workspace_object' for protein in proteins)
    
    def test_parse_workspace_object_feature_data(self, test_config):
        """Test parsing workspace object with feature data."""
        processor = WorkspaceObjectProcessor(test_config)
        
        object_data = {
            'data': {
                'features': [
                    {
                        'id': 'feature1',
                        'protein_translation': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
                    },
                    {
                        'id': 'feature2',
                        'protein_translation': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
                    }
                ]
            },
            'info': [1, 'test_object', 'KBaseGenomes.Genome-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]
        }
        
        proteins = processor._parse_workspace_object(object_data)
        
        assert len(proteins) == 2
        assert proteins[0]['protein_id'] == 'feature1'
        assert proteins[1]['protein_id'] == 'feature2'
        assert all('sequence' in protein for protein in proteins)
    
    def test_parse_workspace_object_empty_data(self, test_config):
        """Test parsing workspace object with empty data."""
        processor = WorkspaceObjectProcessor(test_config)
        
        object_data = {
            'data': {},
            'info': [1, 'test_object', 'KBaseGenomes.Genome-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]
        }
        
        proteins = processor._parse_workspace_object(object_data)
        
        assert len(proteins) == 0
    
    def test_parse_workspace_object_unknown_structure(self, test_config):
        """Test parsing workspace object with unknown structure."""
        processor = WorkspaceObjectProcessor(test_config)
        
        object_data = {
            'data': {
                'unknown_field': [
                    {'id': 'unknown1', 'data': 'some_data'}
                ]
            },
            'info': [1, 'test_object', 'Unknown.Object-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]
        }
        
        proteins = processor._parse_workspace_object(object_data)
        
        # Should handle gracefully
        assert len(proteins) == 0
    
    def test_workspace_data_dataclass(self):
        """Test WorkspaceData dataclass."""
        workspace_data = WorkspaceData(
            object_ref='1/2/3',
            object_type='KBaseGenomes.Genome-1.0',
            protein_records=[
                {'protein_id': 'P12345', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'}
            ],
            metadata={'organism': 'Test organism'}
        )
        
        assert workspace_data.object_ref == '1/2/3'
        assert workspace_data.object_type == 'KBaseGenomes.Genome-1.0'
        assert len(workspace_data.protein_records) == 1
        assert workspace_data.metadata['organism'] == 'Test organism'
    
    def test_processing_time_tracking(self, test_config, mock_workspace_client):
        """Test that processing time is tracked."""
        processor = WorkspaceObjectProcessor(test_config)
        
        with patch.object(processor, '_get_workspace_client', return_value=mock_workspace_client):
            input_data = {
                'input_type': 'workspace_object',
                'workspace_object': '1/2/3'
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert 'processing_time' in result
            assert isinstance(result['processing_time'], float)
            assert result['processing_time'] >= 0
    
    def test_error_handling_malformed_object_data(self, test_config):
        """Test error handling with malformed object data."""
        processor = WorkspaceObjectProcessor(test_config)
        
        # Test with malformed data structure
        object_data = {
            'data': None,  # Malformed
            'info': [1, 'test_object', 'KBaseGenomes.Genome-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]
        }
        
        proteins = processor._parse_workspace_object(object_data)
        
        # Should handle gracefully
        assert len(proteins) == 0
    
    def test_metadata_extraction(self, test_config, mock_workspace_client):
        """Test that metadata is properly extracted from workspace object."""
        processor = WorkspaceObjectProcessor(test_config)
        
        with patch.object(processor, '_get_workspace_client', return_value=mock_workspace_client):
            input_data = {
                'input_type': 'workspace_object',
                'workspace_object': '1/2/3'
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert 'metadata' in result
            assert 'workspace_info' in result['metadata']
    
    def test_batch_object_processing(self, test_config, mock_workspace_client):
        """Test processing multiple workspace objects."""
        processor = WorkspaceObjectProcessor(test_config)
        
        # Mock multiple objects
        mock_client = Mock()
        mock_client.get_objects2.return_value = {
            'data': [
                {
                    'data': {
                        'proteins': [
                            {'id': 'P12345', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'}
                        ]
                    },
                    'info': [1, 'test_object1', 'KBaseGenomes.Genome-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object1', 'test_checksum', 1, {}]
                },
                {
                    'data': {
                        'proteins': [
                            {'id': 'P67890', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'}
                        ]
                    },
                    'info': [1, 'test_object2', 'KBaseGenomes.Genome-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object2', 'test_checksum', 1, {}]
                }
            ]
        }
        
        with patch.object(processor, '_get_workspace_client', return_value=mock_client):
            input_data = {
                'input_type': 'workspace_object',
                'workspace_object': '1/2/3'
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert len(result['proteins']) == 2

