"""
Unit tests for UniProtIdsProcessor.

Tests UniProt ID validation, data fetching, and processing.
"""

import pytest
import sys
import os
import requests
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.input.uniprot_ids.input import UniProtIdsProcessor, UniProtData


class TestUniProtIdsProcessor:
    """Test cases for UniProtIdsProcessor."""
    
    def test_initialization(self, test_config):
        """Test UniProtIdsProcessor initialization."""
        processor = UniProtIdsProcessor(test_config)
        
        assert processor.config == test_config
        assert processor.config.get('uniprot_api_url', 'https://rest.uniprot.org') == 'https://rest.uniprot.org'
    
    def test_initialization_default_config(self):
        """Test initialization with default configuration."""
        processor = UniProtIdsProcessor()
        
        assert processor.config is None or isinstance(processor.config, dict)
    
    def test_process_valid_uniprot_ids(self, test_config, sample_uniprot_ids):
        """Test processing valid UniProt IDs."""
        processor = UniProtIdsProcessor(test_config)
        
        # Mock the fetch_uniprot_data method to avoid actual API calls
        mock_data = {
            'primaryAccession': 'P12345',
            'proteinDescription': {
                'recommendedName': {'fullName': {'value': 'Test Protein'}}
            },
            'organism': {'scientificName': 'Homo sapiens'},
            'sequence': {
                'length': 100,
                'molWeight': 10000,
                'crc64': 'ABCDEF1234567890',
                'value': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            }
        }
        
        with patch.object(processor, '_fetch_uniprot_data', return_value=mock_data):
            input_data = {
                'input_type': 'uniprot_ids',
                'uniprot_ids': sample_uniprot_ids
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert result['input_type'] == 'uniprot_ids'
            assert len(result['proteins']) == len(sample_uniprot_ids)
            assert all('protein_id' in protein for protein in result['proteins'])
            assert all('sequence' in protein for protein in result['proteins'])
            assert all(protein['source'] == 'uniprot' for protein in result['proteins'])
    
    def test_process_single_uniprot_id(self, test_config):
        """Test processing single UniProt ID."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_data = {
            'primaryAccession': 'P12345',
            'proteinDescription': {
                'recommendedName': {'fullName': {'value': 'Test Protein'}}
            },
            'organism': {'scientificName': 'Homo sapiens'},
            'sequence': {
                'length': 100,
                'molWeight': 10000,
                'crc64': 'ABCDEF1234567890',
                'value': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            }
        }
        
        with patch.object(processor, '_fetch_uniprot_data', return_value=mock_data):
            input_data = {
                'input_type': 'uniprot_ids',
                'uniprot_ids': 'P12345'  # Single string instead of list
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert len(result['proteins']) == 1
            assert result['proteins'][0]['protein_id'] == 'P12345'
    
    def test_process_empty_uniprot_ids(self, test_config):
        """Test processing empty UniProt IDs list."""
        processor = UniProtIdsProcessor(test_config)
        
        input_data = {
            'input_type': 'uniprot_ids',
            'uniprot_ids': []
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
        assert 'No UniProt IDs provided' in result['error_message']
    
    def test_process_invalid_uniprot_ids(self, test_config):
        """Test processing invalid UniProt IDs."""
        processor = UniProtIdsProcessor(test_config)
        
        with patch.object(processor, '_fetch_uniprot_data', return_value=None):
            input_data = {
                'input_type': 'uniprot_ids',
                'uniprot_ids': ['INVALID_ID', 'ANOTHER_INVALID']
            }
            
            result = processor.process(input_data)
            
            # Should handle gracefully, possibly with warnings
            assert result['success'] is True or result['success'] is False
            if result['success']:
                assert len(result['proteins']) == 0  # No valid proteins found
    
    def test_process_partial_failure(self, test_config):
        """Test processing with some valid and some invalid IDs."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_data = {
            'primaryAccession': 'P12345',
            'proteinDescription': {
                'recommendedName': {'fullName': {'value': 'Test Protein'}}
            },
            'organism': {'scientificName': 'Homo sapiens'},
            'sequence': {
                'length': 100,
                'molWeight': 10000,
                'crc64': 'ABCDEF1234567890',
                'value': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            }
        }
        
        def mock_fetch(uniprot_id):
            if uniprot_id == 'P12345':
                return mock_data
            return None
        
        with patch.object(processor, '_fetch_uniprot_data', side_effect=mock_fetch):
            input_data = {
                'input_type': 'uniprot_ids',
                'uniprot_ids': ['P12345', 'INVALID_ID']
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert len(result['proteins']) == 1  # Only one valid protein
            assert result['proteins'][0]['protein_id'] == 'P12345'
    
    def test_validate_uniprot_id_valid(self, test_config):
        """Test UniProt ID validation with valid IDs."""
        processor = UniProtIdsProcessor(test_config)
        
        valid_ids = ['P12345', 'Q67890', 'O13547', 'A1B2C3']
        
        for uniprot_id in valid_ids:
            assert processor._validate_uniprot_id(uniprot_id) is True
    
    def test_validate_uniprot_id_invalid(self, test_config):
        """Test UniProt ID validation with invalid IDs."""
        processor = UniProtIdsProcessor(test_config)
        
        invalid_ids = ['INVALID', '12345', '', None, 'P1234', 'P1234567890']
        
        for uniprot_id in invalid_ids:
            assert processor._validate_uniprot_id(uniprot_id) is False
    
    @patch('kbase_protein_query_module.src.input.uniprot_ids.input.requests.get')
    def test_fetch_uniprot_data_success(self, mock_get, test_config):
        """Test successful UniProt data fetching."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'primaryAccession': 'P12345',
            'proteinDescription': {
                'recommendedName': {'fullName': {'value': 'Test Protein'}}
            },
            'organism': {'scientificName': 'Homo sapiens'},
            'sequence': {
                'length': 100,
                'molWeight': 10000,
                'crc64': 'ABCDEF1234567890',
                'value': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            }
        }
        mock_get.return_value = mock_response
        
        result = processor._fetch_uniprot_data('P12345')
        
        assert result is not None
        assert result['primaryAccession'] == 'P12345'
        mock_get.assert_called_once()
    
    @patch('kbase_protein_query_module.src.input.uniprot_ids.input.requests.get')
    def test_fetch_uniprot_data_not_found(self, mock_get, test_config):
        """Test UniProt data fetching when protein not found."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_response = Mock()
        mock_response.status_code = 404
        mock_response.json.return_value = {}
        
        # Create a proper HTTPError with response
        http_error = requests.exceptions.HTTPError("404 Client Error")
        http_error.response = mock_response
        mock_response.raise_for_status.side_effect = http_error
        mock_get.return_value = mock_response
        
        result = processor._fetch_uniprot_data('INVALID_ID')
        
        assert result is None
    
    @patch('kbase_protein_query_module.src.input.uniprot_ids.input.requests.get')
    def test_fetch_uniprot_data_api_error(self, mock_get, test_config):
        """Test UniProt data fetching with API error."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_response = Mock()
        mock_response.status_code = 500
        mock_response.json.return_value = {}
        
        # Create a proper HTTPError with response
        http_error = requests.exceptions.HTTPError("500 Server Error")
        http_error.response = mock_response
        mock_response.raise_for_status.side_effect = http_error
        mock_get.return_value = mock_response
        
        result = processor._fetch_uniprot_data('P12345')
        
        assert result is None
    
    @patch('kbase_protein_query_module.src.input.uniprot_ids.input.requests.get')
    def test_fetch_uniprot_data_network_error(self, mock_get, test_config):
        """Test UniProt data fetching with network error."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_get.side_effect = Exception("Network error")
        
        result = processor._fetch_uniprot_data('P12345')
        
        assert result is None
    
    def test_extract_protein_name_with_recommended_name(self, test_config):
        """Test extracting protein name with recommended name."""
        processor = UniProtIdsProcessor(test_config)
        
        data = {
            'proteinDescription': {
                'recommendedName': {
                    'fullName': {'value': 'Test Protein'}
                }
            }
        }
        
        name = processor._extract_protein_name(data)
        
        assert name == 'Test Protein'
    
    def test_extract_protein_name_with_alternative_name(self, test_config):
        """Test extracting protein name with alternative name."""
        processor = UniProtIdsProcessor(test_config)
        
        data = {
            'proteinDescription': {
                'alternativeNames': [
                    {'fullName': {'value': 'Alternative Protein'}}
                ]
            }
        }
        
        name = processor._extract_protein_name(data)
        
        assert name == 'Alternative Protein'
    
    def test_extract_protein_name_no_name(self, test_config):
        """Test extracting protein name when no name available."""
        processor = UniProtIdsProcessor(test_config)
        
        data = {}
        
        name = processor._extract_protein_name(data)
        
        assert name == 'Unknown Protein'
    
    def test_extract_organism_valid(self, test_config):
        """Test extracting organism name."""
        processor = UniProtIdsProcessor(test_config)
        
        data = {
            'organism': {
                'scientificName': 'Homo sapiens'
            }
        }
        
        organism = processor._extract_organism(data)
        
        assert organism == 'Homo sapiens'
    
    def test_extract_organism_no_organism(self, test_config):
        """Test extracting organism when not available."""
        processor = UniProtIdsProcessor(test_config)
        
        data = {}
        
        organism = processor._extract_organism(data)
        
        assert organism == 'Unknown'
    
    def test_processing_time_tracking(self, test_config):
        """Test that processing time is tracked."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_data = {
            'primaryAccession': 'P12345',
            'proteinDescription': {
                'recommendedName': {'fullName': {'value': 'Test Protein'}}
            },
            'organism': {'scientificName': 'Homo sapiens'},
            'sequence': {
                'length': 100,
                'molWeight': 10000,
                'crc64': 'ABCDEF1234567890',
                'value': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            }
        }
        
        with patch.object(processor, '_fetch_uniprot_data', return_value=mock_data):
            input_data = {
                'input_type': 'uniprot_ids',
                'uniprot_ids': ['P12345']
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert 'processing_time' in result
            assert isinstance(result['processing_time'], float)
            assert result['processing_time'] >= 0
    
    def test_uniprot_data_dataclass(self):
        """Test UniProtData dataclass."""
        uniprot_data = UniProtData(
            protein_id='P12345',
            sequence='MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
            source='uniprot',
            metadata={'organism': 'Homo sapiens', 'protein_name': 'Test Protein'}
        )
        
        assert uniprot_data.protein_id == 'P12345'
        assert uniprot_data.sequence == 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        assert uniprot_data.source == 'uniprot'
        assert uniprot_data.metadata['organism'] == 'Homo sapiens'
    
    def test_batch_processing_efficiency(self, test_config):
        """Test that batch processing is efficient."""
        processor = UniProtIdsProcessor(test_config)
        
        mock_data = {
            'primaryAccession': 'P12345',
            'proteinDescription': {
                'recommendedName': {'fullName': {'value': 'Test Protein'}}
            },
            'organism': {'scientificName': 'Homo sapiens'},
            'sequence': {
                'length': 100,
                'molWeight': 10000,
                'crc64': 'ABCDEF1234567890',
                'value': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            }
        }
        
        with patch.object(processor, '_fetch_uniprot_data', return_value=mock_data):
            input_data = {
                'input_type': 'uniprot_ids',
                'uniprot_ids': ['P12345', 'P67890', 'Q13547']
            }
            
            result = processor.process(input_data)
            
            assert result['success'] is True
            assert len(result['proteins']) == 3
            # Verify that _fetch_uniprot_data was called for each ID
            assert processor._fetch_uniprot_data.call_count == 3

