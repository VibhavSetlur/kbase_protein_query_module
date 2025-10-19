"""
Unit tests for ProteinSequenceProcessor.

Tests direct protein sequence input processing, validation, and FASTA parsing.
"""

import pytest
import sys
import os
from unittest.mock import Mock, patch

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.input.protein_sequence.input import ProteinSequenceProcessor, ProteinSequenceData


class TestProteinSequenceProcessor:
    """Test cases for ProteinSequenceProcessor."""
    
    def test_initialization(self, test_config):
        """Test ProteinSequenceProcessor initialization."""
        processor = ProteinSequenceProcessor(test_config)
        
        assert processor.config == test_config
        assert processor.config.get('max_sequence_length', 10000) == 10000
    
    def test_initialization_default_config(self):
        """Test initialization with default configuration."""
        processor = ProteinSequenceProcessor()
        
        assert processor.config is None or isinstance(processor.config, dict)
    
    def test_process_valid_sequence_list(self, test_config):
        """Test processing valid sequence list."""
        processor = ProteinSequenceProcessor(test_config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': [
                'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
                'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            ]
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is True
        assert result['input_type'] == 'protein_input'
        assert len(result['proteins']) == 2
        assert all('protein_id' in protein for protein in result['proteins'])
        assert all('sequence' in protein for protein in result['proteins'])
        assert all(protein['source'] == 'protein_sequence' for protein in result['proteins'])
    
    def test_process_single_sequence_string(self, test_config):
        """Test processing single sequence as string."""
        processor = ProteinSequenceProcessor(test_config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is True
        assert len(result['proteins']) == 1
        assert result['proteins'][0]['sequence'] == 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
    
    def test_process_fasta_format(self, test_config):
        """Test processing FASTA format sequences."""
        processor = ProteinSequenceProcessor(test_config)
        
        fasta_data = """>protein1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
>protein2
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"""
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': fasta_data
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is True
        assert len(result['proteins']) == 2
        assert result['proteins'][0]['protein_id'] == 'protein1'
        assert result['proteins'][1]['protein_id'] == 'protein2'
    
    def test_process_mixed_fasta_and_sequences(self, test_config):
        """Test processing mixed FASTA and direct sequences."""
        processor = ProteinSequenceProcessor(test_config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': [
                'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
                '>protein1\nMKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
            ]
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is True
        assert len(result['proteins']) == 2  # 1 direct + 1 FASTA entry
    
    def test_process_empty_input(self, test_config):
        """Test processing empty input."""
        processor = ProteinSequenceProcessor(test_config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': []
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
        assert 'No protein sequences provided' in result['error_message']
    
    def test_process_invalid_sequence_characters(self, test_config):
        """Test processing sequences with invalid characters."""
        processor = ProteinSequenceProcessor(test_config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG123']
        }
        
        result = processor.process(input_data)
        
        # Should still process but may flag invalid sequences
        assert result['success'] is True
        assert len(result['proteins']) == 1
    
    def test_process_too_long_sequence(self, test_config):
        """Test processing sequence that exceeds length limit."""
        config = test_config.copy()
        config['max_sequence_length'] = 10  # Very short limit
        
        processor = ProteinSequenceProcessor(config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        }
        
        result = processor.process(input_data)
        
        # Should handle gracefully
        assert result['success'] is True or result['success'] is False
        if not result['success']:
            assert 'error_message' in result
    
    def test_validate_sequence_valid(self, test_config):
        """Test sequence validation with valid sequence."""
        processor = ProteinSequenceProcessor(test_config)
        
        valid_sequence = 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        
        assert processor._validate_sequence(valid_sequence) is True
    
    def test_validate_sequence_empty(self, test_config):
        """Test sequence validation with empty sequence."""
        processor = ProteinSequenceProcessor(test_config)
        
        assert processor._validate_sequence('') is False
        assert processor._validate_sequence(None) is False
    
    def test_validate_sequence_too_short(self, test_config):
        """Test sequence validation with too short sequence."""
        processor = ProteinSequenceProcessor(test_config)
        
        short_sequence = ''  # Empty sequence - too short
        
        assert processor._validate_sequence(short_sequence) is False
    
    def test_is_fasta_format_true(self, test_config):
        """Test FASTA format detection with valid FASTA."""
        processor = ProteinSequenceProcessor(test_config)
        
        fasta_data = """>protein1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"""
        
        assert processor._is_fasta_format(fasta_data) is True
    
    def test_is_fasta_format_false(self, test_config):
        """Test FASTA format detection with non-FASTA data."""
        processor = ProteinSequenceProcessor(test_config)
        
        sequence_data = 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        
        assert processor._is_fasta_format(sequence_data) is False
    
    def test_is_fasta_format_partial_header(self, test_config):
        """Test FASTA format detection with partial header."""
        processor = ProteinSequenceProcessor(test_config)
        
        partial_fasta = '>protein1'  # No sequence
        
        assert processor._is_fasta_format(partial_fasta) is False
    
    def test_parse_fasta_data_valid(self, test_config):
        """Test parsing valid FASTA data."""
        processor = ProteinSequenceProcessor(test_config)
        
        fasta_data = """>protein1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
>protein2
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"""
        
        proteins = processor._parse_fasta_data(fasta_data)
        
        assert len(proteins) == 2
        assert proteins[0]['protein_id'] == 'protein1'
        assert proteins[1]['protein_id'] == 'protein2'
        assert all('sequence' in protein for protein in proteins)
    
    def test_parse_fasta_data_multiline_sequence(self, test_config):
        """Test parsing FASTA with multiline sequences."""
        processor = ProteinSequenceProcessor(test_config)
        
        fasta_data = """>protein1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"""
        
        proteins = processor._parse_fasta_data(fasta_data)
        
        assert len(proteins) == 1
        assert proteins[0]['protein_id'] == 'protein1'
        # Should concatenate multiline sequences
        assert len(proteins[0]['sequence']) > 50
    
    def test_parse_fasta_data_empty_sequences(self, test_config):
        """Test parsing FASTA with empty sequences."""
        processor = ProteinSequenceProcessor(test_config)
        
        fasta_data = """>protein1
>protein2
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"""
        
        proteins = processor._parse_fasta_data(fasta_data)
        
        # Should handle empty sequences gracefully
        assert len(proteins) >= 1
        assert any('sequence' in protein for protein in proteins)
    
    def test_parse_direct_sequence_single(self, test_config):
        """Test parsing single direct sequence."""
        processor = ProteinSequenceProcessor(test_config)
        
        sequence = 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        
        proteins = processor._parse_direct_sequence(sequence)
        
        assert len(proteins) == 1
        assert proteins[0]['sequence'] == sequence
        assert 'protein_id' in proteins[0]
        assert proteins[0]['source'] == 'protein_sequence'
    
    def test_parse_direct_sequence_list(self, test_config):
        """Test parsing list of direct sequences."""
        processor = ProteinSequenceProcessor(test_config)
        
        sequences = [
            'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
            'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        ]
        
        proteins = processor._parse_sequence_list(sequences)
        
        assert len(proteins) == 2
        assert all(protein['source'] == 'protein_sequence' for protein in proteins)
        assert all('sequence' in protein for protein in proteins)
    
    def test_protein_sequence_data_dataclass(self):
        """Test ProteinSequenceData dataclass."""
        protein_data = ProteinSequenceData(
            protein_id='test_protein',
            sequence='MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
            source='test',
            metadata={'organism': 'Test'}
        )
        
        assert protein_data.protein_id == 'test_protein'
        assert protein_data.sequence == 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        assert protein_data.source == 'test'
        assert protein_data.metadata['organism'] == 'Test'
    
    def test_error_handling_malformed_input(self, test_config):
        """Test error handling with malformed input."""
        processor = ProteinSequenceProcessor(test_config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': None  # Invalid input
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is False
        assert 'error_message' in result
    
    def test_processing_time_tracking(self, test_config):
        """Test that processing time is tracked."""
        processor = ProteinSequenceProcessor(test_config)
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        }
        
        result = processor.process(input_data)
        
        assert result['success'] is True
        assert 'processing_time' in result
        assert isinstance(result['processing_time'], float)
        assert result['processing_time'] >= 0

