"""
Unit tests for ProteinEmbeddingGenerator.

Tests embedding generation, batch processing, and model management.
"""

import pytest
import sys
import os
import numpy as np
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.util.embeddings.generator import ProteinEmbeddingGenerator


class TestProteinEmbeddingGenerator:
    """Test cases for ProteinEmbeddingGenerator."""
    
    def test_initialization_default(self):
        """Test default initialization."""
        generator = ProteinEmbeddingGenerator()
        
        assert generator.model_name == "esm2_t6_8M_UR50D"
        assert generator.device == "auto"
        assert generator.model is None
        assert generator.tokenizer is None
        assert generator.embedding_dim is None
    
    def test_initialization_custom(self):
        """Test initialization with custom parameters."""
        generator = ProteinEmbeddingGenerator(
            model_name="custom_model",
            device="cpu"
        )
        
        assert generator.model_name == "custom_model"
        assert generator.device == "cpu"
    
    def test_setup_device_auto_cpu(self):
        """Test automatic device setup for CPU."""
        generator = ProteinEmbeddingGenerator(device="auto")
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.cuda.is_available', return_value=False):
            with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.device') as mock_device:
                generator._setup_device("auto")
                mock_device.assert_called_with("cpu")
    
    def test_setup_device_auto_gpu(self):
        """Test automatic device setup for GPU."""
        generator = ProteinEmbeddingGenerator(device="auto")
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.cuda.is_available', return_value=True):
            with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.device') as mock_device:
                generator._setup_device("auto")
                mock_device.assert_called_with("cuda")
    
    def test_setup_device_explicit(self):
        """Test explicit device setup."""
        generator = ProteinEmbeddingGenerator(device="cpu")
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.device') as mock_device:
            generator._setup_device("cpu")
            mock_device.assert_called_with("cpu")
    
    @patch('kbase_protein_query_module.src.util.embeddings.generator.torch')
    def test_load_model_success(self, mock_torch):
        """Test successful model loading."""
        generator = ProteinEmbeddingGenerator()
        
        # Mock model and tokenizer
        mock_model = Mock()
        mock_tokenizer = Mock()
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.AutoTokenizer.from_pretrained', return_value=mock_tokenizer):
            with patch('kbase_protein_query_module.src.util.embeddings.generator.AutoModel.from_pretrained', return_value=mock_model):
                generator._load_model()
                
                assert generator.model == mock_model
                assert generator.tokenizer == mock_tokenizer
                mock_model.eval.assert_called_once()
    
    @patch('kbase_protein_query_module.src.util.embeddings.generator.torch')
    def test_load_model_import_error(self, mock_torch):
        """Test model loading with import error."""
        generator = ProteinEmbeddingGenerator()
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.AutoTokenizer.from_pretrained', side_effect=ImportError("No transformers")):
            # Should not raise exception, just set model to None
            generator._load_model()
            assert generator.model is None
            assert generator.tokenizer is None
    
    def test_generate_embedding_valid_sequence(self):
        """Test generating embedding for valid sequence."""
        generator = ProteinEmbeddingGenerator()
        
        # Mock the model and tokenizer
        mock_model = Mock()
        mock_tokenizer = Mock()
        mock_model.eval.return_value = None
        mock_model.to.return_value = mock_model
        
        # Mock tokenizer output
        mock_input_ids = Mock()
        mock_input_ids.to.return_value = mock_input_ids
        mock_attention_mask = Mock()
        mock_attention_mask.to.return_value = mock_attention_mask
        
        mock_tokens = {
            'input_ids': mock_input_ids,
            'attention_mask': mock_attention_mask
        }
        mock_tokenizer.return_value = mock_tokens
        
        # Mock model output with proper tensor methods
        mock_tensor = Mock()
        mock_tensor.mean.return_value = Mock()
        mock_tensor.mean.return_value.squeeze.return_value = Mock()
        mock_tensor.mean.return_value.squeeze.return_value.cpu.return_value = Mock()
        mock_tensor.mean.return_value.squeeze.return_value.cpu.return_value.numpy.return_value = np.random.rand(320)
        
        mock_output = Mock()
        mock_output.last_hidden_state = mock_tensor
        mock_model.return_value = mock_output
        
        generator.model = mock_model
        generator.tokenizer = mock_tokenizer
        generator.embedding_dim = 320
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.no_grad'):
            with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.mean') as mock_mean:
                mock_mean.return_value = np.random.rand(320)
                
                embedding = generator.generate_embedding("MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG")
                
                assert embedding is not None
                assert len(embedding.shape) == 1  # Should be 1D
                assert embedding.shape[0] == 320
    
    def test_generate_embedding_invalid_sequence(self):
        """Test generating embedding for invalid sequence."""
        generator = ProteinEmbeddingGenerator()
        
        # Test with empty sequence
        embedding = generator.generate_embedding("")
        assert embedding is None
        
        # Test with None sequence
        embedding = generator.generate_embedding(None)
        assert embedding is None
    
    def test_generate_embedding_no_model(self):
        """Test generating embedding when model is not loaded."""
        generator = ProteinEmbeddingGenerator()
        generator.model = None
        
        embedding = generator.generate_embedding("MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG")
        
        assert embedding is None
    
    def test_generate_embeddings_batch(self):
        """Test generating embeddings in batch."""
        generator = ProteinEmbeddingGenerator()
        
        # Mock the generate_embedding method
        mock_embeddings = {
            'P12345': np.random.rand(320),
            'P67890': np.random.rand(320)
        }
        
        def mock_generate(sequence, protein_id=None):
            return mock_embeddings.get(protein_id, np.random.rand(320))
        
        generator.generate_embedding = mock_generate
        
        sequences = ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'] * 2
        protein_ids = ['P12345', 'P67890']
        
        result = generator.generate_embeddings_batch(sequences, protein_ids)
        
        assert isinstance(result, dict)
        assert len(result) == 2
        assert 'P12345' in result
        assert 'P67890' in result
    
    def test_generate_embeddings_batch_empty(self):
        """Test generating embeddings for empty batch."""
        generator = ProteinEmbeddingGenerator()
        
        result = generator.generate_embeddings_batch([], [])
        
        assert result == {}
    
    def test_generate_embeddings_batch_mismatched_lengths(self):
        """Test generating embeddings with mismatched sequence and ID lengths."""
        generator = ProteinEmbeddingGenerator()
        
        # Mock the generate_embedding method
        mock_embeddings = {
            'P12345': np.random.rand(320)
        }
        
        def mock_generate(sequence, protein_id=None):
            return mock_embeddings.get(protein_id, np.random.rand(320))
        
        generator.generate_embedding = mock_generate
        
        sequences = ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
        protein_ids = ['P12345', 'P67890']  # Mismatched length
        
        result = generator.generate_embeddings_batch(sequences, protein_ids)
        
        # Should handle gracefully, only process matching pairs
        assert len(result) == 1
    
    def test_generate_embeddings_from_dict(self):
        """Test generating embeddings from dictionary."""
        generator = ProteinEmbeddingGenerator()
        
        # Mock the generate_embeddings_batch method
        mock_result = {
            'P12345': np.random.rand(320),
            'P67890': np.random.rand(320)
        }
        
        generator.generate_embeddings_batch = Mock(return_value=mock_result)
        
        inputs = {
            'P12345': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
            'P67890': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        }
        
        result = generator.generate_embeddings(inputs)
        
        assert result == mock_result
        generator.generate_embeddings_batch.assert_called_once()
    
    def test_save_embeddings(self, temp_dir):
        """Test saving embeddings to file."""
        generator = ProteinEmbeddingGenerator()
        
        embeddings = {
            'P12345': np.random.rand(320),
            'P67890': np.random.rand(320)
        }
        
        output_file = os.path.join(temp_dir, 'embeddings.npz')
        
        generator.save_embeddings(embeddings, output_file)
        
        assert os.path.exists(output_file)
        
        # Test loading back
        loaded_embeddings, loaded_sequences, loaded_metadata = generator.load_embeddings(output_file)
        
        assert 'P12345' in loaded_embeddings
        assert 'P67890' in loaded_embeddings
        np.testing.assert_array_equal(loaded_embeddings['P12345'], embeddings['P12345'])
    
    def test_save_embeddings_with_sequences(self, temp_dir):
        """Test saving embeddings with sequences."""
        generator = ProteinEmbeddingGenerator()
        
        embeddings = {
            'P12345': np.random.rand(320),
            'P67890': np.random.rand(320)
        }
        
        sequences = {
            'P12345': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
            'P67890': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        }
        
        output_file = os.path.join(temp_dir, 'embeddings_with_sequences.npz')
        
        generator.save_embeddings(embeddings, output_file, sequences_dict=sequences)
        
        loaded_embeddings, loaded_sequences, loaded_metadata = generator.load_embeddings(output_file)
        
        assert loaded_sequences == sequences
    
    def test_save_embeddings_with_metadata(self, temp_dir):
        """Test saving embeddings with metadata."""
        generator = ProteinEmbeddingGenerator()
        
        embeddings = {
            'P12345': np.random.rand(320),
            'P67890': np.random.rand(320)
        }
        
        metadata = {
            'P12345': {'organism': 'Homo sapiens', 'function': 'Test'},
            'P67890': {'organism': 'Mus musculus', 'function': 'Test'}
        }
        
        output_file = os.path.join(temp_dir, 'embeddings_with_metadata.npz')
        
        generator.save_embeddings(embeddings, output_file, metadata=metadata)
        
        loaded_embeddings, loaded_sequences, loaded_metadata = generator.load_embeddings(output_file)
        
        assert loaded_metadata == metadata
    
    def test_load_embeddings_nonexistent_file(self):
        """Test loading embeddings from non-existent file."""
        generator = ProteinEmbeddingGenerator()
        
        with pytest.raises(FileNotFoundError):
            generator.load_embeddings('nonexistent_file.npz')
    
    def test_normalize_embeddings(self):
        """Test embedding normalization."""
        generator = ProteinEmbeddingGenerator()
        
        # Create test embeddings
        embeddings = np.array([[3.0, 4.0], [1.0, 1.0]])
        
        normalized = generator.normalize_embeddings(embeddings)
        
        # Check that each row has unit norm
        norms = np.linalg.norm(normalized, axis=1)
        np.testing.assert_array_almost_equal(norms, [1.0, 1.0])
    
    def test_validate_sequence_valid(self):
        """Test sequence validation with valid sequences."""
        generator = ProteinEmbeddingGenerator()
        
        valid_sequences = [
            'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
            'A' * 100,  # Long sequence
            'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        ]
        
        for sequence in valid_sequences:
            assert generator._validate_sequence(sequence) is True
    
    def test_validate_sequence_invalid(self):
        """Test sequence validation with invalid sequences."""
        generator = ProteinEmbeddingGenerator()
        
        invalid_sequences = [
            '',  # Empty
            None,  # None
            'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG123',  # Invalid characters
            'MK',  # Too short
            'X' * 10001  # Too long
        ]
        
        for sequence in invalid_sequences:
            assert generator._validate_sequence(sequence) is False
    
    def test_preprocess_sequence(self):
        """Test sequence preprocessing."""
        generator = ProteinEmbeddingGenerator()
        
        # Test removing spaces and converting to uppercase
        sequence = ' mktvrqer lksivr ilersk '
        processed = generator._preprocess_sequence(sequence)
        
        assert processed == 'MKTVRQERLKSIVRILERSK'
    
    def test_get_embedding_dimension(self):
        """Test getting embedding dimension."""
        generator = ProteinEmbeddingGenerator()
        generator.embedding_dim = 320
        
        assert generator.get_embedding_dimension() == 320
        
        # Test default when not set
        generator.embedding_dim = None
        assert generator.get_embedding_dimension() == 320
    
    def test_get_model_info(self):
        """Test getting model information."""
        generator = ProteinEmbeddingGenerator()
        generator.embedding_dim = 320
        
        model_info = generator.get_model_info()
        
        assert isinstance(model_info, dict)
        assert 'model' in model_info
        assert 'dimension' in model_info
        assert model_info['model'] == "esm2_t6_8M_UR50D"
        assert model_info['dimension'] == 320
    
    def test_compute_similarity(self):
        """Test computing similarity between embeddings."""
        generator = ProteinEmbeddingGenerator()
        
        # Create test embeddings
        embedding1 = np.array([1.0, 0.0, 0.0])
        embedding2 = np.array([0.0, 1.0, 0.0])
        embedding3 = np.array([1.0, 0.0, 0.0])  # Same as embedding1
        
        # Test orthogonal vectors
        similarity = generator.compute_similarity(embedding1, embedding2)
        assert abs(similarity) < 0.1  # Should be close to 0
        
        # Test identical vectors
        similarity = generator.compute_similarity(embedding1, embedding3)
        assert abs(similarity - 1.0) < 0.1  # Should be close to 1
    
    def test_error_handling_model_loading(self):
        """Test error handling during model loading."""
        generator = ProteinEmbeddingGenerator()
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.AutoTokenizer.from_pretrained', side_effect=Exception("Model loading error")):
            # Should not raise exception
            generator._load_model()
            assert generator.model is None
            assert generator.tokenizer is None
    
    def test_error_handling_embedding_generation(self):
        """Test error handling during embedding generation."""
        generator = ProteinEmbeddingGenerator()
        
        # Mock model that raises exception
        mock_model = Mock()
        mock_model.side_effect = Exception("Model error")
        
        # Mock tokenizer that returns proper format
        mock_tokenizer = Mock()
        mock_input_ids = Mock()
        mock_input_ids.to.return_value = mock_input_ids
        mock_attention_mask = Mock()
        mock_attention_mask.to.return_value = mock_attention_mask
        
        mock_tokenizer.return_value = {
            'input_ids': mock_input_ids,
            'attention_mask': mock_attention_mask
        }
        
        generator.model = mock_model
        generator.tokenizer = mock_tokenizer
        
        with patch('kbase_protein_query_module.src.util.embeddings.generator.torch.no_grad'):
            embedding = generator.generate_embedding("MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG")
        
        assert embedding is None
    
    def test_batch_size_handling(self):
        """Test handling of different batch sizes."""
        generator = ProteinEmbeddingGenerator()
        
        # Mock generate_embedding to track calls
        generator.generate_embedding = Mock(return_value=np.random.rand(320))
        
        sequences = ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'] * 5
        protein_ids = [f'P{i}' for i in range(5)]
        
        result = generator.generate_embeddings_batch(sequences, protein_ids, batch_size=2)
        
        assert len(result) == 5
        assert generator.generate_embedding.call_count == 5

