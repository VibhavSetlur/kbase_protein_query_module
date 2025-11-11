"""
Protein Embedding Generator Module

Uses ESM-2 models from fair-esm package (official Facebook Research implementation).
Models are automatically downloaded on first use.

Follows official ESM usage patterns from:
https://github.com/facebookresearch/esm
"""

import os
import numpy as np
from typing import Optional, Union, Tuple
import logging

try:
    import torch
    import esm
    HAS_ESM = True
except ImportError:
    HAS_ESM = False
    torch = None
    esm = None

logger = logging.getLogger(__name__)

class ProteinEmbeddingGenerator:
    """
    Generates protein embeddings using ESM-2 models from fair-esm package.
    
    Uses official Facebook Research ESM models.
    Simple, clean, efficient, and functional.
    """

    def __init__(self, model_name: str = "esm2_t6_8M_UR50D", device: str = "auto"):
        """Initialize the embedding generator.
        
        Args:
            model_name: ESM model name (e.g., "esm2_t6_8M_UR50D", "esm2_t33_650M_UR50D")
                       Default is smallest model for fast inference
            device: Device to use ("auto", "cpu", or "cuda")
        """
        if not HAS_ESM:
            raise ImportError("torch and esm (fair-esm) are required for embedding generation. Install with: pip install fair-esm")
        
        self.model_name = model_name
        self.device = self._setup_device(device)
        self.model = None
        self.alphabet = None
        self.batch_converter = None
        self.embedding_dim: Optional[int] = None
        self.repr_layer: Optional[int] = None
        
        self._load_model()
    
    def _setup_device(self, device: str) -> str:
        """Setup the device for model inference."""
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        return device
    
    def _load_model(self):
        """Load ESM-2 model using official fair-esm package."""
        try:
            # Map model names to ESM pretrained functions
            model_map = {
                "esm2_t6_8M_UR50D": esm.pretrained.esm2_t6_8M_UR50D,
                "esm2_t12_35M_UR50D": esm.pretrained.esm2_t12_35M_UR50D,
                "esm2_t30_150M_UR50D": esm.pretrained.esm2_t30_150M_UR50D,
                "esm2_t33_650M_UR50D": esm.pretrained.esm2_t33_650M_UR50D,
            }
            
            if self.model_name not in model_map:
                available = ", ".join(model_map.keys())
                raise ValueError(f"Unknown model name: {self.model_name}. Available models: {available}")
            
            # Load model and alphabet
            logger.info(f"Loading ESM model '{self.model_name}'...")
            self.model, self.alphabet = model_map[self.model_name]()
            self.batch_converter = self.alphabet.get_batch_converter()
            
            # Move model to device and set to eval mode
            self.model = self.model.to(self.device)
            self.model.eval()
            
            # Determine representation layer (last layer) and embedding dimension
            # ESM-2 models: layer count varies by model size
            # Use the last layer (num_layers) for all models - this is the standard approach
            self.repr_layer = self.model.num_layers
            
            # Get embedding dimension from model
            self.embedding_dim = self.model.embed_dim
            
            logger.info(f"Loaded ESM model '{self.model_name}' on {self.device} with dim {self.embedding_dim}, layer {self.repr_layer}")
        except Exception as e:
            raise RuntimeError(f"Failed to load ESM model '{self.model_name}': {e}")
    
    def generate_embedding(self, sequence: str, pooling: str = "mean") -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """Generate embeddings for a single protein sequence.
        
        Args:
            sequence: Protein sequence string
            pooling: Pooling strategy - "mean" for mean pooling, "none" for full residue-level embeddings
                    Returns mean embedding if "mean", or (residue_level, mean) tuple if "none"
        
        Returns:
            If pooling="mean": numpy array of shape [D] (mean-pooled embedding)
            If pooling="none": tuple of (residue_level [L x D], mean [D])
        """
        if not sequence or len(sequence.strip()) == 0:
            raise ValueError("Empty or invalid protein sequence")
        
        sequence = self._preprocess_sequence(sequence)
        
        # Prepare data for batch converter (single sequence)
        data = [("protein", sequence)]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_tokens = batch_tokens.to(self.device)
        
        # Calculate sequence length (excluding padding)
        batch_lens = (batch_tokens != self.alphabet.padding_idx).sum(1)
        tokens_len = batch_lens[0].item()
        
        # Generate embeddings
        with torch.no_grad():
            results = self.model(batch_tokens, repr_layers=[self.repr_layer])
        
        # Extract token representations - shape: [batch_size, seq_len, hidden_size]
        token_representations = results["representations"][self.repr_layer]
        
        # Extract sequence representation (exclude padding and special tokens)
        # Token 0 is beginning-of-sequence token, so first residue is token 1
        # Last token before padding is end-of-sequence
        # We want tokens 1 to tokens_len-1 (exclude BOS and EOS)
        if tokens_len >= 2:
            # Extract residue embeddings: tokens 1 to tokens_len-1
            residue_embeddings = token_representations[0, 1:tokens_len-1].cpu().numpy().astype(np.float32)
        else:
            # Very short sequence - use all tokens except padding
            residue_embeddings = token_representations[0, :tokens_len].cpu().numpy().astype(np.float32)
        
        # Compute mean embedding
        mean_embedding = residue_embeddings.mean(axis=0).astype(np.float32)
        
        if pooling == "mean":
            return mean_embedding
        elif pooling == "none":
            return residue_embeddings, mean_embedding
        else:
            raise ValueError(f"Unknown pooling strategy: {pooling}. Use 'mean' or 'none'")
    
    def _preprocess_sequence(self, sequence: str) -> str:
        """Preprocess protein sequence by cleaning and validating.
        
        Args:
            sequence: Raw protein sequence string
            
        Returns:
            Cleaned and validated sequence string
            
        Raises:
            ValueError: If sequence is invalid
        """
        if not sequence or not isinstance(sequence, str):
            raise ValueError("Sequence must be a non-empty string")
        
        # Clean the sequence
        sequence_clean = sequence.strip().upper()
        
        # Remove whitespace and newlines
        sequence_clean = ''.join(sequence_clean.split())
        
        # Validate minimum length
        if len(sequence_clean) < 3:
            raise ValueError(f"Sequence too short: minimum 3 amino acids required, got {len(sequence_clean)}")
        
        # Validate amino acids (ESM models expect standard amino acids)
        # ESM alphabet includes standard 20 amino acids plus special tokens
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        invalid_chars = set(sequence_clean) - valid_aa
        if invalid_chars:
            logger.warning(f"Sequence contains invalid characters: {invalid_chars}. Proceeding anyway.")
        
        return sequence_clean


def main():
    """Test the embedding generator."""
    ok = True
    try:
        test_sequence = "ACDEFGHIKLMNPQRSTVWY"
        
        # Test with default model (smallest, fastest)
        gen = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D", device="cpu")
        
        # Test preprocessing
        preprocessed = gen._preprocess_sequence(test_sequence)
        if preprocessed != test_sequence:
            raise RuntimeError(f"Preprocessing changed sequence: {test_sequence} -> {preprocessed}")
        
        # Test mean pooling
        mean_embedding = gen.generate_embedding(test_sequence, pooling="mean")
        if mean_embedding.size == 0:
            raise RuntimeError("Generated empty mean embedding")
        if mean_embedding.shape[0] != gen.embedding_dim:
            raise RuntimeError(f"Mean embedding dimension {mean_embedding.shape[0]} != {gen.embedding_dim}")
        if not isinstance(mean_embedding, np.ndarray):
            raise RuntimeError(f"Mean embedding is not a numpy array: {type(mean_embedding)}")
        
        # Test full embeddings (no pooling)
        residue_emb, mean_emb = gen.generate_embedding(test_sequence, pooling="none")
        if residue_emb.size == 0 or mean_emb.size == 0:
            raise RuntimeError("Generated empty embeddings")
        if residue_emb.shape[1] != gen.embedding_dim:
            raise RuntimeError(f"Residue embedding dimension mismatch: {residue_emb.shape[1]} != {gen.embedding_dim}")
        if mean_emb.shape[0] != gen.embedding_dim:
            raise RuntimeError(f"Mean embedding dimension mismatch: {mean_emb.shape[0]} != {gen.embedding_dim}")
        
        # Verify mean_emb matches mean-pooled version
        computed_mean = residue_emb.mean(axis=0)
        if not np.allclose(mean_emb, computed_mean, atol=1e-5):
            raise RuntimeError("Mean embedding from 'none' pooling doesn't match computed mean")
        
        print(f"Embedding generator test: SUCCESS")
        print(f"  Model: {gen.model_name}")
        print(f"  Device: {gen.device}")
        print(f"  Embedding dim: {gen.embedding_dim}")
        print(f"  Representation layer: {gen.repr_layer}")
        print(f"  Sequence length: {len(test_sequence)}")
        print(f"  Mean embedding shape: {mean_embedding.shape}")
        print(f"  Residue embedding shape: {residue_emb.shape}")
    except Exception as e:
        ok = False
        print(f"Embedding generator test: FAILED - {e}")
        import traceback
        traceback.print_exc()
    
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
