"""
Protein Embedding Generator Module

Uses ESM-2 models from HuggingFace for reliable embedding generation.
Models are automatically downloaded on first use.
"""

import os
import numpy as np
from typing import Optional, Union, Tuple
import logging

try:
    import torch
    from transformers import EsmModel, EsmTokenizer
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False
    torch = None
    EsmModel = None
    EsmTokenizer = None

logger = logging.getLogger(__name__)

class ProteinEmbeddingGenerator:
    """
    Generates protein embeddings using ESM-2 models from HuggingFace.
    
    Uses facebook/esm2-* models which are automatically downloaded on first use.
    Simple, reliable, and functional.
    """

    def __init__(self, model_name: str = "facebook/esm2_t6_8M_UR50D", device: str = "auto"):
        """Initialize the embedding generator.
        
        Args:
            model_name: HuggingFace model identifier (e.g., "facebook/esm2_t6_8M_UR50D")
            device: Device to use ("auto", "cpu", or "cuda")
        """
        if not HAS_TORCH:
            raise ImportError("torch and transformers are required for embedding generation")
        
        self.model_name = model_name
        self.device = self._setup_device(device)
        self.model = None
        self.tokenizer = None
        self.embedding_dim: Optional[int] = None
        
        self._load_model()
    
    def _setup_device(self, device: str) -> str:
        """Setup the device for model inference."""
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        return device
    
    def _load_model(self):
        """Load ESM-2 model from HuggingFace."""
        try:
            # Load tokenizer and model from HuggingFace
            # Models are automatically downloaded on first use
            self.tokenizer = EsmTokenizer.from_pretrained(self.model_name)
            self.model = EsmModel.from_pretrained(self.model_name)
            self.model.eval()
            self.model = self.model.to(self.device)
            
            # Get embedding dimension from model config
            self.embedding_dim = self.model.config.hidden_size
            
            logger.info(f"Loaded ESM model '{self.model_name}' on {self.device} with dim {self.embedding_dim}")
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
        
        # Tokenize sequence
        tokens = self.tokenizer(sequence, return_tensors="pt", max_length=1024, truncation=True, padding=False)
        tokens = {k: v.to(self.device) for k, v in tokens.items()}
        
        # Generate embeddings
        with torch.no_grad():
            outputs = self.model(**tokens)
        
        # Extract embeddings - shape: [1, seq_len, hidden_size]
        embeddings = outputs.last_hidden_state.squeeze(0)  # [seq_len, hidden_size]
        
        # Remove special tokens (CLS at start, EOS/SEP at end)
        # ESM-2 tokenizer: <cls> sequence <eos>
        seq_len = embeddings.shape[0]
        if seq_len >= 2:
            # Remove first token (CLS) and last token (EOS)
            residue_embeddings = embeddings[1:seq_len-1].cpu().numpy().astype(np.float32)
        else:
            residue_embeddings = embeddings.cpu().numpy().astype(np.float32)
        
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
        
        # Validate amino acids (allow standard 20 plus common modifications)
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
        gen = ProteinEmbeddingGenerator(model_name="facebook/esm2_t6_8M_UR50D", device="cpu")
        
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
        
        # Test full embeddings (no pooling)
        residue_emb, mean_emb = gen.generate_embedding(test_sequence, pooling="none")
        if residue_emb.size == 0 or mean_emb.size == 0:
            raise RuntimeError("Generated empty embeddings")
        if residue_emb.shape[1] != gen.embedding_dim:
            raise RuntimeError(f"Residue embedding dimension mismatch")
        if mean_emb.shape[0] != gen.embedding_dim:
            raise RuntimeError(f"Mean embedding dimension mismatch")
        
        print(f"Embedding generator test: SUCCESS (model={gen.model_name}, dim={gen.embedding_dim})")
    except Exception as e:
        ok = False
        print(f"Embedding generator test: FAILED - {e}")
        import traceback
        traceback.print_exc()
    
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
