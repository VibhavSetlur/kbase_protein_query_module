"""
Protein Embedding Generator Module

This module handles the generation of protein embeddings from amino acid sequences
using ESM-2 models. It provides efficient batch processing and storage capabilities.
"""

import os
import torch
import numpy as np
import pandas as pd
import h5py
from transformers import AutoTokenizer, EsmModel
from typing import List, Dict, Optional, Union
import logging
from tqdm import tqdm
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*gradient checkpointing.*")
warnings.filterwarnings("ignore", message=".*position_ids.*")
warnings.filterwarnings("ignore", message=".*attention_mask.*")

logger = logging.getLogger(__name__)


class ProteinEmbeddingGenerator:
    """
    Generates protein embeddings using ESM-2 models.
    
    This class handles the conversion of protein sequences to high-dimensional
    embeddings that capture structural and functional information.
    """
    
    def __init__(self, model_name: str = "esm2_t48_15B_UR50D", device: str = "auto"):
        """
        Initialize the embedding generator.
        
        Args:
            model_name: Name of the ESM-2 model to use
            device: Device to run the model on ("auto", "cpu", or "cuda")
        """
        self.model_name = model_name
        self.device = self._setup_device(device)
        self.tokenizer = None
        self.model = None
        self.embedding_dim = None
        self._load_model()
        
    def _setup_device(self, device: str) -> torch.device:
        """Setup the device for model inference."""
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        return torch.device(device)
    
    def _load_model(self):
        """Load the ESM-2 model and tokenizer with proper dtype handling."""
        try:
            logger.info(f"Loading ESM-2 model: {self.model_name}")
            
            # Always use local path for esm2_t6_8M_UR50D
            local_model_dir = "data/esm2_t6_8M_UR50D_local"
            if self.model_name == "esm2_t6_8M_UR50D":
                model_path = local_model_dir
            else:
                model_path = f"facebook/{self.model_name}"

            # Load tokenizer
            self.tokenizer = AutoTokenizer.from_pretrained(model_path)
            
            # Determine appropriate dtype based on device
            model_dtype = torch.float16 if self.device.type == 'cuda' else torch.float32
            
            # Load model with proper dtype for large models
            self.model = EsmModel.from_pretrained(
                model_path, 
                torch_dtype=model_dtype
            )
            self.model = self.model.to(self.device)
            self.model.eval()
            
            # Get embedding dimension
            self.embedding_dim = self.model.config.hidden_size
            
            logger.info(f"Model loaded successfully on {self.device} with dtype {model_dtype}")
            logger.info(f"Embedding dimension: {self.embedding_dim}")
            
        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            raise
    
    def generate_embedding(self, sequence: str, pooling_method: str = "mean") -> np.ndarray:
        """
        Generate embedding for a single protein sequence.
        
        Args:
            sequence: Amino acid sequence string
            pooling_method: Method to pool token embeddings ("mean", "cls", "max")
            
        Returns:
            Protein embedding as numpy array
        """
        with torch.no_grad():
            # Tokenize the sequence with proper handling for large models
            tokenized_inputs = self.tokenizer(
                sequence, 
                return_tensors="pt", 
                truncation=True, 
                max_length=1024,
                add_special_tokens=True
            )
            input_ids = tokenized_inputs['input_ids'].to(self.device)
            attention_mask = tokenized_inputs['attention_mask'].to(self.device)
            
            # Generate embeddings
            outputs = self.model(input_ids, attention_mask=attention_mask)
            embeddings = outputs.last_hidden_state
            
            # Pool the embeddings (excluding special tokens)
            if pooling_method == "cls":
                # Use the [CLS] token embedding (position 0)
                pooled_embedding = embeddings[:, 0, :]
            elif pooling_method == "mean":
                # Mean pooling over all tokens (excluding special tokens)
                # ESM2 adds <cls> at 0 and <eos> at end, so slice 1:original_len+1
                original_len = len(sequence)
                # Use attention mask to determine valid tokens
                valid_tokens = attention_mask.squeeze(0) == 1
                # Exclude special tokens (first and last)
                valid_tokens[0] = False  # Exclude <cls>
                valid_tokens[-1] = False  # Exclude <eos>
                
                if valid_tokens.sum() > 0:
                    masked_embeddings = embeddings.squeeze(0)[valid_tokens, :]
                    pooled_embedding = masked_embeddings.mean(dim=0, keepdim=True)
                else:
                    # Fallback to mean of all tokens if no valid tokens
                    pooled_embedding = embeddings.mean(dim=1)
            elif pooling_method == "max":
                # Max pooling over all tokens (excluding special tokens)
                original_len = len(sequence)
                valid_tokens = attention_mask.squeeze(0) == 1
                valid_tokens[0] = False  # Exclude <cls>
                valid_tokens[-1] = False  # Exclude <eos>
                
                if valid_tokens.sum() > 0:
                    masked_embeddings = embeddings.squeeze(0)[valid_tokens, :]
                    pooled_embedding = masked_embeddings.max(dim=0, keepdim=True)[0]
                else:
                    # Fallback to max of all tokens if no valid tokens
                    pooled_embedding = embeddings.max(dim=1)[0]
            else:
                raise ValueError(f"Unknown pooling method: {pooling_method}")
            
            return pooled_embedding.cpu().numpy().flatten()
    
    def generate_embeddings_batch(self, sequences: List[str], 
                                protein_ids: List[str],
                                pooling_method: str = "mean",
                                batch_size: int = 8) -> Dict[str, np.ndarray]:
        """
        Generate embeddings for a batch of protein sequences.
        
        Args:
            sequences: List of amino acid sequences
            protein_ids: List of protein IDs corresponding to sequences
            pooling_method: Method to pool token embeddings
            batch_size: Batch size for processing (reduced for large models)
            
        Returns:
            Dictionary mapping protein IDs to embeddings
        """
        embeddings_dict = {}
        
        for i in tqdm(range(0, len(sequences), batch_size), desc="Generating embeddings"):
            batch_sequences = sequences[i:i + batch_size]
            batch_ids = protein_ids[i:i + batch_size]
            
            try:
                # Tokenize batch with proper handling for large models
                tokenized_inputs = self.tokenizer(
                    batch_sequences, 
                    return_tensors="pt", 
                    truncation=True, 
                    max_length=1024,
                    padding="longest",
                    add_special_tokens=True
                )
                input_ids = tokenized_inputs['input_ids'].to(self.device)
                attention_mask = tokenized_inputs['attention_mask'].to(self.device)
                
                # Generate embeddings
                with torch.no_grad():
                    outputs = self.model(input_ids, attention_mask=attention_mask)
                    embeddings = outputs.last_hidden_state
                    
                    # Pool embeddings for each sequence in batch
                    batch_embeddings = []
                    for j in range(len(batch_sequences)):
                        seq_embeddings = embeddings[j]  # [seq_len, hidden_dim]
                        seq_attention = attention_mask[j]  # [seq_len]
                        
                        if pooling_method == "cls":
                            # Use the [CLS] token embedding
                            pooled_embedding = seq_embeddings[0]
                        elif pooling_method == "mean":
                            # Mean pooling over valid tokens (excluding special tokens)
                            valid_tokens = seq_attention == 1
                            valid_tokens[0] = False  # Exclude <cls>
                            valid_tokens[-1] = False  # Exclude <eos>
                            
                            if valid_tokens.sum() > 0:
                                valid_embeddings = seq_embeddings[valid_tokens]
                                pooled_embedding = valid_embeddings.mean(dim=0)
                            else:
                                pooled_embedding = seq_embeddings.mean(dim=0)
                        elif pooling_method == "max":
                            # Max pooling over valid tokens (excluding special tokens)
                            valid_tokens = seq_attention == 1
                            valid_tokens[0] = False  # Exclude <cls>
                            valid_tokens[-1] = False  # Exclude <eos>
                            
                            if valid_tokens.sum() > 0:
                                valid_embeddings = seq_embeddings[valid_tokens]
                                pooled_embedding = valid_embeddings.max(dim=0)[0]
                            else:
                                pooled_embedding = seq_embeddings.max(dim=0)[0]
                        else:
                            raise ValueError(f"Unknown pooling method: {pooling_method}")
                        
                        batch_embeddings.append(pooled_embedding.cpu().numpy())
                    
                    # Store embeddings
                    for j, protein_id in enumerate(batch_ids):
                        embeddings_dict[protein_id] = batch_embeddings[j]
                        
            except RuntimeError as e:
                logger.error(f"RuntimeError processing batch starting with {batch_ids[0]}: {e}")
                # Continue with next batch
                continue
            except Exception as e:
                logger.error(f"Error processing batch starting with {batch_ids[0]}: {e}")
                # Continue with next batch
                continue
            
            # Clear GPU cache if using CUDA
            if self.device.type == 'cuda':
                torch.cuda.empty_cache()
        
        return embeddings_dict
    
    def save_embeddings(self, embeddings_dict: Dict[str, np.ndarray], 
                       output_file: str, metadata: Optional[pd.DataFrame] = None):
        """
        Save embeddings to HDF5 file.
        
        Args:
            embeddings_dict: Dictionary mapping protein IDs to embeddings
            output_file: Path to output HDF5 file
            metadata: Optional metadata DataFrame to save alongside embeddings
        """
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with h5py.File(output_file, 'w') as f:
            # Save embeddings
            protein_ids = list(embeddings_dict.keys())
            embeddings = np.array([embeddings_dict[pid] for pid in protein_ids])
            
            f.create_dataset('embeddings', data=embeddings, compression='gzip')
            f.create_dataset('protein_ids', data=protein_ids, dtype=h5py.special_dtype(vlen=str))
            
            # Save metadata if provided
            if metadata is not None:
                metadata.to_csv(output_file.replace('.h5', '_metadata.csv'), index=False)
        
        logger.info(f"Saved {len(embeddings_dict)} embeddings to {output_file}")
    
    def load_embeddings(self, input_file: str) -> tuple:
        """
        Load embeddings from HDF5 file.
        
        Args:
            input_file: Path to input HDF5 file
            
        Returns:
            Tuple of (embeddings_array, protein_ids_list)
        """
        with h5py.File(input_file, 'r') as f:
            embeddings = f['embeddings'][:]
            protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                          for pid in f['protein_ids'][:]]
        
        logger.info(f"Loaded {len(protein_ids)} embeddings from {input_file}")
        return embeddings, protein_ids
    
    def normalize_embeddings(self, embeddings: np.ndarray) -> np.ndarray:
        """
        Normalize embeddings to unit length.
        
        Args:
            embeddings: Embeddings array
            
        Returns:
            Normalized embeddings array
        """
        norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
        norms[norms == 0] = 1  # Avoid division by zero
        return embeddings / norms


def generate_embeddings_from_fasta(fasta_file: str, 
                                 output_file: str,
                                 model_name: str = "esm2_t48_15B_UR50D",
                                 pooling_method: str = "mean",
                                 batch_size: int = 8,
                                 device: str = "auto") -> Dict[str, np.ndarray]:
    """
    Generate embeddings from a FASTA file.
    
    Args:
        fasta_file: Path to input FASTA file
        output_file: Path to output HDF5 file
        model_name: ESM-2 model name
        pooling_method: Embedding pooling method
        batch_size: Batch size for processing
        device: Device to run model on
        
    Returns:
        Dictionary mapping protein IDs to embeddings
    """
    from Bio import SeqIO
    
    # Read sequences from FASTA
    sequences = []
    protein_ids = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        protein_ids.append(record.id)
    
    # Generate embeddings
    generator = ProteinEmbeddingGenerator(model_name=model_name, device=device)
    embeddings_dict = generator.generate_embeddings_batch(
        sequences, protein_ids, pooling_method, batch_size
    )
    
    # Save embeddings
    generator.save_embeddings(embeddings_dict, output_file)
    
    return embeddings_dict 

    if __name__ == "__main__":
        import tempfile
        import os

        # Create a small test FASTA file
        test_fasta_content = """>protein1
        MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ
        >protein2
        GAWKDSYVGDEAQSKRGILTLKYPIEHGIITNWDDMEK
        """
        with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".fasta") as fasta_tmp:
            fasta_tmp.write(test_fasta_content)
            fasta_tmp.flush()
            fasta_path = fasta_tmp.name

        # Output HDF5 file
        with tempfile.NamedTemporaryFile("w+b", delete=False, suffix=".h5") as h5_tmp:
            h5_path = h5_tmp.name

        try:
            print("Running embedding generation test...")
            embeddings = generate_embeddings_from_fasta(
                fasta_file=fasta_path,
                output_file=h5_path,
                model_name="esm2_t6_8M_UR50D",
                pooling_method="mean",
                batch_size=2
                )
            print("Embeddings generated for proteins:", list(embeddings.keys()))
            for pid, emb in embeddings.items():
                print(f"{pid}: shape={emb.shape}, norm={np.linalg.norm(emb):.4f}")
            # Check HDF5 file exists and contains expected keys
            import h5py
            with h5py.File(h5_path, "r") as h5f:
                print("HDF5 keys:", list(h5f.keys()))
        except Exception as e:
            print("Test failed:", e)
        finally:
            os.remove(fasta_path)
            os.remove(h5_path)
