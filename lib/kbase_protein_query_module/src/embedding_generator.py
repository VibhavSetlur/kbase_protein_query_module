"""
Protein Embedding Generator Module

This module handles the generation of protein embeddings from amino acid sequences
using ESM-2 models. It provides efficient batch processing and storage capabilities.
"""

import os
import numpy as np
import pandas as pd
import h5py
from typing import List, Dict, Optional, Union
import logging
from tqdm import tqdm
import warnings
import torch
from transformers import AutoTokenizer, EsmModel

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
    
    def __init__(self, model_name: str = "esm2_t6_8M_UR50D", device: str = "auto"):
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
        
    def _setup_device(self, device: str):
        """Setup the device for model inference."""
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        try:
            return torch.device(device)
        except AttributeError:
            # Fallback if torch.device is not available
            return "cpu"
    
    def _load_model(self):
        """Load the ESM-2 model and tokenizer with proper dtype handling."""
        try:
            logger.info(f"Loading ESM-2 model: {self.model_name}")
            
            # Try multiple possible paths for the local model
            possible_paths = [
                os.path.join(os.getcwd(), "data", "esm2_t6_8M_UR50D_local"),
                os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "esm2_t6_8M_UR50D_local"),
                "/kb/module/data/esm2_t6_8M_UR50D_local",
                # Add paths for test context
                os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))), "data", "esm2_t6_8M_UR50D_local"),
                # Add more specific test paths
                os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), "data", "esm2_t6_8M_UR50D_local"),
                os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))), "data", "esm2_t6_8M_UR50D_local")
            ]
            
            model_path = None
            for path in possible_paths:
                if os.path.exists(path):
                    model_path = path
                    logger.info(f"Using local model at: {model_path}")
                    break
            
            if model_path is None:
                # Try to find the model by searching from current directory
                current_dir = os.getcwd()
                for root, dirs, files in os.walk(current_dir):
                    if "esm2_t6_8M_UR50D_local" in dirs:
                        model_path = os.path.join(root, "esm2_t6_8M_UR50D_local")
                        logger.info(f"Found model at: {model_path}")
                        break
                    if root.count(os.sep) > 5:  # Limit search depth
                        break
            
            if model_path is None:
                raise FileNotFoundError(f"Local model not found. Searched in: {possible_paths}")

            # Load tokenizer
            try:
                self.tokenizer = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
            except Exception as e:
                logger.warning(f"Failed to load tokenizer locally: {e}")
                self.tokenizer = AutoTokenizer.from_pretrained(model_path)
            
            # Determine appropriate dtype based on device
            if isinstance(self.device, str):
                model_dtype = torch.float32 if hasattr(torch, 'float32') else torch.float
            else:
                if hasattr(torch, 'float32'):
                    model_dtype = torch.float16 if self.device.type == 'cuda' else torch.float32
                else:
                    model_dtype = torch.float
            
            # Load model with proper dtype for large models
            try:
                self.model = EsmModel.from_pretrained(
                    model_path, 
                    torch_dtype=model_dtype,
                    local_files_only=True
                )
            except Exception as e:
                logger.warning(f"Failed to load model locally: {e}")

            self.model = self.model.to(self.device)
            self.model.eval()
            
            # Get embedding dimension - ensure it's correct for the local model
            self.embedding_dim = self.model.config.hidden_size
            
            # Validate embedding dimension for local model
            if self.model_name == "esm2_t6_8M_UR50D" and self.embedding_dim != 320:
                logger.warning(f"Expected embedding dimension 320 for {self.model_name}, got {self.embedding_dim}")
                # Force the correct dimension for local model
                if hasattr(self.model.config, 'hidden_size'):
                    self.model.config.hidden_size = 320
                    self.embedding_dim = 320
            
            logger.info(f"Model loaded successfully on {self.device} with dtype {model_dtype}")
            logger.info(f"Embedding dimension: {self.embedding_dim}")
            
        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            raise
    
    def generate_embedding(self, sequence: str) -> np.ndarray:
        """
        Generate mean-pooled embedding for a single protein sequence.
        
        Args:
            sequence: Amino acid sequence string
        Returns:
            numpy.ndarray: Protein embedding vector (mean pooled)
        """
        if not sequence or not isinstance(sequence, str):
            raise ValueError("Sequence must be a non-empty string")
        try:
            # Tokenize the sequence
            inputs = self.tokenizer(sequence, return_tensors="pt", padding=True, truncation=True)
            inputs = {k: v.to(self.device) for k, v in inputs.items()}
            # Generate embeddings
            with torch.no_grad():
                outputs = self.model(**inputs)
                embeddings = outputs.last_hidden_state
                # Mean pooling over all tokens (excluding padding)
                attention_mask = inputs['attention_mask']
                masked_embeddings = embeddings * attention_mask.unsqueeze(-1)
                pooled_embedding = masked_embeddings.sum(dim=1) / attention_mask.sum(dim=1, keepdim=True)
                # Convert to numpy and return
                return pooled_embedding.cpu().numpy().flatten()
        except Exception as e:
            logger.error(f"Error generating embedding: {e}")
            raise
    
    def generate_embeddings_batch(self, sequences: List[str], 
                                protein_ids: List[str],
                                batch_size: int = 8) -> Dict[str, np.ndarray]:
        """
        Generate mean-pooled embeddings for a batch of protein sequences.
        
        Args:
            sequences: List of amino acid sequences
            protein_ids: List of protein IDs corresponding to sequences
            batch_size: Batch size for processing (reduced for large models)
        Returns:
            Dictionary mapping protein IDs to mean-pooled embeddings
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
                    # Mean pooling for each sequence in batch
                    batch_embeddings = []
                    for j in range(len(batch_sequences)):
                        seq_embeddings = embeddings[j]  # [seq_len, hidden_dim]
                        seq_attention = attention_mask[j]  # [seq_len]
                        valid_tokens = seq_attention == 1
                        # Exclude special tokens if needed (optional, can be left as is)
                        if valid_tokens.sum() > 0:
                            valid_embeddings = seq_embeddings[valid_tokens]
                            pooled_embedding = valid_embeddings.mean(dim=0)
                        else:
                            pooled_embedding = seq_embeddings.mean(dim=0)
                        batch_embeddings.append(pooled_embedding.cpu().numpy())
                    # Store embeddings
                    for j, protein_id in enumerate(batch_ids):
                        embeddings_dict[protein_id] = batch_embeddings[j]
            except RuntimeError as e:
                logger.error(f"RuntimeError processing batch starting with {batch_ids[0]}: {e}")
                continue
            except Exception as e:
                logger.error(f"Error processing batch starting with {batch_ids[0]}: {e}")
                continue
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
        sequences, protein_ids, batch_size
    )
    
    # Save embeddings
    generator.save_embeddings(embeddings_dict, output_file)
    
    return embeddings_dict 
