"""
FAISS Index Manager Module

This module provides a utility class `FaissIndexManager` for creating, managing,
and searching FAISS indexes for protein embeddings. It supports both float and
binary indexes and integrates with the ProteinEmbeddingGenerator and UniProt API
to build indexes from various inputs (embeddings, sequences, UniProt IDs).
"""

import os
import logging
import json
import pickle
import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Union, Tuple, Any

# Try importing faiss
try:
    import faiss
    HAS_FAISS = True
except ImportError:
    HAS_FAISS = False

# Import ProteinEmbeddingGenerator
try:
    from ..embeddings.generator import ProteinEmbeddingGenerator
except ImportError:
    # Fallback for when running as a script from different location
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embeddings'))
    from generator import ProteinEmbeddingGenerator

# Import UniProt API
try:
    from ..uniprot import api as uniprot_api
except ImportError:
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'uniprot'))
    import api as uniprot_api

logger = logging.getLogger(__name__)

class FaissIndexManager:
    """
    Manager for FAISS indexes.
    
    Supports:
    - Creating new indexes (Flat, IVF, HNSW, Binary)
    - Adding vectors directly
    - Adding vectors from protein sequences
    - Adding vectors from UniProt IDs
    - Searching indexes
    - Saving and loading indexes
    """

    def __init__(
        self, 
        dimension: int, 
        index_type: str = "Flat", 
        metric: str = "L2", 
        use_gpu: bool = False,
        nlist: int = 100,  # For IVF indexes
        embedding_generator: Optional[ProteinEmbeddingGenerator] = None
    ):
        """
        Initialize the FAISS Index Manager.

        Args:
            dimension: Dimension of the vectors.
            index_type: Type of index ("Flat", "IVF", "HNSW", "BinaryFlat", "BinaryIVF").
            metric: Distance metric ("L2", "IP" (Inner Product)). For Binary, usually Hamming.
            use_gpu: Whether to use GPU resources if available.
            nlist: Number of clusters for IVF indexes.
            embedding_generator: Optional instance of ProteinEmbeddingGenerator.
        """
        if not HAS_FAISS:
            raise ImportError("faiss package is not installed.")

        self.dimension = dimension
        self.index_type = index_type
        self.metric = metric
        self.use_gpu = use_gpu
        self.nlist = nlist
        
        self.index = None
        self.id_map = {}  # Maps internal FAISS ID to external ID (e.g., protein ID)
        self.is_binary = "Binary" in index_type
        
        self.embedding_generator = embedding_generator

        self._create_index()

    def _create_index(self):
        """Creates the FAISS index based on configuration."""
        
        # Metric mapping
        if self.metric == "L2":
            faiss_metric = faiss.METRIC_L2
        elif self.metric == "IP":
            faiss_metric = faiss.METRIC_INNER_PRODUCT
        else:
            # Default for float
            faiss_metric = faiss.METRIC_L2
            
        if self.is_binary:
            # Binary Index Creation
            if self.index_type == "BinaryFlat":
                self.index = faiss.IndexBinaryFlat(self.dimension)
            elif self.index_type == "BinaryIVF":
                quantizer = faiss.IndexBinaryFlat(self.dimension)
                self.index = faiss.IndexBinaryIVF(quantizer, self.dimension, self.nlist)
            else:
                raise ValueError(f"Unsupported binary index type: {self.index_type}")
        else:
            # Float Index Creation
            if self.index_type == "Flat":
                self.index = faiss.IndexFlat(self.dimension, faiss_metric)
            elif self.index_type == "IVF":
                quantizer = faiss.IndexFlat(self.dimension, faiss_metric)
                self.index = faiss.IndexIVFFlat(quantizer, self.dimension, self.nlist, faiss_metric)
            elif self.index_type == "HNSW":
                self.index = faiss.IndexHNSWFlat(self.dimension, 32, faiss_metric) # 32 is typical M
            else:
                # Try using index_factory for complex strings
                try:
                    self.index = faiss.index_factory(self.dimension, self.index_type, faiss_metric)
                except Exception as e:
                    raise ValueError(f"Unsupported index type or factory string: {self.index_type}. Error: {e}")

        if self.use_gpu:
            self._to_gpu()

    def _to_gpu(self):
        """Moves index to GPU if available."""
        try:
            res = faiss.StandardGpuResources()
            self.index = faiss.index_cpu_to_gpu(res, 0, self.index)
        except Exception as e:
            logger.warning(f"Failed to move index to GPU: {e}. Continuing with CPU.")

    def train(self, vectors: np.ndarray):
        """
        Trains the index if required (e.g., for IVF).
        
        Args:
            vectors: Numpy array of vectors for training.
        """
        if not self.index.is_trained:
            logger.info("Training index...")
            if self.is_binary:
                 if vectors.dtype != np.uint8:
                     raise ValueError("Binary index requires uint8 vectors")
            else:
                 if vectors.dtype != np.float32:
                     vectors = vectors.astype(np.float32)
            
            self.index.train(vectors)
            logger.info("Index trained.")

    def add_vectors(self, vectors: np.ndarray, ids: List[str] = None):
        """
        Adds vectors to the index.

        Args:
            vectors: Numpy array of vectors.
            ids: List of external IDs corresponding to the vectors.
                 If provided, must match length of vectors.
        """
        if not self.index.is_trained:
            self.train(vectors)

        num_vectors = vectors.shape[0]
        if ids:
            if len(ids) != num_vectors:
                raise ValueError("Length of ids must match number of vectors.")
            
            start_id = self.index.ntotal
            # Update ID map
            for i, ext_id in enumerate(ids):
                self.id_map[start_id + i] = ext_id
        
        # Ensure correct type
        if self.is_binary:
            if vectors.dtype != np.uint8:
                vectors = vectors.astype(np.uint8)
        else:
            if vectors.dtype != np.float32:
                vectors = vectors.astype(np.float32)

        self.index.add(vectors)
        logger.info(f"Added {num_vectors} vectors. Total: {self.index.ntotal}")

    def add_from_sequences(self, sequences: List[str], ids: List[str]):
        """
        Generates embeddings for sequences and adds them to the index.
        
        Args:
            sequences: List of protein sequences.
            ids: List of corresponding IDs.
        """
        if not self.embedding_generator:
            raise ValueError("Embedding generator not initialized.")
            
        if len(sequences) != len(ids):
            raise ValueError("Sequences and IDs must have same length.")
            
        vectors = []
        valid_ids = []
        
        for seq, pid in zip(sequences, ids):
            try:
                emb = self.embedding_generator.generate_embedding(seq)
                vectors.append(emb)
                valid_ids.append(pid)
            except Exception as e:
                logger.warning(f"Failed to generate embedding for {pid}: {e}")
                
        if vectors:
            vec_array = np.vstack(vectors)
            self.add_vectors(vec_array, valid_ids)

    def add_from_uniprot_ids(self, uniprot_ids: List[str]):
        """
        Fetches sequences for UniProt IDs, generates embeddings, and adds to index.
        
        Args:
            uniprot_ids: List of UniProt IDs.
        """
        logger.info(f"Fetching sequences for {len(uniprot_ids)} UniProt IDs...")
        id_seq_map = uniprot_api.fetch_sequences(uniprot_ids)
        
        if not id_seq_map:
            logger.warning("No sequences found for provided UniProt IDs.")
            return
            
        ids = list(id_seq_map.keys())
        sequences = list(id_seq_map.values())
        
        self.add_from_sequences(sequences, ids)

    def search(self, query_vectors: np.ndarray, k: int = 5) -> Tuple[np.ndarray, np.ndarray, List[List[str]]]:
        """
        Searches the index.

        Args:
            query_vectors: Numpy array of query vectors.
            k: Number of nearest neighbors to retrieve.

        Returns:
            distances: Distances to the nearest neighbors.
            indices: Internal FAISS indices of the nearest neighbors.
            external_ids: List of lists of external IDs (if available).
        """
        if self.is_binary:
            if query_vectors.dtype != np.uint8:
                query_vectors = query_vectors.astype(np.uint8)
        else:
            if query_vectors.dtype != np.float32:
                query_vectors = query_vectors.astype(np.float32)

        distances, indices = self.index.search(query_vectors, k)
        
        # Map back to external IDs
        external_ids = []
        for row_indices in indices:
            row_ext_ids = []
            for idx in row_indices:
                if idx != -1 and idx in self.id_map:
                    row_ext_ids.append(self.id_map[idx])
                else:
                    row_ext_ids.append(None) # Or str(idx)
            external_ids.append(row_ext_ids)

        return distances, indices, external_ids

    def save_index(self, filepath: str):
        """
        Saves the index and ID map to disk.
        
        Args:
            filepath: Base path to save the index. 
                      Will create {filepath}.index and {filepath}.ids.json
        """
        # Save FAISS index
        index_path = f"{filepath}.index"
        if self.use_gpu:
            # Move to CPU before saving
            index_cpu = faiss.index_gpu_to_cpu(self.index)
            if self.is_binary:
                faiss.write_index_binary(index_cpu, index_path)
            else:
                faiss.write_index(index_cpu, index_path)
        else:
            if self.is_binary:
                faiss.write_index_binary(self.index, index_path)
            else:
                faiss.write_index(self.index, index_path)
        
        # Save ID map
        ids_path = f"{filepath}.ids.pickle"
        with open(ids_path, 'wb') as f:
            pickle.dump(self.id_map, f)
            
        # Save metadata
        meta_path = f"{filepath}.meta.json"
        meta = {
            "dimension": self.dimension,
            "index_type": self.index_type,
            "metric": self.metric,
            "is_binary": self.is_binary
        }
        with open(meta_path, 'w') as f:
            json.dump(meta, f, indent=2)
            
        logger.info(f"Index saved to {filepath}.*")

    @classmethod
    def load_index(cls, filepath: str, use_gpu: bool = False, embedding_generator: Optional[ProteinEmbeddingGenerator] = None):
        """
        Loads an index from disk.
        
        Args:
            filepath: Base path of the index (without extensions).
            use_gpu: Whether to load onto GPU.
            embedding_generator: Optional generator to attach.
            
        Returns:
            FaissIndexManager instance.
        """
        meta_path = f"{filepath}.meta.json"
        with open(meta_path, 'r') as f:
            meta = json.load(f)
            
        manager = cls(
            dimension=meta["dimension"],
            index_type=meta["index_type"],
            metric=meta["metric"],
            use_gpu=use_gpu,
            embedding_generator=embedding_generator
        )
        
        # Load FAISS index
        index_path = f"{filepath}.index"
        if manager.is_binary:
            manager.index = faiss.read_index_binary(index_path)
        else:
            manager.index = faiss.read_index(index_path)
            
        if use_gpu:
            manager._to_gpu()
            
        # Load ID map
        ids_path = f"{filepath}.ids.pickle"
        if os.path.exists(ids_path):
            with open(ids_path, 'rb') as f:
                manager.id_map = pickle.load(f)
                
        return manager

    def build_from_files(
        self, 
        file_paths: List[str], 
        id_col: str, 
        seq_col: str, 
        batch_size: int = 32,
        delimiter: str = "\t"
    ):
        """
        Builds the index from protein sequence files.

        Args:
            file_paths: List of paths to CSV/TSV files.
            id_col: Name of the column containing IDs.
            seq_col: Name of the column containing sequences.
            batch_size: Batch size for processing.
            delimiter: Delimiter for the files (default tab).
        """
        if not self.embedding_generator:
            raise ValueError("Embedding generator not initialized.")

        for file_path in file_paths:
            logger.info(f"Processing file: {file_path}")
            try:
                # Determine delimiter if not explicitly set or if it's CSV
                sep = delimiter
                if file_path.endswith('.csv'):
                    sep = ','
                
                # Read file in chunks
                for chunk in pd.read_csv(file_path, sep=sep, chunksize=batch_size, usecols=[id_col, seq_col]):
                    chunk_ids = chunk[id_col].astype(str).tolist()
                    chunk_seqs = chunk[seq_col].astype(str).tolist()
                    
                    self.add_from_sequences(chunk_seqs, chunk_ids)
                        
            except Exception as e:
                logger.error(f"Error processing file {file_path}: {e}")

def main():
    """CLI Entry point for testing or simple usage."""
    import argparse
    parser = argparse.ArgumentParser(description="FAISS Index Manager")
    parser.add_argument("--mode", choices=["create", "search"], required=True)
    parser.add_argument("--index_path", required=True, help="Path to save/load index")
    parser.add_argument("--files", nargs="+", help="Input files for creation")
    parser.add_argument("--id_col", default="id")
    parser.add_argument("--seq_col", default="sequence")
    parser.add_argument("--dim", type=int, default=320, help="Embedding dimension (default 320 for esm2_t6_8M)")
    parser.add_argument("--uniprot_ids", nargs="+", help="UniProt IDs to add")
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    
    if args.mode == "create":
        # Initialize generator
        gen = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D")
        
        manager = FaissIndexManager(dimension=args.dim, index_type="Flat", metric="L2", embedding_generator=gen)
        
        if args.files:
            manager.build_from_files(args.files, args.id_col, args.seq_col)
            
        if args.uniprot_ids:
            manager.add_from_uniprot_ids(args.uniprot_ids)
            
        manager.save_index(args.index_path)
        print(f"Index created at {args.index_path}")
        
    elif args.mode == "search":
        # Placeholder for search CLI
        print("Search CLI not fully implemented yet.")

if __name__ == "__main__":
    main()
