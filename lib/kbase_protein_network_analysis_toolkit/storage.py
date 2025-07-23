"""
Storage Module for Large-Scale Protein Data

This module provides efficient storage solutions for massive protein datasets
(250M+ proteins) with hierarchical organization, chunking, and compression.
"""

import os
import h5py
import numpy as np
import pandas as pd
import pickle
import gzip
import zlib
from typing import Dict, List, Tuple, Optional, Union, Iterator
import logging
from collections import defaultdict
import json
from pathlib import Path

logger = logging.getLogger(__name__)

class ProteinStorage:
    """
    Storage system for massive protein datasets.
    
    Features:
    - Hierarchical HDF5 structure by protein families
    - Chunked storage for memory-efficient access
    - Compression with multiple algorithms
    - Metadata indexing for fast lookups
    - Streaming access for large datasets
    """
    
    def __init__(self, 
                 base_dir: str = "data",
                 chunk_size: int = 1000,
                 compression: str = "gzip",
                 compression_opts: int = 6,
                 max_family_size: int = 100000):
        """
        Initialize storage.
        
        Args:
            base_dir: Base directory for data storage
            chunk_size: Number of proteins per chunk
            compression: Compression algorithm ('gzip', 'lzf', 'szip')
            compression_opts: Compression level (1-9 for gzip)
            max_family_size: Maximum proteins per family before splitting
        """
        self.base_dir = Path(base_dir)
        self.chunk_size = chunk_size
        self.compression = compression
        self.compression_opts = compression_opts
        self.max_family_size = max_family_size
        
        # Create directory structure
        self._create_directory_structure()
        
        # Storage paths
        self.family_dir = self.base_dir / "families"
        self.index_dir = self.base_dir / "indexes"
        self.metadata_dir = self.base_dir / "metadata"
        self.cache_dir = self.base_dir / "cache"
        
    def _create_directory_structure(self):
        """Create optimized directory structure."""
        directories = [
            self.base_dir / "families",
            self.base_dir / "indexes", 
            self.base_dir / "metadata",
            self.base_dir / "cache",
            self.base_dir / "temp"
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
    
    def save_html_file(self, html_content: str, filename: str, subdir: Optional[str] = None) -> str:
        """
        Save an HTML file to the storage-managed directory structure.

        Args:
            html_content: The HTML content to write
            filename: The name of the HTML file (should end with .html)
            subdir: Optional subdirectory under the base_dir to store the file (e.g., 'results' or 'html_reports')
        Returns:
            The full path to the saved HTML file
        """
        if subdir:
            output_dir = self.base_dir / subdir
        else:
            output_dir = self.base_dir / "results"
        output_dir.mkdir(parents=True, exist_ok=True)
        html_path = output_dir / filename
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        logger.info(f"Saved HTML file to {html_path}")
        return str(html_path)

    def load_html_file(self, filename: str, subdir: Optional[str] = None) -> str:
        """
        Load an HTML file from the storage-managed directory structure.
        Args:
            filename: The name of the HTML file
            subdir: Optional subdirectory under the base_dir to load from
        Returns:
            The HTML content as a string
        """
        if subdir:
            html_path = self.base_dir / subdir / filename
        else:
            html_path = self.base_dir / "results" / filename
        if not html_path.exists():
            raise FileNotFoundError(f"HTML file not found: {html_path}")
        with open(html_path, "r", encoding="utf-8") as f:
            return f.read()
    
    def store_family_embeddings(self, 
                               family_id: str,
                               embeddings: np.ndarray,
                               protein_ids: List[str],
                               metadata: Optional[pd.DataFrame] = None) -> str:
        """
        Store embeddings for a protein family with chunking.
        Args:
            family_id: Family identifier
            embeddings: Embedding array (N x D//8, np.uint8)
            protein_ids: List of protein IDs
            metadata: Optional metadata DataFrame
        Returns:
            Path to stored family file
        """
        if embeddings.dtype != np.uint8:
            raise ValueError("Embeddings must be np.uint8 for binary storage.")
        family_file = self.family_dir / f"family_{family_id}.h5"
        embedding_dim = embeddings.shape[1]
        num_proteins = embeddings.shape[0]
        optimal_chunk_size = min(self.chunk_size, max(1, min(num_proteins, 10000 // embedding_dim)))
        with h5py.File(family_file, 'w') as f:
            f.create_dataset(
                'embeddings',
                data=embeddings,
                chunks=(min(optimal_chunk_size, num_proteins), embedding_dim),
                compression=self.compression,
                compression_opts=self.compression_opts,
                shuffle=True
            )
            f.create_dataset(
                'protein_ids',
                data=protein_ids,
                dtype=h5py.special_dtype(vlen=str),
                chunks=(min(optimal_chunk_size, num_proteins),),
                compression=self.compression,
                compression_opts=self.compression_opts
            )
            if metadata is not None:
                metadata_file = self.metadata_dir / f"family_{family_id}_metadata.parquet"
                metadata.to_parquet(metadata_file, compression='gzip')
                f.attrs['metadata_file'] = str(metadata_file)
            f.attrs['num_proteins'] = len(protein_ids)
            f.attrs['embedding_dim'] = embedding_dim
            f.attrs['chunk_size'] = optimal_chunk_size
            f.attrs['compression'] = self.compression
        logger.info(f"Stored family {family_id}: {len(protein_ids)} proteins, {embedding_dim} bytes, {optimal_chunk_size} chunk size")
        return str(family_file)
    
    def load_family_embeddings(self, 
                              family_id: str,
                              start_idx: Optional[int] = None,
                              end_idx: Optional[int] = None) -> Tuple[np.ndarray, List[str]]:
        """
        Load embeddings for a protein family with optional slicing.
        
        Args:
            family_id: Family identifier
            start_idx: Start index for slicing
            end_idx: End index for slicing
            
        Returns:
            Tuple of (embeddings, protein_ids)
        """
        family_file = self.family_dir / f"family_{family_id}.h5"
        
        if not family_file.exists():
            raise FileNotFoundError(f"Family file not found: {family_file}")
        
        with h5py.File(family_file, 'r') as f:
            # Load with slicing if specified
            if start_idx is not None or end_idx is not None:
                embeddings = f['embeddings'][start_idx:end_idx]
                protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][start_idx:end_idx]]
            else:
                embeddings = f['embeddings'][:]
                protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][:]]
        
        return embeddings, protein_ids
    
    def create_family_index(self, family_id: str, **kwargs) -> str:
        """
        Create and store FAISS IVF binary similarity index for a family.
        Args:
            family_id: Family identifier
        Returns:
            Path to stored index file
        """
        import faiss
        embeddings, protein_ids = self.load_family_embeddings(family_id)
        # Convert embeddings to binary (uint8) if not already
        if embeddings.dtype != np.uint8:
            raise ValueError("Embeddings must be np.uint8 for binary FAISS indexing.")
        dimension = embeddings.shape[1] * 8  # D in bits
        nlist = min(100, len(embeddings) // 10) or 1
        quantizer = faiss.IndexBinaryFlat(dimension)
        index = faiss.IndexBinaryIVF(quantizer, dimension, nlist)
        index.train(embeddings)
        index.add(embeddings)
        index_file = self.index_dir / f"family_{family_id}.faissbin"
        faiss.write_index(index, str(index_file))
        metadata = {
            'protein_ids': protein_ids,
            'index_type': 'faiss_binary',
            'dimension': dimension,
            'num_proteins': len(protein_ids)
        }
        metadata_file = self.index_dir / f"family_{family_id}_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Created FAISS IVF binary index for family {family_id}: {len(protein_ids)} proteins, {dimension} bits")
        return str(index_file)
    
    def get_family_list(self) -> List[str]:
        """Get list of all available families."""
        family_files = list(self.family_dir.glob("family_*.h5"))
        return [f.stem.replace("family_", "") for f in family_files]
    
    def get_family_stats(self) -> Dict[str, Dict]:
        """Get statistics for all families."""
        stats = {}
        
        for family_id in self.get_family_list():
            family_file = self.family_dir / f"family_{family_id}.h5"
            
            with h5py.File(family_file, 'r') as f:
                stats[family_id] = {
                    'num_proteins': f.attrs['num_proteins'],
                    'embedding_dim': f.attrs['embedding_dim'],
                    'chunk_size': f.attrs['chunk_size'],
                    'compression': f.attrs['compression'],
                    'file_size_mb': family_file.stat().st_size / (1024 * 1024)
                }
        
        return stats
    
    def stream_family_embeddings(self, 
                                family_id: str,
                                batch_size: int = 1000) -> Iterator[Tuple[np.ndarray, List[str]]]:
        """
        Stream embeddings for a family in batches.
        
        Args:
            family_id: Family identifier
            batch_size: Number of proteins per batch
            
        Yields:
            Tuples of (embeddings_batch, protein_ids_batch)
        """
        family_file = self.family_dir / f"family_{family_id}.h5"
        
        with h5py.File(family_file, 'r') as f:
            num_proteins = f.attrs['num_proteins']
            
            for start_idx in range(0, num_proteins, batch_size):
                end_idx = min(start_idx + batch_size, num_proteins)
                
                embeddings = f['embeddings'][start_idx:end_idx]
                protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][start_idx:end_idx]]
                
                yield embeddings, protein_ids

    def load_metadata(self, family_id: str) -> pd.DataFrame:
        """
        Load metadata for a given family using CompressedMetadataStorage.
        Args:
            family_id: Family identifier
        Returns:
            Metadata DataFrame
        """
        meta_storage = CompressedMetadataStorage(metadata_dir=str(self.metadata_dir))
        return meta_storage.load_metadata(family_id=family_id)


class CompressedMetadataStorage:
    """
    Efficient metadata storage with compression and indexing.
    """
    
    def __init__(self, metadata_dir: str = "data/metadata"):
        self.metadata_dir = Path(metadata_dir)
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
    
    def store_metadata(self, 
                      metadata: pd.DataFrame,
                      family_id: Optional[str] = None,
                      compression: str = "gzip") -> str:
        """
        Store metadata with compression.
        
        Args:
            metadata: Metadata DataFrame
            family_id: Optional family identifier
            compression: Compression method
            
        Returns:
            Path to stored metadata file
        """
        if family_id:
            filename = f"family_{family_id}_metadata.parquet"
        else:
            filename = "global_metadata.parquet"
        
        filepath = self.metadata_dir / filename
        
        # Store as Parquet with compression
        metadata.to_parquet(filepath, compression=compression)
        
        # Create index for fast lookups
        self._create_metadata_index(metadata, family_id)
        
        return str(filepath)
    
    def _create_metadata_index(self, metadata: pd.DataFrame, family_id: Optional[str] = None):
        """Create index for fast metadata lookups."""
        # Create protein ID to row index mapping
        protein_to_idx = {pid: idx for idx, pid in enumerate(metadata.index)}
        
        index_data = {
            'protein_to_idx': protein_to_idx,
            'columns': list(metadata.columns),
            'shape': metadata.shape
        }
        
        if family_id:
            index_file = self.metadata_dir / f"family_{family_id}_index.pkl"
        else:
            index_file = self.metadata_dir / "global_index.pkl"
        
        with gzip.open(index_file, 'wb') as f:
            pickle.dump(index_data, f)
    
    def load_metadata(self, 
                     family_id: Optional[str] = None,
                     protein_ids: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Load metadata with optional filtering.
        
        Args:
            family_id: Optional family identifier
            protein_ids: Optional list of protein IDs to filter
            
        Returns:
            Metadata DataFrame
        """
        if family_id:
            filepath = self.metadata_dir / f"family_{family_id}_metadata.parquet"
        else:
            filepath = self.metadata_dir / "global_metadata.parquet"
        
        if not filepath.exists():
            raise FileNotFoundError(f"Metadata file not found: {filepath}")
        
        metadata = pd.read_parquet(filepath)
        
        if protein_ids is not None:
            metadata = metadata.loc[metadata.index.intersection(protein_ids)]
        
        return metadata


class MemoryEfficientLoader:
    """
    Memory-efficient loading for large datasets.
    """
    
    def __init__(self, storage: ProteinStorage):
        self.storage = storage
    
    def load_families_batch(self, 
                           family_ids: List[str],
                           max_memory_gb: float = 4.0) -> Iterator[Tuple[str, np.ndarray, List[str]]]:
        """
        Load multiple families in memory-efficient batches.
        
        Args:
            family_ids: List of family IDs to load
            max_memory_gb: Maximum memory usage in GB
            
        Yields:
            Tuples of (family_id, embeddings, protein_ids)
        """
        max_memory_bytes = max_memory_gb * 1024**3
        
        for family_id in family_ids:
            # Estimate memory usage
            stats = self.storage.get_family_stats()[family_id]
            estimated_memory = stats['num_proteins'] * stats['embedding_dim'] * 4  # float32
            
            if estimated_memory > max_memory_bytes:
                # Load in chunks
                for embeddings, protein_ids in self.storage.stream_family_embeddings(family_id):
                    yield family_id, embeddings, protein_ids
            else:
                # Load entire family
                embeddings, protein_ids = self.storage.load_family_embeddings(family_id)
                yield family_id, embeddings, protein_ids


def create_storage_structure(embeddings_file: str,
                         metadata_file: str,
                         output_dir: str = "data",
                         family_column: str = "Protein families",
                         max_family_size: int = 100000) -> ProteinStorage:
    """
    Convert existing data to storage structure.
    
    Args:
        embeddings_file: Path to existing embeddings H5 file
        metadata_file: Path to existing metadata CSV file
        output_dir: Output directory for optimized storage
        family_column: Column name for family information
        max_family_size: Maximum proteins per family
        
    Returns:
        ProteinStorage instance
    """
    logger.info("Creating optimized storage structure...")
    
    # Initialize storage
    storage = ProteinStorage(base_dir=output_dir, max_family_size=max_family_size)
    
    # Load existing data
    with h5py.File(embeddings_file, 'r') as f:
        embeddings = f['embeddings'][:]
        protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                      for pid in f['protein_ids'][:]]
    
    metadata = pd.read_csv(metadata_file, index_col=0)
    
    # Group by family
    if family_column in metadata.columns:
        family_groups = metadata.groupby(family_column)
    else:
        # Create artificial families based on size
        family_groups = _create_artificial_families(protein_ids, max_family_size)
    
    # Store each family
    for family_id, group_metadata in family_groups:
        family_protein_ids = list(group_metadata.index)
        
        # Get embeddings for this family
        family_indices = [protein_ids.index(pid) for pid in family_protein_ids if pid in protein_ids]
        family_embeddings = embeddings[family_indices]
        
        # Store family
        storage.store_family_embeddings(
            str(family_id),
            family_embeddings,
            family_protein_ids,
            group_metadata
        )
    
    logger.info(f"Optimized storage created in {output_dir}")
    return storage


def _create_artificial_families(protein_ids: List[str], max_family_size: int) -> Iterator[Tuple[int, pd.DataFrame]]:
    """Create artificial families when no family information is available."""
    for i in range(0, len(protein_ids), max_family_size):
        family_proteins = protein_ids[i:i + max_family_size]
        family_metadata = pd.DataFrame(index=family_proteins)
        yield i // max_family_size, family_metadata 