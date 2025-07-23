"""
Similarity Indexing Module for Large-Scale Protein Search

This module provides efficient indexing and search solutions for large protein datasets
with hierarchical structure, quantization, and memory optimization.
"""

import os
import numpy as np
import pandas as pd
import pickle
import gzip
import json
from typing import Dict, List, Tuple, Optional, Union, Iterator
import logging
from pathlib import Path
import faiss
import h5py
from collections import defaultdict

logger = logging.getLogger(__name__)

class HierarchicalIndex:
    """
    Hierarchical index structure for massive protein datasets.
    
    Features:
    - Multi-level indexing (family -> subfamily -> individual)
    - Quantization for memory efficiency
    - Streaming search capabilities
    - Caching for frequently accessed families
    """
    
    def __init__(self, 
                 base_dir: str = "data/indexes",
                 index_type: str = "faiss",
                 quantization: str = "pq",
                 pq_m: int = 8,
                 pq_bits: int = 8,
                 cache_size: int = 10):
        """
        Initialize hierarchical index.
        
        Args:
            base_dir: Base directory for indexes
            index_type: Type of index ('faiss', 'annoy')
            quantization: Quantization method ('pq', 'sq', 'none')
            pq_m: Number of subvectors for Product Quantization
            pq_bits: Bits per subvector
            cache_size: Number of family indexes to cache
        """
        self.base_dir = Path(base_dir)
        self.index_type = index_type
        self.quantization = quantization
        self.pq_m = pq_m
        self.pq_bits = pq_bits
        self.cache_size = cache_size
        
        # Create directory structure
        self._create_directory_structure()
        
        # Initialize cache
        self._family_cache = {}
        self._cache_order = []
        
        # Load family mapping
        self.family_mapping = self._load_family_mapping()
    
    def _create_directory_structure(self):
        """Create directory structure for indexes."""
        directories = [
            self.base_dir,
            self.base_dir / "families",
            self.base_dir / "quantized",
            self.base_dir / "metadata",
            self.base_dir / "cache"
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
    
    def _load_family_mapping(self) -> Dict[str, str]:
        """Load family to file mapping."""
        mapping_file = self.base_dir / "family_mapping.json"
        
        if mapping_file.exists():
            with open(mapping_file, 'r') as f:
                return json.load(f)
        else:
            return {}
    
    def _save_family_mapping(self):
        """Save family to file mapping."""
        mapping_file = self.base_dir / "family_mapping.json"
        
        with open(mapping_file, 'w') as f:
            json.dump(self.family_mapping, f, indent=2)
    
    def create_family_index(self, family_id: str, embeddings: np.ndarray, protein_ids: List[str], **kwargs) -> str:
        dimension = embeddings.shape[1] * 8  # bits
        num_proteins = len(embeddings)
        if embeddings.dtype != np.uint8:
            raise ValueError("Embeddings must be np.uint8 for binary FAISS indexing.")
        nlist = min(100, num_proteins // 10) or 1
        quantizer = faiss.IndexBinaryFlat(dimension)
        index = faiss.IndexBinaryIVF(quantizer, dimension, nlist)
        index.train(embeddings)
        index.add(embeddings)
        index_file = self.base_dir / "families" / f"family_{family_id}.faissbin"
        faiss.write_index(index, str(index_file))
        metadata = {
            'protein_ids': protein_ids,
            'dimension': dimension,
            'num_proteins': num_proteins,
            'index_type': 'faiss_binary'
        }
        metadata_file = self.base_dir / "metadata" / f"family_{family_id}_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Created FAISS IVF binary index for family {family_id}: {num_proteins} proteins, {dimension} bits")
        self.family_mapping[family_id] = str(index_file)
        self._save_family_mapping()
        return str(index_file)

    def _get_cached_index(self, family_id: str):
        if family_id in self._family_cache:
            self._cache_order.remove(family_id)
            self._cache_order.append(family_id)
            return self._family_cache[family_id]
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in index mapping")
        index_file = self.family_mapping[family_id]
        metadata_file = self.base_dir / "metadata" / f"family_{family_id}_metadata.json"
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        index = faiss.read_index_binary(index_file)
        self._family_cache[family_id] = (index, metadata)
        self._cache_order.append(family_id)
        if len(self._family_cache) > self.cache_size:
            oldest_family = self._cache_order.pop(0)
            del self._family_cache[oldest_family]
        return index, metadata

    def search_family(self, family_id: str, query_embedding: np.ndarray, top_k: int = 50, **kwargs) -> Tuple[np.ndarray, List[str]]:
        index, metadata = self._get_cached_index(family_id)
        if query_embedding.dtype != np.uint8:
            raise ValueError("Query embedding must be np.uint8 for binary FAISS search.")
        D, I = index.search(query_embedding.reshape(1, -1), top_k)
        D = D[0]
        I = I[0]
        protein_ids = [metadata['protein_ids'][i] for i in I if i != -1]
        D = D[:len(protein_ids)]
        return D, protein_ids
    
    def search_all_families(self, 
                           query_embedding: np.ndarray,
                           top_k: int = 50,
                           max_families: Optional[int] = None,
                           **kwargs) -> List[Tuple[str, np.ndarray, List[str]]]:
        """
        Search across all families.
        
        Args:
            query_embedding: Query embedding vector
            top_k: Number of results per family
            max_families: Maximum number of families to search
            **kwargs: Additional search parameters
            
        Returns:
            List of (family_id, similarities, protein_ids) tuples
        """
        results = []
        family_ids = list(self.family_mapping.keys())
        
        if max_families:
            family_ids = family_ids[:max_families]
        
        for family_id in family_ids:
            try:
                similarities, protein_ids = self.search_family(
                    family_id, query_embedding, top_k, **kwargs
                )
                
                if len(protein_ids) > 0:
                    results.append((family_id, similarities, protein_ids))
                    
            except Exception as e:
                logger.warning(f"Error searching family {family_id}: {e}")
                continue
        
        return results
    
    def get_family_stats(self) -> Dict[str, Dict]:
        """Get statistics for all families."""
        stats = {}
        
        for family_id in self.family_mapping.keys():
            metadata_file = self.base_dir / "metadata" / f"family_{family_id}_metadata.json"
            
            if metadata_file.exists():
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                stats[family_id] = {
                    'num_proteins': metadata['num_proteins'],
                    'dimension': metadata['dimension'],
                    'index_type': metadata['index_type'],
                    'quantization': metadata.get('quantization', 'none')
                }
        
        return stats


class StreamingIndex:
    """
    Streaming index for memory-efficient processing of massive datasets.
    """
    
    def __init__(self, 
                 storage_dir: str = "data/streaming",
                 batch_size: int = 10000,
                 max_memory_gb: float = 4.0):
        """
        Initialize streaming index.
        
        Args:
            storage_dir: Directory for streaming data
            batch_size: Number of proteins per batch
            max_memory_gb: Maximum memory usage in GB
        """
        self.storage_dir = Path(storage_dir)
        self.batch_size = batch_size
        self.max_memory_gb = max_memory_gb
        
        self.storage_dir.mkdir(parents=True, exist_ok=True)
    
    def create_streaming_index(self, 
                              embeddings_file: str,
                              protein_ids: List[str],
                              family_assignments: Dict[str, str]) -> str:
        """
        Create streaming index structure.
        
        Args:
            embeddings_file: Path to embeddings H5 file
            protein_ids: List of protein IDs
            family_assignments: Mapping of protein ID to family ID
            
        Returns:
            Path to streaming index file
        """
        # Group proteins by family
        family_proteins = defaultdict(list)
        for protein_id in protein_ids:
            family_id = family_assignments.get(protein_id, "unknown")
            family_proteins[family_id].append(protein_id)
        
        # Create streaming structure
        streaming_file = self.storage_dir / "streaming_index.h5"
        
        with h5py.File(streaming_file, 'w') as f:
            # Store family assignments
            f.create_dataset('protein_ids', data=protein_ids, dtype=h5py.special_dtype(vlen=str))
            f.create_dataset('family_assignments', data=[family_assignments.get(pid, "unknown") for pid in protein_ids], 
                           dtype=h5py.special_dtype(vlen=str))
            
            # Store family statistics
            family_stats = {}
            for family_id, proteins in family_proteins.items():
                family_stats[family_id] = {
                    'num_proteins': len(proteins),
                    'start_idx': protein_ids.index(proteins[0]) if proteins else 0,
                    'end_idx': protein_ids.index(proteins[-1]) + 1 if proteins else 0
                }
            
            # Store family stats as JSON
            f.attrs['family_stats'] = json.dumps(family_stats)
        
        return str(streaming_file)
    
    def stream_search(self, 
                     query_embedding: np.ndarray,
                     embeddings_file: str,
                     streaming_file: str,
                     top_k: int = 50,
                     similarity_threshold: float = 0.1) -> List[Tuple[str, float]]:
        """
        Perform streaming search.
        
        Args:
            query_embedding: Query embedding vector
            embeddings_file: Path to embeddings file
            streaming_file: Path to streaming index file
            top_k: Number of results to return
            similarity_threshold: Minimum similarity threshold
            
        Returns:
            List of (protein_id, similarity) tuples
        """
        # Load streaming index
        with h5py.File(streaming_file, 'r') as f:
            protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                          for pid in f['protein_ids'][:]]
            family_assignments = [fid.decode('utf-8') if isinstance(fid, bytes) else fid 
                                 for fid in f['family_assignments'][:]]
            family_stats = json.loads(f.attrs['family_stats'])
        
        # Normalize query embedding
        query_norm = query_embedding / (np.linalg.norm(query_embedding) + 1e-8)
        
        results = []
        
        # Search each family in batches
        with h5py.File(embeddings_file, 'r') as f:
            embeddings = f['embeddings']
            
            for family_id, stats in family_stats.items():
                start_idx = stats['start_idx']
                end_idx = stats['end_idx']
                
                # Process family in batches
                for batch_start in range(start_idx, end_idx, self.batch_size):
                    batch_end = min(batch_start + self.batch_size, end_idx)
                    
                    # Load batch
                    batch_embeddings = embeddings[batch_start:batch_end]
                    
                    # Compute similarities
                    similarities = np.dot(batch_embeddings, query_norm)
                    
                    # Filter by threshold and add to results
                    for i, similarity in enumerate(similarities):
                        if similarity >= similarity_threshold:
                            protein_idx = batch_start + i
                            protein_id = protein_ids[protein_idx]
                            results.append((protein_id, float(similarity)))
        
        # Sort by similarity and return top_k
        results.sort(key=lambda x: x[1], reverse=True)
        return results[:top_k]


def create_optimized_indexes(embeddings_file: str,
                           metadata_file: str,
                           family_column: str = "Protein families",
                           output_dir: str = "data/indexes",
                           **kwargs) -> HierarchicalIndex:
    """
    Create optimized indexes for massive protein dataset.
    
    Args:
        embeddings_file: Path to embeddings H5 file
        metadata_file: Path to metadata CSV file
        family_column: Column name for family information
        output_dir: Output directory for indexes
        **kwargs: Additional parameters for index creation
        
    Returns:
        HierarchicalIndex instance
    """
    logger.info("Creating optimized indexes...")
    
    # Initialize hierarchical index
    index = HierarchicalIndex(base_dir=output_dir, **kwargs)
    
    # Load data
    with h5py.File(embeddings_file, 'r') as f:
        embeddings = f['embeddings'][:]
        protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                      for pid in f['protein_ids'][:]]
    
    metadata = pd.read_csv(metadata_file, index_col=0)
    
    # Group by family
    if family_column in metadata.columns:
        family_groups = metadata.groupby(family_column)
    else:
        # Create artificial families
        family_groups = _create_artificial_family_groups(protein_ids, 100000)
    
    # Create indexes for each family
    for family_id, group_metadata in family_groups:
        family_protein_ids = list(group_metadata.index)
        
        # Get embeddings for this family
        family_indices = [protein_ids.index(pid) for pid in family_protein_ids if pid in protein_ids]
        family_embeddings = embeddings[family_indices]
        
        # Create family index
        index.create_family_index(
            str(family_id),
            family_embeddings,
            family_protein_ids
        )
    
    logger.info(f"Created optimized indexes in {output_dir}")
    return index


def _create_artificial_family_groups(protein_ids: List[str], max_family_size: int) -> Iterator[Tuple[int, pd.DataFrame]]:
    """Create artificial family groups when no family information is available."""
    for i in range(0, len(protein_ids), max_family_size):
        family_proteins = protein_ids[i:i + max_family_size]
        family_metadata = pd.DataFrame(index=family_proteins)
        yield i // max_family_size, family_metadata 
        