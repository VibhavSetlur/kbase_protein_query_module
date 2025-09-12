"""
Storage Module for Large-Scale Protein Data

This module provides efficient storage solutions for massive protein datasets
(250M+ proteins) with hierarchical organization, chunking, and compression.
Optimized for loading only mapped families with advanced hybrid FAISS indexing.
"""

import os
import h5py
import numpy as np
import pandas as pd
import pickle
import gzip
import zlib
from typing import Dict, List, Tuple, Optional, Union, Iterator, Any
import logging
from collections import defaultdict
import json
from pathlib import Path

logger = logging.getLogger(__name__)

class ProteinStorage:
    """
    Advanced storage system for massive protein datasets with hybrid FAISS indexing.
    
    Features:
    - Selective loading of mapped families only (memory efficient)
    - Hybrid FAISS indexing: Binary for families, Float for within-family
    - HNSW graphs for approximate nearest neighbor search
    - Adaptive centroid selection based on family size
    - Multi-level indexing for maximum efficiency
    - Hierarchical HDF5 structure by protein families
    - Chunked storage for memory-efficient access
    - Compression with multiple algorithms
    - Metadata indexing for fast lookups
    - Streaming access for large datasets
    - Memory usage optimized for 250M protein sequences
    """
    
    def __init__(self, 
                 base_dir: str = "data",
                 chunk_size: int = 1000,
                 compression: str = "gzip",
                 compression_opts: int = 6,
                 max_family_size: int = 100000,
                 memory_limit_gb: float = 4.0,
                 use_hnsw: bool = True,
                 hnsw_m: int = 16,
                 hnsw_ef_construction: int = 200):
        """
        Initialize storage with advanced indexing.
        
        Args:
            base_dir: Base directory for data storage
            chunk_size: Number of proteins per chunk
            compression: Compression algorithm ('gzip', 'lzf', 'szip')
            compression_opts: Compression level (1-9 for gzip)
            max_family_size: Maximum proteins per family before splitting
            memory_limit_gb: Memory limit in GB for loading families
            use_hnsw: Whether to use HNSW graphs for approximate search
            hnsw_m: Number of connections per layer in HNSW
            hnsw_ef_construction: Search depth during HNSW construction
        """
        self.base_dir = Path(base_dir)
        self.chunk_size = chunk_size
        self.compression = compression
        self.compression_opts = compression_opts
        self.max_family_size = max_family_size
        self.memory_limit_gb = memory_limit_gb
        self.memory_limit_bytes = memory_limit_gb * 1024**3
        self.use_hnsw = use_hnsw
        self.hnsw_m = hnsw_m
        self.hnsw_ef_construction = hnsw_ef_construction
        
        # Create directory structure
        self._create_directory_structure()
        
        # Storage paths
        self.family_dir = self.base_dir / "families"
        self.index_dir = self.base_dir / "indexes"
        self.metadata_dir = self.base_dir / "metadata"
        self.cache_dir = self.base_dir / "cache"
        
        # Load family mapping for selective loading
        self._load_family_mapping()
        
    def _load_family_mapping(self):
        """Load family mapping for selective loading."""
        mapping_file = self.base_dir / "family_mapping.json"
        if mapping_file.exists():
            with open(mapping_file, 'r') as f:
                self.family_mapping = json.load(f)
            logger.info(f"Loaded family mapping with {len(self.family_mapping)} families")
            return
        # No mapping file: infer from existing H5 files under data/families to avoid hard failures in tests
        self.family_mapping = {}
        try:
            families_dir = self.base_dir / "families"
            if families_dir.exists():
                for fname in os.listdir(families_dir):
                    if fname.endswith('.h5'):
                        file_safe_id = fname[:-3]
                        family_file = families_dir / fname
                        try:
                            with h5py.File(family_file, 'r') as f:
                                num_proteins = int(f.attrs.get('num_proteins', len(f['protein_ids'][:])))
                                embedding_dim = int(f.attrs.get('embedding_dim', f['embeddings'][:].shape[1]))
                        except Exception:
                            # Best-effort defaults
                            num_proteins = 0
                            embedding_dim = 0
                        # Use original id as human family id if possible
                        family_id = file_safe_id
                        self.family_mapping[family_id] = {
                            'file_safe_id': file_safe_id,
                            'num_proteins': num_proteins,
                            'embedding_dim': embedding_dim,
                            'file_size_mb': family_file.stat().st_size / (1024 * 1024)
                        }
                if self.family_mapping:
                    self._save_family_mapping()
                    logger.info(f"Inferred family mapping from existing files: {len(self.family_mapping)} families")
                else:
                    logger.info("No family mapping found, will create on first use")
            else:
                logger.info("No families directory found, will create on first use")
        except Exception as e:
            logger.warning(f"Failed to infer family mapping: {e}")
            logger.info("No family mapping found, will create on first use")
        
    def _save_family_mapping(self):
        """Save family mapping for future selective loading."""
        mapping_file = self.base_dir / "family_mapping.json"
        with open(mapping_file, 'w') as f:
            json.dump(self.family_mapping, f, indent=2)
        
    def _make_file_safe(self, family_id: str) -> str:
        """Convert any family ID to a file-safe format."""
        import re
        file_safe = re.sub(r'[^a-zA-Z0-9]', '_', family_id)
        file_safe = re.sub(r'_+', '_', file_safe)
        file_safe = file_safe.strip('_')
        if not file_safe:
            file_safe = 'unknown_family'
        return file_safe
        
    def _create_directory_structure(self):
        """Create optimized directory structure."""
        directories = [
            self.base_dir / "families",
            self.base_dir / "indexes", 
            self.base_dir / "indexes" / "families",
            self.base_dir / "indexes" / "centroids",
            self.base_dir / "indexes" / "metadata",
            self.base_dir / "indexes" / "hnsw",
            self.base_dir / "metadata",
            self.base_dir / "cache",
            self.base_dir / "temp"
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
    
    def _calculate_optimal_nlist(self, num_proteins: int, embedding_dim: int) -> int:
        """
        Calculate optimal number of centroids based on family size and dimension.
        Based on FAISS best practices for 250M protein scale.
        """
        # For large families (>100K proteins), use more centroids for diversity
        if num_proteins > 100000:
            # Use 1 centroid per 1000-2000 proteins for large families
            nlist = max(100, min(1000, num_proteins // 1500))
        elif num_proteins > 10000:
            # Use 1 centroid per 500-1000 proteins for medium families
            nlist = max(20, min(100, num_proteins // 750))
        else:
            # Use 1 centroid per 100-500 proteins for small families
            nlist = max(1, min(20, num_proteins // 250))
        
        # Ensure we have enough training points (39 per centroid minimum)
        min_training_points = nlist * 39
        if num_proteins < min_training_points:
            nlist = max(1, num_proteins // 39)
        
        return nlist
    
    def create_hybrid_family_index(self, family_id: str, force: bool = False) -> Dict[str, str]:
        """
        Create hybrid FAISS index for a family: Binary for family-level, Float for within-family.
        
        Args:
            family_id: Family identifier (must be mapped)
            force: Whether to recreate existing indexes
            
        Returns:
            Dict with paths to binary and float indexes
        """
        import faiss
        
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in mapped families")
        
        file_safe_id = self.family_mapping[family_id]['file_safe_id']
        
        # Check if indexes already exist
        binary_index_file = self.index_dir / "families" / f"{file_safe_id}_binary.faiss"
        float_index_file = self.index_dir / "families" / f"{file_safe_id}_float.faiss"
        hnsw_index_file = self.index_dir / "hnsw" / f"{file_safe_id}_hnsw.faiss"
        
        if not force and binary_index_file.exists() and float_index_file.exists():
            logger.info(f"Indexes already exist for family {family_id}")
            return {
                'binary_index': str(binary_index_file),
                'float_index': str(float_index_file),
                'hnsw_index': str(hnsw_index_file) if hnsw_index_file.exists() else None
            }
        
        # Load family embeddings
        embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
        embedding_dim = embeddings.shape[1]
        num_proteins = len(protein_ids)
        
        # Ensure embeddings are float32
        if embeddings.dtype != np.float32:
            embeddings = embeddings.astype(np.float32)
        
        # Calculate optimal parameters
        nlist = self._calculate_optimal_nlist(num_proteins, embedding_dim)
        
        logger.info(f"Creating hybrid indexes for family {family_id}: {num_proteins} proteins, {embedding_dim} dim, nlist={nlist}")
        
        indexes = {}
        
        # 1. Create Binary Index for Family-Level Search (Memory Efficient)
        try:
            # Binarize embeddings: sign(x) > 0 → 1, else 0
            binary_embeddings = (embeddings > 0).astype(np.uint8)
            
            # FAISS binary index expects dimension in bits
            binary_dim_bits = embedding_dim
            
            # Create binary index with correct dimension
            binary_index = faiss.IndexBinaryFlat(binary_dim_bits)
            binary_index.add(binary_embeddings)
            
            # Save binary index
            faiss.write_index(binary_index, str(binary_index_file))
            indexes['binary_index'] = str(binary_index_file)
            
            logger.info(f"Created binary index for family {family_id} with dimension {binary_dim_bits} bits")
            
        except Exception as e:
            logger.warning(f"Failed to create binary index for family {family_id}: {e}")
            # Skip binary index creation if it fails
            pass
        
        # 2. Create Float Index for Within-Family Search (High Accuracy)
        try:
            if num_proteins >= (nlist * 39):
                # Use IVF for large families
                quantizer = faiss.IndexFlatL2(embedding_dim)
                float_index = faiss.IndexIVFFlat(quantizer, embedding_dim, nlist)
                
                # Train and add
                float_index.train(embeddings)
                float_index.add(embeddings)
                index_type = 'IVF_float'
            else:
                # Use flat index for small families
                float_index = faiss.IndexFlatL2(embedding_dim)
                float_index.add(embeddings)
                index_type = 'Flat_float'
            
            # Save float index
            faiss.write_index(float_index, str(float_index_file))
            indexes['float_index'] = str(float_index_file)
            
            logger.info(f"Created {index_type} index for family {family_id}")
            
        except Exception as e:
            logger.warning(f"Failed to create float index for family {family_id}: {e}")
        
        # 3. Create HNSW Index for Approximate Search (Optional)
        if self.use_hnsw and num_proteins > 1000:
            try:
                # Create HNSW index for approximate nearest neighbor search
                hnsw_index = faiss.IndexHNSWFlat(embedding_dim, self.hnsw_m)
                hnsw_index.hnsw.efConstruction = self.hnsw_ef_construction
                hnsw_index.add(embeddings)
                
                # Save HNSW index
                faiss.write_index(hnsw_index, str(hnsw_index_file))
                indexes['hnsw_index'] = str(hnsw_index_file)
                
                logger.info(f"Created HNSW index for family {family_id}")
                
            except Exception as e:
                logger.warning(f"Failed to create HNSW index for family {family_id}: {e}")
        
        # Save metadata
        metadata = {
            'family_id': family_id,
            'num_proteins': num_proteins,
            'embedding_dim': embedding_dim,
            'nlist': nlist,
            'index_types': list(indexes.keys()),
            'protein_ids': protein_ids,
            'creation_timestamp': pd.Timestamp.now().isoformat()
        }
        
        metadata_file = self.index_dir / "metadata" / f"{file_safe_id}_hybrid_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        return indexes
    
    def search_within_family_hybrid(self, 
                                  family_id: str, 
                                  query_vector: np.ndarray, 
                                  k: int = 10,
                                  search_type: str = 'auto') -> List[Dict[str, Any]]:
        """
        Search within a family using hybrid indexing strategy.
        
        Args:
            family_id: Family identifier (must be mapped)
            query_vector: Query embedding (float32)
            k: Number of results to return
            search_type: 'auto', 'binary', 'float', or 'hnsw'
            
        Returns:
            List of search results with protein_id and similarity score
        """
        import faiss
        
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in mapped families")
        
        file_safe_id = self.family_mapping[family_id]['file_safe_id']
        
        # Auto-select best search method based on family size
        if search_type == 'auto':
            family_size = self.family_mapping[family_id]['num_proteins']
            if family_size > 100000:
                search_type = 'hnsw'  # Use HNSW for very large families
            elif family_size > 10000:
                search_type = 'float'  # Use float for medium families
            else:
                search_type = 'binary'  # Use binary for small families
        
        # Load protein IDs for this family
        embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
        
        results = []
        
        try:
            if search_type == 'hnsw':
                # Use HNSW for approximate search (fastest for large families)
                hnsw_index_file = self.index_dir / "hnsw" / f"{file_safe_id}_hnsw.faiss"
                if hnsw_index_file.exists():
                    index = faiss.read_index(str(hnsw_index_file))
                    D, I = index.search(query_vector.reshape(1, -1), k)
                    
                    for i, (dist, idx) in enumerate(zip(D[0], I[0])):
                        if idx < len(protein_ids):
                            results.append({
                                'protein_id': protein_ids[idx],
                                'similarity': float(1.0 / (1.0 + dist)),  # Convert distance to similarity
                                'rank': i + 1,
                                'search_method': 'hnsw'
                            })
            
            elif search_type == 'float':
                # Use float index for high accuracy
                float_index_file = self.index_dir / "families" / f"{file_safe_id}_float.faiss"
                if float_index_file.exists():
                    index = faiss.read_index(str(float_index_file))
                    D, I = index.search(query_vector.reshape(1, -1), k)
                    
                    for i, (dist, idx) in enumerate(zip(D[0], I[0])):
                        if idx < len(protein_ids):
                            results.append({
                                'protein_id': protein_ids[idx],
                                'similarity': float(1.0 / (1.0 + dist)),
                                'rank': i + 1,
                                'search_method': 'float'
                            })
            
            elif search_type == 'binary':
                # Use binary index for memory efficiency
                binary_index_file = self.index_dir / "families" / f"{file_safe_id}_binary.faiss"
                if binary_index_file.exists():
                    index = faiss.read_index(str(binary_index_file))
                    
                    # Binarize query
                    query_binary = (query_vector > 0).astype(np.uint8)
                    
                    # Ensure query dimension matches index dimension
                    if query_binary.shape[0] != index.d:
                        if query_binary.shape[0] < index.d:
                            padding = np.zeros(index.d - query_binary.shape[0], dtype=np.uint8)
                            query_binary = np.hstack([query_binary, padding])
                        else:
                            query_binary = query_binary[:index.d]
                    
                    D, I = index.search(query_binary.reshape(1, -1), k)
                    
                    for i, (dist, idx) in enumerate(zip(D[0], I[0])):
                        if idx < len(protein_ids):
                            # Convert Hamming distance to similarity
                            max_distance = query_binary.shape[0]
                            similarity = float(1.0 - (dist / max_distance))
                            results.append({
                                'protein_id': protein_ids[idx],
                                'similarity': similarity,
                                'rank': i + 1,
                                'search_method': 'binary'
                            })
            
            # Fallback to brute force if no index found
            if not results:
                logger.warning(f"No FAISS index found for family {family_id}, using brute force search")
                similarities = np.dot(embeddings, query_vector) / (np.linalg.norm(embeddings, axis=1) * np.linalg.norm(query_vector))
                top_indices = np.argsort(similarities)[::-1][:k]
                
                for i, idx in enumerate(top_indices):
                    results.append({
                        'protein_id': protein_ids[idx],
                        'similarity': float(similarities[idx]),
                        'rank': i + 1,
                        'search_method': 'brute_force'
                    })
        
        except Exception as e:
            logger.error(f"Search failed for family {family_id}: {e}")
            # Fallback to brute force
            similarities = np.dot(embeddings, query_vector) / (np.linalg.norm(embeddings, axis=1) * np.linalg.norm(query_vector))
            top_indices = np.argsort(similarities)[::-1][:k]
            
            for i, idx in enumerate(top_indices):
                results.append({
                    'protein_id': protein_ids[idx],
                    'similarity': float(similarities[idx]),
                    'rank': i + 1,
                    'search_method': 'brute_force_fallback'
                })
        
        return results
    
    def create_family_centroid_binary_advanced(self, output_npz: str = None) -> str:
        """
        Create advanced binary centroids with adaptive selection for mapped families.
        Optimized for 250M protein scale with diversity preservation.
        """
        centroids = []
        eigenprotein_ids = []
        family_ids = []
        family_metadata = []
        
        for family_id in self.family_mapping.keys():
            try:
                embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
                if embeddings.dtype != np.float32:
                    raise ValueError(f"Family {family_id} embeddings must be float32 for centroid computation.")
                
                # Adaptive centroid selection based on family size
                family_size = len(embeddings)
                if family_size > 100000:
                    # Large families: use multiple centroids for diversity
                    n_centroids = min(10, family_size // 10000)
                    centroids_family = []
                    eigenproteins_family = []
                    
                    # K-means clustering for multiple centroids
                    from sklearn.cluster import KMeans
                    kmeans = KMeans(n_clusters=n_centroids, random_state=42, n_init=10)
                    cluster_labels = kmeans.fit_predict(embeddings)
                    
                    for i in range(n_centroids):
                        cluster_mask = cluster_labels == i
                        cluster_embeddings = embeddings[cluster_mask]
                        cluster_proteins = [protein_ids[j] for j, mask in enumerate(cluster_mask) if mask]
                        
                        if len(cluster_embeddings) > 0:
                            centroid = cluster_embeddings.mean(axis=0)
                            centroid_bin = (centroid > 0).astype(np.uint8)
                            centroid_bin_packed = np.packbits(centroid_bin)
                            
                            # Find eigenprotein (closest to centroid)
                            dists = np.linalg.norm(cluster_embeddings - centroid, axis=1)
                            medoid_idx = int(np.argmin(dists))
                            eigenprotein_ids.append(cluster_proteins[medoid_idx])
                            centroids_family.append(centroid_bin_packed)
                    
                    centroids.extend(centroids_family)
                    family_ids.extend([family_id] * len(centroids_family))
                    family_metadata.extend([{
                        'family_id': family_id,
                        'centroid_type': 'multi_centroid',
                        'n_centroids': len(centroids_family),
                        'family_size': family_size
                    }] * len(centroids_family))
                    
                else:
                    # Small/medium families: use single centroid
                    centroid = embeddings.mean(axis=0)
                    centroid_bin = (centroid > 0).astype(np.uint8)
                    centroid_bin_packed = np.packbits(centroid_bin)
                    
                    # Find eigenprotein (closest to centroid)
                    dists = np.linalg.norm(embeddings - centroid, axis=1)
                    medoid_idx = int(np.argmin(dists))
                    eigenprotein_ids.append(protein_ids[medoid_idx])
                    centroids.append(centroid_bin_packed)
                    family_ids.append(family_id)
                    family_metadata.append({
                        'family_id': family_id,
                        'centroid_type': 'single_centroid',
                        'n_centroids': 1,
                        'family_size': family_size
                    })
                
            except Exception as e:
                logger.warning(f"Failed to process family {family_id}: {e}")
                continue
        
        if not centroids:
            raise ValueError("No valid families found for centroid computation")
        
        family_ids_arr = np.array(family_ids)
        centroids_arr = np.stack(centroids)
        eigenprotein_ids_arr = np.array(eigenprotein_ids)
        
        if output_npz is None:
            output_npz = self.base_dir / "family_centroids_binary_advanced.npz"
        
        np.savez_compressed(
            output_npz, 
            family_ids=family_ids_arr, 
            centroids=centroids_arr, 
            eigenprotein_ids=eigenprotein_ids_arr,
            family_metadata=family_metadata
        )
        
        logger.info(f"Created advanced binary centroids for {len(set(family_ids))} families with {len(centroids)} total centroids")
        return str(output_npz)
    
    def assign_family_advanced(self, query_vector: np.ndarray, centroids_npz: str = None) -> dict:
        """
        Assign query to family using advanced binary centroid search with diversity.
        """
        import faiss
        if query_vector.dtype != np.float32:
            raise ValueError("Query vector must be float32 for assignment.")
        if centroids_npz is None:
            centroids_npz = self.base_dir / "family_centroids_binary_advanced.npz"
        
        if not centroids_npz.exists():
            raise FileNotFoundError(f"Advanced centroids file not found: {centroids_npz}. Run create_family_centroid_binary_advanced first.")
        
        data = np.load(centroids_npz, allow_pickle=True)
        family_ids = data['family_ids']
        centroids = data['centroids']
        eigenprotein_ids = data['eigenprotein_ids']
        family_metadata = data['family_metadata']
        
        # Binarize query
        query_bin = (query_vector > 0).astype(np.uint8)
        
        # Pad query if needed to match centroid dimension
        centroid_dim = centroids.shape[1] * 8
        if query_bin.shape[0] < centroid_dim:
            padding = np.zeros(centroid_dim - query_bin.shape[0], dtype=np.uint8)
            query_bin = np.hstack([query_bin, padding])
        elif query_bin.shape[0] > centroid_dim:
            query_bin = query_bin[:centroid_dim]
        
        query_bin_packed = np.packbits(query_bin)
        
        # Use FAISS IndexBinaryFlat for Hamming search
        d = centroids.shape[1] * 8
        index = faiss.IndexBinaryFlat(d)
        index.add(centroids)
        D, I = index.search(query_bin_packed.reshape(1, -1), 1)
        idx = int(I[0][0])
        
        # Convert Hamming distance to confidence (positive value)
        max_distance = centroid_dim
        confidence = float(1.0 - (abs(D[0][0]) / max_distance))
        
        # Get family metadata
        family_meta = family_metadata[idx] if len(family_metadata) > idx else {}
        
        result = {
            'family_id': str(family_ids[idx]),
            'confidence': confidence,
            'eigenprotein_id': str(eigenprotein_ids[idx]),
            'centroid_type': family_meta.get('centroid_type', 'unknown'),
            'n_centroids': family_meta.get('n_centroids', 1),
            'family_size': family_meta.get('family_size', 0)
        }
        
        logger.info(f"Assigned query to family {result['family_id']} (confidence={confidence}, type={result['centroid_type']})")
        return result
    
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
        Store embeddings for a protein family with chunking and mapping.
        Args:
            family_id: Family identifier
            embeddings: Embedding array (N x D, np.float32 or np.uint8)
            protein_ids: List of protein IDs
            metadata: Optional metadata DataFrame
        Returns:
            Path to stored family file
        """
        if embeddings.dtype not in [np.uint8, np.float32]:
            raise ValueError("Embeddings must be np.uint8 (binary) or np.float32 (float) for storage.")
        
        # Create file-safe family ID
        file_safe_id = self._make_file_safe(family_id)
        family_file = self.family_dir / f"{file_safe_id}.h5"
        
        embedding_dim = embeddings.shape[1]
        num_proteins = embeddings.shape[0]
        optimal_chunk_size = min(self.chunk_size, max(1, min(num_proteins, max(1, 10000 // embedding_dim))))
        
        # Ensure chunk size is at least 1
        optimal_chunk_size = max(1, optimal_chunk_size)
        
        with h5py.File(family_file, 'w') as f:
            f.create_dataset(
                'embeddings',
                data=embeddings,
                chunks=(max(1, min(optimal_chunk_size, num_proteins)), max(1, embedding_dim)),
                compression=self.compression,
                compression_opts=self.compression_opts,
                shuffle=True
            )
            f.create_dataset(
                'protein_ids',
                data=protein_ids,
                dtype=h5py.special_dtype(vlen=str),
                chunks=(max(1, min(optimal_chunk_size, num_proteins)),),
                compression=self.compression,
                compression_opts=self.compression_opts
            )
            if metadata is not None:
                try:
                    metadata_file = self.metadata_dir / f"{file_safe_id}_metadata.parquet"
                    metadata.to_parquet(metadata_file, compression='gzip')
                except Exception:
                    metadata_file = self.metadata_dir / f"{file_safe_id}_metadata.csv"
                    metadata.to_csv(metadata_file)
                f.attrs['metadata_file'] = str(metadata_file)
            f.attrs['num_proteins'] = len(protein_ids)
            f.attrs['embedding_dim'] = embedding_dim
            f.attrs['chunk_size'] = optimal_chunk_size
            f.attrs['compression'] = self.compression
            f.attrs['embedding_dtype'] = str(embeddings.dtype)
        
        # Update family mapping for selective loading
        self.family_mapping[family_id] = {
            'file_safe_id': file_safe_id,
            'num_proteins': num_proteins,
            'embedding_dim': embedding_dim,
            'file_size_mb': family_file.stat().st_size / (1024 * 1024)
        }
        self._save_family_mapping()
        
        logger.info(f"Stored family {family_id}: {len(protein_ids)} proteins, {embedding_dim} dim, {optimal_chunk_size} chunk size, dtype={embeddings.dtype}")
        return str(family_file)
    
    def load_family_embeddings(self, 
                              family_id: str,
                              start_idx: Optional[int] = None,
                              end_idx: Optional[int] = None,
                              check_memory: bool = True) -> Tuple[np.ndarray, List[str]]:
        """
        Load embeddings for a mapped protein family with memory optimization.
        
        Args:
            family_id: Family identifier (must be in mapped families)
            start_idx: Start index for slicing
            end_idx: End index for slicing
            check_memory: Whether to check memory usage before loading
            
        Returns:
            Tuple of (embeddings, protein_ids)
        """
        # Check if family is mapped
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in mapped families. Available: {list(self.family_mapping.keys())}")
        
        file_safe_id = self.family_mapping[family_id]['file_safe_id']
        family_file = self.family_dir / f"{file_safe_id}.h5"
        
        if not family_file.exists():
            raise FileNotFoundError(f"Family file not found: {family_file}")
        
        # Check memory usage if requested
        if check_memory:
            estimated_memory = self.family_mapping[family_id]['num_proteins'] * self.family_mapping[family_id]['embedding_dim'] * 4  # float32
            if estimated_memory > self.memory_limit_bytes:
                logger.warning(f"Family {family_id} requires {estimated_memory / 1024**3:.2f}GB, exceeding limit {self.memory_limit_gb}GB. Use stream_family_embeddings instead.")
        
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
    
    def create_family_index_ivf_pq(self, family_id: str, nlist: int = 100, m: int = 8) -> str:
        """
        Create and store FAISS IVF-PQ index for a mapped family with Product Quantization.
        This provides optimal compression and search speed for large families.
        
        Args:
            family_id: Family identifier (must be mapped)
            nlist: Number of clusters for IVF
            m: Number of subquantizers for PQ
        Returns:
            Path to stored index file
        """
        import faiss
        
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in mapped families")
        
        embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
        
        # Ensure embeddings are float32
        if embeddings.dtype != np.float32:
            embeddings = embeddings.astype(np.float32)
        
        dimension = embeddings.shape[1]
        
        # Ensure we have enough training points for FAISS clustering
        # According to FAISS FAQ: minimum 39 training points per centroid
        min_training_points = nlist * 39
        
        if len(embeddings) < min_training_points:
            # Adjust nlist to meet training requirements
            nlist = max(1, len(embeddings) // 39)
            logger.info(f"Adjusted nlist to {nlist} for family {family_id} to meet training requirements")
        
        # Create FAISS IVF-PQ index for optimal compression
        quantizer = faiss.IndexFlatL2(dimension)
        index = faiss.IndexIVFPQ(quantizer, dimension, nlist, m, 8)  # 8 bits per subquantizer
        
        # Train and add embeddings
        index.train(embeddings)
        index.add(embeddings)
        
        file_safe_id = self.family_mapping[family_id]['file_safe_id']
        index_file = self.index_dir / "families" / f"{file_safe_id}_ivf_pq.faiss"
        
        faiss.write_index(index, str(index_file))
        
        # Create metadata
        metadata = {
            'family_id': family_id,
            'num_proteins': len(protein_ids),
            'embedding_dim': dimension,
            'index_type': 'faiss_ivf_pq',
            'nlist': nlist,
            'm': m,
            'bits_per_subquantizer': 8,
            'protein_ids': protein_ids
        }
        
        metadata_file = self.index_dir / "metadata" / f"{file_safe_id}_ivf_pq_metadata.json"
        
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Created FAISS IVF-PQ index for family {family_id}: {len(protein_ids)} proteins, {dimension} dim, nlist={nlist}, m={m}")
        return str(index_file)
    
    def create_family_index_float(self, family_id: str, nlist: int = 10) -> str:
        """
        Create and store FAISS IVF float32 similarity index for a family.
        This is used for within-family similarity search (more accurate).
        
        Args:
            family_id: Family identifier
            nlist: Number of clusters for IVF
        Returns:
            Path to stored index file
        """
        import faiss
        embeddings, protein_ids = self.load_family_embeddings(family_id)
        
        # Ensure embeddings are float32
        if embeddings.dtype != np.float32:
            embeddings = embeddings.astype(np.float32)
        
        dimension = embeddings.shape[1]
        
        # Ensure we have enough training points for FAISS clustering
        # According to FAISS FAQ: minimum 39 training points per centroid
        min_training_points = nlist * 39
        
        if len(embeddings) < min_training_points:
            # Adjust nlist to meet training requirements
            nlist = max(1, len(embeddings) // 39)
            logger.info(f"Adjusted nlist to {nlist} for family {family_id} to meet training requirements")
        
        # Create FAISS IVF float index
        quantizer = faiss.IndexFlatL2(dimension)
        index = faiss.IndexIVFFlat(quantizer, dimension, nlist)
        
        # Train and add embeddings
        index.train(embeddings)
        index.add(embeddings)
        
        # Create file-safe family ID for any family name format
        file_safe_id = self._make_file_safe(family_id)
        
        index_file = self.base_dir / "indexes" / "families" / f"{file_safe_id}.faiss"
        index_file.parent.mkdir(parents=True, exist_ok=True)
        
        faiss.write_index(index, str(index_file))
        
        # Create metadata
        metadata = {
            'family_id': family_id,
            'num_proteins': len(protein_ids),
            'embedding_dim': dimension,
            'index_type': 'faiss_ivf_float',
            'nlist': nlist,
            'protein_ids': protein_ids
        }
        
        metadata_file = self.base_dir / "indexes" / "metadata" / f"{file_safe_id}_float_metadata.json"
        metadata_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Created FAISS IVF float index for family {family_id}: {len(protein_ids)} proteins, {dimension} dim")
        return str(index_file)
    
    def create_family_index(self, family_id: str, **kwargs) -> str:
        """
        Create and store FAISS IVF binary similarity index for a family.
        NOTE: This is deprecated. Use create_family_index_float for within-family search.
        Binary indexes should only be used for centroids.
        
        Args:
            family_id: Family identifier
        Returns:
            Path to stored index file
        """
        logger.warning("create_family_index is deprecated. Use create_family_index_float for within-family search.")
        return self.create_family_index_float(family_id)
    
    def create_family_centroid_binary(self, output_npz: str = None) -> str:
        """
        Compute and save binary centroids for mapped families from float32 means.
        Only processes mapped families for memory efficiency.
        Args:
            output_npz: Optional output .npz path
        Returns:
            Path to saved .npz file
        """
        centroids = []
        eigenprotein_ids = []
        family_ids = []
        
        for family_id in self.family_mapping.keys():
            try:
                embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
                if embeddings.dtype != np.float32:
                    raise ValueError(f"Family {family_id} embeddings must be float32 for centroid computation.")
                
                centroid = embeddings.mean(axis=0)
                # Binarize: sign(x) > 0 → 1, else 0
                centroid_bin = (centroid > 0).astype(np.uint8)
                # Pack bits to uint8
                centroid_bin_packed = np.packbits(centroid_bin)
                # Find eigenprotein (closest to centroid in L2)
                dists = np.linalg.norm(embeddings - centroid, axis=1)
                medoid_idx = int(np.argmin(dists))
                eigenprotein_ids.append(protein_ids[medoid_idx])
                centroids.append(centroid_bin_packed)
                family_ids.append(family_id)
                
            except Exception as e:
                logger.warning(f"Failed to process family {family_id}: {e}")
                continue
        
        if not centroids:
            raise ValueError("No valid families found for centroid computation")
        
        family_ids_arr = np.array(family_ids)
        centroids_arr = np.stack(centroids)
        eigenprotein_ids_arr = np.array(eigenprotein_ids)
        
        if output_npz is None:
            output_npz = self.base_dir / "family_centroids_binary.npz"
        
        np.savez_compressed(output_npz, family_ids=family_ids_arr, centroids=centroids_arr, eigenprotein_ids=eigenprotein_ids_arr)
        logger.info(f"Saved binary centroids for {len(family_ids)} mapped families to {output_npz}")
        return str(output_npz)
    
    def assign_family(self, query_vector: np.ndarray, centroids_npz: str = None) -> dict:
        """
        Assign a float32 query vector to the closest mapped family using binary centroid Hamming search.
        Args:
            query_vector: Query embedding (float32, shape [D])
            centroids_npz: Optional path to centroids .npz
        Returns:
            Dict with keys: 'family_id', 'confidence', 'eigenprotein_id'
        """
        import faiss
        if query_vector.dtype != np.float32:
            raise ValueError("Query vector must be float32 for assignment.")
        if centroids_npz is None:
            centroids_npz = self.base_dir / "family_centroids_binary.npz"
        
        if not centroids_npz.exists():
            raise FileNotFoundError(f"Centroids file not found: {centroids_npz}. Run create_family_centroid_binary first.")
        
        data = np.load(centroids_npz, allow_pickle=True)
        family_ids = data['family_ids']
        centroids = data['centroids']
        eigenprotein_ids = data['eigenprotein_ids']
        
        # Binarize query: sign(x) > 0 → 1, else 0, then packbits
        query_bin = (query_vector > 0).astype(np.uint8)
        query_bin_packed = np.packbits(query_bin)
        
        # Use FAISS IndexBinaryFlat for Hamming search
        d = centroids.shape[1] * 8
        index = faiss.IndexBinaryFlat(d)
        index.add(centroids)
        D, I = index.search(query_bin_packed.reshape(1, -1), 1)
        idx = int(I[0][0])
        confidence = float(-D[0][0])  # negative Hamming distance
        
        result = {
            'family_id': str(family_ids[idx]),
            'confidence': confidence,
            'eigenprotein_id': str(eigenprotein_ids[idx])
        }
        logger.info(f"Assigned query to mapped family {result['family_id']} (confidence={confidence})")
        return result
    
    def get_family_list(self) -> List[str]:
        """Get list of all mapped families only."""
        return list(self.family_mapping.keys())
    
    def get_family_stats(self) -> Dict[str, Dict]:
        """Get statistics for mapped families only."""
        return self.family_mapping.copy()
    
    def stream_family_embeddings(self, 
                                family_id: str,
                                batch_size: int = 1000) -> Iterator[Tuple[np.ndarray, List[str]]]:
        """
        Stream embeddings for a mapped family in batches.
        
        Args:
            family_id: Family identifier (must be mapped)
            batch_size: Number of proteins per batch
            
        Yields:
            Tuples of (embeddings_batch, protein_ids_batch)
        """
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in mapped families")
        
        file_safe_id = self.family_mapping[family_id]['file_safe_id']
        family_file = self.family_dir / f"{file_safe_id}.h5"
        
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
        Load metadata for a mapped family using CompressedMetadataStorage.
        Args:
            family_id: Family identifier (must be mapped)
        Returns:
            Metadata DataFrame
        """
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in mapped families")
        
        meta_storage = CompressedMetadataStorage(metadata_dir=str(self.metadata_dir))
        return meta_storage.load_metadata(family_id=family_id)
    
    def search_within_family(self, 
                           family_id: str, 
                           query_vector: np.ndarray, 
                           k: int = 10,
                           use_ivf_pq: bool = True) -> List[Dict[str, Any]]:
        """
        Search within a mapped family using FAISS index.
        
        Args:
            family_id: Family identifier (must be mapped)
            query_vector: Query embedding (float32)
            k: Number of results to return
            use_ivf_pq: Whether to use IVF-PQ index if available
            
        Returns:
            List of search results with protein_id and similarity score
        """
        import faiss
        
        if family_id not in self.family_mapping:
            raise ValueError(f"Family {family_id} not found in mapped families")
        
        file_safe_id = self.family_mapping[family_id]['file_safe_id']
        
        # Try to load IVF-PQ index first
        if use_ivf_pq:
            index_file = self.index_dir / "families" / f"{file_safe_id}_ivf_pq.faiss"
            if index_file.exists():
                index = faiss.read_index(str(index_file))
                # IVF-PQ requires normalization
                faiss.normalize_L2(query_vector.reshape(1, -1))
                D, I = index.search(query_vector.reshape(1, -1), k)
                
                # Load protein IDs for this family
                embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
                
                results = []
                for i, (dist, idx) in enumerate(zip(D[0], I[0])):
                    if idx < len(protein_ids):
                        results.append({
                            'protein_id': protein_ids[idx],
                            'similarity': float(1.0 - dist),  # Convert distance to similarity
                            'rank': i + 1
                        })
                return results
        
        # Fallback to IVF float index
        index_file = self.index_dir / "families" / f"{file_safe_id}.faiss"
        if index_file.exists():
            index = faiss.read_index(str(index_file))
            D, I = index.search(query_vector.reshape(1, -1), k)
            
            # Load protein IDs for this family
            embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
            
            results = []
            for i, (dist, idx) in enumerate(zip(D[0], I[0])):
                if idx < len(protein_ids):
                    results.append({
                        'protein_id': protein_ids[idx],
                        'similarity': float(1.0 / (1.0 + dist)),  # Convert distance to similarity
                        'rank': i + 1
                    })
            return results
        
        # Fallback to brute force search
        logger.warning(f"No FAISS index found for family {family_id}, using brute force search")
        embeddings, protein_ids = self.load_family_embeddings(family_id, check_memory=False)
        
        # Compute similarities
        similarities = np.dot(embeddings, query_vector) / (np.linalg.norm(embeddings, axis=1) * np.linalg.norm(query_vector))
        top_indices = np.argsort(similarities)[::-1][:k]
        
        results = []
        for i, idx in enumerate(top_indices):
            results.append({
                'protein_id': protein_ids[idx],
                'similarity': float(similarities[idx]),
                'rank': i + 1
            })
        
        return results


class CompressedMetadataStorage:
    """
    Efficient metadata storage with compression and indexing.
    """
    
    def __init__(self, metadata_dir: str = "data/metadata"):
        self.metadata_dir = Path(metadata_dir)
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
        
    def _make_file_safe(self, family_id: str) -> str:
        """
        Convert any family ID to a file-safe format.
        
        Args:
            family_id: Family ID (can be any string)
            
        Returns:
            File-safe version of the family ID
        """
        import re
        # Replace any non-alphanumeric characters with underscores
        # This ensures the filename is safe for all operating systems
        file_safe = re.sub(r'[^a-zA-Z0-9]', '_', family_id)
        # Remove multiple consecutive underscores
        file_safe = re.sub(r'_+', '_', file_safe)
        # Remove leading/trailing underscores
        file_safe = file_safe.strip('_')
        # Ensure it's not empty
        if not file_safe:
            file_safe = 'unknown_family'
        return file_safe
    
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
            # Create file-safe family ID
            file_safe_id = self._make_file_safe(family_id)
            filename = f"{file_safe_id}_metadata.parquet"
        else:
            filename = "global_metadata.parquet"
        
        filepath = self.metadata_dir / filename
        
        # Store as Parquet with compression when engine available; fallback to CSV
        try:
            metadata.to_parquet(filepath, compression=compression)
        except Exception:
            # Fallback to CSV to avoid optional dependency on pyarrow/fastparquet in test env
            filepath = self.metadata_dir / filename.replace('.parquet', '.csv')
            metadata.to_csv(filepath)
        
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
            file_safe_id = self._make_file_safe(family_id)
            index_file = self.metadata_dir / f"{file_safe_id}_index.pkl"
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
            file_safe_id = self._make_file_safe(family_id)
            filepath = self.metadata_dir / f"{file_safe_id}_metadata.parquet"
        else:
            filepath = self.metadata_dir / "global_metadata.parquet"
        
        if not filepath.exists():
            # Try CSV fallback
            filepath_csv = str(filepath).replace('.parquet', '.csv')
            if not Path(filepath_csv).exists():
                raise FileNotFoundError(f"Metadata file not found for family {family_id}: {filepath}")
            metadata = pd.read_csv(filepath_csv, index_col=0)
        else:
            try:
                metadata = pd.read_parquet(filepath)
            except Exception:
                # Fallback to CSV if parquet engine not available
                filepath_csv = str(filepath).replace('.parquet', '.csv')
                if Path(filepath_csv).exists():
                    metadata = pd.read_csv(filepath_csv, index_col=0)
                else:
                    raise
        
        if protein_ids is not None:
            metadata = metadata.loc[metadata.index.intersection(protein_ids)]
        
        return metadata


class MemoryEfficientLoader:
    """
    Memory-efficient loading for large datasets with mapped families only.
    """
    
    def __init__(self, storage: ProteinStorage):
        self.storage = storage
    
    def load_families_batch(self, 
                           family_ids: List[str],
                           max_memory_gb: float = 4.0) -> Iterator[Tuple[str, np.ndarray, List[str]]]:
        """
        Load multiple mapped families in memory-efficient batches.
        
        Args:
            family_ids: List of mapped family IDs to load
            max_memory_gb: Maximum memory usage in GB
            
        Yields:
            Tuples of (family_id, embeddings, protein_ids)
        """
        max_memory_bytes = max_memory_gb * 1024**3
        
        for family_id in family_ids:
            if family_id not in self.storage.family_mapping:
                logger.warning(f"Family {family_id} not found in mapped families, skipping")
                continue
            
            # Estimate memory usage
            stats = self.storage.family_mapping[family_id]
            estimated_memory = stats['num_proteins'] * stats['embedding_dim'] * 4  # float32
            
            if estimated_memory > max_memory_bytes:
                # Load in chunks
                for embeddings, protein_ids in self.storage.stream_family_embeddings(family_id):
                    yield family_id, embeddings, protein_ids
            else:
                # Load entire family
                embeddings, protein_ids = self.storage.load_family_embeddings(family_id)
                yield family_id, embeddings, protein_ids


class ProteinIDsIndex:
    """
    Efficient protein IDs index for fast searching (exact UniProt ID match only).
    Optimized for mapped families only.
    """
    def __init__(self, base_dir: str = "data"):
        self.base_dir = Path(base_dir)
        self.index_dir = self.base_dir / "indexes" / "protein_ids"
        self.index_file = self.index_dir / "protein_ids_index.json"
        self.mapping_file = self.index_dir / "protein_to_family.json"
        self.protein_ids_index = {}
        self.protein_to_family = {}
        self._load_index()

    def _load_index(self):
        """Load the protein IDs index from disk."""
        try:
            if self.index_file.exists():
                with open(self.index_file, 'r') as f:
                    self.protein_ids_index = json.load(f)
                logger.info(f"Loaded protein IDs index with {len(self.protein_ids_index)} entries")
            if self.mapping_file.exists():
                with open(self.mapping_file, 'r') as f:
                    self.protein_to_family = json.load(f)
                logger.info(f"Loaded protein to family mapping with {len(self.protein_to_family)} entries")
        except Exception as e:
            logger.warning(f"Failed to load protein IDs index: {e}")

    def search_protein(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Search for a protein by UniProt ID (exact match only).
        Args:
            uniprot_id: UniProt ID to search for
        Returns:
            Dict with protein information if found, None otherwise
        """
        if not uniprot_id:
            return None
        uniprot_id = uniprot_id.strip()
        return self.protein_ids_index.get(uniprot_id)

    def get_protein_family(self, protein_id: str) -> Optional[str]:
        """Get the family ID for a protein ID (mapped families only)."""
        return self.protein_to_family.get(protein_id)

    def get_all_proteins(self) -> List[str]:
        """Get all protein IDs in the index (mapped families only)."""
        return list(self.protein_to_family.keys())

    def get_proteins_by_family(self, family_id: str) -> List[str]:
        """Get all protein IDs for a specific mapped family."""
        return [pid for pid, fid in self.protein_to_family.items() if fid == family_id]


def create_storage_structure(embeddings_file: str,
                         metadata_file: str,
                         output_dir: str = "data",
                         family_column: str = "Protein families",
                         max_family_size: int = 100000) -> ProteinStorage:
    """
    Convert existing data to optimized storage structure with mapped families only.
    
    Args:
        embeddings_file: Path to existing embeddings H5 file
        metadata_file: Path to existing metadata CSV file
        output_dir: Output directory for optimized storage
        family_column: Column name for family information
        max_family_size: Maximum proteins per family
        
    Returns:
        ProteinStorage instance with mapped families
    """
    logger.info("Creating optimized storage structure with mapped families only...")
    
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
    
    logger.info(f"Optimized storage created in {output_dir} with {len(storage.family_mapping)} mapped families")
    return storage


def _create_artificial_families(protein_ids: List[str], max_family_size: int) -> Iterator[Tuple[int, pd.DataFrame]]:
    """Create artificial families when no family information is available."""
    for i in range(0, len(protein_ids), max_family_size):
        family_proteins = protein_ids[i:i + max_family_size]
        family_metadata = pd.DataFrame(index=family_proteins)
        yield i // max_family_size, family_metadata 