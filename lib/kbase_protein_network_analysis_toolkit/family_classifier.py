"""
Protein Family Classifier Module

This module handles the classification of proteins into families using either
embedding similarity or MSA-based methods. It provides efficient family indexing
for fast retrieval during query operations.
"""

import numpy as np
import pandas as pd
import h5py
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics.pairwise import cosine_similarity
from typing import List, Dict, Tuple, Optional, Union
import logging
from tqdm import tqdm
import os
from collections import defaultdict
import joblib

logger = logging.getLogger(__name__)


class ProteinFamilyClassifier:
    """
    Classifies proteins into families using embedding similarity or MSA-based methods.
    
    This class provides methods to group proteins into functional families
    and create efficient indexes for fast retrieval during query operations.
    """
    
    def __init__(self, method: str = "embedding_similarity", **kwargs):
        """
        Initialize the family classifier.
        
        Args:
            method: Classification method ("embedding_similarity" or "msa_based")
            **kwargs: Additional parameters for the classification method
        """
        self.method = method
        self.family_centroids = None
        self.family_assignments = None
        self.family_proteins = None
        self.family_embeddings = None
        self.family_metadata = None
        
        # Method-specific parameters
        if method == "embedding_similarity":
            self.similarity_threshold = kwargs.get('similarity_threshold', 0.7)
            self.min_family_size = kwargs.get('min_family_size', 10)
            self.max_families = kwargs.get('max_families', 100)
            self.clustering_method = kwargs.get('clustering_method', 'kmeans')
        elif method == "msa_based":
            self.identity_threshold = kwargs.get('identity_threshold', 0.3)
            self.coverage_threshold = kwargs.get('coverage_threshold', 0.8)
        else:
            raise ValueError(f"Unknown classification method: {method}")
    
    def classify_by_embedding_similarity(self, embeddings: np.ndarray, 
                                       protein_ids: List[str],
                                       metadata: Optional[pd.DataFrame] = None) -> Dict[str, int]:
        """
        Classify proteins into families using embedding similarity.
        
        Args:
            embeddings: Protein embeddings array (N x D)
            protein_ids: List of protein IDs
            metadata: Optional metadata DataFrame
            
        Returns:
            Dictionary mapping protein IDs to family IDs
        """
        logger.info("Classifying proteins by embedding similarity...")
        
        # Normalize embeddings
        embeddings_norm = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        
        if self.clustering_method == 'kmeans':
            # Use K-means clustering
            n_clusters = min(self.max_families, len(embeddings) // self.min_family_size)
            n_clusters = max(2, n_clusters)  # At least 2 clusters
            
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            family_labels = kmeans.fit_predict(embeddings_norm)
            
            # Store centroids
            self.family_centroids = kmeans.cluster_centers_
            
        elif self.clustering_method == 'dbscan':
            # Use DBSCAN clustering
            dbscan = DBSCAN(eps=1-self.similarity_threshold, min_samples=self.min_family_size)
            family_labels = dbscan.fit_predict(embeddings_norm)
            
            # Convert -1 labels (noise) to separate families
            noise_mask = family_labels == -1
            if noise_mask.any():
                max_label = family_labels.max()
                family_labels[noise_mask] = np.arange(max_label + 1, max_label + 1 + noise_mask.sum())
        
        # Create family assignments dictionary
        family_assignments = {pid: int(label) for pid, label in zip(protein_ids, family_labels)}
        
        # Group proteins by family
        self.family_proteins = defaultdict(list)
        self.family_embeddings = defaultdict(list)
        
        for i, protein_id in enumerate(protein_ids):
            family_id = family_labels[i]
            self.family_proteins[family_id].append(protein_id)
            self.family_embeddings[family_id].append(embeddings[i])
        
        # Store metadata by family
        if metadata is not None:
            self.family_metadata = defaultdict(list)
            for i, protein_id in enumerate(protein_ids):
                family_id = family_labels[i]
                if protein_id in metadata.index:
                    self.family_metadata[family_id].append(metadata.loc[protein_id].to_dict())
        
        logger.info(f"Classified {len(protein_ids)} proteins into {len(set(family_labels))} families")
        return family_assignments
    
    def classify_by_msa(self, sequences: List[str], 
                       protein_ids: List[str],
                       metadata: Optional[pd.DataFrame] = None) -> Dict[str, int]:
        """
        Classify proteins into families using MSA-based methods.
        
        Args:
            sequences: List of protein sequences
            protein_ids: List of protein IDs
            metadata: Optional metadata DataFrame
            
        Returns:
            Dictionary mapping protein IDs to family IDs
        """
        logger.info("Classifying proteins by MSA...")
        
        try:
            from Bio import AlignIO, Align
            from Bio.Align.Applications import ClustalOmegaCommandline
        except ImportError:
            raise ImportError("Biopython is required for MSA-based classification")
        
        # This is a simplified implementation
        # In practice, you would use more sophisticated MSA tools
        
        # For now, use a simple approach based on sequence similarity
        family_assignments = {}
        family_id = 0
        assigned_proteins = set()
        
        for i, seq1 in enumerate(sequences):
            if protein_ids[i] in assigned_proteins:
                continue
                
            # Start a new family with this protein
            family_assignments[protein_ids[i]] = family_id
            assigned_proteins.add(protein_ids[i])
            
            # Find similar sequences
            for j, seq2 in enumerate(sequences):
                if i == j or protein_ids[j] in assigned_proteins:
                    continue
                
                # Calculate sequence similarity (simplified)
                similarity = self._calculate_sequence_similarity(seq1, seq2)
                if similarity >= self.identity_threshold:
                    family_assignments[protein_ids[j]] = family_id
                    assigned_proteins.add(protein_ids[j])
            
            family_id += 1
        
        logger.info(f"Classified {len(protein_ids)} proteins into {family_id} families using MSA")
        return family_assignments
    
    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity between two protein sequences."""
        # Simple identity calculation
        min_len = min(len(seq1), len(seq2))
        if min_len == 0:
            return 0.0
        
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / min_len
    
    def classify_proteins(self, embeddings: np.ndarray, 
                         protein_ids: List[str],
                         sequences: Optional[List[str]] = None,
                         metadata: Optional[pd.DataFrame] = None) -> Dict[str, int]:
        """
        Classify proteins into families using the specified method.
        
        Args:
            embeddings: Protein embeddings array
            protein_ids: List of protein IDs
            sequences: List of protein sequences (required for MSA method)
            metadata: Optional metadata DataFrame
            
        Returns:
            Dictionary mapping protein IDs to family IDs
        """
        if self.method == "embedding_similarity":
            return self.classify_by_embedding_similarity(embeddings, protein_ids, metadata)
        elif self.method == "msa_based":
            if sequences is None:
                raise ValueError("Sequences are required for MSA-based classification")
            return self.classify_by_msa(sequences, protein_ids, metadata)
        else:
            raise ValueError(f"Unknown classification method: {self.method}")
    
    def get_family_representatives(self) -> Dict[int, np.ndarray]:
        """
        Get representative embeddings for each family.
        
        Returns:
            Dictionary mapping family IDs to representative embeddings
        """
        if self.family_centroids is not None:
            # Use cluster centroids
            return {i: centroid for i, centroid in enumerate(self.family_centroids)}
        elif self.family_embeddings is not None:
            # Use mean of family embeddings
            representatives = {}
            for family_id, embeddings in self.family_embeddings.items():
                representatives[family_id] = np.mean(embeddings, axis=0)
            return representatives
        else:
            raise ValueError("No family information available")
    
    def predict_family(self, query_embedding: np.ndarray) -> Tuple[int, float]:
        """
        Predict the family for a query protein.
        
        Args:
            query_embedding: Query protein embedding
            
        Returns:
            Tuple of (predicted_family_id, confidence_score)
        """
        if self.family_centroids is None:
            raise ValueError("No family centroids available")
        
        # Check if family centroids are empty or have wrong shape
        if len(self.family_centroids) == 0:
            raise ValueError("Family centroids are empty - no families available for classification")
        
        # Ensure query embedding is 1D
        if query_embedding.ndim > 1:
            query_embedding = query_embedding.flatten()
        
        # Check if query embedding has the right dimension
        if query_embedding.shape[0] != self.family_centroids.shape[1]:
            raise ValueError(f"Query embedding dimension ({query_embedding.shape[0]}) "
                           f"doesn't match family centroids dimension ({self.family_centroids.shape[1]})")
        
        # Normalize query embedding
        query_norm = query_embedding / (np.linalg.norm(query_embedding) + 1e-8)
        
        # Ensure query_norm is 2D for cosine_similarity
        if query_norm.ndim == 1:
            query_norm = query_norm.reshape(1, -1)
        
        # Calculate similarities to all family centroids
        try:
            similarities = cosine_similarity(query_norm, self.family_centroids)[0]
        except Exception as e:
            logger.error(f"Error in cosine similarity calculation: {e}")
            logger.error(f"Query norm shape: {query_norm.shape}")
            logger.error(f"Family centroids shape: {self.family_centroids.shape}")
            raise ValueError(f"Failed to calculate similarities: {e}")
        
        # Check if similarities array is empty
        if len(similarities) == 0:
            raise ValueError("No similarities calculated - family centroids may be corrupted")
        
        # Find the most similar family
        best_family = np.argmax(similarities)
        confidence = similarities[best_family]
        
        return int(best_family), float(confidence)
    
    def save_family_index(self, output_file: str):
        """
        Save family index to HDF5 file.
        
        Args:
            output_file: Path to output HDF5 file
        """
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with h5py.File(output_file, 'w') as f:
            # Save family assignments
            if self.family_assignments is not None:
                protein_ids = list(self.family_assignments.keys())
                family_ids = list(self.family_assignments.values())
                f.create_dataset('protein_ids', data=protein_ids, dtype=h5py.special_dtype(vlen=str))
                f.create_dataset('family_ids', data=family_ids)
            
            # Save family centroids
            if self.family_centroids is not None:
                f.create_dataset('family_centroids', data=self.family_centroids, compression='gzip')
            
            # Save family proteins
            if self.family_proteins is not None:
                for family_id, proteins in self.family_proteins.items():
                    f.create_dataset(f'family_{family_id}_proteins', 
                                   data=proteins, dtype=h5py.special_dtype(vlen=str))
            
            # Save family embeddings
            if self.family_embeddings is not None:
                for family_id, embeddings in self.family_embeddings.items():
                    f.create_dataset(f'family_{family_id}_embeddings', 
                                   data=embeddings, compression='gzip')
        
        logger.info(f"Saved family index to {output_file}")
    
    def load_family_index(self, input_file: str):
        """
        Load family index from HDF5 file.
        
        Args:
            input_file: Path to input HDF5 file
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Family index file not found: {input_file}")
        
        with h5py.File(input_file, 'r') as f:
            # Load family assignments
            if 'protein_ids' in f and 'family_ids' in f:
                protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][:]]
                family_ids = f['family_ids'][:]
                self.family_assignments = {pid: int(fid) for pid, fid in zip(protein_ids, family_ids)}
                logger.info(f"Loaded {len(self.family_assignments)} family assignments")
            else:
                logger.warning("No family assignments found in index file")
                self.family_assignments = {}
            
            # Load family centroids
            if 'family_centroids' in f:
                self.family_centroids = f['family_centroids'][:]
                if len(self.family_centroids) == 0:
                    raise ValueError("Family centroids are empty in the index file")
                logger.info(f"Loaded {len(self.family_centroids)} family centroids with shape {self.family_centroids.shape}")
            else:
                logger.warning("No family centroids found in index file")
                self.family_centroids = None
            
            # Load family proteins
            self.family_proteins = defaultdict(list)
            family_proteins_count = 0
            for key in f.keys():
                if key.startswith('family_') and key.endswith('_proteins'):
                    family_id = int(key.split('_')[1])
                    proteins = [p.decode('utf-8') if isinstance(p, bytes) else p 
                               for p in f[key][:]]
                    self.family_proteins[family_id] = proteins
                    family_proteins_count += len(proteins)
            
            logger.info(f"Loaded {len(self.family_proteins)} families with {family_proteins_count} total proteins")
            
            # Load family embeddings
            self.family_embeddings = defaultdict(list)
            family_embeddings_count = 0
            for key in f.keys():
                if key.startswith('family_') and key.endswith('_embeddings'):
                    family_id = int(key.split('_')[1])
                    embeddings = f[key][:]
                    self.family_embeddings[family_id] = embeddings
                    family_embeddings_count += len(embeddings)
            
            logger.info(f"Loaded {len(self.family_embeddings)} family embedding sets with {family_embeddings_count} total embeddings")
        
        # Validate that we have the minimum required data
        if self.family_centroids is None:
            raise ValueError("Family centroids are required but not found in the index file")
        
        if len(self.family_centroids) == 0:
            raise ValueError("Family centroids are empty - the index file may be corrupted")
        
        logger.info(f"Successfully loaded family index from {input_file}")


def create_family_index(embeddings_file: str,
                       metadata_file: str,
                       output_file: str,
                       method: str = "embedding_similarity",
                       **kwargs) -> ProteinFamilyClassifier:
    """
    Create a family index from embeddings and metadata.
    
    Args:
        embeddings_file: Path to embeddings HDF5 file
        metadata_file: Path to metadata CSV file
        output_file: Path to output family index file
        method: Classification method
        **kwargs: Additional parameters for classification
        
    Returns:
        ProteinFamilyClassifier instance
    """
    # Load embeddings
    with h5py.File(embeddings_file, 'r') as f:
        embeddings = f['embeddings'][:]
        protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                      for pid in f['protein_ids'][:]]
    
    # Load metadata
    metadata = pd.read_csv(metadata_file, index_col=0)
    
    # Create classifier
    classifier = ProteinFamilyClassifier(method=method, **kwargs)
    
    # Classify proteins
    family_assignments = classifier.classify_proteins(embeddings, protein_ids, metadata=metadata)
    classifier.family_assignments = family_assignments
    
    # Save family index
    classifier.save_family_index(output_file)
    
    return classifier 


def validate_family_index(family_index_file: str) -> bool:
    """
    Validate that a family index file is properly formatted and contains required data.
    
    Args:
        family_index_file: Path to family index file
        
    Returns:
        True if valid, False otherwise
    """
    if not os.path.exists(family_index_file):
        logger.error(f"Family index file not found: {family_index_file}")
        return False
    
    try:
        with h5py.File(family_index_file, 'r') as f:
            # Check for required datasets
            required_datasets = ['family_centroids']
            for dataset in required_datasets:
                if dataset not in f:
                    logger.error(f"Required dataset '{dataset}' not found in family index")
                    return False
                if len(f[dataset]) == 0:
                    logger.error(f"Dataset '{dataset}' is empty")
                    return False
            
            # Check for optional but important datasets
            if 'protein_ids' in f and 'family_ids' in f:
                protein_ids = f['protein_ids'][:]
                family_ids = f['family_ids'][:]
                if len(protein_ids) == 0 or len(family_ids) == 0:
                    logger.error("Family assignments are empty")
                    return False
                if len(protein_ids) != len(family_ids):
                    logger.error("Mismatch between protein_ids and family_ids lengths")
                    return False
            
            logger.info(f"Family index file {family_index_file} is valid")
            return True
            
    except Exception as e:
        logger.error(f"Error validating family index: {e}")
        return False


def recreate_family_index_if_needed(embeddings_file: str,
                                   metadata_file: str,
                                   family_index_file: str,
                                   force_recreate: bool = False,
                                   **kwargs) -> bool:
    """
    Recreate family index if it's corrupted or missing.
    
    Args:
        embeddings_file: Path to embeddings file
        metadata_file: Path to metadata file
        family_index_file: Path to family index file
        force_recreate: Whether to force recreation even if file exists
        **kwargs: Additional parameters for family classification
        
    Returns:
        True if successful, False otherwise
    """
    if force_recreate or not os.path.exists(family_index_file) or not validate_family_index(family_index_file):
        logger.info("Recreating family index...")
        try:
            classifier = create_family_index(
                embeddings_file=embeddings_file,
                metadata_file=metadata_file,
                output_file=family_index_file,
                **kwargs
            )
            logger.info(f"Successfully recreated family index: {family_index_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to recreate family index: {e}")
            return False
    else:
        logger.info("Family index is valid, no recreation needed")
        return True 