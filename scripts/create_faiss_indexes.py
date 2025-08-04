#!/usr/bin/env python3
"""
FAISS Index Creation Script for Protein Network Analysis

This script creates FAISS indexes for protein families with proper architecture:
- Binary FAISS for family centroids (mapping query to family)
- Float FAISS for within-family similarity search

Usage:
    python scripts/create_faiss_indexes.py [--family FAMILY_ID] [--all] [--force]
"""

import os
import sys
import argparse
import logging
import h5py
import numpy as np
import json
from pathlib import Path
from typing import List, Dict, Optional
import faiss

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from kbase_protein_query_module.src.storage import ProteinStorage
from kbase_protein_query_module.src.similarity_index import HierarchicalIndex

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FAISSIndexCreator:
    """
    Creates FAISS indexes for protein families with proper architecture.
    
    Architecture:
    - Binary FAISS: Used only for family centroids to map query proteins to families
    - Float FAISS: Used for within-family similarity search (more accurate)
    """
    
    def __init__(self, data_dir: str = "data"):
        self.data_dir = Path(data_dir)
        self.families_dir = self.data_dir / "families"
        self.indexes_dir = self.data_dir / "indexes"
        self.centroids_dir = self.data_dir / "family_centroids"
        
        # Create directories if they don't exist
        self.indexes_dir.mkdir(parents=True, exist_ok=True)
        self.centroids_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize storage
        self.storage = ProteinStorage(base_dir=str(self.data_dir))
    
    def get_available_families(self) -> List[str]:
        """Get list of available family files."""
        if not self.families_dir.exists():
            logger.warning(f"Families directory not found: {self.families_dir}")
            return []
        
        families = []
        for h5_file in self.families_dir.glob("*.h5"):
            family_id = h5_file.stem
            families.append(family_id)
        
        logger.info(f"Found {len(families)} families: {families}")
        return families
    
    def create_family_centroids_binary(self, force: bool = False) -> str:
        """
        Create binary FAISS index for family centroids.
        This is used to map query proteins to their most similar family.
        
        Args:
            force: Force recreation even if file exists
            
        Returns:
            Path to centroids file
        """
        centroids_file = self.centroids_dir / "family_centroids_binary.npz"
        
        if centroids_file.exists() and not force:
            logger.info(f"Centroids file already exists: {centroids_file}")
            return str(centroids_file)
        
        logger.info("Creating family centroids binary index...")
        
        families = self.get_available_families()
        if not families:
            logger.error("No families found to create centroids")
            return None
        
        centroids = []
        family_ids = []
        eigenprotein_ids = []
        
        for family_id in families:
            try:
                family_file = self.families_dir / f"{family_id}.h5"
                with h5py.File(family_file, 'r') as f:
                    embeddings = f['embeddings'][:]
                    protein_ids = f['protein_ids'][:]
                    
                    # Calculate centroid
                    centroid = np.mean(embeddings, axis=0)
                    centroids.append(centroid)
                    family_ids.append(family_id)
                    eigenprotein_ids.extend(protein_ids)
                    
                    logger.info(f"Created centroid for {family_id}: {len(protein_ids)} proteins")
                    
            except Exception as e:
                logger.error(f"Failed to create centroid for {family_id}: {e}")
        
        if not centroids:
            logger.error("No centroids created")
            return None
        
        # Save centroids with binary format for FAISS
        centroids_array = np.array(centroids, dtype=np.float32)
        
        # Create binary centroids (sign function)
        binary_centroids = (centroids_array > 0).astype(np.uint8)
        
        np.savez_compressed(
            centroids_file,
            centroids=binary_centroids,
            family_ids=np.array(family_ids, dtype='<U50'),
            eigenprotein_ids=np.array(eigenprotein_ids, dtype='<U50')
        )
        
        logger.info(f"Created binary centroids file: {centroids_file}")
        logger.info(f"Centroids shape: {binary_centroids.shape}")
        
        return str(centroids_file)
    
    def create_family_float_index(self, family_id: str, force: bool = False) -> Optional[str]:
        """
        Create float FAISS index for within-family similarity search.
        This provides more accurate similarity search within a family.
        
        Args:
            family_id: Family identifier
            force: Force recreation even if index exists
            
        Returns:
            Path to index file or None if failed
        """
        family_file = self.families_dir / f"{family_id}.h5"
        index_file = self.indexes_dir / "families" / f"{family_id}.faiss"
        metadata_file = self.indexes_dir / "metadata" / f"{family_id}_float_metadata.json"
        
        # Create subdirectories
        index_file.parent.mkdir(parents=True, exist_ok=True)
        metadata_file.parent.mkdir(parents=True, exist_ok=True)
        
        if index_file.exists() and not force:
            logger.info(f"Float index already exists for {family_id}: {index_file}")
            return str(index_file)
        
        try:
            # Load family data
            with h5py.File(family_file, 'r') as f:
                embeddings = f['embeddings'][:]
                protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][:]]
            
            logger.info(f"Creating float FAISS index for {family_id}: {len(protein_ids)} proteins, {embeddings.shape[1]} dimensions")
            
            # Ensure embeddings are float32
            embeddings = embeddings.astype(np.float32)
            
            # Create FAISS IVF index for efficient similarity search
            dimension = embeddings.shape[1]
            nlist = min(100, max(1, len(embeddings) // 10))  # Number of clusters
            
            # Ensure we have enough training points
            min_training_points = nlist * 39
            if len(embeddings) < min_training_points:
                nlist = max(1, len(embeddings) // 39)
                logger.info(f"Adjusted nlist to {nlist} for family {family_id} to meet training requirements")
            
            # Create quantizer and index
            quantizer = faiss.IndexFlatL2(dimension)
            index = faiss.IndexIVFFlat(quantizer, dimension, nlist)
            
            # Train the index
            index.train(embeddings)
            index.add(embeddings)
            
            # Save index
            faiss.write_index(index, str(index_file))
            
            # Save metadata
            metadata = {
                'family_id': family_id,
                'num_proteins': len(protein_ids),
                'embedding_dim': dimension,
                'index_type': 'faiss_ivf_float',
                'nlist': nlist,
                'protein_ids': protein_ids
            }
            
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            logger.info(f"Created float FAISS index for {family_id}: {len(protein_ids)} proteins, {dimension} dim")
            return str(index_file)
            
        except Exception as e:
            logger.error(f"Failed to create float index for {family_id}: {e}")
            return None
    
    def create_all_family_indexes(self, force: bool = False) -> Dict[str, str]:
        """
        Create float FAISS indexes for all families.
        
        Args:
            force: Force recreation of existing indexes
            
        Returns:
            Dictionary mapping family_id to index file path
        """
        families = self.get_available_families()
        created_indexes = {}
        
        for family_id in families:
            index_path = self.create_family_float_index(family_id, force=force)
            if index_path:
                created_indexes[family_id] = index_path
        
        logger.info(f"Created {len(created_indexes)} family indexes")
        return created_indexes
    
    def create_family_mapping_file(self) -> str:
        """Create family mapping file for easy lookup."""
        families = self.get_available_families()
        mapping = {}
        
        for family_id in families:
            family_file = self.families_dir / f"{family_id}.h5"
            index_file = self.indexes_dir / "families" / f"{family_id}.faiss"
            
            if family_file.exists():
                mapping[family_id] = str(family_file)
        
        mapping_file = self.indexes_dir / "family_mapping.json"
        with open(mapping_file, 'w') as f:
            json.dump(mapping, f, indent=2)
        
        logger.info(f"Created family mapping file: {mapping_file}")
        return str(mapping_file)
    
    def create_all_indexes(self, force: bool = False) -> bool:
        """
        Create all FAISS indexes (centroids + family indexes).
        
        Args:
            force: Force recreation of existing indexes
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Create family centroids (binary)
            centroids_file = self.create_family_centroids_binary(force=force)
            if not centroids_file:
                logger.error("Failed to create centroids")
                return False
            
            # Create family indexes (float)
            family_indexes = self.create_all_family_indexes(force=force)
            if not family_indexes:
                logger.error("Failed to create family indexes")
                return False
            
            # Create family mapping
            self.create_family_mapping_file()
            
            logger.info("Successfully created all FAISS indexes")
            logger.info(f"Centroids: {centroids_file}")
            logger.info(f"Family indexes: {len(family_indexes)}")
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to create indexes: {e}")
            return False

def main():
    parser = argparse.ArgumentParser(description="Create FAISS indexes for protein families")
    parser.add_argument("--family", help="Create index for specific family")
    parser.add_argument("--all", action="store_true", help="Create indexes for all families")
    parser.add_argument("--force", action="store_true", help="Force recreation of existing indexes")
    parser.add_argument("--data-dir", default="data", help="Data directory")
    
    args = parser.parse_args()
    
    creator = FAISSIndexCreator(data_dir=args.data_dir)
    
    if args.family:
        # Create index for specific family
        index_path = creator.create_family_float_index(args.family, force=args.force)
        if index_path:
            logger.info(f"Successfully created index for {args.family}: {index_path}")
        else:
            logger.error(f"Failed to create index for {args.family}")
            sys.exit(1)
    
    elif args.all:
        # Create all indexes
        success = creator.create_all_indexes(force=args.force)
        if not success:
            logger.error("Failed to create all indexes")
            sys.exit(1)
    
    else:
        # Show available families
        families = creator.get_available_families()
        if families:
            logger.info(f"Available families: {families}")
            logger.info("Use --family FAMILY_ID to create index for specific family")
            logger.info("Use --all to create indexes for all families")
        else:
            logger.error("No families found")

if __name__ == "__main__":
    main() 