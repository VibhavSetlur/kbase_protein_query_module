#!/usr/bin/env python3
"""
Script to create FAISS indexes for all existing family data.

This script will:
1. Load all family data from H5 files
2. Create FAISS IVF float32 indexes for each family
3. Store indexes in the correct directory structure
4. Create metadata files for each index
"""

import os
import sys
import h5py
import numpy as np
import json
import faiss
from pathlib import Path
import logging

# Add the lib directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from kbase_protein_query_module.src.storage import ProteinStorage

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_faiss_indexes():
    """Create FAISS indexes for all existing family data."""
    
    # Initialize storage
    storage = ProteinStorage(base_dir="data")
    
    # Get all family files
    families_dir = Path("data/families")
    if not families_dir.exists():
        logger.error("Families directory not found")
        return
    
    family_files = list(families_dir.glob("*.h5"))
    logger.info(f"Found {len(family_files)} family files")
    
    # Create indexes directory structure
    indexes_dir = Path("data/indexes")
    families_index_dir = indexes_dir / "families"
    metadata_dir = indexes_dir / "metadata"
    
    families_index_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    
    created_count = 0
    failed_count = 0
    
    for family_file in family_files:
        family_id = family_file.stem  # Remove .h5 extension
        
        # Check if index already exists
        index_file = families_index_dir / f"{family_id}.faiss"
        metadata_file = metadata_dir / f"{family_id}_float_metadata.json"
        
        if index_file.exists() and metadata_file.exists():
            logger.info(f"FAISS index already exists for family {family_id}")
            continue
        
        try:
            logger.info(f"Creating FAISS index for family {family_id}")
            
            # Load family data
            with h5py.File(family_file, 'r') as f:
                embeddings = f['embeddings'][:]
                protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][:]]
            
            # Ensure embeddings are float32
            if embeddings.dtype != np.float32:
                embeddings = embeddings.astype(np.float32)
            
            # Create FAISS index
            dimension = embeddings.shape[1]
            
            # Ensure we have enough training points for FAISS clustering
            # According to FAISS FAQ: minimum 39 training points per centroid
            nlist = min(10, max(1, len(embeddings)//10))
            min_training_points = nlist * 39
            
            if len(embeddings) < min_training_points:
                # Adjust nlist to meet training requirements
                nlist = max(1, len(embeddings) // 39)
                if nlist < 1:
                    # If we still don't have enough points, use a flat index instead
                    logger.warning(f"Not enough training points for IVF index in family {family_id}. Using flat index.")
                    index = faiss.IndexFlatL2(dimension)
                    index.add(embeddings)
                else:
                    # Use adjusted nlist
                    logger.info(f"Adjusted nlist to {nlist} for family {family_id} to meet training requirements")
                    quantizer = faiss.IndexFlatL2(dimension)
                    index = faiss.IndexIVFFlat(quantizer, dimension, nlist, faiss.METRIC_L2)
                    index.train(embeddings)
                    index.add(embeddings)
            else:
                # We have enough training points, proceed normally
                quantizer = faiss.IndexFlatL2(dimension)
                index = faiss.IndexIVFFlat(quantizer, dimension, nlist, faiss.METRIC_L2)
                index.train(embeddings)
                index.add(embeddings)
            
            # Save the index
            faiss.write_index(index, str(index_file))
            
            # Create metadata
            metadata = {
                'protein_ids': protein_ids,
                'index_type': 'faiss_float',
                'dimension': dimension,
                'num_proteins': len(protein_ids),
                'nlist': nlist,
                'metric': 'L2'
            }
            
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            logger.info(f"Created FAISS index for family {family_id}: {len(protein_ids)} proteins, {dimension} dim")
            created_count += 1
            
        except Exception as e:
            logger.error(f"Failed to create FAISS index for family {family_id}: {e}")
            failed_count += 1
            continue
    
    logger.info(f"FAISS index creation completed: {created_count} created, {failed_count} failed")
    
    # Create family mapping file
    family_mapping = {}
    for family_file in family_files:
        family_id = family_file.stem
        family_mapping[family_id] = str(family_file)
    
    mapping_file = indexes_dir / "family_mapping.json"
    with open(mapping_file, 'w') as f:
        json.dump(family_mapping, f, indent=2)
    
    logger.info(f"Created family mapping file: {mapping_file}")

if __name__ == "__main__":
    create_faiss_indexes() 