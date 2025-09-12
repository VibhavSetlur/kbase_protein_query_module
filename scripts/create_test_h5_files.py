#!/usr/bin/env python3
"""
Create test H5 files for testing purposes.
This script creates mock H5 files to replace the Git LFS pointer files.
"""

import h5py
import numpy as np
import os
import sys

def create_test_h5_file(filepath, family_name):
    """Create a test H5 file with mock data."""
    try:
        with h5py.File(filepath, 'w') as f:
            # Create mock embeddings (320 dimensions for ESM2)
            num_proteins = 100
            embedding_dim = 320
            
            # Create mock protein embeddings
            embeddings = np.random.randn(num_proteins, embedding_dim).astype(np.float32)
            f.create_dataset('embeddings', data=embeddings)
            
            # Create mock protein IDs
            protein_ids = [f"{family_name}_protein_{i:03d}" for i in range(num_proteins)]
            f.create_dataset('protein_ids', data=[p.encode('utf-8') for p in protein_ids])
            
            # Create mock metadata
            f.attrs['family_name'] = family_name
            f.attrs['num_proteins'] = num_proteins
            f.attrs['embedding_dim'] = embedding_dim
            f.attrs['model_name'] = 'esm2_t6_8M_UR50D'
            
        print(f"Created test H5 file: {filepath}")
        return True
    except Exception as e:
        print(f"Error creating {filepath}: {e}")
        return False

def main():
    """Create test H5 files for all family files."""
    families_dir = "/scratch/vsetlur/kbase_protein_query_module/data/families"
    
    if not os.path.exists(families_dir):
        print(f"Families directory not found: {families_dir}")
        return 1
    
    # Get list of H5 files (currently LFS pointer files)
    h5_files = [f for f in os.listdir(families_dir) if f.endswith('.h5')]
    
    if not h5_files:
        print("No H5 files found in families directory")
        return 1
    
    print(f"Found {len(h5_files)} H5 files to create")
    
    success_count = 0
    for h5_file in h5_files:
        filepath = os.path.join(families_dir, h5_file)
        family_name = os.path.splitext(h5_file)[0]
        
        # Backup the original LFS pointer file
        backup_path = filepath + '.lfs_backup'
        if not os.path.exists(backup_path):
            os.rename(filepath, backup_path)
            print(f"Backed up LFS pointer file: {backup_path}")
        
        # Create the test H5 file
        if create_test_h5_file(filepath, family_name):
            success_count += 1
    
    print(f"\nSuccessfully created {success_count}/{len(h5_files)} test H5 files")
    return 0 if success_count == len(h5_files) else 1

if __name__ == "__main__":
    sys.exit(main())
