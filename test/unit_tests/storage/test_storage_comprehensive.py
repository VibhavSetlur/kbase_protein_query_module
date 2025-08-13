"""
Comprehensive Storage Module Tests

This module tests all storage functionality including:
- Protein storage and retrieval
- Family assignment and management
- Protein existence checking
- Metadata storage and compression
- Memory-efficient loading
- Storage structure creation
"""

import unittest
import numpy as np
import os
import sys
import tempfile
import h5py
import pandas as pd
from pathlib import Path

from kbase_protein_query_module.src.storage import (
    ProteinStorage,
    ProteinFamilyAssigner,
    ProteinExistenceChecker,
    CompressedMetadataStorage,
    MemoryEfficientLoader,
    ProteinIDsIndex,
    create_storage_structure
)

class TestStorageComprehensive(unittest.TestCase):
    """Comprehensive tests for all storage functionality."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test environment with real data."""
        cls.temp_dir = tempfile.mkdtemp()
        
        # Use real data from data/families/
        cls.families_dir = "/home/vibhav/Downloads/Work/ANL/Research/kbase_protein_query_module/data/families"
        
        # Load real family data
        family_id = 'FAM0'
        family_file = os.path.join(cls.families_dir, f'{family_id}.h5')
        
        if not os.path.exists(family_file):
            cls.fail(f"Real family data not found: {family_file}. Tests must use actual data from data/families/")
        
        # Load real embeddings and protein IDs
        with h5py.File(family_file, 'r') as f:
            cls.embeddings = f['embeddings'][:300]  # Use first 300 for testing
            cls.protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                              for pid in f['protein_ids'][:300]]
        
        # Initialize storage components
        cls.storage = ProteinStorage(base_dir=cls.temp_dir)
        cls.family_assigner = ProteinFamilyAssigner()
        cls.existence_checker = ProteinExistenceChecker(base_dir=cls.temp_dir)
        cls.metadata_storage = CompressedMetadataStorage(metadata_dir=cls.temp_dir)
        cls.memory_loader = MemoryEfficientLoader(cls.storage)
        
        # Create test metadata
        cls.metadata = pd.DataFrame({
            'protein_id': cls.protein_ids[:100],
            'organism': ['E. coli'] * 100,
            'family': ['family_0'] * 100,
            'function': ['Unknown'] * 100,
            'length': [len(pid) for pid in cls.protein_ids[:100]]
        }).set_index('protein_id')
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(cls.temp_dir)
    
    def test_protein_storage_basic_operations(self):
        """Test basic protein storage operations."""
        family_id = 'test_family_basic'
        
        # Store family embeddings
        self.storage.store_family_embeddings(
            family_id, self.embeddings[:50], self.protein_ids[:50]
        )
        
        # Load family embeddings
        loaded_emb, loaded_ids = self.storage.load_family_embeddings(family_id)
        
        # Verify data integrity
        np.testing.assert_array_equal(loaded_emb, self.embeddings[:50])
        self.assertEqual(loaded_ids, self.protein_ids[:50])
        
        # Test partial loading
        partial_emb, partial_ids = self.storage.load_family_embeddings(
            family_id, start_idx=10, end_idx=30
        )
        self.assertEqual(len(partial_emb), 20)
        self.assertEqual(len(partial_ids), 20)
    
    def test_protein_storage_streaming(self):
        """Test streaming operations for large datasets."""
        family_id = 'test_family_streaming'
        
        # Store larger dataset
        self.storage.store_family_embeddings(
            family_id, self.embeddings[:200], self.protein_ids[:200]
        )
        
        # Test streaming with different batch sizes
        batch_sizes = [10, 50, 100]
        for batch_size in batch_sizes:
            batches = list(self.storage.stream_family_embeddings(family_id, batch_size=batch_size))
            
            total_proteins = sum(len(batch_ids) for _, batch_ids in batches)
            self.assertEqual(total_proteins, 200)
            
            # Verify first batch
            first_batch_emb, first_batch_ids = batches[0]
            self.assertLessEqual(len(first_batch_emb), batch_size)
            self.assertLessEqual(len(first_batch_ids), batch_size)
    
    def test_family_assignment(self):
        """Test protein family assignment functionality."""
        # Test with real embeddings
        test_embeddings = self.embeddings[:10]
        
        # Test family assignment for each embedding
        for i, embedding in enumerate(test_embeddings):
            result = self.family_assigner.assign_family(embedding)
            
            self.assertIn('family_id', result)
            self.assertIn('confidence', result)
            self.assertIn('similarity_score', result)
            self.assertIsInstance(result['family_id'], str)
            self.assertIsInstance(result['confidence'], float)
            self.assertGreaterEqual(result['confidence'], 0.0)
            self.assertLessEqual(result['confidence'], 1.0)
    
    def test_family_assignment_with_monitoring(self):
        """Test family assignment with monitoring."""
        test_embedding = self.embeddings[0]
        
        result = self.family_assigner.assign_family_with_monitoring(test_embedding)
        
        self.assertIn('family_id', result)
        self.assertIn('confidence', result)
        self.assertIn('processing_time', result)
        self.assertIn('memory_usage', result)
        self.assertIsInstance(result['processing_time'], float)
        self.assertIsInstance(result['memory_usage'], dict)
    
    def test_protein_existence_checking(self):
        """Test protein existence checking functionality."""
        # Test with real protein IDs
        test_protein_ids = self.protein_ids[:10]
        
        for protein_id in test_protein_ids:
            result = self.existence_checker.check_protein_existence(protein_id)
            
            self.assertIn('exists', result)
            self.assertIn('protein_id', result)
            self.assertIn('check_time', result)
            self.assertEqual(result['protein_id'], protein_id)
            self.assertIsInstance(result['exists'], bool)
    
    def test_metadata_storage_and_retrieval(self):
        """Test metadata storage and retrieval operations."""
        family_id = 'test_metadata_family'
        
        # Store metadata
        metadata_path = self.metadata_storage.store_metadata(
            self.metadata, family_id=family_id
        )
        self.assertTrue(os.path.exists(metadata_path))
        
        # Load all metadata
        loaded_metadata = self.metadata_storage.load_metadata(family_id=family_id)
        self.assertIsNotNone(loaded_metadata)
        self.assertEqual(len(loaded_metadata), 100)
        
        # Load specific protein metadata
        test_protein_ids = self.protein_ids[:5]
        specific_metadata = self.metadata_storage.load_metadata(
            family_id=family_id, protein_ids=test_protein_ids
        )
        self.assertEqual(len(specific_metadata), 5)
        
        # Verify data integrity
        for protein_id in test_protein_ids:
            self.assertIn(protein_id, specific_metadata.index)
    
    def test_memory_efficient_loading(self):
        """Test memory-efficient loading functionality."""
        # Store multiple families
        family_ids = ['family_1', 'family_2', 'family_3']
        for i, family_id in enumerate(family_ids):
            start_idx = i * 50
            end_idx = start_idx + 50
            self.storage.store_family_embeddings(
                family_id, 
                self.embeddings[start_idx:end_idx], 
                self.protein_ids[start_idx:end_idx]
            )
        
        # Test batch loading with memory constraints
        batches = list(self.memory_loader.load_families_batch(
            family_ids, max_memory_gb=0.1  # Very low memory limit
        ))
        
        self.assertEqual(len(batches), 3)
        for family_id, embeddings, protein_ids in batches:
            self.assertIn(family_id, family_ids)
            self.assertEqual(len(embeddings), 50)
            self.assertEqual(len(protein_ids), 50)
    
    def test_protein_ids_index(self):
        """Test protein IDs indexing functionality."""
        # Create index
        index = ProteinIDsIndex(base_dir=self.temp_dir)
        
        # Test protein search
        test_protein_id = self.protein_ids[0]
        result = index.search_protein(test_protein_id)
        
        if result is not None:
            self.assertIn('protein_id', result)
            self.assertIn('family_id', result)
            self.assertEqual(result['protein_id'], test_protein_id)
        
        # Test family retrieval
        family_proteins = index.get_proteins_by_family('family_0')
        self.assertIsInstance(family_proteins, list)
    
    def test_storage_structure_creation(self):
        """Test complete storage structure creation."""
        # Create test embeddings and metadata files
        test_embeddings_file = os.path.join(self.temp_dir, 'test_embeddings.h5')
        test_metadata_file = os.path.join(self.temp_dir, 'test_metadata.csv')
        
        # Save test embeddings
        with h5py.File(test_embeddings_file, 'w') as f:
            f.create_dataset('embeddings', data=self.embeddings[:100])
            f.create_dataset('protein_ids', data=np.array(self.protein_ids[:100], dtype='S'))
        
        # Save test metadata
        test_metadata = pd.DataFrame({
            'protein_id': self.protein_ids[:100],
            'organism': ['E. coli'] * 100,
            'family': ['family_0'] * 100,
            'function': ['Unknown'] * 100
        })
        test_metadata.to_csv(test_metadata_file, index=False)
        
        # Create storage structure
        storage = create_storage_structure(
            test_embeddings_file,
            test_metadata_file,
            output_dir=os.path.join(self.temp_dir, 'structured_storage'),
            family_column='family',
            max_family_size=50
        )
        
        # Verify structure was created
        self.assertIsInstance(storage, ProteinStorage)
        
        # Test that families were created
        family_list = storage.get_family_list()
        self.assertGreater(len(family_list), 0)
        
        # Test family statistics
        family_stats = storage.get_family_stats()
        self.assertIsInstance(family_stats, dict)
        self.assertGreater(len(family_stats), 0)
    
    def test_storage_error_handling(self):
        """Test error handling in storage operations."""
        # Test loading non-existent family
        with self.assertRaises(FileNotFoundError):
            self.storage.load_family_embeddings('non_existent_family')
        
        # Test with invalid embeddings
        invalid_embeddings = np.random.rand(5, 100)  # Wrong dimension
        invalid_ids = ["INVALID1", "INVALID2", "INVALID3", "INVALID4", "INVALID5"]
        
        # Should handle gracefully
        try:
            self.storage.store_family_embeddings(
                'invalid_family', invalid_embeddings, invalid_ids
            )
        except Exception as e:
            # Expected to fail due to dimension mismatch
            self.assertIsInstance(e, Exception)
    
    def test_storage_performance(self):
        """Test storage performance with larger datasets."""
        if len(self.protein_ids) >= 200:
            large_embeddings = self.embeddings[:200]
            large_protein_ids = self.protein_ids[:200]
            
            # Test storage performance
            import time
            start_time = time.time()
            
            self.storage.store_family_embeddings(
                'large_family', large_embeddings, large_protein_ids
            )
            
            storage_time = time.time() - start_time
            self.assertLess(storage_time, 10)  # Should complete within 10 seconds
            
            # Test loading performance
            start_time = time.time()
            
            loaded_emb, loaded_ids = self.storage.load_family_embeddings('large_family')
            
            load_time = time.time() - start_time
            self.assertLess(load_time, 5)  # Should complete within 5 seconds
    
    def test_storage_integration(self):
        """Test integration between different storage components."""
        family_id = 'integration_test_family'
        
        # Store embeddings
        self.storage.store_family_embeddings(
            family_id, self.embeddings[:50], self.protein_ids[:50]
        )
        
        # Store metadata
        family_metadata = self.metadata.loc[self.protein_ids[:50]]
        self.metadata_storage.store_metadata(family_metadata, family_id=family_id)
        
        # Test family assignment
        test_embedding = self.embeddings[0]
        assignment_result = self.family_assigner.assign_family(test_embedding)
        
        # Test existence checking
        test_protein_id = self.protein_ids[0]
        existence_result = self.existence_checker.check_protein_existence(test_protein_id)
        
        # All operations should work together
        self.assertIsNotNone(assignment_result)
        self.assertIsNotNone(existence_result)
        
        # Verify data consistency
        loaded_emb, loaded_ids = self.storage.load_family_embeddings(family_id)
        loaded_metadata = self.metadata_storage.load_metadata(family_id=family_id)
        
        self.assertEqual(len(loaded_ids), len(loaded_metadata))
        self.assertEqual(len(loaded_emb), len(loaded_metadata))

if __name__ == '__main__':
    unittest.main()
