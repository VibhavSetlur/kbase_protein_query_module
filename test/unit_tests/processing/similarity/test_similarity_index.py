import unittest
import numpy as np
import os
import sys
import tempfile
import h5py
from kbase_protein_query_module.src.processing.similarity.hierarchical_index import HierarchicalIndex

class TestSimilarityIndex(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.family_id = 'FAM0'  # Use real family ID (actual file name)
        
        # Use real data from data/families/ instead of synthetic data
        self.families_dir = "/home/vibhav/Downloads/Work/ANL/Research/kbase_protein_query_module/data/families"
        
        # Load real family data
        family_file = os.path.join(self.families_dir, f'{self.family_id}.h5')
        
        if not os.path.exists(family_file):
            self.fail(f"Real family data not found: {family_file}. Tests must use actual data from data/families/")
        
        # Load real embeddings and protein IDs
        with h5py.File(family_file, 'r') as f:
            self.embeddings = f['embeddings'][:400]  # Use first 400 for testing
            self.protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                               for pid in f['protein_ids'][:400]]
        
        self.index = HierarchicalIndex(base_dir=self.temp_dir)
        
    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)
        
    def test_search_family_float(self):
        # Test with proper data types using real data dimensions
        query = np.random.normal(0, 1, size=(self.embeddings.shape[1],)).astype(np.float32)
        
        # This should work without type errors
        try:
            # Create index first
            self.index.create_family_index_float(self.family_id, self.embeddings, self.protein_ids)
            # Then search
            distances, indices = self.index.search_family_float(self.family_id, query, top_k=5)
            self.assertGreater(len(distances), 0)
            self.assertEqual(len(distances), len(indices))
            self.assertLessEqual(len(distances), 5)
        except Exception as e:
            self.fail(f"Search failed with error: {e}")
            
    def test_search_family_binary(self):
        # Test binary search with real data converted to binary
        embeddings_uint8 = (self.embeddings > 0).astype(np.uint8)
        
        # This should work without clustering warnings
        try:
            # Create binary index first
            self.index.create_family_index(self.family_id, embeddings_uint8, self.protein_ids)
            # Then search
            query = np.random.randint(0, 2, size=(self.embeddings.shape[1],), dtype=np.uint8)
            distances, indices = self.index.search_family(self.family_id, query, top_k=5)
            self.assertGreater(len(distances), 0)
            self.assertEqual(len(distances), len(indices))
            self.assertLessEqual(len(distances), 5)
        except Exception as e:
            self.fail(f"Binary search failed with error: {e}")
            
    def test_search_with_large_dataset(self):
        # Test with real data that meets FAISS requirements
        # Use more proteins from the real data if available
        if len(self.protein_ids) >= 500:
            large_embeddings = self.embeddings[:500]
            large_protein_ids = self.protein_ids[:500]
        else:
            # Use all available real data
            large_embeddings = self.embeddings
            large_protein_ids = self.protein_ids
        
        # Create index with real dataset
        self.index.create_family_index_float(self.family_id, large_embeddings, large_protein_ids)
        
        # Test search
        query = np.random.normal(0, 1, size=(large_embeddings.shape[1],)).astype(np.float32)
        try:
            distances, indices = self.index.search_family_float(self.family_id, query, top_k=10)
            self.assertGreater(len(distances), 0)
            self.assertEqual(len(distances), len(indices))
            self.assertLessEqual(len(distances), 10)
        except Exception as e:
            self.fail(f"Large dataset search failed with error: {e}")
            
    def test_binary_search_with_large_dataset(self):
        # Test binary search with real data converted to binary
        # Use more proteins from the real data if available
        if len(self.protein_ids) >= 500:
            large_embeddings = self.embeddings[:500]
            large_protein_ids = self.protein_ids[:500]
        else:
            # Use all available real data
            large_embeddings = self.embeddings
            large_protein_ids = self.protein_ids
        
        # Convert real embeddings to binary
        large_binary_embeddings = (large_embeddings > 0).astype(np.uint8)
        
        # Create binary index with real dataset
        self.index.create_family_index(self.family_id, large_binary_embeddings, large_protein_ids)
        
        # Test binary search
        query = np.random.randint(0, 2, size=(large_binary_embeddings.shape[1],), dtype=np.uint8)
        try:
            distances, indices = self.index.search_family(self.family_id, query, top_k=10)
            self.assertGreater(len(distances), 0)
            self.assertEqual(len(distances), len(indices))
            self.assertLessEqual(len(distances), 10)
        except Exception as e:
            self.fail(f"Large binary dataset search failed with error: {e}")

if __name__ == '__main__':
    unittest.main() 