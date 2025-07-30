import unittest
import numpy as np
import tempfile
import os
from kbase_protein_query_module.src.similarity_index import HierarchicalIndex

class TestSimilarityIndex(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.family_id = 'test_family'
        # Increase dataset size to meet FAISS clustering requirements (minimum 400 proteins)
        self.embeddings = np.random.normal(0, 1, size=(400, 8)).astype(np.float32)
        self.protein_ids = [f'P{i:06d}' for i in range(400)]
        self.index = HierarchicalIndex(base_dir=self.temp_dir)
        
    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)
        
    def test_search_family_float(self):
        # Test with proper data types
        query = np.random.normal(0, 1, size=(8,)).astype(np.float32)
        
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
        # Test binary search with larger dataset
        embeddings_uint8 = np.random.randint(0, 2, size=(400, 8), dtype=np.uint8)
        
        # This should work without clustering warnings
        try:
            # Create binary index first
            self.index.create_family_index(self.family_id, embeddings_uint8, self.protein_ids)
            # Then search
            query = np.random.randint(0, 2, size=(8,), dtype=np.uint8)
            distances, indices = self.index.search_family(self.family_id, query, top_k=5)
            self.assertGreater(len(distances), 0)
            self.assertEqual(len(distances), len(indices))
            self.assertLessEqual(len(distances), 5)
        except Exception as e:
            self.fail(f"Binary search failed with error: {e}")
            
    def test_search_with_large_dataset(self):
        # Test with large dataset that meets FAISS requirements
        large_embeddings = np.random.normal(0, 1, size=(500, 320)).astype(np.float32)
        large_protein_ids = [f'P{i:06d}' for i in range(500)]
        
        # Create index with large dataset
        self.index.create_family_index_float(self.family_id, large_embeddings, large_protein_ids)
        
        # Test search
        query = np.random.normal(0, 1, size=(320,)).astype(np.float32)
        try:
            distances, indices = self.index.search_family_float(self.family_id, query, top_k=10)
            self.assertGreater(len(distances), 0)
            self.assertEqual(len(distances), len(indices))
            self.assertLessEqual(len(distances), 10)
        except Exception as e:
            self.fail(f"Large dataset search failed with error: {e}")
            
    def test_binary_search_with_large_dataset(self):
        # Test binary search with large dataset
        large_binary_embeddings = np.random.randint(0, 256, size=(500, 40), dtype=np.uint8)
        large_protein_ids = [f'P{i:06d}' for i in range(500)]
        
        # Create binary index with large dataset
        self.index.create_family_index(self.family_id, large_binary_embeddings, large_protein_ids)
        
        # Test binary search
        query = np.random.randint(0, 256, size=(40,), dtype=np.uint8)
        try:
            distances, indices = self.index.search_family(self.family_id, query, top_k=10)
            self.assertGreater(len(distances), 0)
            self.assertEqual(len(distances), len(indices))
            self.assertLessEqual(len(distances), 10)
        except Exception as e:
            self.fail(f"Large binary dataset search failed with error: {e}")

if __name__ == '__main__':
    unittest.main() 