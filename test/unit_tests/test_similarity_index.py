import unittest
import tempfile
import shutil
import numpy as np
import os
from kbase_protein_network_analysis_toolkit.similarity_index import HierarchicalIndex, StreamingIndex

class TestSimilarityIndex(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.embeddings = np.random.normal(0, 1, size=(10, 8)).astype(np.float32)
        self.protein_ids = [f"P{i:03d}" for i in range(10)]
        self.index = HierarchicalIndex(base_dir=self.temp_dir, index_type='faiss', quantization='none', cache_size=2)
        self.family_id = 'FAMX'
        self.index.create_family_index_float(self.family_id, self.embeddings, self.protein_ids)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_search_family(self):
        query = np.random.normal(0, 1, size=(8,)).astype(np.float32)
        sims, ids = self.index.search_family_float(self.family_id, query, top_k=5)
        self.assertEqual(len(ids), 5)
        self.assertEqual(len(sims), 5)
        # Error on wrong dtype
        with self.assertRaises(ValueError):
            self.index.search_family_float(self.family_id, np.random.randint(0, 256, size=(8,), dtype=np.uint8), top_k=5)

    def test_cache_eviction(self):
        # Fill cache and check eviction
        # Use uint8 embeddings for binary index
        embeddings_uint8 = np.random.randint(0, 2, size=(10, 8), dtype=np.uint8)
        for i in range(5):
            fam = f'FAM{i}'
            self.index.create_family_index(fam, embeddings_uint8, self.protein_ids)
            self.index.search_family(fam, np.random.randint(0, 2, size=(8,), dtype=np.uint8), top_k=1)
        self.assertLessEqual(len(self.index._family_cache), 2)

    def test_streaming_index(self):
        streaming = StreamingIndex(storage_dir=self.temp_dir, batch_size=5)
        # Simulate family assignments
        fam_assign = {pid: self.family_id for pid in self.protein_ids}
        emb_file = os.path.join(self.temp_dir, 'embeddings.h5')
        import h5py
        with h5py.File(emb_file, 'w') as f:
            f.create_dataset('embeddings', data=self.embeddings)
            f.create_dataset('protein_ids', data=self.protein_ids, dtype=h5py.special_dtype(vlen=str))
        streaming_file = streaming.create_streaming_index(emb_file, self.protein_ids, fam_assign)
        query = np.random.randint(0, 256, size=(8,), dtype=np.uint8)
        results = streaming.stream_search(query, emb_file, streaming_file, top_k=3)
        self.assertLessEqual(len(results), 3)

if __name__ == '__main__':
    unittest.main() 