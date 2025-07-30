import unittest
import numpy as np
import tempfile
import os
import h5py
import pandas as pd
from kbase_protein_query_module.src.storage import ProteinStorage, CompressedMetadataStorage

class TestProteinStorage(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.storage = ProteinStorage(base_dir=self.temp_dir, chunk_size=2)
        # Increase dataset size to meet FAISS clustering requirements (minimum 400 proteins)
        self.embeddings = np.random.randint(0, 256, size=(400, 8), dtype=np.uint8)
        self.protein_ids = [f'P{i:06d}' for i in range(400)]
        
    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)
        
    def test_store_and_load_family(self):
        family_id = 'test_family'
        self.storage.store_family_embeddings(family_id, self.embeddings, self.protein_ids)
        
        loaded_emb, loaded_ids = self.storage.load_family_embeddings(family_id)
        np.testing.assert_array_equal(loaded_emb, self.embeddings)
        self.assertEqual(loaded_ids, self.protein_ids)
        
    def test_stream_family_embeddings(self):
        family_id = 'test_family'
        self.storage.store_family_embeddings(family_id, self.embeddings, self.protein_ids)
        
        batches = list(self.storage.stream_family_embeddings(family_id, batch_size=100))
        self.assertEqual(len(batches), 4)  # 400 proteins / 100 batch size = 4 batches
        
        # Check first batch
        first_batch_emb, first_batch_ids = batches[0]
        self.assertEqual(first_batch_emb.shape[0], 100)
        self.assertEqual(len(first_batch_ids), 100)
        
    def test_metadata_storage(self):
        # Create larger metadata for testing
        N = 400  # Increased from 10 to meet FAISS requirements
        D = 64   # Binary dimension
        embeddings = np.random.randint(0, 256, size=(N, D // 8), dtype=np.uint8)
        protein_ids = [f'P{i:06d}' for i in range(N)]
        
        metadata = pd.DataFrame({
            'protein_id': protein_ids,
            'organism': ['Homo sapiens'] * N,
            'family': ['Kinase'] * N,
            'function': ['Signal Transduction'] * N
        }).set_index('protein_id')
        
        # Use the correct constructor signature
        storage = CompressedMetadataStorage(metadata_dir=self.temp_dir)
        storage.store_metadata(metadata, family_id='test_family')
        
        loaded_metadata = storage.load_metadata(family_id='test_family', protein_ids=protein_ids[:10])
        self.assertIsNotNone(loaded_metadata)
        self.assertGreater(len(loaded_metadata), 0)

    def test_artificial_families(self):
        from kbase_protein_query_module.src.storage import _create_artificial_families
        ids = [f'P{i:05d}' for i in range(7)]
        fams = list(_create_artificial_families(ids, max_family_size=3))
        self.assertEqual(len(fams), 3)
        self.assertEqual(len(fams[0][1]), 3)

    def test_missing_family(self):
        with self.assertRaises(FileNotFoundError):
            self.storage.load_family_embeddings('NOFAM')

if __name__ == '__main__':
    unittest.main() 