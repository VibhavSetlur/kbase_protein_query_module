import unittest
import tempfile
import shutil
import os
import pandas as pd
import numpy as np
from kbase_protein_network_analysis_toolkit.storage import ProteinStorage, CompressedMetadataStorage

class TestProteinStorage(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.storage = ProteinStorage(base_dir=self.temp_dir, chunk_size=2)
        self.family_id = 'FAM2'
        self.protein_ids = ['A', 'B', 'C', 'D']
        self.embeddings = np.random.randint(0, 256, size=(4, 4), dtype=np.uint8)
        self.metadata = pd.DataFrame({'desc': ['a', 'b', 'c', 'd']}, index=self.protein_ids)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_store_and_load_family(self):
        N = 10
        D = 320  # must be multiple of 8
        embeddings = np.random.randint(0, 256, size=(N, D // 8), dtype=np.uint8)
        protein_ids = [f'prot_{i}' for i in range(N)]
        self.storage.store_family_embeddings('test_family', embeddings, protein_ids)
        loaded_embeddings, loaded_ids = self.storage.load_family_embeddings('test_family')
        assert np.array_equal(embeddings, loaded_embeddings)
        assert protein_ids == loaded_ids

    def test_chunking(self):
        self.storage.store_family_embeddings(self.family_id, self.embeddings, self.protein_ids)
        # Should be able to stream in chunks
        batches = list(self.storage.stream_family_embeddings(self.family_id, batch_size=2))
        self.assertEqual(len(batches), 2)
        self.assertEqual(batches[0][0].shape[0], 2)

    def test_metadata_storage_and_index(self):
        meta_storage = CompressedMetadataStorage(metadata_dir=str(self.storage.metadata_dir))
        meta_storage.store_metadata(self.metadata, family_id=self.family_id)
        loaded = meta_storage.load_metadata(family_id=self.family_id, protein_ids=['A', 'C'])
        self.assertEqual(list(loaded.index), ['A', 'C'])

    def test_artificial_families(self):
        from kbase_protein_network_analysis_toolkit.storage import _create_artificial_families
        ids = [f'P{i:05d}' for i in range(7)]
        fams = list(_create_artificial_families(ids, max_family_size=3))
        self.assertEqual(len(fams), 3)
        self.assertEqual(len(fams[0][1]), 3)

    def test_missing_family(self):
        with self.assertRaises(FileNotFoundError):
            self.storage.load_family_embeddings('NOFAM')

if __name__ == '__main__':
    unittest.main() 