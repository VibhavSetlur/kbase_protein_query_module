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
        
        # Use real data from data/families/ instead of synthetic data
        self.data_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
        self.families_dir = os.path.join(self.data_dir, 'families')
        
        # Load real family data
        family_id = 'FAM0'  # Use a real family (actual file name)
        family_file = os.path.join(self.families_dir, f'{family_id}.h5')
        
        if not os.path.exists(family_file):
            self.fail(f"Real family data not found: {family_file}. Tests must use actual data from data/families/")
        
        # Load real embeddings and protein IDs
        with h5py.File(family_file, 'r') as f:
            self.embeddings = f['embeddings'][:400]  # Use first 400 for testing
            self.protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                               for pid in f['protein_ids'][:400]]
        
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
        # Use real data for metadata testing
        N = len(self.protein_ids)  # Use actual number of proteins from real data
        D = self.embeddings.shape[1]  # Use actual embedding dimension from real data
        
        metadata = pd.DataFrame({
            'protein_id': self.protein_ids,
            'organism': ['E. coli'] * N,  # Use realistic organism data
            'family': ['family_0'] * N,   # Use actual family ID
            'function': ['Unknown'] * N    # Use realistic function data
        }).set_index('protein_id')
        
        # Use the correct constructor signature
        storage = CompressedMetadataStorage(metadata_dir=self.temp_dir)
        storage.store_metadata(metadata, family_id='test_family')
        
        loaded_metadata = storage.load_metadata(family_id='test_family', protein_ids=self.protein_ids[:10])
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