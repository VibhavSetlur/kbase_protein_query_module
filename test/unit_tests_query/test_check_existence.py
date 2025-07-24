import unittest
import tempfile
import shutil
import os
import pandas as pd
import numpy as np
from kbase_protein_query_module_src.check_existence import ProteinExistenceChecker
from kbase_protein_query_module_src.storage import ProteinStorage, CompressedMetadataStorage

class TestProteinExistenceChecker(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for storage
        self.temp_dir = tempfile.mkdtemp()
        self.storage = ProteinStorage(base_dir=self.temp_dir)
        # Create a small family with embeddings and metadata
        self.family_id = 'FAM1'
        self.protein_ids = ['P00001', 'P00002', 'P00003']
        self.embeddings = np.random.randint(0, 256, size=(3, 8), dtype=np.uint8)
        self.metadata = pd.DataFrame({
            'Protein names': ['Alpha', 'Beta', 'Gamma'],
            'Organism': ['E. coli', 'E. coli', 'E. coli']
        }, index=self.protein_ids)
        self.storage.store_family_embeddings(self.family_id, self.embeddings, self.protein_ids, self.metadata)
        # Store compressed metadata
        self.metadata_storage = CompressedMetadataStorage(metadata_dir=str(self.storage.metadata_dir))
        self.metadata_storage.store_metadata(self.metadata, family_id=self.family_id)
        self.checker = ProteinExistenceChecker(storage=self.storage)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_protein_exists(self):
        result = self.checker.check_protein_existence('P00001')
        self.assertTrue(result['exists'])
        self.assertEqual(result['family_id'], self.family_id)
        self.assertIsInstance(result['metadata'], dict)
        self.assertEqual(result['metadata']['Protein names'], 'Alpha')

    def test_protein_not_exists(self):
        result = self.checker.check_protein_existence('P99999')
        self.assertFalse(result['exists'])
        self.assertIsNone(result['family_id'])
        self.assertIsNone(result['metadata'])

    def test_empty_protein_id(self):
        with self.assertRaises(ValueError):
            self.checker.check_protein_existence('')

    def test_corrupted_metadata(self):
        # Remove metadata file to simulate corruption
        meta_file = os.path.join(self.storage.metadata_dir, f'family_{self.family_id}_metadata.parquet')
        os.remove(meta_file)
        # Should still find by fallback
        result = self.checker.check_protein_existence('P00002')
        self.assertTrue(result['exists'])

    def test_performance_batch(self):
        import time
        ids = self.protein_ids * 100
        start = time.time()
        for pid in ids:
            self.checker.check_protein_existence(pid)
        elapsed = time.time() - start
        self.assertLess(elapsed, 2.0)  # Should be fast for small data

if __name__ == '__main__':
    unittest.main() 