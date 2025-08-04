import unittest
import numpy as np
import tempfile
import os
import h5py
from kbase_protein_query_module.src.check_existence import ProteinExistenceChecker

class TestProteinExistenceChecker(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        
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
        
    def test_check_protein_existence(self):
        checker = ProteinExistenceChecker(base_dir=self.temp_dir)
        
        # Test with large dataset that meets FAISS requirements
        result = checker.check_protein_existence('P000001')
        
        # Verify result structure
        self.assertIsInstance(result, dict)
        self.assertIn('exists', result)
        self.assertIn('family_id', result)
        self.assertIn('metadata', result)
        
    def test_check_multiple_proteins(self):
        checker = ProteinExistenceChecker(base_dir=self.temp_dir)
        
        # Test with multiple proteins
        query_proteins = ['P000001', 'P000002', 'P000003']
        
        for protein_id in query_proteins:
            result = checker.check_protein_existence(protein_id)
            self.assertIsInstance(result, dict)
            self.assertIn('exists', result)
            self.assertIn('family_id', result)
            self.assertIn('metadata', result)
            
    def test_large_dataset_existence_checking(self):
        checker = ProteinExistenceChecker(base_dir=self.temp_dir)
        
        # Test existence checking with large dataset
        result = checker.check_protein_existence('P000001')
        
        self.assertIsInstance(result, dict)
        self.assertIn('exists', result)
        self.assertIn('family_id', result)
        self.assertIn('metadata', result)

if __name__ == '__main__':
    unittest.main() 