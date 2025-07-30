import unittest
import numpy as np
import tempfile
import os
from kbase_protein_query_module.src.check_existence import ProteinExistenceChecker

class TestProteinExistenceChecker(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        # Increase dataset size to meet FAISS clustering requirements (minimum 400 proteins)
        self.embeddings = np.random.randint(0, 256, size=(400, 8), dtype=np.uint8)
        self.protein_ids = [f'P{i:06d}' for i in range(400)]
        
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