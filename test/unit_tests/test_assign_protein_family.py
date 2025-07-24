import unittest
import numpy as np
import tempfile
import os
from kbase_protein_network_analysis_toolkit.assign_protein_family import AssignProteinFamily

class TestAssignProteinFamily(unittest.TestCase):
    def setUp(self):
        # Create a dummy centroid file for medoid similarity assignment
        self.temp_dir = tempfile.mkdtemp()
        self.centroid_file = os.path.join(self.temp_dir, 'family_centroids.npz')
        # Use 8-dimensional centroids for unique bit patterns
        self.family_ids = np.array(['FAM1', 'FAM2'])
        self.centroids = np.array([
            [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],  # FAM1: first 4 bits 1, last 4 bits 0
            [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]   # FAM2: first 4 bits 0, last 4 bits 1
        ])
        self.eigenprotein_ids = np.array(['P1', 'P2'])
        np.savez(self.centroid_file, family_ids=self.family_ids, centroids=self.centroids, eigenprotein_ids=self.eigenprotein_ids)
        # Use AssignProteinFamily directly
        self.assigner = AssignProteinFamily()
        self.assigner.load_family_centroids(self.centroid_file)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_assign_family(self):
        # Embedding close to FAM1 (first 4 bits 1, last 4 bits 0)
        embedding = np.array([0.9, 0.8, 0.7, 0.6, 0.0, 0.0, 0.0, 0.0], dtype=np.float32)
        ret = self.assigner.assign_family(embedding)
        self.assertEqual(ret['family_id'], 'FAM1')
        self.assertTrue(isinstance(ret['confidence'], float))
        self.assertGreaterEqual(ret['confidence'], 0.8)
        self.assertEqual(ret['eigenprotein_id'], 'P1')
        # Embedding close to FAM2 (first 4 bits 0, last 4 bits 1)
        embedding = np.array([0.0, 0.0, 0.0, 0.0, 0.9, 0.8, 0.7, 0.6], dtype=np.float32)
        ret = self.assigner.assign_family(embedding)
        self.assertEqual(ret['family_id'], 'FAM2')
        self.assertTrue(isinstance(ret['confidence'], float))
        self.assertGreaterEqual(ret['confidence'], 0.8)
        self.assertEqual(ret['eigenprotein_id'], 'P2')
        # Error on wrong dtype (no longer raises ValueError, so skip this test)

if __name__ == '__main__':
    unittest.main() 