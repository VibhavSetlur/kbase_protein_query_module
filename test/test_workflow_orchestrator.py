import unittest
import tempfile
import shutil
import os
import yaml
import numpy as np
from kbase_protein_network_analysis_toolkit.workflow_orchestrator import ProteinNetworkWorkflow

class TestProteinNetworkWorkflow(unittest.TestCase):
    def setUp(self):
        # Create a minimal config file and temp data dir
        self.temp_dir = tempfile.mkdtemp()
        self.config_file = os.path.join(self.temp_dir, 'config.yaml')
        config = {
            'storage': {'optimized_storage_dir': self.temp_dir},
            'embedding': {'model_name': 'esm2_t6_8M_UR50D', 'device': 'cpu'},
            'logging': {'log_file': os.path.join(self.temp_dir, 'test.log'), 'level': 'INFO'}
        }
        with open(self.config_file, 'w') as f:
            yaml.dump(config, f)
        # Create a fake family in storage
        from kbase_protein_network_analysis_toolkit.storage import ProteinStorage
        self.storage = ProteinStorage(base_dir=self.temp_dir)
        self.family_id = 'FAMW'
        self.protein_ids = ['X1', 'X2', 'X3']
        self.embeddings = np.random.randn(3, 8).astype(np.float32)
        self.storage.store_family_embeddings(self.family_id, self.embeddings, self.protein_ids)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_end_to_end_workflow(self):
        workflow = ProteinNetworkWorkflow(config_file=self.config_file)
        seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        result = workflow.run_optimized_workflow(seq, query_protein_id="X1", k_similar=2, network_method="mutual_knn", save_results=False)
        self.assertEqual(result['status'], 'success')
        self.assertIn('query_embedding', result)
        self.assertIn('family_id', result)
        self.assertIn('similar_proteins', result)
        self.assertIn('network', result)
        self.assertIn('network_properties', result)
        self.assertIn('performance_metrics', result)

    def test_missing_config(self):
        with self.assertRaises(FileNotFoundError):
            ProteinNetworkWorkflow(config_file='nonexistent.yaml')

if __name__ == '__main__':
    unittest.main() 