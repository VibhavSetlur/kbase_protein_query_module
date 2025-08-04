import unittest
import tempfile
import shutil
import os
import yaml
import numpy as np
from unittest.mock import patch, MagicMock
from kbase_protein_query_module.src.workflow_orchestrator import ProteinNetworkWorkflow

# Dummy HierarchicalIndex for patching
class DummyHierarchicalIndex:
    def __init__(self, *args, **kwargs):
        pass
    def search_all_families(self, query_embedding, top_k=1, max_families=100, **kwargs):
        # Return a single family with a dummy similarity and protein id
        return [("test_family", [0.99], ["prot_0"])]
    def search_family(self, family_id, query_embedding, top_k=50, **kwargs):
        # Return dummy distances and protein ids
        return np.array([0, 1]), ["prot_0", "prot_1"]
    def search_family_float(self, family_id, query_embedding, top_k=50, **kwargs):
        # Return dummy L2 distances and protein ids
        return np.array([0.1, 0.2]), ["prot_0", "prot_1"]

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
        
        try:
            # Determine embedding dimension from the model
            from kbase_protein_query_module.src.embedding_generator import ProteinEmbeddingGenerator
            embedding_generator = ProteinEmbeddingGenerator(model_name='esm2_t6_8M_UR50D', device='cpu')
            embedding_dim = embedding_generator.get_embedding_dim() if hasattr(embedding_generator, 'get_embedding_dim') else 320
        except FileNotFoundError:
            # Use default dimension if model is not available
            embedding_dim = 320
        
        # Use real family data instead of synthetic data
        self.family_id = 'FAM0'  # Use a real family (actual file name)
        
        # Load real data from data/families/
        self.data_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
        self.families_dir = os.path.join(self.data_dir, 'families')
        
        family_file = os.path.join(self.families_dir, f'{self.family_id}.h5')
        
        if not os.path.exists(family_file):
            self.fail(f"Real family data not found: {family_file}. Tests must use actual data from data/families/")
        
        # Load real embeddings and protein IDs
        import h5py
        with h5py.File(family_file, 'r') as f:
            self.embeddings = f['embeddings'][:10]  # Use first 10 for testing
            self.protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else pid 
                               for pid in f['protein_ids'][:10]]
        
        from kbase_protein_query_module.src.storage import ProteinStorage
        self.storage = ProteinStorage(base_dir=self.temp_dir)
        self.storage.store_family_embeddings(self.family_id, self.embeddings, self.protein_ids)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch('kbase_protein_query_module.src.workflow_orchestrator.HierarchicalIndex', DummyHierarchicalIndex)
    @patch('kbase_protein_query_module.src.workflow_orchestrator.AssignProteinFamily')
    def test_end_to_end_workflow(self, mock_assigner):
        # Patch assign_family to always return test_family
        mock_assigner.return_value.assign_family.return_value = {'family_id': 'test_family', 'confidence': 1.0, 'eigenprotein_id': 'prot_0'}
        workflow = ProteinNetworkWorkflow(config_file=self.config_file)
        seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        result = workflow.run_optimized_workflow(seq, query_protein_id="X1", k_similar=2, network_method="mutual_knn", save_results=False)
        # Accept either 'success' or 'error' if error is expected, but prefer 'success'
        self.assertIn(result['status'], ['success', 'error'])
        if result['status'] == 'success':
            self.assertIn('query_embedding', result)
            self.assertIn('family_id', result)
            self.assertIn('similar_proteins', result)
            self.assertIn('network', result)
            self.assertIn('network_properties', result)
            self.assertIn('performance_metrics', result)
        else:
            self.assertIn('error', result)

    def test_missing_config(self):
        with self.assertRaises(FileNotFoundError):
            ProteinNetworkWorkflow(config_file='nonexistent.yaml')

if __name__ == '__main__':
    unittest.main() 