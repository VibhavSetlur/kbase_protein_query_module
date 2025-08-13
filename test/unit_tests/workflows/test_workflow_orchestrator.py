import unittest
import tempfile
import shutil
import os
import yaml
import numpy as np
from unittest.mock import patch, MagicMock
from kbase_protein_query_module.src.workflows import ProteinQueryWorkflow as ProteinNetworkWorkflow, WorkflowResult

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
        # Create a minimal config and temp data dir
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a simple PipelineConfig for testing
        from kbase_protein_query_module.src.core import PipelineConfig
        self.config = PipelineConfig(
            input_proteins=[{'protein_id': 'test_protein', 'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'}],
            storage_config={'base_dir': self.temp_dir},
            embedding_model='esm2_t6_8M_UR50D',
            embedding_device='cpu'
        )
        
        try:
            # Determine embedding dimension from the model
            from kbase_protein_query_module.src.processing.embeddings.generator import ProteinEmbeddingGenerator
            embedding_generator = ProteinEmbeddingGenerator(model_name='esm2_t6_8M_UR50D', device='cpu')
            embedding_dim = embedding_generator.embedding_dim if hasattr(embedding_generator, 'embedding_dim') else 320
        except FileNotFoundError:
            # Use default dimension if model is not available
            embedding_dim = 320
        
        # Use real family data instead of synthetic data
        self.family_id = 'FAM0'  # Use a real family (actual file name)
        
        # Load real data from data/families/
        self.families_dir = "/home/vibhav/Downloads/Work/ANL/Research/kbase_protein_query_module/data/families"
        
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

    def test_end_to_end_workflow(self):
        # Test basic workflow execution
        workflow = ProteinNetworkWorkflow(config=self.config)
        
        # Test with simple input data
        input_data = {
            'protein_id': 'test_protein',
            'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'
        }
        
        result = workflow.execute(input_data)
        
        # Check that we got a WorkflowResult
        self.assertIsInstance(result, WorkflowResult)
        self.assertIsInstance(result.success, bool)
        self.assertIsInstance(result.run_id, str)
        self.assertIsInstance(result.stages_completed, list)

    def test_missing_config(self):
        # Test that workflow works with a minimal config
        from kbase_protein_query_module.src.core import PipelineConfig
        minimal_config = PipelineConfig(
            input_proteins=[{'protein_id': 'test_protein', 'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'}]
        )
        workflow = ProteinNetworkWorkflow(config=minimal_config)
        self.assertIsNotNone(workflow.config)

if __name__ == '__main__':
    unittest.main() 