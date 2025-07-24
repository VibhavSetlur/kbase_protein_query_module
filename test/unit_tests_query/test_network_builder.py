import unittest
import numpy as np
import pandas as pd
import os
from unittest.mock import MagicMock
from kbase_protein_query_module_src.network_builder import DynamicNetworkBuilder
from kbase_protein_query_module_src.similarity_index import HierarchicalIndex

class TestNetworkBuilder(unittest.TestCase):
    def setUp(self):
        self.embeddings = np.random.normal(0, 1, size=(10, 8)).astype(np.float32)
        self.protein_ids = [f"P{i:03d}" for i in range(10)]
        self.metadata = pd.DataFrame({
            'Protein names': [f'Prot{i}' for i in range(10)],
            'Organism': ['E. coli']*10
        }, index=self.protein_ids)
        self.query_emb = np.random.normal(0, 1, size=(1, 8)).astype(np.float32)
        self.query_id = "QUERY"
        self.builder = DynamicNetworkBuilder(k_neighbors=3, similarity_threshold=0.1)
        # Mock HierarchicalIndex
        self.mock_index = MagicMock()
        # Always return 3 neighbors with increasing distance
        self.mock_index.search_family_float.side_effect = lambda fam, q, top_k: (np.arange(1, top_k+1, dtype=np.float32), np.arange(top_k))
        self.family_id = 'FAMX'

    def test_mutual_knn_network(self):
        # Provide mock index and family_id as required
        G = self.builder.build_mutual_knn_network(self.embeddings, self.protein_ids, family_id=self.family_id, index=self.mock_index)
        self.assertGreaterEqual(len(G.nodes), 10)
        self.assertGreater(len(G.edges), 0)

    def test_threshold_network(self):
        G = self.builder.build_threshold_network(self.embeddings, self.protein_ids, family_id=self.family_id, index=self.mock_index)
        self.assertGreaterEqual(len(G.nodes), 10)

    def test_hybrid_network(self):
        G = self.builder.build_hybrid_network(self.embeddings, self.protein_ids, family_id=self.family_id, index=self.mock_index)
        self.assertGreaterEqual(len(G.nodes), 10)

    def test_network_properties(self):
        G = self.builder.build_mutual_knn_network(self.embeddings, self.protein_ids, family_id=self.family_id, index=self.mock_index)
        props = self.builder.analyze_network_properties(G)
        self.assertIn('density', props)
        self.assertIn('average_degree', props)

    def test_visualization(self):
        # Only test that it runs and produces a file
        out = 'test_network.html'
        from kbase_protein_query_module_src.network_builder import visualize_interactive_protein_network
        fig, G = visualize_interactive_protein_network(self.embeddings, self.protein_ids, self.metadata, query_embedding=self.query_emb, query_protein_id=self.query_id, output_file=out)
        self.assertIsNotNone(fig)
        self.assertTrue(os.path.exists(out))
        os.remove(out)

if __name__ == '__main__':
    unittest.main() 