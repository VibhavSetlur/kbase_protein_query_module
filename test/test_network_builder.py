import unittest
import numpy as np
import pandas as pd
import os
from kbase_protein_network_analysis_toolkit.network_builder import DynamicNetworkBuilder, visualize_interactive_protein_network

class TestNetworkBuilder(unittest.TestCase):
    def setUp(self):
        self.embeddings = np.random.randn(10, 8).astype(np.float32)
        self.protein_ids = [f"P{i:03d}" for i in range(10)]
        self.metadata = pd.DataFrame({
            'Protein names': [f'Prot{i}' for i in range(10)],
            'Organism': ['E. coli']*10
        }, index=self.protein_ids)
        self.query_emb = np.random.randn(1, 8).astype(np.float32)
        self.query_id = "QUERY"
        self.builder = DynamicNetworkBuilder(k_neighbors=3, similarity_threshold=0.1)

    def test_mutual_knn_network(self):
        G = self.builder.build_mutual_knn_network(self.embeddings, self.protein_ids, self.query_emb, self.query_id)
        self.assertGreaterEqual(len(G.nodes), 11)
        self.assertGreater(len(G.edges), 0)

    def test_threshold_network(self):
        G = self.builder.build_threshold_network(self.embeddings, self.protein_ids, self.query_emb, self.query_id)
        self.assertGreaterEqual(len(G.nodes), 11)

    def test_hybrid_network(self):
        G = self.builder.build_hybrid_network(self.embeddings, self.protein_ids, self.query_emb, self.query_id)
        self.assertGreaterEqual(len(G.nodes), 11)

    def test_network_properties(self):
        G = self.builder.build_mutual_knn_network(self.embeddings, self.protein_ids)
        props = self.builder.analyze_network_properties(G)
        self.assertIn('density', props)
        self.assertIn('average_degree', props)

    def test_visualization(self):
        # Only test that it runs and produces a file
        out = 'test_network.html'
        fig, G = visualize_interactive_protein_network(self.embeddings, self.protein_ids, self.metadata, query_embedding=self.query_emb, query_protein_id=self.query_id, output_file=out)
        self.assertIsNotNone(fig)
        self.assertTrue(os.path.exists(out))
        os.remove(out)

if __name__ == '__main__':
    unittest.main() 