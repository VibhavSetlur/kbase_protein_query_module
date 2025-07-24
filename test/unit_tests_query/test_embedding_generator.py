import unittest
import numpy as np
from kbase_protein_query_module_src.embedding_generator import ProteinEmbeddingGenerator

class TestProteinEmbeddingGenerator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Use a small ESM2 model for speed
        cls.generator = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D", device="cpu")
        cls.seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"

    def test_single_embedding_mean(self):
        emb = self.generator.generate_embedding(self.seq, pooling_method="mean")
        self.assertEqual(emb.shape[0], self.generator.embedding_dim)

    def test_single_embedding_cls(self):
        emb = self.generator.generate_embedding(self.seq, pooling_method="cls")
        self.assertEqual(emb.shape[0], self.generator.embedding_dim)

    def test_single_embedding_max(self):
        emb = self.generator.generate_embedding(self.seq, pooling_method="max")
        self.assertEqual(emb.shape[0], self.generator.embedding_dim)

    def test_batch_embeddings(self):
        seqs = [self.seq, self.seq[:20], self.seq[:10]]
        ids = ["A", "B", "C"]
        embs = self.generator.generate_embeddings_batch(seqs, ids, pooling_method="mean", batch_size=2)
        self.assertEqual(len(embs), 3)
        for e in embs.values():
            self.assertEqual(e.shape[0], self.generator.embedding_dim)

    def test_invalid_pooling(self):
        with self.assertRaises(ValueError):
            self.generator.generate_embedding(self.seq, pooling_method="invalid")

if __name__ == '__main__':
    unittest.main() 