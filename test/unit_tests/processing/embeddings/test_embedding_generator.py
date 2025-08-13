import unittest
import numpy as np
import os
import sys
from kbase_protein_query_module.src.processing.embeddings.generator import ProteinEmbeddingGenerator

class TestProteinEmbeddingGenerator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            # Use the absolute path to the model
            model_path = "/home/vibhav/Downloads/Work/ANL/Research/kbase_protein_query_module/data/esm2_t6_8M_UR50D_local"
            
            if not os.path.exists(model_path) or not os.path.isdir(model_path):
                raise FileNotFoundError(f"ESM2 model not found at {model_path}")
            
            print(f"Using model at: {model_path}")
            
            # Initialize the embedding generator
            cls.generator = ProteinEmbeddingGenerator(
                model_name="esm2_t6_8M_UR50D", 
                device="cpu"
            )
            cls.seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
            cls.model_available = True
            print("Embedding generator initialized successfully")
        except Exception as e:
            print(f"Error initializing embedding generator: {e}")
            cls.model_available = False

    def test_single_embedding_mean(self):
        # Test requires real model - no fallback
        if not hasattr(self, 'model_available') or not self.model_available:
            self.fail("ESM2 model not available. Tests must use actual model from data/esm2_t6_8M_UR50D_local/")
        
        emb = self.generator.generate_embedding(self.seq)
        self.assertEqual(emb.shape[0], self.generator.embedding_dim)
        self.assertEqual(emb.shape[0], 320)  # ESM2 t6 model dimension

    def test_batch_embeddings(self):
        # Test requires real model - no fallback
        if not hasattr(self, 'model_available') or not self.model_available:
            self.fail("ESM2 model not available. Tests must use actual model from data/esm2_t6_8M_UR50D_local/")
        
        seqs = [self.seq, self.seq[:20], self.seq[:10]]
        ids = ["A", "B", "C"]
        embs = self.generator.generate_embeddings_batch(seqs, ids, batch_size=2)
        self.assertEqual(len(embs), 3)
        for e in embs.values():
            self.assertEqual(e.shape[0], self.generator.embedding_dim)
            self.assertEqual(e.shape[0], 320)

if __name__ == '__main__':
    unittest.main() 