"""
Minimal binary SimilaritySearch using Hamming distance.

Data expectations (produced by test/generate_test_data.py):
- data/indexes/id_to_row.json: mapping uniprot_id -> row index (ordering of TSV)
- data/indexes/binary_embeddings.npy: uint8 matrix [N, D] where value>0 -> 1 else 0

If binary file missing, falls back to computing binary from storage embeddings in-memory.
"""

# Handle both script execution and module import - MUST BE FIRST
import sys
import os
if __name__ == "__main__":
    # Add parent directories to path for script execution
    current_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(os.path.dirname(current_dir))
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

import json
from typing import List, Tuple, Optional

import numpy as np

# Import storage - try relative first, then absolute
try:
    from .storage import ProteinStorage
except (ImportError, ValueError):
    from util.storage.storage import ProteinStorage


class SimilaritySearch:
    def __init__(self, index_dir: str, storage: ProteinStorage) -> None:
        self.index_dir = index_dir
        self.storage = storage
        self._binary_matrix: Optional[np.ndarray] = None
        self._ids: Optional[List[str]] = None
        self._load_binary()

    def _load_binary(self) -> None:
        bin_path = os.path.join(self.index_dir, "binary_embeddings.npy")
        if os.path.exists(bin_path):
            self._binary_matrix = np.load(bin_path).astype(np.uint8)
            self._ids = self.storage.list_proteins()
        else:
            # Compute on the fly
            ids = self.storage.list_proteins()
            mat = np.vstack([self.storage.get_embedding(pid) for pid in ids]).astype(np.float32)
            self._binary_matrix = (mat > 0).astype(np.uint8)
            self._ids = ids

    def _hamming_similarity(self, A_bin: np.ndarray, b_bin: np.ndarray) -> np.ndarray:
        # Similarity = 1 - (Hamming distance / D)
        D = float(A_bin.shape[1])
        # XOR then sum
        diffs = np.bitwise_xor(A_bin, b_bin).sum(axis=1).astype(np.float32)
        return 1.0 - (diffs / D)

    def query(self, embedding: np.ndarray, top_k: int = 5) -> List[Tuple[str, float]]:
        if embedding.dtype != np.float32:
            embedding = embedding.astype(np.float32)
        assert self._binary_matrix is not None and self._ids is not None
        q_bin = (embedding > 0).astype(np.uint8)
        sims = self._hamming_similarity(self._binary_matrix, q_bin)
        order = np.argsort(sims)[::-1][: max(1, int(top_k))]
        return [(self._ids[i], float(sims[i])) for i in order]

    def batch_query(self, embeddings: List[np.ndarray], top_k: int = 5) -> List[List[Tuple[str, float]]]:
        return [self.query(emb, top_k=top_k) for emb in embeddings]

    # No rebuild step: binary matrix is derived directly from embeddings (or precomputed file)


def main() -> int:
    # Simple self-test: ensure identical vector returns itself as top-1
    try:
        storage = ProteinStorage("data/embeddings/embeddings.tsv")
        search = SimilaritySearch(index_dir="data/indexes", storage=storage)
        pid0 = storage.list_proteins()[0]
        q = storage.get_embedding(pid0)
        res = search.query(q, top_k=10)
        if not res or res[0][0] != pid0:
            raise RuntimeError("Top-1 is not the exact id for identical query")
        print("SIMILARITY_OK")
        return 0
    except Exception as e:
        print(f"SIMILARITY_FAIL: {e}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


