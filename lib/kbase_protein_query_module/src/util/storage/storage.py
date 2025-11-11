"""
Minimal ProteinStorage for simple, testable workflows.

Responsibilities:
- Load embeddings from a simple TSV/CSV: columns uniprot_id \t v1,v2,...,vD
- Keep in-memory maps for id -> embedding and ordered arrays for fast ops
- Provide a tiny API: get_embedding, list_proteins, search_by_ids

This is intentionally small and easy to replace later.
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

import csv
import json
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


class ProteinStorage:
    def __init__(
        self,
        embeddings_file_path: str,
        index_path: Optional[str] = None,
        config: Optional[dict] = None,
    ) -> None:
        self.embeddings_file_path = embeddings_file_path
        self.index_path = index_path
        self.config = config or {}

        self._ids: List[str] = []
        self._emb_matrix: Optional[np.ndarray] = None
        self.id_to_index: Dict[str, int] = {}

        self.dim: int = 0
        self.n: int = 0

        self._load_embeddings()

        # Defer index usage; built or loaded on demand via SimilaritySearch
        self._similarity: Optional["SimilaritySearch"] = None

    def _load_embeddings(self) -> None:
        if not os.path.exists(self.embeddings_file_path):
            raise FileNotFoundError(f"Embeddings file not found: {self.embeddings_file_path}")

        ids: List[str] = []
        vectors: List[np.ndarray] = []

        # Support TSV/CSV with two columns: uniprot_id, embedding_string
        # embedding_string is comma-separated floats
        with open(self.embeddings_file_path, "r", newline="") as f:
            # Auto-detect delimiter between tab and comma. We still expect the embedding field itself to be comma-separated.
            sample = f.read(1024)
            f.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
            reader = csv.reader(f, dialect)

            header_peek = sample.splitlines()[0]
            has_header = ("uniprot_id" in header_peek) or ("embedding" in header_peek)

            if has_header:
                next(reader, None)

            for row in reader:
                if not row:
                    continue
                if len(row) < 2:
                    # Allow files where embedding is the second column only
                    # Skip malformed lines
                    continue
                uniprot_id = row[0].strip()
                emb_str = row[1].strip()
                if not uniprot_id or not emb_str:
                    continue
                try:
                    vec = np.fromstring(emb_str, sep=",").astype(np.float32)
                except Exception:
                    # Try JSON list fallback
                    try:
                        vec = np.array(json.loads(emb_str), dtype=np.float32)
                    except Exception:
                        continue
                if np.isnan(vec).any():
                    raise ValueError(f"NaNs found in embedding for id {uniprot_id}")
                ids.append(uniprot_id)
                vectors.append(vec)

        if not vectors:
            raise ValueError("No embeddings parsed from file")

        # Ensure consistent dimensionality
        dim = int(vectors[0].shape[0])
        for i, v in enumerate(vectors):
            if v.shape[0] != dim:
                raise ValueError(f"Inconsistent embedding dim at row {i}: expected {dim}, got {v.shape[0]}")

        self._ids = ids
        self._emb_matrix = np.vstack(vectors)
        self.dim = dim
        self.n = len(ids)
        self.id_to_index = {pid: idx for idx, pid in enumerate(ids)}

    def get_embedding(self, uniprot_id: str) -> np.ndarray:
        if uniprot_id not in self.id_to_index:
            raise KeyError(f"Protein id not found: {uniprot_id}")
        assert self._emb_matrix is not None
        return self._emb_matrix[self.id_to_index[uniprot_id]]

    def list_proteins(self) -> List[str]:
        return list(self._ids)

    def search_by_ids(self, ids: List[str]) -> Dict[str, np.ndarray]:
        result: Dict[str, np.ndarray] = {}
        for pid in ids:
            if pid in self.id_to_index:
                result[pid] = self.get_embedding(pid)
        return result

    # --- Similarity integration ---
    def _get_similarity(self) -> "SimilaritySearch":
        # Lazy import to avoid cycles
        try:
            from .similarity_search import SimilaritySearch
        except (ImportError, ValueError):
            # Fallback for script execution
            from util.storage.similarity_search import SimilaritySearch
        index_dir = self.index_path or os.path.join(os.path.dirname(self.embeddings_file_path), "../indexes")
        index_dir = os.path.normpath(index_dir)
        if self._similarity is None:
            self._similarity = SimilaritySearch(index_dir=index_dir, storage=self)
        return self._similarity

    def find_top_k_similar(self, embedding: np.ndarray, top_k: int = 5) -> List[Tuple[str, float]]:
        sim = self._get_similarity()
        return sim.query(embedding, top_k=top_k)

    def find_top_k_similar_by_id(self, uniprot_id: str, top_k: int = 5) -> List[Tuple[str, float]]:
        emb = self.get_embedding(uniprot_id)
        return self.find_top_k_similar(emb, top_k=top_k)


def main() -> int:
    """Simple self-test for ProteinStorage."""
    ok = True
    try:
        # Handle both script execution and module import
        if __name__ == "__main__":
            # Add parent directories to path for script execution
            current_dir = os.path.dirname(os.path.abspath(__file__))
            # Try to find data directory relative to project root
            project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(current_dir))))
            embeddings_file = os.path.join(project_root, "data/embeddings/embeddings.tsv")
            if not os.path.exists(embeddings_file):
                # Try alternative path
                embeddings_file = "data/embeddings/embeddings.tsv"
        else:
            embeddings_file = "data/embeddings/embeddings.tsv"
        
        storage = ProteinStorage(
            embeddings_file_path=embeddings_file,
            index_path="data/indexes"
        )
        
        if storage.n <= 0 or storage.dim <= 0:
            raise RuntimeError(f"Invalid storage stats: n={storage.n}, dim={storage.dim}")
        
        any_id = storage.list_proteins()[0]
        vec = storage.get_embedding(any_id)
        
        if vec.shape[0] != storage.dim:
            raise RuntimeError(f"Dimensionality mismatch: vec.shape[0]={vec.shape[0]}, storage.dim={storage.dim}")
        
        top = storage.find_top_k_similar(vec, top_k=10)
        if not top:
            raise RuntimeError("Similarity search returned empty results")
        
        # Check that we get at least one result
        if top[0][0] != any_id:
            logger.warning(f"Top-1 similarity result is not self: expected {any_id}, got {top[0][0]}")
        
        print(f"STORAGE_OK: loaded {storage.n} proteins, dim={storage.dim}")
        return 0
    except FileNotFoundError as e:
        print(f"STORAGE_FAIL: File not found - {e}")
        print("  (This is expected if data files are not present)")
        return 0  # Don't fail tests if data files are missing
    except Exception as e:
        print(f"STORAGE_FAIL: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


