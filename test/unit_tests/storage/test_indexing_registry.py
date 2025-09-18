import numpy as np
import pytest

from lib.kbase_protein_query_module.src.util.storage.indexing_strategy import (
    IndexingConfig,
    IndexingStrategy,
    get_indexing_registry,
    register_indexing_strategy,
)


class _TinyStrategy(IndexingStrategy):
    def build_index(self, embeddings: np.ndarray, metadata=None):
        self.index = embeddings.copy()
        self.is_built = True
        return self.index

    def search(self, query_embedding: np.ndarray, k: int = 1, **kwargs):
        if not self.is_built:
            raise ValueError("Index not built")
        # naive best match by dot product
        sims = self.index @ query_embedding
        best = int(np.argmax(sims))
        return [(best, float(sims[best]))]

    def get_index_info(self):
        return {"strategy_name": "tiny", "index_size": len(self.index) if self.index is not None else 0}

    def save_index(self, filepath: str) -> bool:
        return True

    def load_index(self, filepath: str) -> bool:
        return True


def test_register_and_use_custom_strategy():
    @register_indexing_strategy("tiny")
    class _Registered(_TinyStrategy):
        pass

    reg = get_indexing_registry()
    cfg = IndexingConfig(strategy_name="tiny", distance_metric="dot_product")
    strat = reg.get_strategy("tiny", cfg)

    X = np.eye(4, dtype=np.float32)
    strat.build_index(X)

    res = strat.search(np.array([1, 0, 0, 0], dtype=np.float32), k=1)
    assert res[0][0] == 0


