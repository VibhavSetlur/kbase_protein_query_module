# Indexing and Similarity (Concise)

Where
- `processing/similarity/hierarchical_index.py`: in-memory/placeholder APIs
- `storage/protein_storage.py`: persistent family/embedding storage, IVF/PQ helpers

Patterns
- Family assignment: centroid-based (binary & advanced methods)
- Similarity search: FAISS float IVF, HNSW or streaming fallbacks

Add/adjust
- New index strategy: see `storage/indexing_strategy.py` and `register_indexing_strategy`
- Tune nlist/m/pq bits conservatively for server memory footprints

Outputs
- Top matches → CSV/TSV
- Network building consumes embeddings + matches
