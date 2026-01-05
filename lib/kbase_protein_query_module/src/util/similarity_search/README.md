# Similarity Search Utility

This module provides general-purpose similarity search capabilities using cosine similarity on embeddings.

## Components

*   **`SimilaritySearch`**: The main class for computing similarities and finding top-K matches.

## Usage

```python
from kbase_protein_query_module.src.util.similarity_search.similarity_search import SimilaritySearch
import numpy as np

# Initialize
search = SimilaritySearch()

# Compute similarity
sim = search.cosine_similarity(emb1, emb2)

# Find top-K
results = search.find_top_k_similar(query_emb, target_embeddings, top_k=5)
```
