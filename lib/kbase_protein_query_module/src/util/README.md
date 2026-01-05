# Utility Modules

This directory contains shared utility modules used across the KBase Protein Query Module.

## Components

*   **`storage/`**: Handles interactions with the KBase Protein Language Model (PLM) API, including querying for homologs and managing protein data.
*   **`similarity_search/`**: Provides general-purpose similarity search capabilities using cosine similarity on embeddings.
*   **`faiss/`**: (Optional) Wrappers for FAISS-based similarity search for efficient large-scale retrieval.
*   **`embeddings/`**: Utilities for handling protein embeddings.
*   **`uniprot/`**: Utilities for interacting with the UniProt API.

## Usage

These utilities are typically imported by higher-level managers (e.g., `AnalysisManager`, `InputManager`) but can be used independently.

```python
from kbase_protein_query_module.src.util.storage.storage import ProteinStorage
from kbase_protein_query_module.src.util.similarity_search.similarity_search import SimilaritySearch

# Initialize storage
storage = ProteinStorage()

# Initialize similarity search
search = SimilaritySearch()
```
