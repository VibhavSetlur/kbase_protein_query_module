# Storage Utility

This module handles interactions with the KBase Protein Language Model (PLM) API and manages protein data storage.

## Components

*   **`ProteinStorage`**: The main class for querying homologs and managing protein data.

## Usage

```python
from kbase_protein_query_module.src.util.storage.storage import ProteinStorage

# Initialize
storage = ProteinStorage()

# Query similar proteins
results = storage.query_similar_proteins(
    proteins={
        "p1": {"sequence": "MKLL...", "original_id": "p1"}
    },
    max_hits=100
)
```
