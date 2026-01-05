# UniProt API Utility

This module provides lightweight functions to fetch sequences and metadata for UniProt IDs.

## Components

*   **`api.py`**: Contains functions like `fetch_sequences`, `fetch_metadata`, etc.

## Usage

```python
from kbase_protein_query_module.src.util.uniprot import api

# Fetch sequences
sequences = api.fetch_sequences(["P12345", "Q67890"])

# Fetch metadata
metadata = api.fetch_metadata(["P12345"])
```
