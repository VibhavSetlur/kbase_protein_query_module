# Embedding Generator Utility

This module provides utilities for generating protein embeddings using ESM-2 models.

## Components

*   **`ProteinEmbeddingGenerator`**: Generates embeddings using local ESM models or the KBase PLM API.

## Usage

```python
from kbase_protein_query_module.src.util.embeddings.generator import ProteinEmbeddingGenerator

# Initialize (uses server by default if configured, or local if available)
generator = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D")

# Generate embedding
embedding = generator.generate_embedding("MKLL...")
```
