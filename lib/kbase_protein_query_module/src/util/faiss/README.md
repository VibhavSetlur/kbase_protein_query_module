# FAISS Index Manager Utility

This utility provides a high-level interface for managing FAISS indexes for protein embeddings. It simplifies the creation, training, searching, and persistence of indexes, and integrates with the `ProteinEmbeddingGenerator` and `UniProt API` to build indexes from various inputs.

## Features

*   **Supported Index Types**: Flat, IVF, HNSW, BinaryFlat, BinaryIVF.
*   **Persistence**: Save and load indexes with associated metadata and ID mappings.
*   **Flexible Inputs**:
    *   Direct vectors (numpy arrays).
    *   Protein sequences (automatically embedded).
    *   UniProt IDs (automatically fetched and embedded).
    *   CSV/TSV files.
*   **GPU Support**: Optional GPU acceleration (requires `faiss-gpu`).

## Python API Usage

### 1. Initialization

```python
from lib.kbase_protein_query_module.src.util.faiss.faiss import FaissIndexManager
from lib.kbase_protein_query_module.src.util.embeddings.generator import ProteinEmbeddingGenerator

# Initialize embedding generator (optional, but needed for sequence/ID inputs)
generator = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D")

# Initialize manager
manager = FaissIndexManager(
    dimension=320, 
    index_type="Flat", 
    metric="L2",
    embedding_generator=generator
)
```

### 2. Adding Data

#### From Vectors
```python
import numpy as np
vectors = np.random.rand(100, 320).astype(np.float32)
ids = [f"protein_{i}" for i in range(100)]
manager.add_vectors(vectors, ids)
```

#### From Sequences
```python
sequences = ["MKLL...", "MVTV..."]
ids = ["seq1", "seq2"]
manager.add_from_sequences(sequences, ids)
```

#### From UniProt IDs
```python
uniprot_ids = ["P12345", "Q67890"]
manager.add_from_uniprot_ids(uniprot_ids)
```

#### From Files
```python
manager.build_from_files(
    file_paths=["/path/to/proteins.csv"],
    id_col="protein_id",
    seq_col="sequence"
)
```

### 3. Saving and Loading

```python
# Save
manager.save_index("/path/to/output/my_protein_index")

# Load
loaded_manager = FaissIndexManager.load_index(
    "/path/to/output/my_protein_index",
    embedding_generator=generator # Re-attach generator if needed
)
```

### 4. Searching

```python
query_vector = np.random.rand(1, 320).astype(np.float32)
distances, indices, external_ids = manager.search(query_vector, k=5)
```

## CLI Usage

```bash
python3 -m lib.kbase_protein_query_module.src.util.faiss.faiss \
    --mode create \
    --index_path /path/to/output_index \
    --files /path/to/proteins.csv \
    --uniprot_ids P12345 Q67890 \
    --dim 320
```

> **Note**: The file is named `faiss.py`. When importing it, ensure you use the full package path (e.g., `from lib.kbase_protein_query_module.src.util.faiss.faiss import ...`) or ensure the directory containing `faiss.py` is NOT in your `sys.path` to avoid conflicts with the `faiss` library.
