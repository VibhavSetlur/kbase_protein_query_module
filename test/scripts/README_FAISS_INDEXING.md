# FAISS Indexing Guide for Protein Embeddings

This guide explains how to create FAISS indexes for protein embeddings databases, with support for family-based subsetting for efficient hierarchical similarity search.

## Overview

FAISS (Facebook AI Similarity Search) is a library for efficient similarity search and clustering of dense vectors. For protein embeddings, FAISS indexes enable fast similarity search across large databases.

### Key Concepts

1. **Flat Index**: Simple brute-force search, exact results
   - Best for: Small datasets (<10K proteins)
   - Speed: O(n) where n is number of proteins
   - Accuracy: 100% (exact)

2. **IVF Index**: Inverted File Index, approximate search
   - Best for: Medium-large datasets (10K-1M proteins)
   - Speed: O(n/k) where k is number of clusters (nlist)
   - Accuracy: ~95-99% (approximate, depends on nprobe)

3. **Family-based Indexing**: Hierarchical search strategy
   - Level 1: Map query to protein family using family centroids
   - Level 2: Search within selected family using family-specific index
   - Benefits: Faster search, better organization, scalable to millions of proteins

## Installation

```bash
# Install FAISS (CPU version)
pip install faiss-cpu

# Or GPU version (if CUDA available)
pip install faiss-gpu
```

## Usage

### Basic Usage: Create Indexes for All Proteins

```bash
python test/scripts/create_faiss_indexing_guide.py \
    --embeddings data/embeddings/embeddings.tsv \
    --output data/indexes
```

This creates a single index for all proteins in the embeddings file.

### Family-based Indexing

If you have protein families, create indexes for each family:

```bash
# Step 1: Create family mapping file (JSON format)
# family_mapping.json:
# {
#   "P12345": "kinase_enzyme",
#   "P67890": "kinase_enzyme",
#   "P11111": "G_protein_coupled_receptor",
#   ...
# }

# Step 2: Create indexes with family grouping
python test/scripts/create_faiss_indexing_guide.py \
    --embeddings data/embeddings/embeddings.tsv \
    --output data/indexes \
    --families data/families/family_mapping.json
```

This creates:
- One index per family in `data/indexes/families/`
- Family centroids index in `data/indexes/centroids/`
- Metadata files in `data/indexes/metadata/`

### Create Index for Specific Family

```bash
python test/scripts/create_faiss_indexing_guide.py \
    --embeddings data/embeddings/embeddings.tsv \
    --output data/indexes \
    --families data/families/family_mapping.json \
    --family "kinase_enzyme"
```

### Choose Index Type

```bash
# Force flat index (exact search)
python test/scripts/create_faiss_indexing_guide.py \
    --embeddings data/embeddings/embeddings.tsv \
    --output data/indexes \
    --index-type flat

# Force IVF index (approximate search, faster for large datasets)
python test/scripts/create_faiss_indexing_guide.py \
    --embeddings data/embeddings/embeddings.tsv \
    --output data/indexes \
    --index-type ivf

# Auto-select based on dataset size (default)
python test/scripts/create_faiss_indexing_guide.py \
    --embeddings data/embeddings/embeddings.tsv \
    --output data/indexes \
    --index-type auto
```

## Input Format

### Embeddings File (TSV)

The embeddings file should be in TSV format with two columns:

```
uniprot_id	embedding
P12345	0.1,0.2,0.3,...,0.9
P67890	0.2,0.3,0.4,...,1.0
...
```

- First column: Protein ID (e.g., UniProt ID)
- Second column: Comma-separated embedding values (floats)
- Header row is optional (will be auto-detected)

### Family Mapping File (JSON, Optional)

If using family-based indexing, provide a JSON file mapping protein IDs to family IDs:

```json
{
  "P12345": "kinase_enzyme",
  "P67890": "kinase_enzyme",
  "P11111": "G_protein_coupled_receptor",
  "P22222": "G_protein_coupled_receptor"
}
```

## Output Structure

After running the script, you'll have the following structure:

```
data/indexes/
├── families/
│   ├── kinase_enzyme.faiss
│   ├── G_protein_coupled_receptor.faiss
│   └── ...
├── metadata/
│   ├── kinase_enzyme_metadata.json
│   ├── G_protein_coupled_receptor_metadata.json
│   └── ...
├── centroids/
│   ├── family_centroids.faiss
│   └── family_centroids_metadata.json
├── family_mapping.json
└── indexing_summary.json
```

### Index Files (.faiss)

Binary FAISS index files that can be loaded for similarity search.

### Metadata Files (.json)

JSON files containing:
- Protein IDs in the index
- Number of proteins
- Embedding dimension
- Index type and parameters
- File paths

### Family Mapping (family_mapping.json)

JSON file mapping each protein ID to its family ID.

### Summary (indexing_summary.json)

Summary of the indexing process:
- Total number of proteins
- Total number of families
- List of families
- Number of indexes created
- Embedding dimension

## Using the Indexes

### Load and Search a Family Index

```python
import faiss
import numpy as np
import json

# Load index
index = faiss.read_index("data/indexes/families/kinase_enzyme.faiss")

# Load metadata
with open("data/indexes/metadata/kinase_enzyme_metadata.json", 'r') as f:
    metadata = json.load(f)

# Prepare query embedding (normalized)
query_embedding = np.array([...], dtype=np.float32)
faiss.normalize_L2(query_embedding.reshape(1, -1))

# Search
k = 10  # Number of results
distances, indices = index.search(query_embedding.reshape(1, -1), k)

# Get protein IDs
protein_ids = [metadata['protein_ids'][idx] for idx in indices[0]]
```

### Hierarchical Search (Family-based)

```python
import faiss
import numpy as np
import json

# Step 1: Load centroids index
centroids_index = faiss.read_index("data/indexes/centroids/family_centroids.faiss")
with open("data/indexes/centroids/family_centroids_metadata.json", 'r') as f:
    centroids_metadata = json.load(f)

# Step 2: Find most similar family
query_embedding = np.array([...], dtype=np.float32)
faiss.normalize_L2(query_embedding.reshape(1, -1))
distances, indices = centroids_index.search(query_embedding.reshape(1, -1), 1)
family_id = centroids_metadata['family_ids'][indices[0][0]]

# Step 3: Search within family
family_index = faiss.read_index(f"data/indexes/families/{family_id}.faiss")
with open(f"data/indexes/metadata/{family_id}_metadata.json", 'r') as f:
    family_metadata = json.load(f)

distances, indices = family_index.search(query_embedding.reshape(1, -1), k=10)
protein_ids = [family_metadata['protein_ids'][idx] for idx in indices[0]]
```

## Best Practices

### 1. Index Type Selection

- **< 10K proteins**: Use flat index (exact, fast enough)
- **10K - 100K proteins**: Use IVF index with nlist = sqrt(n)
- **> 100K proteins**: Use IVF index, consider increasing nlist

### 2. Family Organization

- Group related proteins together (e.g., by function, structure, or annotation)
- Aim for 100-10K proteins per family for optimal performance
- Too many small families: Overhead from centroids search
- Too few large families: Less benefit from hierarchical search

### 3. Normalization

- Always normalize embeddings (L2 normalization) for cosine similarity
- The script automatically normalizes embeddings before indexing
- Normalize query embeddings before search

### 4. IVF Parameters

- **nlist**: Number of clusters (cells)
  - Rule of thumb: nlist = sqrt(n_proteins)
  - Minimum: 1, Maximum: 4096
  - Requires at least 39 * nlist training points

- **nprobe**: Number of cells to search
  - Higher nprobe = more accurate but slower
  - Default: min(10, nlist)
  - For high accuracy: nprobe = nlist (exhaustive search)

### 5. Memory Considerations

- Flat index: ~4 bytes * dim * n_proteins
- IVF index: ~4 bytes * dim * n_proteins + overhead
- For 1M proteins with 320-dim embeddings: ~1.3 GB

## Troubleshooting

### Error: "faiss-cpu not installed"

```bash
pip install faiss-cpu
```

### Error: "Not enough training points for IVF"

IVF index requires at least 39 * nlist training points. Solutions:
- Use flat index for small families
- Reduce nlist (script does this automatically)
- Combine small families

### Error: "Inconsistent embedding dimension"

All embeddings must have the same dimension. Check your embeddings file for:
- Missing values
- Malformed rows
- Different embedding models

### Error: "NaN values in embedding"

Embeddings should not contain NaN values. Check your embedding generation process.

## Performance Tips

1. **Use GPU**: If available, use `faiss-gpu` for faster indexing and search
2. **Batch Processing**: Process multiple queries in batches
3. **Pre-compute**: Create indexes once, reuse for multiple queries
4. **Cache**: Keep indexes in memory for repeated searches
5. **Parallelize**: Create indexes for different families in parallel

## Examples

See `test/scripts/create_faiss_indexing_guide.py` for complete implementation and examples.

## References

- [FAISS Documentation](https://github.com/facebookresearch/faiss)
- [FAISS Wiki](https://github.com/facebookresearch/faiss/wiki)
- [Protein Embeddings](https://github.com/facebookresearch/esm)

