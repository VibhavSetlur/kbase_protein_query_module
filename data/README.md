# Data Directory

This directory contains reference data for the kbase_protein_query_module.

## Structure

```
data/
├── embeddings/
│   └── embeddings.tsv          # Protein embeddings (TSV format)
└── indexes/
    ├── id_to_row.json          # UniProt ID to row index mapping
    └── binary_embeddings.npy   # Binary embedding matrix (optional)
```

## Test Data vs Reference Data

### Test Data (Small, in Docker Image)

- **Location**: `/kb/module/data/` (inside Docker image)
- **Size**: < 100 MB
- **Purpose**: Testing and development
- **Generation**: Automatically generated during Docker build
- **Status**: Included in Docker image

### Reference Data (Large, Mounted at Runtime)

- **Location**: Mounted at runtime (e.g., `/data`, `/kb/data`, `/refdata`)
- **Size**: TB scale
- **Purpose**: Production use
- **Setup**: See `docs/REFERENCE_DATA_SETUP.md`
- **Status**: Not in Docker image (too large)

## Setup for Reference Data

When using large reference data, you have three options:

### Option 1: Environment Variable (Recommended)

```bash
export KB_DATA_DIR=/path/to/reference/data
```

### Option 2: Standard Mount Point

Mount reference data at:
- `/data/embeddings/embeddings.tsv`
- `/kb/data/embeddings/embeddings.tsv`
- `/refdata/embeddings/embeddings.tsv`

### Option 3: Custom Path

Set in module configuration:
```python
config = {
    'embeddings_file': '/custom/path/embeddings.tsv',
    'index_path': '/custom/path/indexes'
}
```

## File Formats

See `docs/REFERENCE_DATA_SETUP.md` for detailed file format specifications.

## Troubleshooting

If data files are not found:

1. Check environment variables: `KB_DATA_DIR`, `KB_REFDATA_DIR`
2. Check mount points: `/data`, `/kb/data`, `/refdata`
3. Check file permissions
4. See logs for path resolution details
5. See `docs/REFERENCE_DATA_SETUP.md` for detailed troubleshooting

