# Generating Real UniProt Test Data

This script generates real test data from UniProt using actual protein sequences and embeddings.

## Overview

The `generate_test_data.py` script:
1. Queries UniProt for proteins in a specified PFAM family
2. Fetches real protein sequences
3. Generates real embeddings using ESM-2 models
4. Creates the data files needed for testing (`embeddings.tsv`, indexes)

## Why Real Data?

- **Metadata works**: Real UniProt IDs allow metadata fetching to work in tests
- **Realistic testing**: Tests run against actual protein data, not synthetic data
- **Small size**: Uses a small PFAM family to keep data size manageable for GitHub
- **Reproducible**: Uses deterministic parameters (reviewed proteins, sorted by accession)

## Recommended PFAM Families

Small, well-characterized families that work well for testing:

| PFAM ID | Name | Reviewed Proteins | Notes |
|---------|------|-------------------|-------|
| PF00096 | Zinc finger, C2H2 type | ~1000+ | Default, very common domain |
| PF00004 | ATP synthase | ~500+ | Well-characterized |
| PF00005 | ABC transporter | ~400+ | Membrane proteins |
| PF00106 | Short chain dehydrogenase | ~300+ | Small domain |
| PF00076 | RNA recognition motif | ~200+ | RNA-binding |

## Usage

### Basic Usage

```bash
# Generate data with default settings (PF00096, 50 proteins)
python test/scripts/generate_test_data.py
```

### Custom PFAM Family

```bash
# Use a different PFAM family with more proteins
python test/scripts/generate_test_data.py --pfam_id PF00004 --max_proteins 100
```

### Skip Embedding Generation

```bash
# Skip embedding generation (faster, for script testing only)
python test/scripts/generate_test_data.py --skip_embeddings
```

### All Options

```bash
python test/scripts/generate_test_data.py --help
```

## Options

- `--pfam_id`: PFAM family ID (default: PF00096)
- `--max_proteins`: Maximum number of proteins to fetch (default: 50)
- `--reviewed_only`: Only fetch reviewed (Swiss-Prot) proteins (default: True)
- `--model`: ESM model name (default: esm2_t6_8M_UR50D)
- `--out_dir`: Output directory (default: data)
- `--skip_embeddings`: Skip embedding generation

## Output Files

The script creates:

1. **data/embeddings/embeddings.tsv**: TSV file with UniProt IDs and embeddings
   - Format: `uniprot_id\tembedding_vector`
   - Embedding vector: comma-separated floats

2. **data/indexes/id_to_row.json**: Mapping of UniProt ID to row index
   - JSON format: `{"uniprot_id": row_index}`

3. **data/indexes/binary_embeddings.npy**: Binary embedding matrix
   - NumPy array: `(n_proteins, embedding_dim)`
   - Binary encoding: value > 0 -> 1, else 0

## Performance

- **Sequence fetching**: ~0.1s per protein (UniProt API rate limiting)
- **Embedding generation**: ~5-10s per protein (depends on sequence length and hardware)
- **Total time for 50 proteins**: ~5-10 minutes

## Data Size

For 50 proteins with ESM-2 8M model (320 dimensions):
- **embeddings.tsv**: ~500 KB
- **id_to_row.json**: ~5 KB
- **binary_embeddings.npy**: ~16 KB
- **Total**: ~500 KB (fits easily in GitHub)

## Requirements

- Python 3.8+
- numpy
- requests
- torch
- fair-esm (for ESM-2 models)
- Internet connection (for UniProt API and model download)

## Integration with Docker

The Dockerfile automatically generates test data if it doesn't exist:

```dockerfile
RUN if [ ! -f /kb/module/data/embeddings/embeddings.tsv ]; then \
    python test/scripts/generate_test_data.py --pfam_id PF00096 --max_proteins 50; \
    fi
```

## Troubleshooting

### Model Download Fails

The ESM-2 model is downloaded on first use (~50MB). If download fails:
- Check internet connection
- Verify PyTorch is installed correctly
- Check disk space

### UniProt API Rate Limiting

The script includes small delays between requests. If you hit rate limits:
- Reduce `--max_proteins`
- Add longer delays in the script
- Use UniProt's batch API (future enhancement)

### Memory Issues

If you run out of memory:
- Reduce `--max_proteins`
- Use a smaller ESM model (esm2_t6_8M_UR50D is the smallest)
- Process proteins in batches (future enhancement)

## Future Enhancements

- Batch processing for large protein sets
- Caching of embeddings to avoid regeneration
- Support for other protein family databases (InterPro, SMART)
- Parallel embedding generation
- Progress bars for long-running operations

