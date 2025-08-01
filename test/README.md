# KBase Protein Query Module - Testing Documentation

This directory contains comprehensive testing infrastructure for the KBase Protein Query Module, including unit tests, integration tests, and test data.

## Testing Overview

The testing infrastructure follows KBase SDK standards and includes:

1. **Unit Tests** (`unit_tests_query/`) - Individual component testing
2. **Integration Tests** (`kbase_protein_query_module_query_server_test.py`) - Full module testing
3. **Local Test Runner** (`test_local/`) - Local development testing
4. **Actual Data Testing** - Tests using real data from the `data/` folder

## Test Data

The module uses actual test data generated by `scripts/generate_dummy_data.py`:

- **50 protein families** with 500 proteins each
- **ESM2 model** (`esm2_t6_8M_UR50D`) for embedding generation
- **FAISS indexes** for similarity search
- **Metadata** in Parquet format
- **Centroids** for family assignment

## Running Tests

### 1. KBase SDK Tests (Recommended)

```bash
# Run KBase SDK tests (proper way)
make test-kbase

# Or directly
kb-sdk test
```

### 2. Local Tests with Actual Data

```bash
# Set up test data and run tests
make test-with-data

# Or step by step
make setup-test-data
make test-with-actual-data
```

### 3. Unit Tests Only

```bash
# Run unit tests
cd test/unit_tests_query
python -m unittest discover

# Run specific test
python -m unittest test_embedding_generator.TestProteinEmbeddingGenerator
```

### 4. Integration Tests

```bash
# Run integration tests
cd test
python -m unittest kbase_protein_query_module_query_server_test
```

## Test Structure

### Unit Tests (`unit_tests_query/`)

- `test_assign_protein_family_query.py` - Family assignment testing
- `test_check_existence.py` - Protein existence checking
- `test_embedding_generator.py` - Embedding generation
- `test_network_builder.py` - Network construction
- `test_similarity_index.py` - Similarity search
- `test_storage.py` - Data storage
- `test_workflow_orchestrator.py` - End-to-end workflow

### Integration Tests

- `kbase_protein_query_module_query_server_test.py` - Full module testing with KBase SDK

### Local Test Runner

- `test_local/test_with_data.py` - Tests using actual data
- `test_local/test_runner.py` - General test runner
- `test_local/generate_test_data.py` - Test data generation

## Test Data Format

The test data follows this structure:

```
data/
├── esm2_t6_8M_UR50D_local/     # ESM2 model files
├── families/                     # HDF5 files with embeddings
│   ├── family_0.h5
│   ├── family_1.h5
│   └── ...
├── metadata/                     # Parquet metadata files
│   ├── family_0_metadata.parquet
│   ├── family_1_metadata.parquet
│   └── ...
├── indexes/                      # FAISS indexes
│   ├── families/
│   │   ├── family_0.faiss
│   │   ├── family_1.faiss
│   │   └── ...
│   └── family_mapping.json
└── family_centroids_binary.npz  # Family centroids
```

## Protein ID Formats

The test data includes multiple protein ID formats:

- **UniProt-like**: `P00000001`, `P00000002`, etc.
- **Dummy format**: `family_0_prot_1`, `family_25_prot_100`, etc.
- **Generic**: 6-10 character alphanumeric IDs

## Model Requirements

Tests require the ESM2 model (`esm2_t6_8M_UR50D`) to be available in:
- `data/esm2_t6_8M_UR50D_local/` (relative to project root)
- `/kb/module/data/esm2_t6_8M_UR50D_local/` (in Docker container)

## Test Configuration

### Environment Variables

```bash
export KB_AUTH_TOKEN=dummy_token
export SDK_CALLBACK_URL=http://localhost:5000/callback
export KB_DEPLOYMENT_CONFIG=test_local/test.cfg
```

### Test Configuration File

`test_local/test.cfg` contains KBase service endpoints and authentication settings.

## Troubleshooting

### Model Not Found

If tests fail with "Local model not found":

1. Ensure the ESM2 model is downloaded to `data/esm2_t6_8M_UR50D_local/`
2. Run `make setup-test-data` to generate test data
3. Check model file permissions

### Data Not Found

If tests fail with "Data directory not found":

1. Run `make setup-test-data` to generate test data
2. Ensure the `data/` directory exists with proper structure
3. Check file permissions

### KBase SDK Issues

If `kb-sdk test` fails:

1. Ensure KBase SDK is properly installed
2. Check network connectivity to KBase services
3. Verify authentication token in `test_local/test.cfg`

## Coverage

Tests cover:

- ✅ Protein existence checking
- ✅ Embedding generation
- ✅ Family assignment
- ✅ Similarity search
- ✅ Network construction
- ✅ Data storage and retrieval
- ✅ End-to-end workflows
- ✅ Error handling
- ✅ Parameter validation

## Continuous Integration

The module includes:

- **Travis CI** configuration (`.travis.yml`)
- **Docker** testing support
- **Coverage** reporting
- **Automated** test execution

## Best Practices

1. **Use actual data** when possible for realistic testing
2. **Mock external dependencies** for unit tests
3. **Test error conditions** and edge cases
4. **Follow KBase SDK standards** for integration tests
5. **Maintain test data** consistency across environments
 
