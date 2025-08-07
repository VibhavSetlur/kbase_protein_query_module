# KBase Protein Query Analysis Module

A comprehensive protein query analysis system for KBase that provides advanced capabilities for protein network analysis, similarity search, and bioinformatics integration.

## Overview

This module provides comprehensive protein query analysis capabilities using UniProt IDs as the canonical identifier. It enables researchers to perform advanced protein analysis including existence checking, embedding generation, family assignment, similarity search, and network visualization.

## Comprehensive Analysis Workflow

### 1. Protein Existence Check
- **Input**: UniProt ID (e.g., P00001, P12345)
- **Output**: Existence status, family assignment, metadata, optional embedding
- **Features**: Fast index-based search with comprehensive metadata retrieval

### 2. Protein Embedding Generation
- **Input**: Protein sequence or existing protein check results
- **Output**: High-dimensional protein embeddings using ESM-2 model
- **Features**: Two workflow options (direct sequence or from protein check)

### 3. Family Assignment
- **Input**: Protein embedding
- **Output**: Family classification with confidence scores
- **Features**: Binary centroid similarity for fast family assignment

### 4. Similarity Search
- **Input**: Protein embedding and family context
- **Output**: Top similar proteins within family
- **Features**: FAISS-based efficient similarity search

### 5. Network Analysis & Visualization
- **Input**: Analysis results from previous steps
- **Output**: Comprehensive HTML reports with network visualization
- **Features**: Advanced network analysis and protein relationship mapping

## Advanced Capabilities

### Core Analysis Features
- **UniProt ID Canonical System**: Exact match identifier system for reliable protein tracking
- **ESM-2 Language Model**: State-of-the-art protein embedding generation
- **FAISS Similarity Search**: Efficient vector similarity search and clustering
- **Family Assignment**: Binary centroid similarity for protein family classification
- **Metadata Management**: Comprehensive protein information storage and retrieval

### Bioinformatics Integration
- **Protein Database Integration**: Seamless connection with major protein databases
- **Network Analysis**: Advanced protein relationship mapping and analysis
- **Statistical Analysis**: Comprehensive similarity metrics and statistical evaluation
- **HTML Report Generation**: Rich visualization and reporting capabilities

### Workspace Management
- **Object Persistence**: Proper saving and referencing of workspace objects
- **Downstream Analysis**: Seamless integration with other KBase tools
- **Metadata Tracking**: Complete audit trail of analysis steps

## Installation

```bash
# Clone the repository
git clone https://github.com/kbaseapps/kbase_protein_query_module.git

# Install dependencies
pip install -r requirements.txt

# Setup test data
python scripts/setup_test_data.py
```

## Usage

### Basic Protein Analysis

```python
from lib.kbase_protein_query_module.src import check_existence, embedding_generator

# Check protein existence
checker = check_existence.ProteinExistenceChecker()
result = checker.check_protein_existence("P00001")

# Generate embedding
generator = embedding_generator.ProteinEmbeddingGenerator()
embedding = generator.generate_embedding("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ")
```

### Advanced Network Analysis

```python
from lib.kbase_protein_query_module.src import similarity_index, network_builder

# Perform similarity search
index = similarity_index.HierarchicalIndex()
matches = index.search_family_float("ABC_transporter", query_embedding, top_k=50)

# Build protein network
builder = network_builder.ProteinNetworkBuilder()
network = builder.build_network(matches, metadata)
```

## Configuration

The module supports various configuration options:

- **Embedding Models**: ESM-2 variants for different use cases
- **Similarity Metrics**: Cosine, Euclidean, and custom similarity functions
- **FAISS Index Types**: IVF, Flat, and Product Quantization options
- **Metadata Storage**: Parquet, HDF5, and JSON formats

## Testing

```bash
# Run unit tests
python -m pytest test/unit_tests/

# Run integration tests
python -m pytest test/integration_tests/

# Run with conda environment
conda activate kbase-test
python test/run_tests.py
```

## Documentation

- **API Reference**: Complete documentation of all functions and classes
- **Workflow Guide**: Step-by-step analysis workflows
- **Configuration Guide**: Detailed configuration options
- **Troubleshooting**: Common issues and solutions

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this module in your research, please cite:

```
KBase Protein Query Analysis Module
Vibhav Setlur
https://github.com/kbaseapps/kbase_protein_query_module
```

## Support

For questions and support:
- **Documentation**: [KBase Documentation](https://docs.kbase.us/)
- **Community**: [KBase Community Forum](https://community.kbase.us/)
- **Contact**: [KBase Support](https://kbase.us/contact-us/)
