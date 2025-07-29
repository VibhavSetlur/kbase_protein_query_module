# KBase Protein Query Module

A comprehensive KBase module for large-scale protein network analysis, embedding generation, family assignment, and efficient storage/indexing. This toolkit provides scalable solutions for protein similarity search, network construction, and visualization within the KBase platform.

## Overview

The KBase Protein Query Module enables researchers to:
- **Generate protein embeddings** using state-of-the-art ESM-2 models
- **Check protein existence** in large-scale protein databases
- **Assign proteins to families** using fast centroid-based classification
- **Find similar proteins** using efficient FAISS indexing
- **Build protein networks** with interactive visualization
- **Summarize and visualize** analysis results

## Module Structure

```
kbase_protein_network_analysis_toolkit/
├── kbase.yml                    # Module metadata and configuration
├── kbase_protein_query_module.spec # KIDL interface specification
├── Dockerfile                   # Container configuration
├── Makefile                     # Build system
├── deploy.cfg                   # Deployment configuration
├── requirements.txt             # Python dependencies
├── pyproject.toml              # Project configuration
├── lib/                        # Source code implementation
│   ├── kbase_protein_query_module/
│   │   ├── src/               # Core modules
│   │   ├── kbase_protein_query_moduleImpl.py
│   │   └── kbase_protein_query_moduleServer.py
│   └── installed_clients/      # KBase service clients
├── scripts/                    # Utility and deployment scripts
├── test/                       # Test framework and unit tests
├── ui/narrative/methods/       # KBase app UI specifications
└── data/                       # Reference data and examples
```

## Installation

### Prerequisites

1. **KBase SDK**: Install the [KBase SDK](https://kbase.github.io/kb_sdk_docs/tutorial/2_install.html)
2. **Docker**: Required for containerized execution
3. **Python 3.8+**: For local development and testing

### Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/kbase_protein_network_analysis_toolkit.git
   cd kbase_protein_network_analysis_toolkit
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Build the module**:
   ```bash
   make compile
   ```

4. **Configure testing**:
   - Add your KBase developer token to `test_local/test.cfg`
   - Ensure Docker is running

## Usage

### KBase Apps

The module provides five main KBase apps:

1. **Check Protein Existence** (`CheckProteinExistence`)
   - Input: UniProt protein ID
   - Output: Existence status, family ID, and metadata

2. **Generate Protein Embeddings** (`GenerateProteinEmbeddings`)
   - Input: Protein sequence or UniProt ID
   - Output: High-dimensional protein embedding

3. **Assign Protein Family** (`AssignProteinFamilyWithEmbedding`)
   - Input: Protein embedding
   - Output: Family ID, confidence score, and representative protein

4. **Find Top Matches** (`FindTopMatchesFromEmbedding`)
   - Input: Protein embedding and family ID
   - Output: Ranked list of similar proteins

5. **Summarize and Visualize Results** (`SummarizeAndVisualizeResults`)
   - Input: Search results
   - Output: Interactive visualizations and reports

### API Reference

#### Core Functions

- `check_protein_existence(params)` - Check if protein exists in database
- `generate_protein_embedding(params)` - Generate embedding from sequence
- `assign_family_fast(params)` - Assign embedding to protein family
- `find_top_matches_from_embedding(params)` - Find similar proteins
- `summarize_and_visualize_results(params)` - Create visualizations

#### Data Types

- `ProteinEmbeddingResult` - Protein embedding with metadata
- `ProteinFamilyAssignmentResult` - Family assignment results
- `FindTopMatchesFromEmbeddingResults` - Similarity search results
- `CheckProteinExistenceResults` - Protein existence check results

### Local Development

1. **Run tests**:
   ```bash
   kb-sdk test
   ```

2. **Build Docker image**:
   ```bash
   make build
   ```

3. **Run locally**:
   ```bash
   make test
   ```

## Architecture

### Core Modules

- **`embedding_generator.py`**: ESM-2 model integration for protein embeddings
- **`assign_protein_family.py`**: Fast family assignment using centroid similarity
- **`similarity_index.py`**: Hierarchical indexing for large-scale search
- **`network_builder.py`**: Protein network construction and visualization
- **`storage.py`**: Efficient storage for massive protein datasets
- **`workflow_orchestrator.py`**: Complete workflow orchestration

### Key Features

- **Scalability**: Handles datasets with 250M+ proteins
- **Efficiency**: FAISS indexing for fast similarity search
- **Modularity**: Clean separation of concerns
- **Integration**: Seamless KBase platform integration
- **Visualization**: Interactive network plots and reports

## Testing

### Test Structure

- **Unit Tests**: `test/unit_tests_query/` - Individual module testing
- **Integration Tests**: `test/kbase_protein_query_module_query_server_test.py` - Full workflow testing
- **Test Data**: `test/data/` - Sample data for testing

### Running Tests

```bash
# Run all tests
kb-sdk test

# Run specific test file
python -m pytest test/unit_tests_query/test_embedding_generator.py

# Run with coverage
make test
```

## Deployment

### KBase Registration

1. **Push to GitHub**: Ensure all changes are committed and pushed
2. **Register on AppDev**: Use the [KBase AppDev](https://appdev.kbase.us) environment
3. **Test in Development**: Switch to 'D' (develop) mode in the Apps Panel
4. **Move to Beta**: Migrate to beta for user testing
5. **Request Release**: Submit for production release

### Registration Steps

1. Go to [AppDev Module Registration](https://appdev.kbase.us/#appcatalog/module/kbase_protein_query_module)
2. Enter your GitHub repository URL
3. Click "Register" and wait for Docker build
4. Test your app in the development environment
5. Follow the [publishing guide](https://kbase.github.io/kb_sdk_docs/tutorial/8_publish.html) for release

## Documentation

### Code Documentation

- **Docstrings**: All functions and classes are documented
- **Type Hints**: Comprehensive type annotations
- **Comments**: Clear inline documentation
- **Examples**: Usage examples in docstrings

### User Documentation

- **UI Documentation**: Each app has detailed display.yaml files
- **Parameter Descriptions**: Clear input/output specifications
- **Error Handling**: Comprehensive error messages
- **Troubleshooting**: Common issues and solutions

### Developer Documentation

- **API Reference**: Complete function documentation
- **Architecture**: System design and component interactions
- **Testing**: Test coverage and methodology
- **Contributing**: Development guidelines and standards

## Contributing

### Development Guidelines

1. **Code Style**: Follow PEP 8 and KBase conventions
2. **Documentation**: Update docstrings and README files
3. **Testing**: Add tests for new functionality
4. **Version Control**: Use descriptive commit messages

### Pull Request Process

1. Fork the repository
2. Create a feature branch
3. Make your changes with tests
4. Update documentation
5. Submit a pull request

### Release Process

1. Update `RELEASE_NOTES.md` with changes
2. Increment version in `kbase.yml`
3. Update version in UI method specifications
4. Test thoroughly before release
5. Follow KBase release procedures

## Support

### Resources

- **[KBase SDK Documentation](https://kbase.github.io/kb_sdk_docs/)**
- **[KBase Developer Guidelines](https://kbase.github.io/kb_sdk_docs/references/developer_guidelines.html)**
- **[Troubleshooting Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html)**
- **[FAQ](https://kbase.github.io/kb_sdk_docs/references/questions_and_answers.html)**

### Getting Help

- **KBase Community**: [KBase Forum](https://kbase.us/community/)
- **Developer Support**: Contact KBase team for developer account
- **Bug Reports**: Use GitHub Issues
- **Feature Requests**: Submit via GitHub Issues

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this module in your research, please cite:

```
KBase Protein Query Module (2024). 
https://github.com/yourusername/kbase_protein_network_analysis_toolkit
```

## Acknowledgments

- **KBase Team**: For the SDK and platform infrastructure
- **ESM-2 Authors**: For the protein language models
- **FAISS Team**: For the efficient similarity search library
- **Open Source Community**: For the various dependencies and tools

## Changelog

See [RELEASE_NOTES.md](RELEASE_NOTES.md) for detailed version history and changes.
