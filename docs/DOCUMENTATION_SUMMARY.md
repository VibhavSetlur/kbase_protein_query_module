# KBase Protein Query Analysis Module - Documentation Summary

## Overview

The KBase Protein Query Analysis Module is a comprehensive protein analysis system that provides advanced capabilities for protein network analysis, similarity search, and bioinformatics integration. This module enables researchers to perform sophisticated protein analysis workflows using UniProt IDs as the canonical identifier system.

## Comprehensive Analysis Capabilities

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

## Analysis Workflow

### 1. Protein Existence Check
**Purpose**: Verify protein exists in database and retrieve metadata
- **Input**: UniProt ID (e.g., P00001, P12345)
- **Output**: Existence status, family assignment, metadata, optional embedding
- **Features**: Fast index-based search with comprehensive metadata retrieval

### 2. Protein Embedding Generation
**Purpose**: Generate high-dimensional protein representations
- **Input**: Protein sequence or existing protein check results
- **Output**: High-dimensional protein embeddings using ESM-2 model
- **Features**: Two workflow options (direct sequence or from protein check)

### 3. Family Assignment
**Purpose**: Classify proteins into functional families
- **Input**: Protein embedding
- **Output**: Family classification with confidence scores
- **Features**: Binary centroid similarity for fast family assignment

### 4. Similarity Search
**Purpose**: Find similar proteins within families
- **Input**: Protein embedding and family context
- **Output**: Top similar proteins within family
- **Features**: FAISS-based efficient similarity search

### 5. Network Analysis & Visualization
**Purpose**: Generate comprehensive analysis reports
- **Input**: Analysis results from previous steps
- **Output**: Comprehensive HTML reports with network visualization
- **Features**: Advanced network analysis and protein relationship mapping

## Technical Architecture

### Core Modules
- **`check_existence.py`**: Protein existence verification and metadata retrieval
- **`embedding_generator.py`**: ESM-2 model integration for protein embeddings
- **`assign_protein_family.py`**: Fast family assignment using centroid similarity
- **`similarity_index.py`**: Hierarchical indexing for large-scale search
- **`network_builder.py`**: Protein network construction and visualization
- **`storage.py`**: Efficient storage for massive protein datasets
- **`workflow_orchestrator.py`**: Complete workflow orchestration
- **`sequence_analyzer.py`**: Bioinformatics sequence analysis
- **`html_report_generator.py`**: Comprehensive report generation

### Key Features
- **Scalability**: Handles datasets with 250M+ proteins
- **Efficiency**: FAISS indexing for fast similarity search
- **Modularity**: Clean separation of concerns
- **Integration**: Seamless KBase platform integration
- **Visualization**: Interactive network plots and reports

## Data Structures

### ProteinExistenceResult
```json
{
  "protein_id": "P00001",
  "exists": true,
  "family_id": "ABC_transporter",
  "metadata": {...},
  "embedding_ref": "workspace_object_ref",
  "embedding": [0.1, 0.2, ...],
  "model_name": "esm2_t6_8M_UR50D",
  "search_timestamp": 1234567890.0,
  "summary": "Protein found in database"
}
```

### ProteinEmbeddingResult
```json
{
  "input_id": "P00001",
  "input_source": "workspace_object",
  "embedding_ref": "workspace_object_ref",
  "embedding": [0.1, 0.2, ...],
  "model_name": "esm2_t6_8M_UR50D",
  "pooling_method": "mean",
  "metadata": {...},
  "sequence_length": 500,
  "embedding_norm": 1.0,
  "embedding_dim": 320,
  "protein_id": "P00001",
  "family_id": "ABC_transporter"
}
```

### ProteinFamilyAssignmentResult
```json
{
  "input_id": "P00001",
  "input_type": "embedding",
  "embedding_ref": "workspace_object_ref",
  "assigned_family_id": "ABC_transporter",
  "similarity_score": 0.85,
  "metadata": {...},
  "eigenprotein_id": "P00002",
  "confidence": 0.92
}
```

## Configuration Options

### Embedding Models
- **ESM-2 Variants**: Different model sizes for various use cases
- **Pooling Methods**: Mean, max, attention-based pooling
- **Normalization**: L2 normalization for similarity search

### Similarity Metrics
- **Cosine Similarity**: Standard for normalized embeddings
- **Euclidean Distance**: For raw embedding comparison
- **Hamming Distance**: For binary centroid comparison

### FAISS Index Types
- **IVF (Inverted File)**: For large-scale similarity search
- **Flat**: For exact nearest neighbor search
- **Product Quantization**: For memory-efficient search

### Metadata Storage
- **Parquet**: Efficient columnar storage for metadata
- **HDF5**: Hierarchical storage for embeddings
- **JSON**: Human-readable configuration and results

## Performance Characteristics

### Scalability
- **Dataset Size**: Supports 250M+ proteins
- **Memory Usage**: Efficient chunking and streaming
- **Search Speed**: Sub-second similarity search
- **Storage**: Compressed storage with indexing

### Accuracy
- **Embedding Quality**: State-of-the-art ESM-2 embeddings
- **Family Assignment**: High-confidence family classification
- **Similarity Search**: Precise similarity ranking
- **Metadata**: Comprehensive protein information

## Integration Points

### KBase Platform
- **Workspace Integration**: Seamless object management
- **App Catalog**: Standard KBase app interface
- **SDK Compliance**: Full KBase SDK compatibility
- **Deployment**: Docker container deployment

### External Databases
- **UniProt**: Primary protein database integration
- **PDB**: Structural information integration
- **InterPro**: Functional annotation integration
- **Custom Databases**: Extensible database support

## Testing Framework

### Test Categories
- **Unit Tests**: Individual module functionality
- **Integration Tests**: End-to-end workflow testing
- **Performance Tests**: Scalability and speed testing
- **UI Tests**: User interface functionality

### Test Data
- **Real Protein Data**: Actual UniProt proteins
- **Synthetic Data**: Generated test sequences
- **Edge Cases**: Boundary condition testing
- **Error Conditions**: Exception handling testing

## Deployment

### KBase Registration
1. **Module Registration**: Register on KBase AppDev
2. **Testing**: Comprehensive testing in development environment
3. **Beta Release**: User testing in beta environment
4. **Production**: Full production release

### Configuration
- **Environment Variables**: Runtime configuration
- **Docker Configuration**: Container setup
- **Resource Limits**: Memory and CPU allocation
- **Network Access**: Database connectivity

## Support and Documentation

### Resources
- **API Documentation**: Complete function documentation
- **User Guide**: Step-by-step usage instructions
- **Developer Guide**: Technical implementation details
- **Troubleshooting**: Common issues and solutions

### Community
- **KBase Forum**: Community support and discussion
- **GitHub Issues**: Bug reports and feature requests
- **Documentation**: Comprehensive online documentation
- **Training**: User training and workshops

## Future Development

### Planned Features
- **Additional Models**: Support for other protein language models
- **Enhanced Visualization**: Advanced network visualization tools
- **Machine Learning**: Integration with ML workflows
- **API Expansion**: Extended API capabilities

### Research Integration
- **Academic Collaboration**: Research partnership opportunities
- **Publication Support**: Citation and attribution tools
- **Data Sharing**: Standardized data export formats
- **Reproducibility**: Complete workflow reproducibility

## Citation

When using this module in research, please cite:

```
KBase Protein Query Analysis Module
Vibhav Setlur
https://github.com/kbaseapps/kbase_protein_query_module
```

## Contact

For questions, support, and collaboration:
- **Documentation**: [KBase Documentation](https://docs.kbase.us/)
- **Community**: [KBase Community Forum](https://community.kbase.us/)
- **Support**: [KBase Support](https://kbase.us/contact-us/)
- **Development**: [GitHub Repository](https://github.com/kbaseapps/kbase_protein_query_module) 