# KBase Protein Query Module - Library Documentation

This directory contains the core implementation of the KBase Protein Query Module, including all source code, installed clients, and supporting infrastructure.

## Directory Structure

```
lib/
├── kbase_protein_query_module/           # Main module implementation
│   ├── src/                             # Core functional modules
│   │   ├── assign_protein_family.py     # Protein family assignment
│   │   ├── check_existence.py           # Protein existence checking
│   │   ├── embedding_generator.py       # ESM-2 embedding generation
│   │   ├── network_builder.py           # Network construction & visualization
│   │   ├── similarity_index.py          # FAISS-based similarity search
│   │   ├── storage.py                   # Hierarchical data storage
│   │   └── workflow_orchestrator.py     # Complete workflow orchestration
│   ├── kbase_protein_query_moduleImpl.py # Main KBase implementation
│   ├── kbase_protein_query_moduleServer.py # KBase server interface
│   └── __init__.py                      # Package initialization
├── installed_clients/                    # KBase service clients
│   ├── KBaseReportClient.py             # Report generation client
│   ├── WorkspaceClient.py               # Workspace data client
│   ├── authclient.py                    # Authentication client
│   └── baseclient.py                    # Base client utilities
└── README.md                            # This documentation
```

## Core Modules

### `src/assign_protein_family.py`
**Purpose**: Fast protein family assignment using centroid similarity

**Key Features**:
- Binary FAISS indexing for efficient similarity search
- Centroid-based family classification
- Confidence scoring for assignments
- Support for large-scale protein families

**Main Classes**:
- `AssignProteinFamily`: Main family assignment class

**Usage**:
```python
from kbase_protein_query_module.src.assign_protein_family import AssignProteinFamily

assigner = AssignProteinFamily()
assigner.load_family_centroids('path/to/centroids.npz')
result = assigner.assign_family(embedding)
```

### `src/check_existence.py`
**Purpose**: Check protein existence in large-scale databases

**Key Features**:
- Fast lookup in hierarchical storage
- Metadata retrieval for existing proteins
- Family information extraction
- Compressed metadata storage support

**Main Classes**:
- `ProteinExistenceChecker`: Main existence checking class

**Usage**:
```python
from kbase_protein_query_module.src.check_existence import ProteinExistenceChecker

checker = ProteinExistenceChecker()
result = checker.check_protein_existence("P12345")
```

### `src/embedding_generator.py`
**Purpose**: Generate protein embeddings using ESM-2 models

**Key Features**:
- ESM-2 model integration (esm2_t6_8M_UR50D)
- Multiple pooling methods (mean, cls, max)
- Batch processing capabilities
- GPU acceleration support
- Memory-efficient processing

**Main Classes**:
- `ProteinEmbeddingGenerator`: Main embedding generation class

**Usage**:
```python
from kbase_protein_query_module.src.embedding_generator import ProteinEmbeddingGenerator

generator = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D")
embedding = generator.generate_embedding("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ")
```

### `src/similarity_index.py`
**Purpose**: Efficient similarity search for large protein datasets

**Key Features**:
- Hierarchical indexing structure
- FAISS-based similarity search
- Streaming capabilities for massive datasets
- Quantization for memory efficiency
- Caching for frequently accessed families

**Main Classes**:
- `HierarchicalIndex`: Main indexing class
- `StreamingIndex`: Streaming search capabilities

**Usage**:
```python
from kbase_protein_query_module.src.similarity_index import HierarchicalIndex

index = HierarchicalIndex(base_dir="data/indexes")
results = index.search_family("family_1", query_embedding, top_k=50)
```

### `src/network_builder.py`
**Purpose**: Protein network construction and visualization

**Key Features**:
- Multiple network construction methods (mutual_knn, threshold, hybrid)
- Interactive visualization with Plotly
- Network analysis and statistics
- Dynamic network building
- Export capabilities

**Main Classes**:
- `DynamicNetworkBuilder`: Main network construction class

**Usage**:
```python
from kbase_protein_query_module.src.network_builder import DynamicNetworkBuilder

builder = DynamicNetworkBuilder(k_neighbors=8)
network = builder.build_mutual_knn_network(embeddings, protein_ids)
```

### `src/storage.py`
**Purpose**: Efficient storage for massive protein datasets

**Key Features**:
- Hierarchical storage structure
- Chunked data access
- Compressed metadata storage
- Memory-efficient loading
- Streaming capabilities

**Main Classes**:
- `ProteinStorage`: Main storage class
- `CompressedMetadataStorage`: Metadata storage
- `MemoryEfficientLoader`: Memory-optimized loading

**Usage**:
```python
from kbase_protein_query_module.src.storage import ProteinStorage

storage = ProteinStorage(base_dir="data")
embeddings, protein_ids = storage.load_family_embeddings("family_1")
```

### `src/workflow_orchestrator.py`
**Purpose**: Complete workflow orchestration

**Key Features**:
- End-to-end workflow management
- Performance monitoring
- Memory optimization
- Configuration management
- Legacy compatibility

**Main Classes**:
- `ProteinNetworkWorkflow`: Main workflow class
- `LegacyProteinNetworkWorkflow`: Legacy compatibility

**Usage**:
```python
from kbase_protein_query_module.src.workflow_orchestrator import ProteinNetworkWorkflow

workflow = ProteinNetworkWorkflow(config_file="config.yaml")
results = workflow.run_optimized_workflow(query_sequence, query_protein_id)
```

## KBase Integration

### `kbase_protein_query_moduleImpl.py`
**Purpose**: Main KBase implementation interface

**Key Features**:
- KBase function implementations
- Parameter validation
- Error handling
- Report generation
- Workspace integration

**Main Functions**:
- `check_protein_existence()`: Check protein existence
- `generate_protein_embedding()`: Generate embeddings
- `assign_family_fast()`: Assign protein families
- `find_top_matches_from_embedding()`: Find similar proteins
- `summarize_and_visualize_results()`: Create visualizations

### `kbase_protein_query_moduleServer.py`
**Purpose**: KBase server interface and API

**Key Features**:
- JSON-RPC service implementation
- Authentication handling
- Request/response processing
- Error handling
- Logging and monitoring

## Installed Clients

### `installed_clients/KBaseReportClient.py`
**Purpose**: Generate KBase reports

**Usage**:
```python
from installed_clients.KBaseReportClient import KBaseReport

report = KBaseReport(callback_url)
report_info = report.create_extended_report({
    'message': 'Analysis complete',
    'objects_created': [],
    'workspace_name': workspace_name
})
```

### `installed_clients/WorkspaceClient.py`
**Purpose**: Interact with KBase workspace

**Usage**:
```python
from installed_clients.WorkspaceClient import Workspace

ws = Workspace(workspace_url)
objects = ws.get_objects([{'ref': workspace_name + '/' + object_id}])
```

## Development

### Code Generation
To regenerate templated code after spec file changes:
```bash
make compile
```

### Testing
Run tests for the library components:
```bash
python -m pytest test/unit_tests_query/
```

### Documentation
All modules include comprehensive docstrings following Google style:
```python
def function_name(param1: str, param2: int) -> Dict[str, Any]:
    """
    Brief description of function.
    
    Args:
        param1: Description of parameter 1
        param2: Description of parameter 2
        
    Returns:
        Description of return value
        
    Raises:
        ValueError: Description of when this is raised
    """
```

## Performance Considerations

### Memory Usage
- Large protein datasets may require significant memory
- Use streaming capabilities for massive datasets
- Enable GPU acceleration for embedding generation

### Scalability
- Supports datasets with 250M+ proteins
- Hierarchical indexing for efficient search
- Chunked storage for memory efficiency

### Optimization
- FAISS indexing for fast similarity search
- Compressed storage for metadata
- Batch processing for embedding generation

## Troubleshooting

### Common Issues
1. **Import Errors**: Ensure PYTHONPATH includes lib directory
2. **Memory Issues**: Use streaming capabilities for large datasets
3. **GPU Issues**: Check CUDA installation for GPU acceleration
4. **Dependency Issues**: Verify all requirements are installed

### Debugging
- Enable debug logging in workflow orchestrator
- Use memory monitoring for large datasets
- Check FAISS index training status

## Contributing

See the main [README.md](../README.md) for contribution guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.
