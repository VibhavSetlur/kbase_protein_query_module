
# API Reference

## Core Classes

### BaseAnalysis
**Location**: `lib/kbase_protein_query_module/src/core/analysis_registry.py:45`

Abstract base class for all protein analyses.

**Required Methods**:
- `analyze(proteins, **kwargs) -> Dict[str, Any]`
- `get_output_files(output_dir) -> List[str]`
- `get_metadata() -> AnalysisMetadata`

**Optional Methods**:
- `validate_input(proteins) -> bool`
- `estimate_resources(proteins) -> Dict[str, Any]`

### IndexingStrategy
**Location**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:65`

Abstract base class for indexing implementations.

**Required Methods**:
- `build_index(embeddings, metadata) -> Any`
- `search(query_embedding, k) -> List[Tuple[int, float]]`
- `get_index_info() -> Dict[str, Any]`
- `save_index(filepath) -> bool`
- `load_index(filepath) -> bool`

**Extension Points**:
- `_preprocess_embeddings(embeddings) -> np.ndarray`
- `_calculate_distance(query, candidates) -> np.ndarray`
- `_apply_quantization(embeddings) -> np.ndarray`

### ResourceManager
**Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:65`

Server-aware resource management for DOE environments.

**Key Methods**:
- `check_resource_availability() -> bool`
- `optimize_batch_size(base_size, data_size) -> int`
- `resource_context(operation_name)` (context manager)

## Configuration Classes

### PipelineConfig
**Location**: `lib/kbase_protein_query_module/src/core/pipeline_config.py:17`

Main configuration for pipeline execution.

### ResourceLimits  
**Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:20`

Server-aware resource limits using percentages.

## Extension Patterns

### Registry Pattern
Used for analyses and indexing strategies:
```python
@register_analysis("my_analysis")
class MyAnalysis(BaseAnalysis): ...

@register_indexing_strategy("my_index")  
class MyIndex(IndexingStrategy): ...
```

### Stage Pattern
Used for pipeline stages:
```python
class MyStage(BaseStage):
    def run(self, input_data): ...
```
