# KBase Protein Query Module - Developer Guide

## Overview

This guide helps developers understand, extend, and contribute to the KBase Protein Query Module. The module is designed with extensibility in mind, allowing easy addition of new analysis types, indexing strategies, and output formats.

## Architecture Overview

### Core Components

```
lib/kbase_protein_query_module/
├── src/
│   ├── core/                    # Core framework and registries
│   │   ├── analysis_registry.py # Framework for adding new analyses
│   │   ├── resource_manager.py  # Server-aware resource management
│   │   ├── parallel_processor.py # Parallel processing framework
│   │   └── pipeline_config.py   # Configuration management
│   ├── stages/                  # Pipeline stages (modular)
│   │   ├── input/              # Input processing stages
│   │   ├── processing/         # Analysis processing stages
│   │   ├── analysis/           # Analysis implementation stages
│   │   └── output/             # Output generation stages
│   ├── storage/                # Data storage and indexing
│   │   ├── indexing_strategy.py # Extensible indexing framework
│   │   └── protein_storage.py  # Main storage implementation
│   ├── workflows/              # Workflow orchestration
│   └── reports/                # Report generation
```

## Adding New Analysis Types

### Step 1: Create Your Analysis Class

Create a new file in `lib/kbase_protein_query_module/src/stages/analysis/`:

```python
# lib/kbase_protein_query_module/src/stages/analysis/my_new_analysis.py

from typing import Dict, Any, List
from ...core.analysis_registry import BaseAnalysis, AnalysisMetadata, register_analysis

@register_analysis("my_new_analysis")
class MyNewAnalysis(BaseAnalysis):
    """
    Custom analysis implementation.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/stages/analysis/my_new_analysis.py:8
    EXTENDS: BaseAnalysis
    USED BY: WorkflowOrchestrator when "my_new_analysis" is requested
    """
    
    def get_metadata(self) -> AnalysisMetadata:
        return AnalysisMetadata(
            name="My New Analysis",
            description="Description of what this analysis does",
            version="1.0.0",
            author="Your Name",
            output_files=["my_analysis.html", "my_results.csv"],
            dependencies=["numpy", "pandas"],
            category="custom",
            computational_complexity="medium",
            memory_requirements="low"
        )
    
    def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
        """Implement your analysis logic here."""
        results = {}
        
        for protein in proteins:
            # Your analysis logic
            analysis_result = self._perform_analysis(protein)
            results[protein.get('id', 'unknown')] = analysis_result
        
        return {
            'analysis_type': 'my_new_analysis',
            'results': results,
            'summary': f"Analyzed {len(proteins)} proteins"
        }
    
    def get_output_files(self, output_dir: str) -> List[str]:
        """Return list of files this analysis will create."""
        return [
            os.path.join(output_dir, "my_analysis.html"),
            os.path.join(output_dir, "my_results.csv")
        ]
    
    def _perform_analysis(self, protein: Any) -> Dict[str, Any]:
        """Your specific analysis implementation."""
        # Implement your analysis logic here
        return {"result": "example"}
```

### Step 2: Create HTML Template

Create `lib/kbase_protein_query_module/src/reports/html/templates/my_new_analysis.html`:

```html
<!DOCTYPE html>
<html>
<head>
    <title>My New Analysis Results</title>
    <style>
        /* Your custom styles */
    </style>
</head>
<body>
    <h1>My New Analysis Results</h1>
    <!-- Your HTML template using {{variable}} placeholders -->
</body>
</html>
```

### Step 3: Update Pipeline Configuration

Your analysis will be automatically available once registered. Users can include it by adding `"my_new_analysis"` to the `analysis_stages` parameter.

## Adding New Indexing Strategies

### Step 1: Implement IndexingStrategy

```python
# lib/kbase_protein_query_module/src/storage/my_custom_index.py

from typing import Dict, Any, List, Tuple
import numpy as np
from .indexing_strategy import IndexingStrategy, IndexingConfig, register_indexing_strategy

@register_indexing_strategy("my_custom_index")
class MyCustomIndexingStrategy(IndexingStrategy):
    """
    Custom indexing strategy implementation.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/my_custom_index.py:8
    EXTENDS: IndexingStrategy
    USED BY: ProteinStorage when strategy="my_custom_index"
    
    CUSTOMIZATION POINTS:
    - Override _build_custom_index() for your indexing algorithm
    - Override _search_custom_index() for your search implementation
    - Override _optimize_for_dataset_size() for size-specific optimizations
    """
    
    def build_index(self, embeddings: np.ndarray, metadata: Dict[str, Any] = None) -> Any:
        """Build your custom index."""
        start_time = time.time()
        
        # Preprocess embeddings using base class method
        processed_embeddings = self._preprocess_embeddings(embeddings)
        
        # Your custom index building logic
        custom_index = self._build_custom_index(processed_embeddings)
        
        # Store results
        self.index = custom_index
        self.metadata = metadata or {}
        self.build_time = time.time() - start_time
        self.is_built = True
        
        return custom_index
    
    def search(self, query_embedding: np.ndarray, k: int = 10, **kwargs) -> List[Tuple[int, float]]:
        """Search using your custom index."""
        # Your search implementation
        return self._search_custom_index(query_embedding, k)
    
    def _build_custom_index(self, embeddings: np.ndarray) -> Any:
        """Implement your custom indexing algorithm here."""
        # Your indexing logic
        pass
    
    def _search_custom_index(self, query: np.ndarray, k: int) -> List[Tuple[int, float]]:
        """Implement your custom search algorithm here."""
        # Your search logic
        pass
```

### Step 2: Integration Points

To use your custom indexing strategy, update the configuration:

```python
# In your pipeline configuration
storage_config = {
    'indexing_strategy': 'my_custom_index',
    'custom_params': {
        'your_param': 'your_value'
    }
}
```

## Class Tracking and Dependencies

### Key Classes and Their Locations

#### Core Framework Classes:
- **`BaseAnalysis`**: `lib/kbase_protein_query_module/src/core/analysis_registry.py:45`
  - **Purpose**: Base class for all analysis implementations
  - **Used by**: All analysis classes in `src/stages/analysis/`
  - **Extends**: `ABC` (Abstract Base Class)

- **`AnalysisRegistry`**: `lib/kbase_protein_query_module/src/core/analysis_registry.py:120`
  - **Purpose**: Registry for managing available analysis types
  - **Used by**: `WorkflowOrchestrator`, `ReportGenerationStage`
  - **Pattern**: Singleton (use `get_registry()`)

- **`ResourceManager`**: `lib/kbase_protein_query_module/src/core/resource_manager.py:65`
  - **Purpose**: Server-aware resource management
  - **Used by**: `ParallelProcessor`, all pipeline stages
  - **Features**: Percentage-based limits for DOE servers

- **`IndexingStrategy`**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:65`
  - **Purpose**: Base class for indexing implementations
  - **Used by**: `ProteinStorage`, `HierarchicalIndex`
  - **Extends**: `ABC`

#### Pipeline Stages:
- **`BaseStage`**: `lib/kbase_protein_query_module/src/stages/base_stage.py:15`
  - **Purpose**: Base class for all pipeline stages
  - **Used by**: All stage implementations
  - **Pattern**: Template method pattern

- **`WorkflowOrchestrator`**: `lib/kbase_protein_query_module/src/workflows/workflow_orchestrator.py:54`
  - **Purpose**: Main workflow coordination
  - **Used by**: Main implementation class
  - **Features**: Stage dependency management, resource coordination

#### Storage and Data:
- **`ProteinStorage`**: `lib/kbase_protein_query_module/src/storage/protein_storage.py:24`
  - **Purpose**: Main storage system for protein data
  - **Features**: Hierarchical storage, memory-efficient loading
  - **Extensible**: Via indexing strategies

### Dependency Tracking

To track where classes are used and how to modify them:

1. **Search for class usage**:
   ```bash
   grep -r "ClassName" lib/kbase_protein_query_module/src/
   ```

2. **Find import statements**:
   ```bash
   grep -r "from.*import.*ClassName" lib/kbase_protein_query_module/
   ```

3. **Check inheritance**:
   ```bash
   grep -r "class.*ClassName" lib/kbase_protein_query_module/
   ```

## Server-Aware Resource Management

### Resource Limits for DOE Servers

The module automatically detects server environments and applies conservative resource limits:

```python
# Default server-aware limits
ResourceLimits(
    max_memory_percent=60.0,    # Use max 60% of server memory
    max_cpu_percent=70.0,       # Use max 70% of server CPU
    max_disk_percent=80.0,      # Use max 80% of server disk
    batch_size_proteins=500,    # Conservative batch sizes
    max_concurrent_tasks=2,     # Limited concurrency for shared servers
    server_safety_margin=0.2    # 20% safety margin
)
```

### Customizing Resource Management

To customize for your environment:

```python
# In your configuration
custom_limits = ResourceLimits(
    max_memory_percent=40.0,  # Even more conservative
    max_cpu_percent=50.0,     # Lower CPU usage
    auto_detect_server_env=True
)
```

## Output Directory Structure

The module creates organized output directories that reflect all analyses performed:

```
output_directory/
├── index.html                 # Main dashboard with links to all analyses
├── sequence_analysis.html     # Sequence analysis results
├── network_analysis.html      # Network analysis results  
├── family_assignment.html     # Family assignment results
├── similarity_search.html     # Similarity search results
├── multi_protein_analysis.html # Multi-protein specific analysis
├── protein_metadata.csv       # Comprehensive protein metadata
├── top_matches.csv            # Top similarity matches
├── pipeline_log.json          # Complete pipeline execution log
└── performance_report.html    # Performance and resource usage report
```

### Adding New Output Types

To add a new output type to the directory structure:

1. **Update ReportGenerationStage**: Add your output generation in `_generate_stage_html_files()`
2. **Update DataExportStage**: Add your data export in `_generate_csv_files()`
3. **Update pipeline log**: Add your analysis info in `_create_pipeline_log()`

## Testing New Extensions

### Unit Testing New Analyses

Create test files following the pattern:

```python
# test/unit_tests/analysis/test_my_new_analysis.py

import unittest
from lib.kbase_protein_query_module.src.core.analysis_registry import get_registry

class TestMyNewAnalysis(unittest.TestCase):
    def setUp(self):
        self.registry = get_registry()
        self.analysis = self.registry.get_analysis("my_new_analysis")
    
    def test_analysis_execution(self):
        """Test that the analysis runs without errors."""
        proteins = [{"id": "test", "sequence": "MKTEST"}]
        results = self.analysis.analyze(proteins)
        self.assertIn('results', results)
```

### Integration Testing

Add your analysis to the integration tests in `test/integration_tests/`.

## Performance Considerations

### Memory Optimization

- Use `np.float32` instead of `np.float64` for embeddings
- Implement streaming for large datasets
- Use batch processing with appropriate sizes
- Monitor memory usage with ResourceManager

### CPU Optimization

- Implement parallel processing where possible
- Use efficient algorithms (O(n log n) or better)
- Leverage NumPy vectorization
- Monitor CPU usage to stay within server limits

### Disk Optimization

- Use compression for large files
- Implement memory mapping for large indexes
- Clean up temporary files
- Monitor disk usage

## Best Practices

### Code Organization

1. **One class per file** for clarity
2. **Clear naming conventions** (AnalysisName + "Analysis")
3. **Comprehensive docstrings** with class locations
4. **Type hints** for all public methods
5. **Logging** for debugging and monitoring

### Error Handling

1. **Graceful degradation** when optional components fail
2. **Clear error messages** with actionable information
3. **Resource cleanup** in finally blocks
4. **Fallback strategies** for critical operations

### Documentation

1. **Class location comments** in docstrings
2. **Usage examples** for complex classes
3. **Extension points** clearly marked
4. **Dependencies** explicitly listed

## Common Extension Patterns

### Adding a New Analysis Type

1. Inherit from `BaseAnalysis`
2. Use `@register_analysis("name")` decorator
3. Implement required abstract methods
4. Create HTML template for results
5. Add unit tests

### Adding a New Indexing Strategy

1. Inherit from `IndexingStrategy`
2. Use `@register_indexing_strategy("name")` decorator
3. Implement required abstract methods
4. Handle serialization/deserialization
5. Add performance benchmarks

### Adding a New Output Format

1. Create new output stage inheriting from `BaseStage`
2. Implement in `src/stages/output/`
3. Update `WorkflowOrchestrator` to include your stage
4. Add configuration options
5. Test with various input sizes

## Debugging and Troubleshooting

### Class Tracking

To find where a class is defined and used:

```bash
# Find class definition
grep -r "class ClassName" lib/kbase_protein_query_module/

# Find class usage
grep -r "ClassName" lib/kbase_protein_query_module/ | grep -v "__pycache__"

# Find imports
grep -r "from.*import.*ClassName" lib/kbase_protein_query_module/
```

### Resource Monitoring

Monitor resource usage during development:

```python
from lib.kbase_protein_query_module.src.core.resource_manager import get_resource_manager

rm = get_resource_manager()
with rm.resource_context("my_operation"):
    # Your code here
    pass

# Get performance summary
summary = rm.get_performance_summary()
print(summary)
```

### Performance Profiling

Enable detailed profiling:

```python
from lib.kbase_protein_query_module.src.core.performance_monitor import get_performance_profiler

profiler = get_performance_profiler()
with profiler.profile_operation("my_analysis"):
    # Your analysis code
    pass

# Generate performance report
profiler.generate_performance_report("performance_output/")
```

## Configuration Management

### Adding New Configuration Options

Update `PipelineConfig` in `src/core/pipeline_config.py`:

```python
@dataclass
class PipelineConfig:
    # ... existing options ...
    
    # Your new configuration options
    my_analysis_enabled: bool = True
    my_analysis_threshold: float = 0.8
    my_custom_params: Dict[str, Any] = field(default_factory=dict)
```

### Environment Variables

The module respects these environment variables:

- `HTML_REPORTS_DIR`: Output directory for HTML reports
- `EXPORTS_DIR`: Output directory for data exports
- `SCRATCH_DIR`: Scratch space for temporary files
- `KB_AUTH_TOKEN`: KBase authentication token
- `SDK_CALLBACK_URL`: KBase SDK callback URL

## KBase Integration

### Workspace Operations

Use KBUtilLib for all workspace operations:

```python
# Get workspace client
workspace_client = self.kb_util.get_workspace_client()

# Save objects
object_ref = self.kb_util.save_workspace_object(
    workspace_name, object_name, object_type, data
)

# Create reports
report_info = self.kb_util.create_report(workspace_name, report_params)
```

### Report Generation

Follow KBase report standards:

```python
report_params = {
    'message': 'Analysis completed successfully',
    'objects_created': [{'ref': object_ref, 'description': 'Analysis result'}],
    'workspace_name': workspace_name,
    'report_object_name': f"{analysis_name}_report",
    'html_links': [{'path': html_file, 'name': 'Analysis Report'}],
    'file_links': [{'path': csv_file, 'name': 'Data Export'}]
}
```

## Testing Guidelines

### Unit Tests

- Test each analysis type independently
- Mock external dependencies (workspace, file system)
- Test edge cases and error conditions
- Verify resource usage stays within limits

### Integration Tests

- Test complete pipeline workflows
- Test with real data from `data/` directory
- Verify output directory structure
- Test resource management under load

### Performance Tests

- Benchmark with different dataset sizes
- Monitor memory usage patterns
- Test parallel processing efficiency
- Verify server resource limits are respected

## Contributing Guidelines

### Code Style

- Follow PEP 8 Python style guidelines
- Use type hints for all public methods
- Include comprehensive docstrings
- Add class location comments

### Documentation

- Update this guide when adding new extension points
- Include examples for complex features
- Document any breaking changes
- Maintain class tracking information

### Testing

- All new features must include tests
- Maintain test coverage above 40%
- Test on both local and server environments
- Include performance benchmarks for new algorithms

## Troubleshooting Common Issues

### Import Errors

If you get import errors for new modules:
1. Check that `__init__.py` files exist in all directories
2. Verify the module path is correct
3. Check for circular imports
4. Ensure dependencies are installed

### Resource Limit Errors

If you hit resource limits:
1. Check current usage with `ResourceManager.get_current_metrics()`
2. Reduce batch sizes or concurrent tasks
3. Implement streaming for large datasets
4. Add memory cleanup in your code

### Registration Errors

If your analysis/indexing strategy isn't found:
1. Ensure the decorator is applied correctly
2. Check that the module is imported somewhere
3. Verify the registry is initialized
4. Check for naming conflicts

## Contact and Support

For questions about extending this module:
1. Check existing implementations for examples
2. Review the class tracking information in docstrings
3. Use the debugging commands to understand data flow
4. Follow the established patterns for consistency

The module is designed to be extensible - if you need to add something that doesn't fit the current patterns, consider whether it should be a new pattern that others can follow!
