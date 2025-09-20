# Developer Guide

## Architecture

Modular backend with thin KBase facade:

```
KBase Facade (Impl) → WorkflowOrchestrator → Modular Backend
                                           ├── Input Processing
                                           ├── Analysis System
                                           ├── Output Generation
                                           └── Utilities
```

## Directory Structure

```
lib/kbase_protein_query_module/src/
├── core/                    # Core orchestration
│   ├── workflow_orchestrator.py
│   ├── pipeline_config.py
│   ├── parallel_processor.py
│   ├── performance_monitor.py
│   └── resource_manager.py
├── analysis/               # Analysis modules
│   ├── analysis_manager.py
│   ├── config.py
│   └── */                  # Individual analysis modules
├── input/                  # Input processing
│   ├── input_validation.py
│   ├── data_extraction.py
│   └── workspace_object.py
├── outputs/                # Output generation
│   ├── output_manager.py
│   ├── config.py
│   └── analysis/           # Analysis-specific outputs
└── util/                   # Utilities
    ├── storage/            # Data storage
    ├── embeddings/         # Embedding generation
    ├── similarity/         # Similarity search
    └── visualization/      # Visualization tools
```

## Key Components

### WorkflowOrchestrator

Central coordinator managing the analysis pipeline:

```python
class WorkflowOrchestrator:
    def execute(self, input_data: Union[Dict, PipelineConfig]) -> WorkflowResult:
        """Main entry point for workflow execution"""
        
    def run_workflow(self, input_data: Dict, output_dir: str) -> WorkflowResult:
        """Execute the complete analysis pipeline"""
```

### AnalysisManager

Manages analysis modules and execution:

```python
class AnalysisManager:
    def run_analysis(self, analysis_name: str, proteins: List, output_dir: str) -> Dict:
        """Run a single analysis"""
        
    def register_analysis(self, analysis_name: str, analysis_module: Any):
        """Register a new analysis module"""
```

### PipelineConfig

Configuration management:

```python
@dataclass
class PipelineConfig:
    input_type: str
    input_data: Union[str, List[str]]
    selected_analyses: List[str]
    output_dir: str
    max_workers: int
    memory_limit_gb: float
```

## Adding Components

### New Analysis

1. Create module in `src/analysis/your_analysis/`
2. Register in `src/analysis/config.py`
3. Add output handler in `src/outputs/analysis/your_analysis/`

### New Input Type

1. Extend `src/input/input_validation.py`
2. Add extraction logic in `src/input/data_extraction.py`

### New Output Format

1. Extend `src/outputs/output_manager.py`
2. Register in `src/outputs/config.py`

## Configuration

Hierarchical configuration system:

1. Runtime environment variables
2. PipelineConfig parameters
3. Default configuration values

## Error Handling

- Graceful degradation for non-critical failures
- Comprehensive logging for debugging
- Automatic retry for transient failures
- User-friendly error messages

## Testing

```bash
# Run all tests
pytest

# Run specific categories
pytest test/unit_tests/core/
pytest test/unit_tests/analysis/

# Run with coverage
pytest --cov=src --cov-report=html
```

## Best Practices

### Code Organization

- Single responsibility per module
- Dependency injection for testability
- Clear interfaces between components
- Configuration-driven behavior

### Error Handling

- Validate inputs early
- Continue execution when possible
- Provide clear error descriptions
- Log all errors with context

### Testing

- Unit tests for individual components
- Integration tests for component interactions
- Mock external dependencies
- Maintain high test coverage

## Migration

### From Legacy Architecture

1. Map old functions to new modules
2. Update import paths
3. Convert configuration to PipelineConfig
4. Update tests for new architecture

### Backward Compatibility

- Legacy KBase methods still work
- Automatic configuration migration
- Maintained API contracts
- Incremental migration support