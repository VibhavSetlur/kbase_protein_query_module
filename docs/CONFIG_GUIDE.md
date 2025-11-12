# Configuration Guide

This guide explains how the configuration system works and how to add new configurations.

## Overview

The configuration system manages:
- Analysis registration and enablement
- Module settings
- Runtime parameters

## Analysis Configuration

Analyses are registered in `lib/kbase_protein_query_module/src/analysis/config.py`.

### Adding an Analysis to Config

Edit the `_ANALYSIS_BASE` dictionary:

```python
# lib/kbase_protein_query_module/src/analysis/config.py

_ANALYSIS_BASE: Dict[str, Dict[str, Any]] = {
    "network_analysis": {
        "name": "Network Analysis",
        "description": "Protein similarity network analysis",
        "category": "network",
        "module_path": "kbase_protein_query_module.src.analysis.network_analysis.network_analysis",
        "class_name": "NetworkAnalysis"
    },
    # Add your analysis here:
    "your_analysis": {
        "name": "Your Analysis",
        "description": "What your analysis does",
        "category": "your_category",
        "module_path": "kbase_protein_query_module.src.analysis.your_analysis.your_analysis",
        "class_name": "YourAnalysis",
        "requires_deps": ["optional_package"]  # Optional
    }
}
```

### Configuration Fields

- **name**: Display name for the analysis
- **description**: What the analysis does
- **category**: Category (e.g., "network", "sequence", "structure")
- **module_path**: Python import path (use dots, not slashes)
- **class_name**: Name of the class to instantiate
- **requires_deps**: Optional list of Python package names to check

### Module Path Format

The `module_path` must match the Python import path:

```
kbase_protein_query_module.src.analysis.your_analysis.your_analysis
```

This corresponds to:
```
lib/kbase_protein_query_module/src/analysis/your_analysis/your_analysis.py
```

### Dependencies

If your analysis requires specific Python packages, list them:

```python
"requires_deps": ["networkx", "sklearn", "pandas"]
```

The system will check if these packages are available. If not, the analysis will be disabled.

## Runtime Configuration

Configuration is passed to analyses and utilities at runtime.

### Workflow Configuration

Configuration is set in `kbase_protein_query_moduleImpl.py`:

```python
workflow_config = {
    'output_dir': output_dir,
    'workspace_name': workspace_name,
    'selected_analyses': selected_analyses,
    'module_root': module_root,
    'embeddings_file': 'data/embeddings/embeddings.tsv',
    'index_path': 'data/indexes',
    # Add your config here
    'your_config_param': 'value'
}
```

### Accessing Configuration

Analyses and utilities receive config in their `__init__` method:

```python
class YourAnalysis:
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.param = self.config.get('your_config_param', 'default')
```

### Configuration Hierarchy

Configuration is merged in this order (later overrides earlier):

1. Global workflow config
2. Analysis-specific config (from `_ANALYSIS_BASE`)
3. Runtime parameters

## Environment Variables

Some configuration comes from environment variables:

- `KB_MODULE_DIR`: Module root directory (usually `/kb/module`)
- `KB_DATA_DIR`: Reference data directory
- `KB_REFDATA_DIR`: Alternative reference data directory
- `SDK_CALLBACK_URL`: KBase callback URL

## Example: Complete Configuration

Here's how to add a new analysis with configuration:

### 1. Add to config.py

```python
_ANALYSIS_BASE: Dict[str, Dict[str, Any]] = {
    "my_analysis": {
        "name": "My Analysis",
        "description": "Does something useful",
        "category": "custom",
        "module_path": "kbase_protein_query_module.src.analysis.my_analysis.my_analysis",
        "class_name": "MyAnalysis",
        "requires_deps": ["numpy"],
        "config": {
            "threshold": 0.5,
            "max_results": 100
        }
    }
}
```

### 2. Use in Analysis

```python
class MyAnalysis:
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        # Get from analysis-specific config or use default
        self.threshold = self.config.get('threshold', 0.5)
        self.max_results = self.config.get('max_results', 100)
```

## Testing Configuration

Test that your analysis is registered:

```python
from kbase_protein_query_module.src.analysis.config import get_enabled_analyses

analyses = get_enabled_analyses()
assert 'my_analysis' in analyses
assert analyses['my_analysis']['enabled'] is True
```

## Common Issues

### Analysis Not Found

- Check `module_path` matches your file structure
- Verify `class_name` matches your class name
- Make sure the module can be imported

### Dependencies Not Found

- Check package names in `requires_deps` are correct
- Verify packages are installed
- Check import names match package names

### Configuration Not Passed

- Verify config is in workflow_config
- Check analysis receives config in `__init__`
- Use `self.config.get('key', 'default')` for safe access

## See Also

- `docs/CREATING_ANALYSES.md` - How to create analyses
- `docs/WORKFLOW_GUIDE.md` - How the workflow uses configuration

