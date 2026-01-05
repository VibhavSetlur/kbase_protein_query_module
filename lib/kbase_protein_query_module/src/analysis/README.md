# Analysis Module

This module manages the execution of various protein analysis workflows.

## Components

*   **`AnalysisManager`**: The main class that orchestrates the loading and execution of analysis modules.
*   **`config.py`**: Configuration for available analyses and dependency checking.

## Subdirectories

*   **`network_analysis/`**: Contains the logic for Protein Similarity Network Analysis.

## Usage

The `AnalysisManager` can be initialized and used to run specific analyses:

```python
from kbase_protein_query_module.src.analysis import AnalysisManager

manager = AnalysisManager()
results = manager.run_analyses("network_analysis", input_data)
```
