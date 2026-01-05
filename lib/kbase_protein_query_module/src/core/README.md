# Core Module

This module contains the core logic for the KBase Protein Query Module.

## Components

*   **`WorkflowOrchestrator`**: The central class that manages the entire workflow. It coordinates:
    1.  **Input Processing**: Using `InputManager`.
    2.  **Analysis Execution**: Using `AnalysisManager`.
    3.  **Output Generation**: Using `OutputManager`.

*   **`pipeline_config.yaml`**: Configuration file for the pipeline (e.g., data directories, logging levels).

## Usage

The `WorkflowOrchestrator` is the main entry point for running the module.

```python
from kbase_protein_query_module.src.core import WorkflowOrchestrator

orchestrator = WorkflowOrchestrator(config={...})
result = orchestrator.run_workflow(input_data, output_dir="/path/to/output")
```
