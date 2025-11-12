# Workflow Guide

This guide explains how the workflow system works and how components interact.

## Overview

The workflow system processes protein queries through these stages:

1. **Input Processing** - Parse and validate input
2. **Analysis Execution** - Run selected analyses
3. **Output Generation** - Create results and reports

## Architecture

```
User Input
    ↓
InputManager (processes input)
    ↓
WorkflowOrchestrator (coordinates workflow)
    ↓
AnalysisManager (runs analyses)
    ↓
OutputManager (saves results)
    ↓
KBase Report (returns to user)
```

## Components

### InputManager

Handles different input types:
- Protein sequences
- UniProt IDs
- Workspace objects

Location: `src/input/input_manager.py`

### WorkflowOrchestrator

Coordinates the entire workflow:
- Initializes components
- Runs analyses sequentially
- Handles errors
- Returns results

Location: `src/core/workflow_orchestrator.py`

### AnalysisManager

Manages analysis execution:
- Loads analyses from config
- Runs analyses with input data
- Handles analysis results

Location: `src/analysis/analysis_manager.py`

### OutputManager

Manages output files:
- Creates output directory structure
- Saves analysis results
- Organizes files by analysis

Location: `src/output/output_manager.py`

## Workflow Flow

### 1. Initialization

```python
orchestrator = WorkflowOrchestrator(config=workflow_config)
orchestrator.initialize_components(output_dir, workspace_name)
```

This creates:
- InputManager
- AnalysisManager
- OutputManager

### 2. Input Processing

```python
processed_data = input_manager.process_input(input_data)
```

Returns:
- `proteins`: List of protein dictionaries
- `input_type`: Type of input
- `success`: Whether processing succeeded

### 3. Analysis Execution

For each selected analysis:

```python
analysis_input_data = processed_data.copy()
analysis_input_data['output_dir'] = analysis_output_dir
result = analysis_manager.run_analyses(analysis_name, analysis_input_data)
```

Each analysis receives:
- `proteins`: List of proteins to analyze
- `output_dir`: Where to write output files
- Other fields from processed input

### 4. Output Saving

```python
output_manager.save_analysis_output(analysis_name, result, root_dir)
```

Files are organized as:
```
outputs/
├── analysis/
│   ├── network_analysis/
│   │   ├── results.json
│   │   ├── file1.tsv
│   │   └── file2.html
│   └── your_analysis/
│       └── results.json
└── metadata/
    ├── run_metadata.json
    └── process_info.json
```

## Input Data Structure

Input data flows through the system:

### User Input

```python
{
    'input_type': 'uniprot_id',
    'uniprot_id': ['P12345', 'P67890'],
    'workspace_name': 'my_workspace'
}
```

### After Input Processing

```python
{
    'input_type': 'uniprot_id',
    'proteins': [
        {
            'protein_id': 'P12345',
            'sequence': 'MKFLVNVALV...',
            'source': 'uniprot'
        },
        {
            'protein_id': 'P67890',
            'sequence': 'MKFLVNVALV...',
            'source': 'uniprot'
        }
    ],
    'success': True,
    'processing_time': 0.5
}
```

### To Analysis

```python
{
    'proteins': [...],  # Same as above
    'input_type': 'uniprot_id',
    'output_dir': '/path/to/analysis/output',
    'processing_time': 0.5
}
```

## Analysis Results

Analyses return results in this format:

```python
{
    'success': True,  # or False
    'output_files': ['/path/to/file1.tsv', '/path/to/file2.html'],
    'processing_time': 1.23,
    # Optional additional fields
    'results': {...},
    'statistics': {...}
}
```

## Error Handling

The workflow handles errors at multiple levels:

### Analysis Failure

If an analysis fails:
- Error is logged
- Analysis is marked as failed
- Other analyses continue
- Workflow returns `success: False` if any analysis fails

### Workflow Failure

If the workflow fails:
- Error message is included in result
- Partial results may be saved
- KBase job fails (raises exception)

## Output Directory Structure

```
outputs/
├── analysis/
│   └── {analysis_name}/
│       ├── results.json          # Analysis results
│       ├── summary.txt           # Summary (optional)
│       └── {output_files}         # Analysis output files
└── metadata/
    ├── run_metadata.json          # Run metadata
    └── process_info.json          # Process information
```

## Adding to Workflow

### Adding a New Analysis

1. Create analysis (see `docs/CREATING_ANALYSES.md`)
2. Register in config (see `docs/CONFIG_GUIDE.md`)
3. Workflow automatically discovers and runs it

### Adding a New Input Type

1. Create input processor in `src/input/your_type/`
2. Register in `InputManager.process_input()`
3. Add to `standardize_input()` method

### Modifying Workflow Behavior

Edit `WorkflowOrchestrator.run_workflow()`:
- Change analysis execution order
- Add preprocessing steps
- Modify error handling
- Add custom output processing

## Configuration Flow

Configuration flows through the system:

```
kbase_protein_query_moduleImpl
    ↓ (creates workflow_config)
WorkflowOrchestrator
    ↓ (passes to components)
AnalysisManager
    ↓ (passes to analyses)
YourAnalysis.__init__(config)
```

## Logging

The workflow uses Python logging:

```python
import logging
logger = logging.getLogger(__name__)

logger.debug("Detailed debug info")
logger.info("General information")
logger.warning("Warning message")
logger.error("Error message")
```

## Testing Workflow

Test the workflow:

```python
from kbase_protein_query_module.src.core.workflow_orchestrator import WorkflowOrchestrator

config = {
    'output_dir': '/tmp/test_output',
    'workspace_name': 'test_workspace'
}

orchestrator = WorkflowOrchestrator(config=config)

input_data = {
    'input_type': 'protein_sequence',
    'protein_sequence': 'ACDEFGHIKLMNPQRSTVWY'
}

result = orchestrator.run_workflow(
    input_data=input_data,
    selected_analyses=['network_analysis']
)

assert result['success'] is True
```

## Common Patterns

### Processing Multiple Proteins

Analyses receive a list of proteins and process each:

```python
def run_network_analysis(self, input_data):
    proteins = input_data.get('proteins', [])
    for protein in proteins:
        # Process each protein
        result = self._process_protein(protein)
```

### Saving Output Files

Always save to `output_dir` from `input_data`:

```python
output_dir = input_data.get('output_dir', 'outputs')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'results.tsv')
# Save file
return {'output_files': [output_file]}
```

### Error Handling

Return error information in result:

```python
try:
    # Do work
    return {'success': True, ...}
except Exception as e:
    logger.error(f"Error: {e}", exc_info=True)
    return {
        'success': False,
        'error_message': str(e),
        'output_files': []
    }
```

## See Also

- `docs/CREATING_ANALYSES.md` - How analyses fit into workflow
- `docs/CONFIG_GUIDE.md` - Configuration system
- `docs/CREATING_UTILITIES.md` - Utilities used in workflow

