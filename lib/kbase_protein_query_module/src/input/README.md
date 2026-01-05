# Input Module

This module handles the processing and normalization of various input types for the protein query module.

## Components

*   **`InputManager`**: The central class that routes input to the appropriate processor based on `input_type`.
*   **`protein_sequence/`**: Handler for raw protein sequence input.
*   **`uniprot_id/`**: Handler for UniProt ID input.
*   **`workspace_object/`**: Handler for KBase Workspace Object input (e.g., KBaseGenomes.Genome).

## Usage

The `InputManager` takes a dictionary of input parameters and returns a standardized dictionary containing protein data.

```python
from kbase_protein_query_module.src.input import InputManager

mgr = InputManager(config={})
processed_data = mgr.process_input({
    'input_type': 'protein_sequence',
    'protein_sequence': 'ACDEFGHIKLMNPQRSTVWY'
})
```
