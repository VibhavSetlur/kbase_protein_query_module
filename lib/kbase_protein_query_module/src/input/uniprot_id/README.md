# UniProt ID Input Handler

This module processes UniProt ID input by fetching protein data from the UniProt API.

## Features

*   **Validation**: Validates UniProt ID format using regex.
*   **Fetching**: Retrieves protein sequence and metadata from the UniProt REST API.
*   **Retries**: Implements retry logic for API requests to handle transient failures.
*   **Standardization**: Converts fetched data into the standardized dictionary format.

## Usage

This module is typically used by the `InputManager`.

```python
from kbase_protein_query_module.src.input.uniprot_id import UniProtIdProcessor

processor = UniProtIdProcessor()
result = processor.process({'uniprot_id': ['P12345', 'Q67890']})
```
