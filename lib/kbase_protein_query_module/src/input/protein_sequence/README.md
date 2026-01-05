# Protein Sequence Input Handler

This module processes raw protein sequence input.

## Features

*   **Format Detection**: Automatically detects and handles:
    *   FASTA format (single or multiple sequences).
    *   Raw sequences separated by newlines, commas, or tabs.
*   **Validation**: Uses `Bio.Seq` to validate protein sequences, ensuring they contain only valid IUPAC protein characters.
*   **Standardization**: Converts all input into a standardized dictionary format used by the rest of the pipeline.

## Usage

This module is typically used by the `InputManager`.

```python
from kbase_protein_query_module.src.input.protein_sequence import ProteinSequenceProcessor

processor = ProteinSequenceProcessor()
result = processor.process({'protein_sequence': 'ACDEFGHIKLMNPQRSTVWY'})
```
