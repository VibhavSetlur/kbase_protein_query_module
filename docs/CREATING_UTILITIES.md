# Creating Utilities

This guide explains how to create utility modules that can be used across analyses and other parts of the system.

## Overview

Utilities are reusable modules that provide common functionality. They live in `lib/kbase_protein_query_module/src/util/`.

## Directory Structure

Utilities are organized by category:

```
src/util/
├── embeddings/      # Embedding generation utilities
├── storage/         # Data storage and retrieval
├── uniprot/         # UniProt API interactions
└── your_category/   # Your new utility category
```

## Step 1: Create Utility Directory

Create a directory for your utility:

```bash
mkdir -p lib/kbase_protein_query_module/src/util/your_utility_name
```

## Step 2: Create Utility Module

Create the main utility file:

```python
# lib/kbase_protein_query_module/src/util/your_utility_name/your_utility_name.py

import logging
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)

class YourUtilityName:
    """
    Description of what this utility does.
    
    This utility provides functionality for X, Y, and Z.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the utility.
        
        Args:
            config: Configuration dictionary (optional)
        """
        self.config = config or {}
        self._initialize()
    
    def _initialize(self):
        """Initialize any required resources."""
        # Setup code here
        pass
    
    def do_something(self, input_data: str) -> Dict[str, Any]:
        """
        Main utility method.
        
        Args:
            input_data: Input data to process
        
        Returns:
            Dictionary with results
        """
        # Your utility logic here
        return {'result': 'value'}
```

## Step 3: Create __init__.py

Create an `__init__.py` file:

```python
# lib/kbase_protein_query_module/src/util/your_utility_name/__init__.py

from .your_utility_name import YourUtilityName

__all__ = ['YourUtilityName']
```

## Step 4: Using Your Utility

Import and use your utility in analyses or other code:

```python
# In an analysis or other module
from ...util.your_utility_name import YourUtilityName

# Or if running as script:
from util.your_utility_name import YourUtilityName

# Use it
utility = YourUtilityName(config={})
result = utility.do_something("input")
```

## Example: Simple Utility

Here's a complete example of a simple utility:

```python
# lib/kbase_protein_query_module/src/util/sequence/sequence_utils.py

import logging
from typing import Dict, Any

logger = logging.getLogger(__name__)

class SequenceUtils:
    """Utilities for protein sequence operations."""
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    
    def validate_sequence(self, sequence: str) -> bool:
        """
        Validate that a sequence contains only valid amino acids.
        
        Args:
            sequence: Protein sequence string
        
        Returns:
            True if valid, False otherwise
        """
        if not sequence:
            return False
        return all(aa in self.valid_amino_acids for aa in sequence.upper())
    
    def calculate_molecular_weight(self, sequence: str) -> float:
        """
        Calculate approximate molecular weight.
        
        Args:
            sequence: Protein sequence
        
        Returns:
            Molecular weight in Daltons
        """
        # Simple calculation (not exact, for example only)
        aa_weights = {
            'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13,
            'F': 165.19, 'G': 75.07, 'H': 155.16, 'I': 131.17,
            'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
            'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
            'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
        }
        
        weight = 0.0
        for aa in sequence.upper():
            weight += aa_weights.get(aa, 0.0)
        
        # Subtract weight of water molecules (n-1 for n amino acids)
        weight -= (len(sequence) - 1) * 18.015
        
        return weight
```

## Utility Categories

### Embeddings

Utilities for generating protein embeddings:
- `ProteinEmbeddingGenerator` - Generates embeddings using ESM models

### Storage

Utilities for data storage and retrieval:
- `ProteinStorage` - Loads and manages protein embeddings
- `SimilaritySearch` - Performs similarity searches

### UniProt

Utilities for UniProt API interactions:
- `fetch_metadata` - Fetch protein metadata
- `fetch_protein_sequence` - Fetch protein sequences

## Best Practices

1. **Keep utilities focused**: Each utility should do one thing well
2. **Use logging**: Log important operations and errors
3. **Handle errors gracefully**: Return None or raise clear exceptions
4. **Document well**: Include docstrings for all public methods
5. **Make it testable**: Write tests for your utilities

## Testing Utilities

Create a test file for your utility:

```python
# test/unit_tests/util/test_your_utility.py

import unittest
from kbase_protein_query_module.src.util.your_utility_name import YourUtilityName

class TestYourUtility(unittest.TestCase):
    def setUp(self):
        self.utility = YourUtilityName()
    
    def test_basic_functionality(self):
        result = self.utility.do_something("test")
        self.assertIsNotNone(result)
```

## Import Paths

When importing utilities, use relative imports within the module:

```python
# From analysis to util
from ...util.your_utility_name import YourUtilityName

# From util to util
from ..other_utility import OtherUtility

# When running as script
from util.your_utility_name import YourUtilityName
```

## Common Patterns

### Configuration

Utilities can accept configuration:

```python
class YourUtility:
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.param = self.config.get('param', 'default')
```

### Lazy Initialization

Initialize expensive resources only when needed:

```python
class YourUtility:
    def __init__(self):
        self._resource = None
    
    def _ensure_resource(self):
        if self._resource is None:
            self._resource = self._load_resource()
    
    def use_resource(self):
        self._ensure_resource()
        return self._resource.do_work()
```

### Error Handling

Handle errors gracefully:

```python
def do_work(self, input_data):
    try:
        # Work here
        return result
    except SpecificError as e:
        logger.error(f"Specific error: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        raise
```

## See Also

- `docs/CREATING_ANALYSES.md` - How to create analyses that use utilities
- `docs/CONFIG_GUIDE.md` - Configuration system
- Existing utilities in `src/util/` for examples

