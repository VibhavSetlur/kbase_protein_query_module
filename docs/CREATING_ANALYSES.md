# Creating New Analyses

This guide explains how to create a new analysis for the kbase_protein_query_module. Follow these steps to add your own analysis.

## Overview

An analysis is a module that processes protein data and produces results. The system automatically discovers and runs analyses based on configuration.

## Step 1: Create the Analysis Directory

Create a new directory for your analysis:

```bash
mkdir -p lib/kbase_protein_query_module/src/analysis/your_analysis_name
```

## Step 2: Create the Analysis Class

Create the main analysis file:

```python
# lib/kbase_protein_query_module/src/analysis/your_analysis_name/your_analysis_name.py

import os
import logging
from typing import Dict, Any, List, Optional
import time

logger = logging.getLogger(__name__)

class YourAnalysisName:
    """
    Your analysis description here.
    
    This analysis does X, Y, and Z.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the analysis.
        
        Args:
            config: Configuration dictionary (optional)
        """
        self.config = config or {}
        self.analysis_name = "your_analysis_name"
        
        # Set any configuration parameters
        self.param1 = self.config.get('param1', 'default_value')
        self.param2 = self.config.get('param2', 10)
        
        logger.debug(f"Initialized {self.analysis_name}")
    
    def run_network_analysis(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Main method that runs the analysis.
        
        This method is called by the workflow orchestrator. The method name
        can be 'run_network_analysis', 'run', or 'analyze' - the system
        will find it automatically.
        
        Args:
            input_data: Dictionary containing:
                - 'proteins': List of protein dictionaries with 'protein_id', 'sequence', 'source'
                - 'input_type': Type of input ('protein_sequence', 'uniprot_id', etc.)
                - 'output_dir': Directory where output files should be written
        
        Returns:
            Dictionary with:
                - 'success': Boolean indicating if analysis succeeded
                - 'output_files': List of file paths created
                - 'processing_time': Time taken in seconds
                - 'error_message': Error message if failed (optional)
        """
        start_time = time.time()
        
        try:
            # Get input data
            proteins = input_data.get('proteins', [])
            output_dir = input_data.get('output_dir', 'outputs')
            
            # Ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            
            logger.info(f"Running {self.analysis_name} on {len(proteins)} proteins")
            
            # Process each protein
            results = []
            for protein in proteins:
                protein_id = protein.get('protein_id', 'unknown')
                sequence = protein.get('sequence', '')
                
                # Do your analysis here
                result = self._analyze_protein(protein_id, sequence)
                results.append(result)
            
            # Save results to files
            output_files = self._save_results(results, output_dir)
            
            execution_time = time.time() - start_time
            
            return {
                'success': True,
                'output_files': output_files,
                'processing_time': execution_time,
                'results': results
            }
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"{self.analysis_name} failed: {e}", exc_info=True)
            return {
                'success': False,
                'output_files': [],
                'processing_time': execution_time,
                'error_message': str(e)
            }
    
    def _analyze_protein(self, protein_id: str, sequence: str) -> Dict[str, Any]:
        """
        Analyze a single protein.
        
        Args:
            protein_id: UniProt ID or protein identifier
            sequence: Protein sequence
        
        Returns:
            Dictionary with analysis results
        """
        # Your analysis logic here
        # Example:
        length = len(sequence)
        molecular_weight = self._calculate_molecular_weight(sequence)
        
        return {
            'protein_id': protein_id,
            'length': length,
            'molecular_weight': molecular_weight
        }
    
    def _calculate_molecular_weight(self, sequence: str) -> float:
        """Calculate molecular weight (example method)."""
        # Your calculation here
        return 0.0
    
    def _save_results(self, results: List[Dict[str, Any]], output_dir: str) -> List[str]:
        """
        Save analysis results to files.
        
        Args:
            results: List of result dictionaries
            output_dir: Directory to save files
        
        Returns:
            List of file paths created
        """
        output_files = []
        
        # Example: Save as JSON
        import json
        json_path = os.path.join(output_dir, f"{self.analysis_name}_results.json")
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2)
        output_files.append(json_path)
        
        # Example: Save as TSV
        import csv
        tsv_path = os.path.join(output_dir, f"{self.analysis_name}_results.tsv")
        with open(tsv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys() if results else [])
            writer.writeheader()
            writer.writerows(results)
        output_files.append(tsv_path)
        
        logger.info(f"Saved {len(output_files)} output files to {output_dir}")
        
        return output_files
```

## Step 3: Create __init__.py

Create an `__init__.py` file to make it a Python package:

```python
# lib/kbase_protein_query_module/src/analysis/your_analysis_name/__init__.py

from .your_analysis_name import YourAnalysisName

__all__ = ['YourAnalysisName']
```

## Step 4: Register in Config

Add your analysis to the configuration file:

```python
# lib/kbase_protein_query_module/src/analysis/config.py

# In the _ANALYSIS_BASE dictionary, add:

"your_analysis_name": {
    "name": "Your Analysis Display Name",
    "description": "What your analysis does",
    "category": "your_category",  # e.g., "sequence", "network", "structure"
    "module_path": "kbase_protein_query_module.src.analysis.your_analysis_name.your_analysis_name",
    "class_name": "YourAnalysisName",
    "requires_deps": ["optional_dependency"]  # Optional: list of Python packages needed
}
```

**Important**: The `module_path` must match the Python import path. Use dots (.) not slashes.

## Step 5: Add Self-Test (REQUIRED)

**All analyses must include a self-test.** This is how the system verifies your analysis works. The test must:

1. Print "ANALYSIS_OK" on success
2. Print "ANALYSIS_FAIL: <error>" on failure  
3. Return 0 on success, 1 on failure

Add a self-test method to your analysis class:

```python
# Add this method to your analysis class

def main():
    """Self-test for the analysis."""
    import shutil
    import tempfile
    
    ok = True
    output_dir = tempfile.mkdtemp(prefix='test_your_analysis_')
    
    try:
        # Create test input
        input_data = {
            'input_type': 'protein_sequence',
            'proteins': [
                {
                    'protein_id': 'TEST_PROTEIN',
                    'sequence': 'ACDEFGHIKLMNPQRSTVWY',
                    'source': 'test'
                }
            ],
            'output_dir': output_dir
        }
        
        # Run analysis
        analysis = YourAnalysisName(config={})
        result = analysis.run_network_analysis(input_data)
        
        # Check result
        if not isinstance(result, dict) or result.get('success') is not True:
            raise RuntimeError(f"Analysis failed: {result}")
        
        # Check output files
        output_files = result.get('output_files', [])
        if not output_files:
            raise RuntimeError("No output files created")
        
        for file_path in output_files:
            if not os.path.exists(file_path):
                raise RuntimeError(f"Output file not found: {file_path}")
        
        print("ANALYSIS_OK")
        
    except Exception as e:
        ok = False
        print(f"ANALYSIS_FAIL: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Cleanup
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir, ignore_errors=True)
    
    return 0 if ok else 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
```

## Step 6: Test Your Analysis

Run the self-test:

```bash
cd lib/kbase_protein_query_module/src/analysis/your_analysis_name
python your_analysis_name.py
```

You should see `ANALYSIS_OK` if it works.

## Method Name Options

The workflow orchestrator will look for these method names (in order):
1. `run_network_analysis(input_data)` - Used by network_analysis
2. `run(input_data, **kwargs)` - Generic run method
3. `analyze(input_data, **kwargs)` - Alternative analyze method

Use whichever name makes sense for your analysis. The system will find it automatically.

## Input Data Structure

Your analysis will receive `input_data` with this structure:

```python
{
    'proteins': [
        {
            'protein_id': 'P12345',  # UniProt ID or identifier
            'sequence': 'MKFLVNVALV...',  # Protein sequence
            'source': 'uniprot'  # Source: 'uniprot', 'protein_sequence', etc.
        },
        # ... more proteins
    ],
    'input_type': 'uniprot_id',  # or 'protein_sequence'
    'output_dir': '/path/to/output/directory',
    # ... other fields from input processing
}
```

## Output Requirements

Your analysis must return a dictionary with:

- `success`: Boolean (True if successful, False if failed)
- `output_files`: List of file paths created (can be empty)
- `processing_time`: Time taken in seconds (float)
- `error_message`: Error message string (only if success is False)

Example:

```python
{
    'success': True,
    'output_files': ['/path/to/file1.tsv', '/path/to/file2.json'],
    'processing_time': 1.23,
    # Optional: any additional fields
    'results': {...},
    'statistics': {...}
}
```

## Using Utilities

You can use existing utilities in your analysis:

```python
# Import utilities
from ...util.embeddings.generator import ProteinEmbeddingGenerator
from ...util.storage.storage import ProteinStorage
from ...util.uniprot.api import fetch_metadata, fetch_protein_sequence

# Use in your analysis
embedding_gen = ProteinEmbeddingGenerator()
embedding = embedding_gen.generate_embedding(sequence)
```

## Example: Complete Simple Analysis

Here's a complete example of a simple analysis:

```python
import os
import logging
import time
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)

class SequenceLengthAnalysis:
    """Simple analysis that calculates sequence lengths."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        self.config = config or {}
        self.analysis_name = "sequence_length_analysis"
    
    def run_network_analysis(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        start_time = time.time()
        
        try:
            proteins = input_data.get('proteins', [])
            output_dir = input_data.get('output_dir', 'outputs')
            os.makedirs(output_dir, exist_ok=True)
            
            results = []
            for protein in proteins:
                protein_id = protein.get('protein_id', 'unknown')
                sequence = protein.get('sequence', '')
                results.append({
                    'protein_id': protein_id,
                    'length': len(sequence)
                })
            
            # Save results
            import json
            output_file = os.path.join(output_dir, 'sequence_lengths.json')
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            return {
                'success': True,
                'output_files': [output_file],
                'processing_time': time.time() - start_time
            }
        except Exception as e:
            return {
                'success': False,
                'output_files': [],
                'processing_time': time.time() - start_time,
                'error_message': str(e)
            }
```

## Testing

All analyses must include self-tests. The test should:
1. Create test input data
2. Run the analysis
3. Verify output files are created
4. Print "ANALYSIS_OK" on success or "ANALYSIS_FAIL: <error>" on failure
5. Return 0 on success, 1 on failure

## Common Issues

### Analysis Not Found

- Check that `module_path` in config matches your file structure
- Check that `class_name` matches your class name exactly
- Make sure `__init__.py` exists and imports your class

### Method Not Found

- Make sure your class has one of: `run_network_analysis`, `run`, or `analyze`
- Check method signature matches: `def method_name(self, input_data: Dict[str, Any])`

### Output Files Not Found

- Make sure you're writing files to `output_dir` from `input_data`
- Check file paths are absolute or relative to `output_dir`
- Verify files are actually created before returning

## Next Steps

1. Test your analysis with the self-test
2. Add it to the config
3. Test it in the full workflow
4. Add any required dependencies to `requirements.txt`
5. Update documentation

## See Also

- `docs/CREATING_UTILITIES.md` - How to create utility modules
- `docs/CONFIG_GUIDE.md` - Understanding configuration
- `docs/WORKFLOW_GUIDE.md` - Understanding the workflow system

