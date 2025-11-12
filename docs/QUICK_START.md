# Quick Start Guide

This guide helps you get started quickly with the kbase_protein_query_module.

## What You Need to Know

The module processes protein queries through a workflow:
1. Input is processed (sequences, UniProt IDs, etc.)
2. Analyses are run on the proteins
3. Results are saved and returned

## Creating Your First Analysis

### Step 1: Create the Files

```bash
# Create directory
mkdir -p lib/kbase_protein_query_module/src/analysis/my_first_analysis

# Create main file
touch lib/kbase_protein_query_module/src/analysis/my_first_analysis/my_first_analysis.py
touch lib/kbase_protein_query_module/src/analysis/my_first_analysis/__init__.py
```

### Step 2: Write the Analysis

Copy this template to `my_first_analysis.py`:

```python
import os
import logging
import time
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)

class MyFirstAnalysis:
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        self.config = config or {}
        self.analysis_name = "my_first_analysis"
    
    def run_network_analysis(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        start_time = time.time()
        
        try:
            proteins = input_data.get('proteins', [])
            output_dir = input_data.get('output_dir', 'outputs')
            os.makedirs(output_dir, exist_ok=True)
            
            logger.info(f"Processing {len(proteins)} proteins")
            
            # Your analysis code here
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
            output_file = os.path.join(output_dir, 'my_first_analysis_results.json')
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

def main():
    """Self-test - required for all analyses."""
    import shutil
    import tempfile
    
    ok = True
    output_dir = tempfile.mkdtemp(prefix='test_my_first_analysis_')
    
    try:
        input_data = {
            'input_type': 'protein_sequence',
            'proteins': [{
                'protein_id': 'TEST',
                'sequence': 'ACDEFGHIKLMNPQRSTVWY',
                'source': 'test'
            }],
            'output_dir': output_dir
        }
        
        analysis = MyFirstAnalysis()
        result = analysis.run_network_analysis(input_data)
        
        if not result.get('success'):
            raise RuntimeError(f"Analysis failed: {result.get('error_message')}")
        
        print("ANALYSIS_OK")
    except Exception as e:
        ok = False
        print(f"ANALYSIS_FAIL: {e}")
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)
    
    return 0 if ok else 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
```

### Step 3: Create __init__.py

```python
# lib/kbase_protein_query_module/src/analysis/my_first_analysis/__init__.py
from .my_first_analysis import MyFirstAnalysis
__all__ = ['MyFirstAnalysis']
```

### Step 4: Register in Config

Edit `lib/kbase_protein_query_module/src/analysis/config.py`:

```python
_ANALYSIS_BASE: Dict[str, Dict[str, Any]] = {
    "network_analysis": { ... },
    # Add this:
    "my_first_analysis": {
        "name": "My First Analysis",
        "description": "A simple example analysis",
        "category": "example",
        "module_path": "kbase_protein_query_module.src.analysis.my_first_analysis.my_first_analysis",
        "class_name": "MyFirstAnalysis"
    }
}
```

### Step 5: Test

```bash
cd lib/kbase_protein_query_module/src/analysis/my_first_analysis
python my_first_analysis.py
```

You should see: `ANALYSIS_OK`

## Next Steps

- Read [CREATING_ANALYSES.md](CREATING_ANALYSES.md) for detailed guide
- Read [CONFIG_GUIDE.md](CONFIG_GUIDE.md) for configuration options
- Read [WORKFLOW_GUIDE.md](WORKFLOW_GUIDE.md) to understand the system

## Common Questions

**Q: What method name should I use?**
A: Use `run_network_analysis`, `run`, or `analyze` - the system finds it automatically.

**Q: Where do I save output files?**
A: Save to `output_dir` from `input_data`. The system handles the rest.

**Q: How do I test my analysis?**
A: Add a `main()` function that prints "ANALYSIS_OK" or "ANALYSIS_FAIL: <error>".

**Q: What if my analysis needs dependencies?**
A: Add them to `requires_deps` in config and to `requirements.txt`.

