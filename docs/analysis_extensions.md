# Adding New Analyses

## Overview

Add new analysis types by creating a module, registering it, adding output handlers, and writing tests. Keep tests concise, professional, and aligned with `kb-sdk test`.

## Step-by-Step Guide

### 1. Create Analysis Module

Create directory and main file:

```bash
mkdir -p lib/kbase_protein_query_module/src/analysis/your_analysis
```

```python
# lib/kbase_protein_query_module/src/analysis/your_analysis/your_analysis.py
from typing import List, Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)

class YourAnalysis:
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        self.config = config or {}
        self.name = "your_analysis"
        self.version = "1.0.0"
        
    def run(self, proteins: List[Dict[str, Any]], output_dir: str, **kwargs) -> Dict[str, Any]:
        try:
            logger.info(f"Starting {self.name} analysis on {len(proteins)} proteins")
            results = self._perform_analysis(proteins, **kwargs)
            self._save_results(results, output_dir)
            
            return {
                "status": "completed",
                "analysis_name": self.name,
                "proteins_analyzed": len(proteins),
                "results": results,
                "output_files": self._get_output_files(output_dir)
            }
        except Exception as e:
            logger.error(f"Analysis {self.name} failed: {e}")
            return {"status": "failed", "analysis_name": self.name, "error": str(e)}
    
    def _perform_analysis(self, proteins: List[Dict[str, Any]], **kwargs) -> Dict[str, Any]:
        results = {}
        for protein in proteins:
            protein_id = protein.get('id', 'unknown')
            results[protein_id] = {
                "score": 0.95,
                "confidence": 0.87,
                "details": "Analysis details"
            }
        return results
    
    def _save_results(self, results: Dict[str, Any], output_dir: str) -> None:
        import json
        import os
        os.makedirs(output_dir, exist_ok=True)
        results_file = os.path.join(output_dir, f"{self.name}_results.json")
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
    
    def _get_output_files(self, output_dir: str) -> List[str]:
        import os
        files = []
        for file in os.listdir(output_dir):
            if file.startswith(self.name):
                files.append(os.path.join(output_dir, file))
        return files
    
    def get_metadata(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "version": self.version,
            "description": "Description of your analysis",
            "category": "custom",
            "dependencies": [],
            "output_type": "your_output_type"
        }
    
    def validate_input(self, proteins: List[Dict[str, Any]]) -> bool:
        if not proteins:
            return False
        for protein in proteins:
            if 'id' not in protein:
                return False
        return True
```

Create `__init__.py`:

```python
# lib/kbase_protein_query_module/src/analysis/your_analysis/__init__.py
from .your_analysis import YourAnalysis
__all__ = ['YourAnalysis']
```

### 2. Register Analysis

Add to configuration:

```python
# lib/kbase_protein_query_module/src/analysis/config.py
ANALYSIS_CONFIG = {
    "your_analysis": {
        "enabled": True,
        "name": "Your Custom Analysis",
        "description": "A custom analysis for protein data",
        "category": "custom",
        "dependencies": [],
        "output_type": "your_output_type",
        "module_path": "analysis.your_analysis.your_analysis",
        "class_name": "YourAnalysis"
    }
}
```

### 3. Add Output Handler

Create output handling:

```bash
mkdir -p lib/kbase_protein_query_module/src/outputs/analysis/your_analysis
```

```python
# lib/kbase_protein_query_module/src/outputs/analysis/your_analysis/output.py
from typing import Dict, Any, List
import os
import json

class YourAnalysisOutput:
    def __init__(self, output_manager):
        self.output_manager = output_manager
    
    def generate_output(self, results: Dict[str, Any], output_dir: str, 
                       format_type: str = "json") -> Dict[str, Any]:
        output_files = []
        
        if format_type == "json":
            output_files.extend(self._generate_json_output(results, output_dir))
        elif format_type == "csv":
            output_files.extend(self._generate_csv_output(results, output_dir))
        
        return {
            "output_files": output_files,
            "format": format_type,
            "analysis": "your_analysis"
        }
    
    def _generate_json_output(self, results: Dict[str, Any], output_dir: str) -> List[str]:
        output_file = os.path.join(output_dir, "your_analysis_results.json")
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        return [output_file]
    
    def _generate_csv_output(self, results: Dict[str, Any], output_dir: str) -> List[str]:
        import csv
        output_file = os.path.join(output_dir, "your_analysis_results.csv")
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['protein_id', 'score', 'confidence', 'details'])
            for protein_id, result in results.get('results', {}).items():
                writer.writerow([
                    protein_id,
                    result.get('score', ''),
                    result.get('confidence', ''),
                    result.get('details', '')
                ])
        return [output_file]
```

### 4. Write Tests

Create minimal and focused tests:

```python
# test/unit_tests/analysis/test_your_analysis.py
import pytest
import tempfile
from kbase_protein_query_module.src.analysis.your_analysis.your_analysis import YourAnalysis

class TestYourAnalysis:
    def setup_method(self):
        self.analysis = YourAnalysis()
        self.sample_proteins = [
            {"id": "P12345", "sequence": "MKFLVNVALVFMVVYISYIYA"},
            {"id": "P67890", "sequence": "MKFLVNVALVFMVVYISYIYA"}
        ]
    
    def test_analysis_initialization(self):
        assert self.analysis.name == "your_analysis"
        assert self.analysis.version == "1.0.0"
    
    def test_run_analysis_success(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            result = self.analysis.run(self.sample_proteins, temp_dir)
            assert result["status"] == "completed"
            assert result["analysis_name"] == "your_analysis"
            assert result["proteins_analyzed"] == 2
            assert "results" in result
            assert "output_files" in result
### 5. Add a server-side analysis test method

Add a minimal per-analysis test in `test/kbase_protein_query_module_server_test.py`:

```python
def test_your_analysis(self):
    params = {
        'workspace_name': self.wsName,
        'input_type': 'uniprot_ids',
        'uniprot_ids': ['P12345'],
        'analysis_name': 'your_analysis_test',
        'analysis_stages': ['your_analysis'],
        'output_config': {'output_dir': self.test_local_tmp}
    }
    _ = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
    self.workflow_mock.run_workflow.assert_called_once()
    self.assertIn('your_analysis', self.workflow_mock.run_workflow.return_value.analyses_completed)
```

    
    def test_validate_input(self):
        assert self.analysis.validate_input(self.sample_proteins) is True
        assert self.analysis.validate_input([]) is False
```

## Advanced Features

### Dependencies

Specify analysis dependencies:

```python
"your_analysis": {
    "dependencies": ["embeddings", "family_assignment"],
    # ... other config
}
```

### Parameters

Define configurable parameters:

```python
"your_analysis": {
    "parameters": {
        "threshold": {
            "type": "float",
            "default": 0.5,
            "description": "Analysis threshold"
        }
    }
}
```

## Best Practices

1. **Error Handling**: Implement comprehensive error handling
2. **Logging**: Use appropriate logging levels
3. **Resource Management**: Be mindful of memory and CPU usage
4. **Testing**: Write comprehensive unit and integration tests
5. **Documentation**: Document your analysis thoroughly

## Examples

See existing analyses for reference:
- `src/analysis/sequence_analysis/` - Sequence-based analysis
- `src/analysis/family_assignment/` - Family classification
- `src/analysis/similarity_search/` - Similarity search