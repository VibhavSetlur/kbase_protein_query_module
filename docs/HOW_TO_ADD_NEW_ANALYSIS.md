# How to Add a New Analysis to the KBase Protein Query Module

## Quick Start Guide

Adding a new analysis to this module is **extremely simple** thanks to the modular extensibility framework. Follow these **4 easy steps**:

## Step 1: Create Your Analysis Class (2 minutes)

Create a new file: `lib/kbase_protein_query_module/src/stages/analysis/your_analysis_name.py`

```python
"""
Your New Analysis Implementation

CLASS LOCATION: lib/kbase_protein_query_module/src/stages/analysis/your_analysis_name.py:10
EXTENDS: BaseAnalysis
USED BY: WorkflowOrchestrator when "your_analysis_name" is requested
"""

import os
import time
import logging
from typing import Dict, Any, List
import numpy as np
import pandas as pd

from ...core.analysis_registry import BaseAnalysis, AnalysisMetadata, register_analysis

@register_analysis("your_analysis_name")  # This makes it available immediately!
class YourAnalysisName(BaseAnalysis):
    """
    Your custom protein analysis implementation.
    
    This class will be automatically discovered and integrated into the pipeline.
    """
    
    def get_metadata(self) -> AnalysisMetadata:
        """Define your analysis metadata."""
        return AnalysisMetadata(
            name="Your Analysis Name",
            description="Brief description of what your analysis does",
            version="1.0.0",
            author="Your Name",
            output_files=["your_analysis.html", "your_data.csv"],
            dependencies=["numpy", "pandas"],  # Add any special dependencies
            category="custom",  # or "sequence", "structure", "functional", etc.
            computational_complexity="medium",  # "low", "medium", "high"
            memory_requirements="low"  # "low", "medium", "high"
        )
    
    def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
        """
        Implement your analysis logic here.
        
        Args:
            proteins: List of protein records with 'id' and 'sequence' fields
            **kwargs: Additional parameters from pipeline configuration
            
        Returns:
            Dictionary with your analysis results
        """
        results = {}
        
        self.logger.info(f"Starting {self.metadata.name} for {len(proteins)} proteins")
        
        for protein in proteins:
            protein_id = protein.get('id', 'unknown')
            sequence = protein.get('sequence', '')
            
            if sequence:
                # YOUR ANALYSIS LOGIC HERE
                analysis_result = self._perform_your_analysis(sequence)
                results[protein_id] = analysis_result
        
        return {
            'analysis_type': 'your_analysis_name',
            'results': results,
            'summary': f"Completed {self.metadata.name} for {len(results)} proteins",
            'metadata': {
                'total_proteins': len(proteins),
                'successful_analyses': len(results),
                'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
            }
        }
    
    def get_output_files(self, output_dir: str) -> List[str]:
        """Define what files your analysis will create."""
        return [
            os.path.join(output_dir, "your_analysis.html"),
            os.path.join(output_dir, "your_data.csv"),
            os.path.join(output_dir, "your_summary.json")
        ]
    
    def _perform_your_analysis(self, sequence: str) -> Dict[str, Any]:
        """
        Your specific analysis implementation.
        
        Example: Calculate sequence properties, predict structure, etc.
        """
        # Example analysis - replace with your logic
        return {
            'length': len(sequence),
            'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) if sequence else 0,
            'your_custom_metric': 42.0,  # Replace with your calculations
            'confidence': 0.95
        }
```

## Step 2: Create HTML Template (1 minute)

Create: `lib/kbase_protein_query_module/src/reports/html/templates/your_analysis_name.html`

```html
<!DOCTYPE html>
<html>
<head>
    <title>{{analysis_name}} Results</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .protein-result { border: 1px solid #ddd; margin: 10px; padding: 15px; }
        .metric { margin: 5px 0; }
        .confidence { font-weight: bold; color: #2e7d32; }
    </style>
</head>
<body>
    <h1>{{analysis_name}} Results</h1>
    <p><strong>Analysis completed:</strong> {{timestamp}}</p>
    <p><strong>Proteins analyzed:</strong> {{protein_count}}</p>
    
    <h2>Results Summary</h2>
    {{#each results}}
    <div class="protein-result">
        <h3>Protein: {{@key}}</h3>
        {{#each this}}
        <div class="metric">
            <strong>{{@key}}:</strong> {{this}}
        </div>
        {{/each}}
    </div>
    {{/each}}
    
    <h2>Download Data</h2>
    <p><a href="your_data.csv" download>Download CSV Data</a></p>
    <p><a href="your_summary.json" download>Download JSON Summary</a></p>
</body>
</html>
```

## Step 3: Test Your Analysis (30 seconds)

```python
# Test your new analysis
python -c "
import sys
sys.path.insert(0, 'lib')
from kbase_protein_query_module.src.core.analysis_registry import get_registry

# Your analysis is automatically available!
registry = get_registry()
analysis = registry.get_analysis('your_analysis_name')

# Test with sample data
proteins = [{'id': 'test_protein', 'sequence': 'MKVLLTLLCLAVAAL'}]
results = analysis.analyze(proteins)
print('✅ Your analysis works!', results)
"
```

## Step 4: Use in Pipeline (immediate)

```python
# Add to any pipeline run
params = {
    'workspace_name': 'your_workspace',
    'input_proteins': ['PROTEIN_SEQUENCE_HERE'],
    'analysis_stages': [
        'sequence_analysis',
        'your_analysis_name',  # Your new analysis is now available!
        'family_assignment'
    ]
}

result = impl.run_protein_query_analysis(ctx, params)
```

## **That's It! 🎉**

Your new analysis is now **fully integrated** and will:
- ✅ Appear automatically in the pipeline
- ✅ Generate its own HTML report section  
- ✅ Export data to CSV and JSON
- ✅ Be included in the comprehensive pipeline log
- ✅ Respect server resource limits
- ✅ Work with parallel processing

## Advanced Customization

### Custom Resource Requirements

```python
def estimate_resources(self, proteins: List[Any]) -> Dict[str, Any]:
    """Override to provide custom resource estimates."""
    return {
        "estimated_time_seconds": len(proteins) * 2.5,  # Your analysis timing
        "estimated_memory_mb": len(proteins) * 50,      # Your memory needs
        "cpu_cores_recommended": 2,                     # Your CPU needs
        "supports_parallel": True                       # Can run in parallel?
    }
```

### Custom Input Validation

```python
def validate_input(self, proteins: List[Any]) -> bool:
    """Override to add custom validation."""
    # Check if your analysis has special requirements
    for protein in proteins:
        if len(protein.get('sequence', '')) < 50:
            return False  # Need longer sequences
    return True
```

### Integration with Existing Data

```python
def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
    """Access results from other analyses."""
    
    # Get results from previous stages
    pipeline_results = kwargs.get('pipeline_results', {})
    
    # Use sequence analysis results if available
    if 'sequence_analysis' in pipeline_results:
        seq_results = pipeline_results['sequence_analysis']['results']
        # Use existing sequence analysis data in your analysis
    
    # Your analysis logic here...
```

## Module Organization

Your analysis will be automatically integrated into the module's organized structure:

```
lib/kbase_protein_query_module/
├── src/
│   ├── core/                           # ✅ Framework (extensibility)
│   │   ├── analysis_registry.py        # ✅ Your analysis registers here
│   │   ├── resource_manager.py         # ✅ Manages your analysis resources  
│   │   └── performance_monitor.py      # ✅ Monitors your analysis performance
│   ├── stages/
│   │   ├── analysis/
│   │   │   ├── sequence_analysis.py    # ✅ Existing analysis
│   │   │   ├── network_analysis.py     # ✅ Existing analysis  
│   │   │   └── your_analysis_name.py   # 🆕 YOUR NEW ANALYSIS HERE
│   │   ├── output/
│   │   │   ├── report_generation.py    # ✅ Will include your analysis
│   │   │   └── data_export.py          # ✅ Will export your data
│   └── reports/html/templates/
│       └── your_analysis_name.html     # 🆕 YOUR TEMPLATE HERE
```

## Output Integration

Your analysis will automatically be included in:

### 1. **Comprehensive HTML Dashboard** (`index.html`)
- Link to your analysis results
- Summary of your findings
- Performance metrics

### 2. **Individual Analysis Report** (`your_analysis_name.html`)  
- Detailed results from your analysis
- Custom visualizations
- Downloadable data files

### 3. **Data Exports**
- `your_data.csv`: Structured data from your analysis
- `your_summary.json`: Metadata and summary information
- `pipeline_log.json`: Your analysis included in complete pipeline log

### 4. **Pipeline Integration**
- Automatic resource management
- Progress tracking and logging
- Error handling and recovery
- Performance monitoring

## Testing Your Analysis

The module includes comprehensive testing infrastructure:

```python
# Create a unit test
# test/unit_tests/analysis/test_your_analysis_name.py

import unittest
from lib.kbase_protein_query_module.src.core.analysis_registry import get_registry

class TestYourAnalysis(unittest.TestCase):
    def setUp(self):
        self.registry = get_registry()
        self.analysis = self.registry.get_analysis("your_analysis_name")
    
    def test_analysis_execution(self):
        """Test that your analysis runs correctly."""
        proteins = [{"id": "test", "sequence": "MKVLLTLLCLAVAAL"}]
        results = self.analysis.analyze(proteins)
        
        self.assertIn('results', results)
        self.assertIn('test', results['results'])
        self.assertEqual(results['analysis_type'], 'your_analysis_name')

# Run tests
kb-sdk test  # Your analysis will be automatically tested!
```

## Real-World Example: Adding Hydrophobicity Analysis

Here's a complete real example:

```python
# lib/kbase_protein_query_module/src/stages/analysis/hydrophobicity_analysis.py

@register_analysis("hydrophobicity_analysis")
class HydrophobicityAnalysis(BaseAnalysis):
    def get_metadata(self) -> AnalysisMetadata:
        return AnalysisMetadata(
            name="Hydrophobicity Analysis",
            description="Analyzes protein hydrophobicity patterns using Kyte-Doolittle scale",
            version="1.0.0",
            author="Research Team", 
            output_files=["hydrophobicity.html", "hydrophobicity_data.csv"],
            dependencies=["numpy", "matplotlib"],
            category="sequence_properties",
            computational_complexity="low",
            memory_requirements="low"
        )
    
    def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
        # Kyte-Doolittle hydrophobicity scale
        hydro_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        results = {}
        for protein in proteins:
            sequence = protein.get('sequence', '')
            protein_id = protein.get('id', 'unknown')
            
            if sequence:
                hydro_values = [hydro_scale.get(aa, 0) for aa in sequence.upper()]
                results[protein_id] = {
                    'avg_hydrophobicity': np.mean(hydro_values),
                    'max_hydrophobicity': max(hydro_values),
                    'min_hydrophobicity': min(hydro_values),
                    'hydrophobic_regions': self._find_hydrophobic_regions(hydro_values)
                }
        
        return {
            'analysis_type': 'hydrophobicity_analysis',
            'results': results,
            'summary': f"Analyzed hydrophobicity for {len(results)} proteins"
        }

# Use immediately:
# params['analysis_stages'] = ['sequence_analysis', 'hydrophobicity_analysis']
```

## Module Integration Benefits

When you add a new analysis, it **automatically gets**:

### ✅ **Resource Management**
- Memory usage monitoring and optimization
- CPU usage limits (70% max on DOE servers)
- Automatic batch size optimization
- Garbage collection when needed

### ✅ **Parallel Processing** 
- Multi-core utilization when beneficial
- Thread pool management
- Task queue optimization
- Load balancing

### ✅ **Performance Monitoring**
- Execution time tracking
- Memory usage profiling  
- Performance bottleneck identification
- Optimization recommendations

### ✅ **Output Integration**
- Automatic HTML report generation
- CSV data export
- JSON metadata export
- Pipeline logging and tracking

### ✅ **Error Handling**
- Graceful failure recovery
- Detailed error logging
- Fallback strategies
- User-friendly error messages

## Advanced Integration Examples

### Accessing Previous Analysis Results

```python
def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
    # Access pipeline results from previous stages
    pipeline_results = kwargs.get('pipeline_results', {})
    
    # Use sequence analysis if available
    if 'sequence_analysis' in pipeline_results:
        seq_data = pipeline_results['sequence_analysis']['results']
        # Enhance your analysis with existing sequence data
    
    # Use family assignment if available  
    if 'family_assignment' in pipeline_results:
        family_data = pipeline_results['family_assignment']['results']
        # Use family information in your analysis
```

### Custom Configuration Options

```python
def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
    # Access custom configuration
    config = kwargs.get('analysis_config', {})
    your_threshold = config.get('your_threshold', 0.5)
    your_method = config.get('your_method', 'default')
    
    # Use configuration in your analysis
    if your_method == 'advanced':
        # Use advanced algorithm
        pass
```

### Resource-Aware Processing

```python
def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
    # Get resource manager for efficient processing
    resource_manager = kwargs.get('resource_manager')
    
    if resource_manager and len(proteins) > 100:
        # Use batch processing for large datasets
        batch_size = resource_manager.optimize_batch_size(50, len(proteins))
        
        results = {}
        for i in range(0, len(proteins), batch_size):
            batch = proteins[i:i+batch_size]
            batch_results = self._analyze_batch(batch)
            results.update(batch_results)
        
        return {'results': results}
    else:
        # Process normally for small datasets
        return self._analyze_all(proteins)
```

## Testing Your New Analysis

### Unit Test Template

```python
# test/unit_tests/analysis/test_your_analysis_name.py

import unittest
import os
import sys

# Add the lib directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..', 'lib'))

from kbase_protein_query_module.src.core.analysis_registry import get_registry

class TestYourAnalysis(unittest.TestCase):
    def setUp(self):
        self.registry = get_registry()
        self.analysis = self.registry.get_analysis("your_analysis_name")
    
    def test_analysis_basic(self):
        """Test basic analysis functionality."""
        proteins = [
            {"id": "protein1", "sequence": "MKVLLTLLCLAVAAL"},
            {"id": "protein2", "sequence": "MATLLTLLCLAVAAL"}
        ]
        
        results = self.analysis.analyze(proteins)
        
        # Verify structure
        self.assertIn('analysis_type', results)
        self.assertEqual(results['analysis_type'], 'your_analysis_name')
        self.assertIn('results', results)
        self.assertIn('summary', results)
        
        # Verify results for each protein
        for protein in proteins:
            self.assertIn(protein['id'], results['results'])
    
    def test_analysis_metadata(self):
        """Test analysis metadata."""
        metadata = self.analysis.get_metadata()
        
        self.assertEqual(metadata.name, "Your Analysis Name")
        self.assertIn("your_analysis.html", metadata.output_files)
        self.assertEqual(metadata.category, "custom")
    
    def test_output_files(self):
        """Test output file generation."""
        output_dir = "/tmp/test_output"
        files = self.analysis.get_output_files(output_dir)
        
        self.assertTrue(any("your_analysis.html" in f for f in files))
        self.assertTrue(any("your_data.csv" in f for f in files))

if __name__ == '__main__':
    unittest.main()
```

### Integration Test

```python
# Add to existing integration test
def test_your_new_analysis_integration(self):
    """Test your analysis in the full pipeline."""
    params = {
        'workspace_name': self.wsName,
        'input_proteins': ['MKVLLTLLCLAVAAL'],
        'analysis_stages': ['sequence_analysis', 'your_analysis_name']
    }
    
    result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
    
    # Verify your analysis was included
    self.assertIn('your_analysis_name', result[0]['stages_completed'])
```

## Running and Verifying

### Test Your Analysis
```bash
# Run all tests (including your new one)
kb-sdk test

# Test just your analysis
python -m pytest test/unit_tests/analysis/test_your_analysis_name.py -v
```

### Generate Documentation
```bash
# Update documentation with your new analysis
python -c "
from lib.kbase_protein_query_module.src.utils.documentation_generator import generate_full_documentation
generate_full_documentation('docs/')
"
```

### Test in Pipeline
```bash
# Test your analysis in the full pipeline
python -c "
import sys, os
sys.path.insert(0, 'lib')
os.makedirs('test/outputs', exist_ok=True)
os.environ['HTML_REPORTS_DIR'] = 'test/outputs'

from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module

impl = kbase_protein_query_module({})
ctx = {'provenance': [{'ws_name': 'demo'}]}
params = {
    'workspace_name': 'demo',
    'input_proteins': ['MKVLLTLLCLAVAAL'],
    'analysis_stages': ['sequence_analysis', 'your_analysis_name']
}

result = impl.run_protein_query_analysis(ctx, params)
print('✅ Pipeline with your analysis completed!')
print('📊 Check test/outputs/ for results')
"
```

## Key Integration Points

### 1. **Automatic Registry Integration**
- Your analysis is automatically discovered via `@register_analysis` decorator
- No manual registration or configuration files needed
- Available immediately in all pipeline runs

### 2. **Resource Management Integration**  
- Your analysis automatically benefits from server-aware resource limits
- Memory and CPU usage are monitored and optimized
- Batch processing is applied when beneficial

### 3. **Output System Integration**
- Your HTML template is automatically used for report generation
- Data exports are automatically created and included
- Output files are properly organized in timestamped directories

### 4. **Pipeline Logging Integration**
- Your analysis results are automatically included in pipeline logs
- Performance metrics are tracked and reported
- Success/failure status is monitored and logged

## Why This Architecture is Highly Scalable

### 1. **No Code Changes Required**
- Adding new analyses requires **zero changes** to existing code
- Registry pattern automatically discovers and integrates new components
- Existing pipeline logic handles new analyses transparently

### 2. **Server-Safe Resource Management**
- Percentage-based limits ensure respectful resource usage on shared servers
- Dynamic scaling adapts to available resources automatically
- Conservative defaults protect other server processes

### 3. **Modular and Maintainable**
- Each analysis is completely independent
- Clear separation of concerns
- Easy to test, debug, and maintain individual components

### 4. **Future-Proof Architecture**
- Extension points are clearly documented
- Base classes provide stable interfaces
- Backward compatibility is maintained automatically

The module is now **perfectly organized** and **ready for any new additions** while maintaining **full KBase SDK compliance** and **production-ready scalability**!
