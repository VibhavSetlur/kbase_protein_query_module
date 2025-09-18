"""
Documentation Generator for Module Extensions

This utility automatically generates documentation for developers showing
class locations, dependencies, and extension points in the module.

USAGE:
------
python -c "
from lib.kbase_protein_query_module.src.utils.documentation_generator import generate_full_documentation
generate_full_documentation('docs/')
"

CLASS LOCATION: lib/kbase_protein_query_module/src/utils/documentation_generator.py:15
USED BY: Developers for understanding module structure
"""

import os
import ast
import inspect
import importlib
from typing import Dict, Any, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)

class ClassTracker:
    """
    Tracks class definitions, usage, and dependencies across the module.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/utils/documentation_generator.py:25
    USED BY: DocumentationGenerator for creating class maps
    """
    
    def __init__(self, base_path: str = "lib/kbase_protein_query_module"):
        self.base_path = base_path
        self.classes: Dict[str, Dict[str, Any]] = {}
        self.imports: Dict[str, List[str]] = {}
        self.inheritance: Dict[str, List[str]] = {}
        
    def scan_codebase(self):
        """Scan the entire codebase for class definitions and usage."""
        for root, dirs, files in os.walk(self.base_path):
            # Skip __pycache__ and other non-source directories
            dirs[:] = [d for d in dirs if not d.startswith('__pycache__')]
            
            for file in files:
                if file.endswith('.py'):
                    filepath = os.path.join(root, file)
                    self._analyze_file(filepath)
    
    def _analyze_file(self, filepath: str):
        """Analyze a single Python file for class information."""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            tree = ast.parse(content)
            
            # Extract class definitions
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    self._extract_class_info(node, filepath)
                elif isinstance(node, ast.Import) or isinstance(node, ast.ImportFrom):
                    self._extract_import_info(node, filepath)
                    
        except Exception as e:
            logger.debug(f"Could not analyze {filepath}: {e}")
    
    def _extract_class_info(self, node: ast.ClassDef, filepath: str):
        """Extract information about a class definition."""
        class_name = node.name
        
        # Extract base classes
        base_classes = []
        for base in node.bases:
            if isinstance(base, ast.Name):
                base_classes.append(base.id)
            elif isinstance(base, ast.Attribute):
                base_classes.append(f"{base.value.id}.{base.attr}")
        
        # Extract docstring
        docstring = ""
        if (node.body and isinstance(node.body[0], ast.Expr) and 
            isinstance(node.body[0].value, ast.Constant)):
            docstring = node.body[0].value.value
        
        # Extract methods
        methods = []
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                methods.append(item.name)
        
        self.classes[class_name] = {
            'file': filepath,
            'line': node.lineno,
            'bases': base_classes,
            'methods': methods,
            'docstring': docstring[:200] + "..." if len(docstring) > 200 else docstring,
            'is_abstract': any('ABC' in base or 'abstract' in base.lower() for base in base_classes)
        }
        
        # Track inheritance
        for base in base_classes:
            if base not in self.inheritance:
                self.inheritance[base] = []
            self.inheritance[base].append(class_name)
    
    def _extract_import_info(self, node: ast.AST, filepath: str):
        """Extract import information."""
        if filepath not in self.imports:
            self.imports[filepath] = []
        
        if isinstance(node, ast.Import):
            for alias in node.names:
                self.imports[filepath].append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            for alias in node.names:
                self.imports[filepath].append(f"{module}.{alias.name}")
    
    def get_class_info(self, class_name: str) -> Optional[Dict[str, Any]]:
        """Get detailed information about a specific class."""
        return self.classes.get(class_name)
    
    def find_class_usage(self, class_name: str) -> List[str]:
        """Find files that import or use a specific class."""
        usage_files = []
        
        for filepath, imports in self.imports.items():
            if any(class_name in imp for imp in imports):
                usage_files.append(filepath)
        
        return usage_files
    
    def generate_class_map(self) -> str:
        """Generate a comprehensive class map document."""
        doc = "# Class Map and Dependencies\n\n"
        
        # Group classes by category
        categories = {
            'Core Framework': [],
            'Pipeline Stages': [],
            'Storage & Indexing': [],
            'Analysis Types': [],
            'Utilities': [],
            'Other': []
        }
        
        for class_name, info in self.classes.items():
            filepath = info['file']
            
            if 'core' in filepath:
                categories['Core Framework'].append((class_name, info))
            elif 'stages' in filepath:
                categories['Pipeline Stages'].append((class_name, info))
            elif 'storage' in filepath:
                categories['Storage & Indexing'].append((class_name, info))
            elif 'analysis' in filepath:
                categories['Analysis Types'].append((class_name, info))
            elif 'utils' in filepath:
                categories['Utilities'].append((class_name, info))
            else:
                categories['Other'].append((class_name, info))
        
        # Generate documentation for each category
        for category, classes in categories.items():
            if classes:
                doc += f"\n## {category}\n\n"
                
                for class_name, info in sorted(classes):
                    doc += f"### `{class_name}`\n"
                    doc += f"- **Location**: `{info['file']}:{info['line']}`\n"
                    
                    if info['bases']:
                        doc += f"- **Extends**: {', '.join(info['bases'])}\n"
                    
                    if class_name in self.inheritance:
                        doc += f"- **Extended by**: {', '.join(self.inheritance[class_name])}\n"
                    
                    usage_files = self.find_class_usage(class_name)
                    if usage_files:
                        doc += f"- **Used in**: {len(usage_files)} files\n"
                    
                    if info['docstring']:
                        doc += f"- **Description**: {info['docstring']}\n"
                    
                    if info['methods']:
                        key_methods = [m for m in info['methods'] if not m.startswith('_')][:5]
                        if key_methods:
                            doc += f"- **Key Methods**: {', '.join(key_methods)}\n"
                    
                    doc += "\n"
        
        return doc

class DocumentationGenerator:
    """
    Generates comprehensive documentation for module extensions.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/utils/documentation_generator.py:120
    USED BY: Developers for understanding and extending the module
    """
    
    def __init__(self, base_path: str = "lib/kbase_protein_query_module"):
        self.base_path = base_path
        self.class_tracker = ClassTracker(base_path)
        
    def generate_extension_guide(self, output_dir: str):
        """Generate complete extension guide for developers."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Scan codebase
        self.class_tracker.scan_codebase()
        
        # Generate class map
        class_map = self.class_tracker.generate_class_map()
        with open(os.path.join(output_dir, 'CLASS_MAP.md'), 'w') as f:
            f.write(class_map)
        
        # Generate extension examples
        extension_examples = self._generate_extension_examples()
        with open(os.path.join(output_dir, 'EXTENSION_EXAMPLES.md'), 'w') as f:
            f.write(extension_examples)
        
        # Generate API reference
        api_reference = self._generate_api_reference()
        with open(os.path.join(output_dir, 'API_REFERENCE.md'), 'w') as f:
            f.write(api_reference)
        
        logger.info(f"Documentation generated in {output_dir}")
    
    def _generate_extension_examples(self) -> str:
        """Generate practical examples for extending the module."""
        return """
# Extension Examples

## Adding a New Protein Analysis

### Example: Hydrophobicity Analysis

```python
# File: lib/kbase_protein_query_module/src/stages/analysis/hydrophobicity_analysis.py

from typing import Dict, Any, List
from ...core.analysis_registry import BaseAnalysis, AnalysisMetadata, register_analysis
import numpy as np

@register_analysis("hydrophobicity_analysis")
class HydrophobicityAnalysis(BaseAnalysis):
    \"\"\"
    Analyzes protein hydrophobicity patterns.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/stages/analysis/hydrophobicity_analysis.py:8
    EXTENDS: BaseAnalysis
    USED BY: WorkflowOrchestrator when "hydrophobicity_analysis" is requested
    \"\"\"
    
    def get_metadata(self) -> AnalysisMetadata:
        return AnalysisMetadata(
            name="Hydrophobicity Analysis",
            description="Analyzes hydrophobicity patterns in protein sequences",
            version="1.0.0",
            author="Research Team",
            output_files=["hydrophobicity.html", "hydrophobicity_data.csv"],
            dependencies=["numpy", "matplotlib"],
            category="sequence_properties",
            computational_complexity="low",
            memory_requirements="low"
        )
    
    def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
        \"\"\"Perform hydrophobicity analysis.\"\"\"
        results = {}
        
        # Hydrophobicity scale (Kyte-Doolittle)
        hydro_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        for protein in proteins:
            sequence = protein.get('sequence', '')
            protein_id = protein.get('id', 'unknown')
            
            if sequence:
                # Calculate hydrophobicity profile
                hydro_values = [hydro_scale.get(aa, 0) for aa in sequence.upper()]
                
                results[protein_id] = {
                    'avg_hydrophobicity': np.mean(hydro_values),
                    'hydrophobic_regions': self._find_hydrophobic_regions(hydro_values),
                    'hydrophobicity_profile': hydro_values
                }
        
        return {
            'analysis_type': 'hydrophobicity_analysis',
            'results': results,
            'summary': f"Analyzed hydrophobicity for {len(results)} proteins"
        }
    
    def get_output_files(self, output_dir: str) -> List[str]:
        return [
            os.path.join(output_dir, "hydrophobicity_analysis.html"),
            os.path.join(output_dir, "hydrophobicity_data.csv")
        ]
    
    def _find_hydrophobic_regions(self, values: List[float], window_size: int = 7) -> List[Tuple[int, int]]:
        \"\"\"Find hydrophobic regions using sliding window.\"\"\"
        regions = []
        for i in range(len(values) - window_size + 1):
            window_avg = np.mean(values[i:i+window_size])
            if window_avg > 1.0:  # Hydrophobic threshold
                regions.append((i, i+window_size))
        return regions
```

## Adding a New Indexing Strategy

### Example: LSH (Locality Sensitive Hashing) Index

```python
# File: lib/kbase_protein_query_module/src/storage/lsh_indexing.py

import numpy as np
from typing import Dict, Any, List, Tuple
from .indexing_strategy import IndexingStrategy, IndexingConfig, register_indexing_strategy
import hashlib

@register_indexing_strategy("lsh")
class LSHIndexingStrategy(IndexingStrategy):
    \"\"\"
    Locality Sensitive Hashing indexing strategy for approximate similarity search.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/lsh_indexing.py:10
    EXTENDS: IndexingStrategy  
    USED BY: ProteinStorage when strategy="lsh"
    
    CUSTOMIZATION POINTS:
    - Override _generate_hash_functions() for different hash families
    - Override _create_buckets() for custom bucketing strategies
    - Override _tune_parameters() for dataset-specific optimization
    \"\"\"
    
    def __init__(self, config: IndexingConfig):
        super().__init__(config)
        self.num_hash_functions = config.custom_params.get('num_hash_functions', 10) if config.custom_params else 10
        self.num_buckets = config.custom_params.get('num_buckets', 100) if config.custom_params else 100
        self.hash_functions = []
        self.buckets = {}
    
    def build_index(self, embeddings: np.ndarray, metadata: Dict[str, Any] = None) -> Any:
        \"\"\"Build LSH index with random projection hash functions.\"\"\"
        start_time = time.time()
        
        processed_embeddings = self._preprocess_embeddings(embeddings)
        dimension = processed_embeddings.shape[1]
        
        # Generate hash functions
        self.hash_functions = self._generate_hash_functions(dimension)
        
        # Create buckets and assign embeddings
        self.buckets = {}
        for idx, embedding in enumerate(processed_embeddings):
            bucket_key = self._hash_embedding(embedding)
            if bucket_key not in self.buckets:
                self.buckets[bucket_key] = []
            self.buckets[bucket_key].append(idx)
        
        self.index = {
            'hash_functions': self.hash_functions,
            'buckets': self.buckets,
            'embeddings': processed_embeddings
        }
        
        self.build_time = time.time() - start_time
        self.is_built = True
        
        return self.index
    
    def _generate_hash_functions(self, dimension: int) -> List[np.ndarray]:
        \"\"\"Generate random projection hash functions.\"\"\"
        hash_functions = []
        for _ in range(self.num_hash_functions):
            # Random projection vector
            random_vector = np.random.normal(0, 1, dimension)
            hash_functions.append(random_vector)
        return hash_functions
```

## Adding Custom Output Formats

### Example: JSON Export Stage

```python
# File: lib/kbase_protein_query_module/src/stages/output/json_export.py

from typing import Dict, Any, List
import json
import os
from ..base_stage import BaseStage

class JSONExportStage(BaseStage):
    \"\"\"
    Export analysis results in structured JSON format.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/stages/output/json_export.py:10
    EXTENDS: BaseStage
    USED BY: WorkflowOrchestrator for JSON output generation
    \"\"\"
    
    def run(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        \"\"\"Export pipeline results to structured JSON files.\"\"\"
        
        pipeline_results = input_data.get('pipeline_results', {})
        output_dir = input_data.get('output_directory')
        
        if not output_dir:
            raise ValueError("Output directory not specified")
        
        json_files = []
        
        # Export each analysis as separate JSON
        for analysis_name, results in pipeline_results.items():
            json_file = os.path.join(output_dir, f"{analysis_name}.json")
            
            with open(json_file, 'w') as f:
                json.dump({
                    'analysis_type': analysis_name,
                    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'results': results
                }, f, indent=2, default=str)
            
            json_files.append(json_file)
        
        return {
            'success': True,
            'output_data': {
                'json_files': json_files,
                'export_count': len(json_files)
            }
        }
```
"""
        
    def _generate_api_reference(self) -> str:
        """Generate API reference documentation."""
        return """
# API Reference

## Core Classes

### BaseAnalysis
**Location**: `lib/kbase_protein_query_module/src/core/analysis_registry.py:45`

Abstract base class for all protein analyses.

**Required Methods**:
- `analyze(proteins, **kwargs) -> Dict[str, Any]`
- `get_output_files(output_dir) -> List[str]`
- `get_metadata() -> AnalysisMetadata`

**Optional Methods**:
- `validate_input(proteins) -> bool`
- `estimate_resources(proteins) -> Dict[str, Any]`

### IndexingStrategy
**Location**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:65`

Abstract base class for indexing implementations.

**Required Methods**:
- `build_index(embeddings, metadata) -> Any`
- `search(query_embedding, k) -> List[Tuple[int, float]]`
- `get_index_info() -> Dict[str, Any]`
- `save_index(filepath) -> bool`
- `load_index(filepath) -> bool`

**Extension Points**:
- `_preprocess_embeddings(embeddings) -> np.ndarray`
- `_calculate_distance(query, candidates) -> np.ndarray`
- `_apply_quantization(embeddings) -> np.ndarray`

### ResourceManager
**Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:65`

Server-aware resource management for DOE environments.

**Key Methods**:
- `check_resource_availability() -> bool`
- `optimize_batch_size(base_size, data_size) -> int`
- `resource_context(operation_name)` (context manager)

## Configuration Classes

### PipelineConfig
**Location**: `lib/kbase_protein_query_module/src/core/pipeline_config.py:17`

Main configuration for pipeline execution.

### ResourceLimits  
**Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:20`

Server-aware resource limits using percentages.

## Extension Patterns

### Registry Pattern
Used for analyses and indexing strategies:
```python
@register_analysis("my_analysis")
class MyAnalysis(BaseAnalysis): ...

@register_indexing_strategy("my_index")  
class MyIndex(IndexingStrategy): ...
```

### Stage Pattern
Used for pipeline stages:
```python
class MyStage(BaseStage):
    def run(self, input_data): ...
```
"""

def generate_full_documentation(output_dir: str = "docs/"):
    """
    Generate complete documentation package for module developers.
    
    Creates:
    - Class map with locations and dependencies
    - Extension examples and patterns
    - API reference documentation
    - Performance guidelines
    """
    generator = DocumentationGenerator()
    generator.generate_extension_guide(output_dir)
    
    print(f"✅ Complete developer documentation generated in {output_dir}")
    print("📚 Files created:")
    print("  - CLASS_MAP.md: Complete class locations and dependencies")
    print("  - EXTENSION_EXAMPLES.md: Practical examples for adding new features")
    print("  - API_REFERENCE.md: Detailed API documentation")

if __name__ == "__main__":
    generate_full_documentation()
