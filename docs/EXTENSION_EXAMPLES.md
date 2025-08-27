
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
    """
    Analyzes protein hydrophobicity patterns.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/stages/analysis/hydrophobicity_analysis.py:8
    EXTENDS: BaseAnalysis
    USED BY: WorkflowOrchestrator when "hydrophobicity_analysis" is requested
    """
    
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
        """Perform hydrophobicity analysis."""
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
        """Find hydrophobic regions using sliding window."""
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
    """
    Locality Sensitive Hashing indexing strategy for approximate similarity search.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/lsh_indexing.py:10
    EXTENDS: IndexingStrategy  
    USED BY: ProteinStorage when strategy="lsh"
    
    CUSTOMIZATION POINTS:
    - Override _generate_hash_functions() for different hash families
    - Override _create_buckets() for custom bucketing strategies
    - Override _tune_parameters() for dataset-specific optimization
    """
    
    def __init__(self, config: IndexingConfig):
        super().__init__(config)
        self.num_hash_functions = config.custom_params.get('num_hash_functions', 10) if config.custom_params else 10
        self.num_buckets = config.custom_params.get('num_buckets', 100) if config.custom_params else 100
        self.hash_functions = []
        self.buckets = {}
    
    def build_index(self, embeddings: np.ndarray, metadata: Dict[str, Any] = None) -> Any:
        """Build LSH index with random projection hash functions."""
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
        """Generate random projection hash functions."""
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
    """
    Export analysis results in structured JSON format.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/stages/output/json_export.py:10
    EXTENDS: BaseStage
    USED BY: WorkflowOrchestrator for JSON output generation
    """
    
    def run(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Export pipeline results to structured JSON files."""
        
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
