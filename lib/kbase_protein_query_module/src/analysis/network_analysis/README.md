# Network Analysis

This module implements Protein Similarity Network Analysis. It constructs a network where nodes are proteins and edges represent similarity (based on embeddings).

## Components

*   **`network_analysis.py`**: The main class `NetworkAnalysis` that orchestrates the analysis. It performs similarity search, retrieves metadata, and prepares data for visualization.
*   **`network_visualizer.py`**: The `NetworkVisualizer` class that uses Plotly and NetworkX to create interactive 2D network visualizations.

## Features

*   **Similarity Search**: Uses FAISS (via `ProteinStorage`) to find similar proteins.
*   **Network Construction**: Builds a graph where edges connect similar proteins (above a threshold).
*   **Visualization**: Interactive HTML plots with:
    *   Kamada-Kawai layout.
    *   Node coloring based on query status.
    *   Rich hover information (metadata, network stats).
*   **Outputs**:
    *   Interactive HTML visualization.
    *   TSV files with network statistics and edge lists.
    *   TSV file with top-k matches and metadata.

## Usage

This module is typically invoked by the `AnalysisManager`.

```python
from kbase_protein_query_module.src.analysis.network_analysis import NetworkAnalysis

analysis = NetworkAnalysis(config={'k_neighbors': 20, 'similarity_threshold': 0.5})
results = analysis.run_network_analysis(input_data)
```
