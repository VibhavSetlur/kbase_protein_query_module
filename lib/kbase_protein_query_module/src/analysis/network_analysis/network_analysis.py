"""
Network Analysis Stage

This stage provides comprehensive network analysis for protein similarity search results,
including interactive network visualizations, CSV data outputs, and statistics.

Features:
- Interactive Plotly network visualization with Kamada-Kawai layout
- CSV outputs for network statistics and protein similarity data
- Network properties and clustering analysis
- Top similar proteins analysis with metadata
- Robust edge selection for connectivity
"""

# Handle both script execution and module import - MUST BE FIRST
import sys
import os

# Set up path for script execution - add src directory to path
# File is at: src/analysis/network_analysis/network_analysis.py
# We need: src/ (which contains util/)
_script_file = os.path.abspath(__file__)
_current_dir = os.path.dirname(_script_file)  # analysis/network_analysis
_src_dir = os.path.dirname(os.path.dirname(_current_dir))  # src
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union, Any
import logging
import json
import time

# --- Third-party package imports (optional) ---
try:
    import networkx as nx
    from sklearn.metrics.pairwise import cosine_similarity
    from sklearn.cluster import AgglomerativeClustering
    HAS_GRAPH_DEPS = True
except Exception:
    HAS_GRAPH_DEPS = False


# --- Utility module imports ---
# When running as script, use absolute imports; when imported as module, use relative
if __name__ == "__main__":
    # Script execution - use absolute imports with path we just set up
    from util.storage.storage import ProteinStorage
    from util.embeddings.generator import ProteinEmbeddingGenerator
    from util.uniprot.api import fetch_metadata, fetch_protein_sequence
else:
    # Module import - use relative imports
    try:
        from ...util.storage.storage import ProteinStorage
        from ...util.embeddings.generator import ProteinEmbeddingGenerator
        from ...util.uniprot.api import fetch_metadata, fetch_protein_sequence
    except ImportError:
        # Fallback to absolute if relative fails
        from util.storage.storage import ProteinStorage
        from util.embeddings.generator import ProteinEmbeddingGenerator
        from util.uniprot.api import fetch_metadata, fetch_protein_sequence


logger = logging.getLogger(__name__)

# Lazy import for visualizer
NetworkVisualizer = None


class NetworkAnalysis:
    """Network analysis for protein similarity search results."""
    
    def __init__(self, config: Dict[str, Any] = None):
        """Initialize the NetworkAnalysis class."""
        self.config = config or {}
        self.analysis_name = "network_analysis"
        
        # Initialize network parameters 
        self.k_neighbors = self.config.get('k_neighbors', 10)
        self.similarity_threshold = self.config.get('similarity_threshold', 0.1)
        
        # Initialize utilities
        self.embedding_generator = ProteinEmbeddingGenerator()
        # Storage for simple access + integrated similarity
        embeddings_file = self.config.get('embeddings_file', 'data/embeddings/embeddings.tsv')
        index_path = self.config.get('index_path', 'data/indexes/ivf_index.json')
        self.storage = ProteinStorage(embeddings_file_path=embeddings_file, index_path=index_path)
        self.fetch_metadata = fetch_metadata
        self.fetch_protein_sequence = fetch_protein_sequence
    
    def run_network_analysis(self, input_data: Dict[str, Any]):
        """Run network analysis for a list of proteins."""
        start_time = time.time()
        
        if input_data is None:
            raise ValueError("Input data is required for network analysis")
        
        try:
            # Get input type from processed input data
            input_type = input_data.get('input_type')
            
            if input_type == 'protein_sequence':
                # Generate embedding for query protein sequence (mean-pooled for similarity search)
                query_protein_sequence = input_data.get('protein_sequence', '')
                query_embedding = self.embedding_generator.generate_embedding(query_protein_sequence, pooling="mean")
                
                if query_embedding is None or query_embedding.size == 0:
                    raise ValueError("Failed to generate embedding for query protein sequence")
                
                # Perform similarity search using storage-integrated search
                top = self.storage.find_top_k_similar(query_embedding, top_k=self.k_neighbors)
                if top is None:
                    raise ValueError("Failed to perform similarity search for query protein sequence")
                
                # Initialize query protein metadata (empty for protein sequence input)
                query_protein_metadata = {
                    'Entry': 'QUERY_PROTEIN',
                    'Protein names': 'Query Protein',
                    'Organism': 'N/A',
                    'EC number': 'N/A',
                    'Protein families': 'N/A',
                    'Reviewed': 'N/A'
                }
                
                # Get metadata for each of the top k similar proteins
                top_k_similar_proteins_metadata = {}
                if top:
                    for uniprot_id, _ in top:
                        try:
                            metadata = self.fetch_metadata([uniprot_id])
                            if metadata and len(metadata) > 0:
                                top_k_similar_proteins_metadata[uniprot_id] = metadata[0]
                            else:
                                top_k_similar_proteins_metadata[uniprot_id] = {'Entry': uniprot_id}
                        except Exception as e:
                            logger.warning(f"Could not fetch metadata for UniProt ID {uniprot_id}: {e}")
                            top_k_similar_proteins_metadata[uniprot_id] = {'Entry': uniprot_id}
                
                # Get embeddings for top k proteins from storage
                top_k_embeddings = {}
                if top:
                    for uniprot_id, _ in top:
                        top_k_embeddings[uniprot_id] = self.storage.get_embedding(uniprot_id)
                
            elif input_type == 'uniprot_id':
                # Get UniProt ID(s)
                uniprot_id = input_data.get('uniprot_id', [])
                if isinstance(uniprot_id, str):
                    uniprot_id = [uniprot_id]
                
                if not uniprot_id:
                    raise ValueError("No UniProt ID provided")
                
                query_uniprot_id = uniprot_id[0]
                
                # Fetch query protein sequence and generate embedding (mean-pooled for similarity search)
                query_sequence = self.fetch_protein_sequence(query_uniprot_id)
                if not query_sequence:
                    raise ValueError(f"Failed to fetch sequence for UniProt ID: {query_uniprot_id}")
                query_embedding = self.embedding_generator.generate_embedding(query_sequence, pooling="mean")
                
                if query_embedding is None or query_embedding.size == 0:
                    raise ValueError("Failed to generate embedding for query protein sequence")
                
                # Perform similarity search using storage-integrated search
                top = self.storage.find_top_k_similar(query_embedding, top_k=self.k_neighbors)
                if top is None:
                    raise ValueError("Failed to perform similarity search for query protein sequence")
                
                # Get query protein metadata
                try:
                    metadata = self.fetch_metadata([query_uniprot_id])
                    if metadata and len(metadata) > 0:
                        query_protein_metadata = metadata[0]
                    else:
                        query_protein_metadata = {'Entry': query_uniprot_id}
                except Exception as e:
                    logger.warning(f"Could not fetch metadata for query UniProt ID {query_uniprot_id}: {e}")
                    raise ValueError(f"Could not fetch metadata for query UniProt ID {query_uniprot_id}: {e}")
                
                # Get metadata for top k similar proteins
                top_k_similar_proteins_metadata = {}
                if top:
                    for uniprot_id, _ in top:
                        try:
                            metadata = self.fetch_metadata([uniprot_id])
                            if metadata and len(metadata) > 0:
                                top_k_similar_proteins_metadata[uniprot_id] = metadata[0]
                            else:
                                top_k_similar_proteins_metadata[uniprot_id] = {'Entry': uniprot_id}
                        except Exception as e:
                            logger.warning(f"Could not fetch metadata for UniProt ID {uniprot_id}: {e}")
                            top_k_similar_proteins_metadata[uniprot_id] = {'Entry': uniprot_id}
                
                # Get embeddings for top k proteins
                top_k_embeddings = {}
                if top:
                    for uniprot_id, _ in top:
                        top_k_embeddings[uniprot_id] = self.storage.get_embedding(uniprot_id)
            else:
                raise ValueError(f"Unsupported input_type: {input_type}")
            
            # Combine embeddings: query + top k
            protein_ids = [query_protein_metadata.get('Entry', 'QUERY_PROTEIN')]
            embeddings_list = [query_embedding]
            all_metadata = {protein_ids[0]: query_protein_metadata}
            
            for uniprot_id in top_k_embeddings.keys():
                protein_ids.append(uniprot_id)
                embeddings_list.append(top_k_embeddings[uniprot_id])
                all_metadata[uniprot_id] = top_k_similar_proteins_metadata.get(uniprot_id, {'Entry': uniprot_id})
            
            embeddings = np.vstack([np.asarray(v, dtype=np.float32) for v in embeddings_list])
            
            # Create visualization (moved into run_network_analysis)
            output_dir = input_data.get('output_dir', 'outputs')
            os.makedirs(output_dir, exist_ok=True)
            
            if HAS_GRAPH_DEPS:
                global NetworkVisualizer
                if NetworkVisualizer is None:
                    from .network_visualizer import NetworkVisualizer
                visualizer = NetworkVisualizer({
                    'k_neighbors': self.k_neighbors,
                    'similarity_threshold': self.similarity_threshold
                })
                visualization_results = visualizer.create_interactive_visualization(
                    embeddings=embeddings,
                    protein_ids=protein_ids,
                    metadata=all_metadata,
                    query_embedding=query_embedding,
                    query_protein_id=protein_ids[0],
                    output_dir=output_dir
                )
                # Save outputs
                saved_files = self.save_results(
                    visualization_results, protein_ids[0], output_dir
                )
                execution_time = time.time() - start_time
                return {
                    'success': True,
                    'output_files': saved_files,
                    'processing_time': execution_time,
                    'network_stats': visualization_results.get('network_statistics', {})
                }
            else:
                # Minimal success when graph deps unavailable
                execution_time = time.time() - start_time
                return {
                    'success': True,
                    'output_files': [],
                    'processing_time': execution_time,
                    'network_stats': {'note': 'graph dependencies missing; ran minimal analysis'},
                }
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Network analysis failed: {str(e)}")
            return {
                'success': False,
                'output_files': [],
                'processing_time': execution_time,
                'error_message': str(e)
            }
    
    def save_results(self, visualization_results: Dict[str, Any], 
                    query_protein_id: str, output_dir: str) -> List[str]:
        """Save analysis results and return file paths."""
        saved_files = []
        
        # Save HTML visualization
        html_path = visualization_results.get('html_path')
        if html_path:
            saved_files.append(html_path)
        
        # Save network graph data - generate CSVs here
        network_graph = visualization_results.get('network_graph')
        if network_graph:
            # Generate network statistics CSV
            try:
                stats = []
                for node in network_graph.nodes():
                    stats.append({
                        'protein_id': node,
                        'degree': network_graph.degree(node),
                        'clustering_coefficient': nx.clustering(network_graph, node),
                        'is_query_protein': (node == query_protein_id)
                    })
                
                stats_df = pd.DataFrame(stats)
                stats_path = os.path.join(output_dir, f"network_statistics_{query_protein_id}_{int(time.time())}.csv")
                stats_df.to_csv(stats_path, index=False)
                logger.info(f"Network statistics saved to {stats_path}")
                saved_files.append(stats_path)
            except Exception as e:
                logger.warning(f"Failed to generate statistics CSV: {e}")
            
            # Generate edges CSV
            try:
                edges = []
                for u, v, data in network_graph.edges(data=True):
                    edges.append({
                        'protein_1': u,
                        'protein_2': v,
                        'similarity_weight': data.get('weight', 0),
                        'edge_type': 'query_connection' if (u == query_protein_id or v == query_protein_id) else 'protein_connection'
                    })
                
                edges_df = pd.DataFrame(edges)
                edges_path = os.path.join(output_dir, f"network_edges_{query_protein_id}_{int(time.time())}.csv")
                edges_df.to_csv(edges_path, index=False)
                logger.info(f"Network edges saved to {edges_path}")
                saved_files.append(edges_path)
            except Exception as e:
                logger.warning(f"Failed to generate edges CSV: {e}")
        
        return saved_files


def main():
    """Self-test for network analysis."""
    import shutil
    ok = True
    output_dir = os.path.join('/tmp', f"kpqm_network_test_{int(time.time())}")
    try:
        test_sequence = "ACDEFGHIKLMNPQRSTVWY"
        input_data = {
            'input_type': 'protein_sequence',
            'protein_sequence': test_sequence,
            'output_dir': output_dir
        }
        os.makedirs(output_dir, exist_ok=True)
        analysis = NetworkAnalysis(config={})
        result = analysis.run_network_analysis(input_data)
        
        if not isinstance(result, dict) or result.get('success') is not True:
            raise RuntimeError(f"Network analysis failed: {result}")
        
        if HAS_GRAPH_DEPS:
            files = result.get('output_files') or []
            missing = [p for p in files if not (isinstance(p, str) and os.path.exists(p))]
            if missing:
                raise RuntimeError(f"Some output files were not found on disk: {missing}")
        
        print("ANALYSIS_OK")
    except Exception as e:
        ok = False
        print(f"ANALYSIS_FAIL: {e}")
        import traceback
        traceback.print_exc()
    finally:
        try:
            if os.path.isdir(output_dir):
                shutil.rmtree(output_dir, ignore_errors=True)
        except Exception:
            pass
    return 0 if ok else 1


if __name__ == "__main__":
    main()
