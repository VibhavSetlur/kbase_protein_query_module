"""
Network Analysis Stage

This stage provides comprehensive network analysis for protein similarity search results,
including interactive network visualizations, TSV data outputs, and statistics.

Features:
- Interactive Plotly network visualization with Kamada-Kawai layout
- TSV outputs for network statistics and protein similarity data
- Network properties and clustering analysis
- Top similar proteins analysis with metadata
- Edge selection for connectivity
- Interactive hover metadata display
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


def _resolve_module_root() -> str:
    """
    Resolve the module root directory.
    In Docker: /kb/module
    In local dev: project root (parent of lib/)
    """
    # Check environment variable first (set by KBase SDK)
    module_root = os.environ.get('KB_MODULE_DIR')
    if module_root and os.path.exists(module_root):
        if os.path.exists(os.path.join(module_root, 'data', 'embeddings')):
            return module_root
    
    # Try standard Docker path
    docker_path = '/kb/module'
    if os.path.exists(docker_path):
        if os.path.exists(os.path.join(docker_path, 'data', 'embeddings')):
            return docker_path
    
    # Try to find module root by looking for common markers
    current_file = os.path.abspath(__file__)
    
    # Navigate up from: lib/kbase_protein_query_module/src/analysis/network_analysis/network_analysis.py
    # To module root (parent of lib/)
    # Current file is at: .../lib/kbase_protein_query_module/src/analysis/network_analysis/network_analysis.py
    # Module root should be 5 levels up
    lib_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(current_file))))  # lib/
    possible_roots = [
        os.path.dirname(lib_dir),  # Parent of lib/ (most likely module root)
        os.path.join(lib_dir, '..'),  # Same, but using join
        os.getcwd(),  # Current working directory
    ]
    
    # Normalize and check each possible root
    for root in possible_roots:
        root = os.path.abspath(os.path.normpath(root))
        # Check for markers that indicate module root
        data_embeddings = os.path.join(root, 'data', 'embeddings')
        if os.path.exists(data_embeddings):
            return root
        # Also check for lib directory and Makefile
        if os.path.exists(os.path.join(root, 'lib', 'kbase_protein_query_module')) and os.path.exists(os.path.join(root, 'Makefile')):
            return root
    
    # Last resort: return parent of lib directory
    return os.path.dirname(lib_dir) if os.path.exists(os.path.dirname(lib_dir)) else os.getcwd()


def _resolve_data_path(relative_path: str, module_root: str = None) -> str:
    """
    Resolve a data file path relative to module root.
    
    Args:
        relative_path: Relative path like 'data/embeddings/embeddings.tsv'
        module_root: Module root directory (auto-detected if None)
    
    Returns:
        Absolute path to the data file
    """
    if module_root is None:
        module_root = _resolve_module_root()
    
    # If path is already absolute, return as-is
    if os.path.isabs(relative_path):
        return relative_path
    
    # Resolve relative to module root
    abs_path = os.path.join(module_root, relative_path)
    abs_path = os.path.normpath(abs_path)
    
    return abs_path


class NetworkAnalysis:
    """Network analysis for protein similarity search results."""
    
    def __init__(self, config: Dict[str, Any] = None):
        """Initialize the NetworkAnalysis class."""
        self.config = config or {}
        self.analysis_name = "network_analysis"
        
        # Initialize network parameters 
        self.k_neighbors = self.config.get('k_neighbors', 10)
        self.similarity_threshold = self.config.get('similarity_threshold', 0.1)
        
        # Resolve module root for data paths
        self.module_root = self.config.get('module_root') or _resolve_module_root()
        logger.debug(f"NetworkAnalysis using module_root: {self.module_root}")
        
        # Initialize utilities
        self.embedding_generator = ProteinEmbeddingGenerator()
        
        # Storage for simple access + integrated similarity
        # Resolve embeddings file path - try config first, then resolve relative to module root
        embeddings_file = self.config.get('embeddings_file')
        if not embeddings_file:
            # Default relative path
            embeddings_file = 'data/embeddings/embeddings.tsv'
        
        # Resolve to absolute path
        embeddings_file = _resolve_data_path(embeddings_file, self.module_root)
        logger.debug(f"NetworkAnalysis using embeddings_file: {embeddings_file}")
        
        # Resolve index path
        index_path = self.config.get('index_path', 'data/indexes')
        if not os.path.isabs(index_path):
            index_path = _resolve_data_path(index_path, self.module_root)
        logger.debug(f"NetworkAnalysis using index_path: {index_path}")
        
        # Store paths for lazy initialization
        self.embeddings_file = embeddings_file
        self.index_path = index_path
        self.storage = None  # Will be initialized on first use
        self.fetch_metadata = fetch_metadata
        self.fetch_protein_sequence = fetch_protein_sequence
        
        # Check if embeddings file exists and log warning if not
        if not os.path.exists(embeddings_file):
            logger.warning(f"Embeddings file not found: {embeddings_file}")
            logger.warning(f"Module root: {self.module_root}")
            logger.warning(f"Tried paths: {embeddings_file}")
            if os.path.exists(self.module_root):
                try:
                    contents = os.listdir(self.module_root)[:10]
                    logger.warning(f"Contents of module root: {contents}")
                except Exception:
                    pass
            # Don't raise error here - allow analysis to be loaded
            # Storage will be initialized lazily when needed, and will raise then if file still missing
        else:
            logger.info(f"Found embeddings file: {embeddings_file}")
    
    def _ensure_storage(self):
        """Ensure storage is initialized. Initialize lazily if not already done."""
        if self.storage is None:
            if not os.path.exists(self.embeddings_file):
                raise FileNotFoundError(
                    f"Embeddings file not found: {self.embeddings_file}\n"
                    f"Module root: {self.module_root}\n"
                    f"Please ensure the embeddings file exists at the specified path."
                )
            self.storage = ProteinStorage(embeddings_file_path=self.embeddings_file, index_path=self.index_path)
            logger.info(f"Initialized ProteinStorage with {self.storage.n} proteins")
    
    def run_network_analysis(self, input_data: Dict[str, Any]):
        """Run network analysis for a list of proteins."""
        start_time = time.time()
        
        if input_data is None:
            raise ValueError("Input data is required for network analysis")
        
        try:
            # Ensure storage is initialized
            self._ensure_storage()
            
            # Get output directory - this should be analysis/network_analysis
            output_dir = input_data.get('output_dir', 'outputs')
            try:
                os.makedirs(output_dir, exist_ok=True)
                # Verify directory is writable
                if not os.path.isdir(output_dir):
                    raise ValueError(f"Output directory is not a directory: {output_dir}")
                if not os.access(output_dir, os.W_OK):
                    raise ValueError(f"Output directory is not writable: {output_dir}")
                logger.info(f"Network analysis output directory: {output_dir} (exists and writable)")
            except Exception as e:
                logger.error(f"Failed to create or verify output directory {output_dir}: {e}")
                raise
            
            # Check if we have a proteins list (from input processing)
            proteins = input_data.get('proteins', [])
            
            if proteins and isinstance(proteins, list) and len(proteins) > 0:
                # Process multiple proteins - run analysis for each protein separately
                logger.info(f"Processing {len(proteins)} proteins for network analysis")
                all_saved_files = []
                all_network_stats = []
                protein_results = []
                
                for idx, protein in enumerate(proteins):
                    try:
                        logger.info(f"Processing protein {idx + 1}/{len(proteins)}: {protein.get('protein_id', 'unknown')}")
                        result = self._process_single_protein(protein, output_dir, idx)
                        if result and result.get('success'):
                            all_saved_files.extend(result.get('output_files', []))
                            all_network_stats.append(result.get('network_stats', {}))
                            protein_results.append({
                                'protein_id': protein.get('protein_id', 'unknown'),
                                'success': True,
                                'output_files': result.get('output_files', [])
                            })
                        else:
                            logger.warning(f"Failed to process protein {protein.get('protein_id', 'unknown')}: {result.get('error_message', 'Unknown error')}")
                            protein_results.append({
                                'protein_id': protein.get('protein_id', 'unknown'),
                                'success': False,
                                'error_message': result.get('error_message', 'Unknown error')
                            })
                    except Exception as e:
                        logger.error(f"Error processing protein {protein.get('protein_id', 'unknown')}: {e}", exc_info=True)
                        protein_results.append({
                            'protein_id': protein.get('protein_id', 'unknown'),
                            'success': False,
                            'error_message': str(e)
                        })
                
                execution_time = time.time() - start_time
                return {
                    'success': True,
                    'output_files': all_saved_files,
                    'processing_time': execution_time,
                    'network_stats': all_network_stats,
                    'proteins_processed': len(proteins),
                    'protein_results': protein_results
                }
            else:
                # Legacy single protein processing (for backward compatibility)
                logger.info("Processing single protein (legacy mode)")
                input_type = input_data.get('input_type')
                
                if input_type == 'protein_sequence':
                    query_protein_sequence = input_data.get('protein_sequence', '')
                    if not query_protein_sequence:
                        raise ValueError("No protein sequence provided")
                    protein = {
                        'protein_id': 'QUERY_PROTEIN',
                        'sequence': query_protein_sequence,
                        'source': 'protein_sequence'
                    }
                elif input_type == 'uniprot_id':
                    uniprot_id = input_data.get('uniprot_id', [])
                    if isinstance(uniprot_id, str):
                        uniprot_id = [uniprot_id]
                    if not uniprot_id:
                        raise ValueError("No UniProt ID provided")
                    query_uniprot_id = uniprot_id[0]
                    
                    # Fetch sequence
                    query_sequence = self.fetch_protein_sequence(query_uniprot_id)
                    if not query_sequence:
                        raise ValueError(f"Failed to fetch sequence for UniProt ID: {query_uniprot_id}")
                    protein = {
                        'protein_id': query_uniprot_id,
                        'sequence': query_sequence,
                        'source': 'uniprot'
                    }
                else:
                    raise ValueError(f"Unsupported input_type: {input_type}")
                
                result = self._process_single_protein(protein, output_dir, 0)
                execution_time = time.time() - start_time
                if result:
                    result['processing_time'] = execution_time
                return result
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Network analysis failed: {str(e)}", exc_info=True)
            return {
                'success': False,
                'output_files': [],
                'processing_time': execution_time,
                'error_message': str(e)
            }
    
    def _process_single_protein(self, protein: Dict[str, Any], output_dir: str, protein_index: int) -> Dict[str, Any]:
        """Process a single protein for network analysis."""
        try:
            # Ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            
            protein_id = protein.get('protein_id', f'PROTEIN_{protein_index}')
            protein_sequence = protein.get('sequence', '')
            protein_source = protein.get('source', 'unknown')
            
            logger.debug(f"Processing protein {protein_index}: {protein_id} (source: {protein_source}, sequence length: {len(protein_sequence) if protein_sequence else 0})")
            
            if not protein_sequence:
                raise ValueError(f"No sequence provided for protein {protein_id}")
            
            # Generate embedding for query protein sequence
            query_embedding = self.embedding_generator.generate_embedding(protein_sequence, pooling="mean")
            if query_embedding is None or query_embedding.size == 0:
                raise ValueError(f"Failed to generate embedding for protein {protein_id}")
            
            # Perform similarity search
            top = self.storage.find_top_k_similar(query_embedding, top_k=self.k_neighbors)
            if top is None or len(top) == 0:
                raise ValueError(f"Failed to perform similarity search for protein {protein_id}")
            
            # Get query protein metadata
            if protein_source == 'uniprot' and protein_id:
                try:
                    metadata_list = self.fetch_metadata([protein_id])
                    if metadata_list and len(metadata_list) > 0:
                        query_protein_metadata = metadata_list[0]
                    else:
                        query_protein_metadata = {'Entry': protein_id}
                except Exception as e:
                    logger.warning(f"Could not fetch metadata for query protein {protein_id}: {e}")
                    query_protein_metadata = {'Entry': protein_id}
            else:
                # For protein_sequence or unknown source, create minimal metadata
                query_protein_metadata = {
                    'Entry': protein_id,
                    'Protein names': f'Query Protein {protein_id}',
                    'Organism': 'N/A',
                    'EC number': 'N/A',
                    'Protein families': 'N/A',
                    'Reviewed': 'N/A'
                }
            
            # Get metadata and sequences for top k similar proteins
            top_k_similar_proteins_metadata = {}
            top_k_sequences = {}
            top_k_similarity_scores = {}
            
            # Extract similarity scores from top results
            for uniprot_id, similarity_score in top:
                top_k_similarity_scores[uniprot_id] = float(similarity_score)
                try:
                    metadata_list = self.fetch_metadata([uniprot_id])
                    if metadata_list and len(metadata_list) > 0:
                        top_k_similar_proteins_metadata[uniprot_id] = metadata_list[0]
                    else:
                        top_k_similar_proteins_metadata[uniprot_id] = {'Entry': uniprot_id}
                except Exception as e:
                    logger.warning(f"Could not fetch metadata for UniProt ID {uniprot_id}: {e}")
                    top_k_similar_proteins_metadata[uniprot_id] = {'Entry': uniprot_id}
                
                # Fetch sequence for each match
                try:
                    seq = self.fetch_protein_sequence(uniprot_id)
                    if seq:
                        top_k_sequences[uniprot_id] = seq
                    else:
                        top_k_sequences[uniprot_id] = ''
                except Exception as e:
                    logger.warning(f"Could not fetch sequence for UniProt ID {uniprot_id}: {e}")
                    top_k_sequences[uniprot_id] = ''
            
            # Get embeddings for top k proteins
            top_k_embeddings = {}
            for uniprot_id, _ in top:
                try:
                    embedding = self.storage.get_embedding(uniprot_id)
                    if embedding is not None:
                        top_k_embeddings[uniprot_id] = embedding
                except Exception as e:
                    logger.warning(f"Could not get embedding for UniProt ID {uniprot_id}: {e}")
            
            # Combine embeddings: query + top k
            query_protein_id_key = query_protein_metadata.get('Entry', protein_id)
            protein_ids = [query_protein_id_key]
            embeddings_list = [query_embedding]
            all_metadata = {query_protein_id_key: query_protein_metadata}
            
            for uniprot_id in top_k_embeddings.keys():
                protein_ids.append(uniprot_id)
                embeddings_list.append(top_k_embeddings[uniprot_id])
                all_metadata[uniprot_id] = top_k_similar_proteins_metadata.get(uniprot_id, {'Entry': uniprot_id})
            
            if len(embeddings_list) == 0:
                raise ValueError(f"No embeddings available for protein {protein_id}")
            
            embeddings = np.vstack([np.asarray(v, dtype=np.float32) for v in embeddings_list])
            
            # Create visualization and save results
            saved_files = []
            network_stats = {}
            
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
                    query_protein_id=query_protein_id_key,
                    output_dir=output_dir
                )
                
                # Save outputs (HTML, network stats, edges)
                saved_files = self.save_results(
                    visualization_results, query_protein_id_key, output_dir, protein_index
                )
                network_stats = visualization_results.get('network_statistics', {})
            else:
                logger.warning("Graph dependencies not available, skipping visualization")
            
            # Create TSV file with top k matches and UniProt data (including query)
            top_k_tsv_path = self._create_top_k_matches_tsv(
                query_protein_id_key,
                protein_sequence,
                query_protein_metadata,
                top_k_similar_proteins_metadata,
                top_k_sequences,
                top_k_similarity_scores,
                output_dir,
                protein_index
            )
            if top_k_tsv_path:
                saved_files.append(top_k_tsv_path)
                logger.info(f"Created top k matches TSV for protein {protein_id}: {top_k_tsv_path}")
            else:
                logger.warning(f"Failed to create top k matches TSV for protein {protein_id}")
            
            # Log all saved files
            logger.info(f"Saved {len(saved_files)} output files for protein {protein_id} to {output_dir}")
            for file_path in saved_files:
                if os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    logger.debug(f"  - {os.path.basename(file_path)} ({file_size} bytes)")
                else:
                    logger.warning(f"  - {os.path.basename(file_path)} (FILE NOT FOUND)")
            
            return {
                'success': True,
                'output_files': saved_files,
                'network_stats': network_stats,
                'protein_id': protein_id
            }
            
        except Exception as e:
            logger.error(f"Error processing single protein: {e}", exc_info=True)
            return {
                'success': False,
                'output_files': [],
                'error_message': str(e),
                'protein_id': protein.get('protein_id', 'unknown')
            }
    
    def save_results(self, visualization_results: Dict[str, Any], 
                    query_protein_id: str, output_dir: str, protein_index: int = 0) -> List[str]:
        """Save analysis results and return file paths."""
        saved_files = []
        
        # Create safe filename from protein_id (remove special characters)
        safe_protein_id = query_protein_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
        
        # Save HTML visualization - always ensure unique filename with protein_index
        html_path = visualization_results.get('html_path')
        if html_path and os.path.exists(html_path):
            # Create a new HTML file in the output directory with unique name based on protein_index
            new_html_path = os.path.join(output_dir, f"network_visualization_{safe_protein_id}_{protein_index}.html")
            
            # Only copy/rename if the paths are different
            if os.path.abspath(os.path.normpath(html_path)) != os.path.abspath(os.path.normpath(new_html_path)):
                try:
                    import shutil
                    shutil.copy2(html_path, new_html_path)
                    logger.info(f"Saved HTML visualization to {new_html_path}")
                    saved_files.append(new_html_path)
                    # Optionally remove the original file if it's in a temp location
                    # (but keep it if it's already in output_dir with a different name)
                except Exception as e:
                    logger.warning(f"Could not copy/rename HTML file: {e}")
                    # Fall back to original path if copy fails
                    if os.path.exists(html_path):
                        saved_files.append(html_path)
            else:
                # File is already at the correct location with correct name
                saved_files.append(html_path)
                logger.debug(f"HTML visualization already at correct location: {html_path}")
        
        # Save network graph data - generate TSV files
        network_graph = visualization_results.get('network_graph')
        if network_graph:
            # Generate network statistics TSV
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
                stats_path = os.path.join(output_dir, f"network_statistics_{safe_protein_id}_{protein_index}.tsv")
                stats_df.to_csv(stats_path, index=False, sep='\t')
                logger.info(f"Network statistics saved to {stats_path}")
                saved_files.append(stats_path)
            except Exception as e:
                logger.warning(f"Failed to generate statistics TSV: {e}")
            
            # Generate edges TSV
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
                edges_path = os.path.join(output_dir, f"network_edges_{safe_protein_id}_{protein_index}.tsv")
                edges_df.to_csv(edges_path, index=False, sep='\t')
                logger.info(f"Network edges saved to {edges_path}")
                saved_files.append(edges_path)
            except Exception as e:
                logger.warning(f"Failed to generate edges TSV: {e}")
        
        return saved_files
    
    def _create_top_k_matches_tsv(self, query_protein_id: str, query_sequence: str,
                                  query_metadata: Dict[str, Any],
                                  top_k_metadata: Dict[str, Dict[str, Any]],
                                  top_k_sequences: Dict[str, str],
                                  top_k_similarity_scores: Dict[str, float],
                                  output_dir: str, protein_index: int) -> Optional[str]:
        """Create a TSV file with top k matches and UniProt data for a protein."""
        try:
            # Create safe filename
            safe_protein_id = query_protein_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
            
            # Build rows for TSV - include query protein first, then top k matches
            rows = []
            
            # Add query protein as first row
            query_row = {
                'uniprot_id': query_protein_id,
                'protein_name': query_metadata.get('Protein names', ''),
                'organism': query_metadata.get('Organism', ''),
                'ec_number': query_metadata.get('EC number', ''),
                'protein_families': query_metadata.get('Protein families', ''),
                'reviewed': query_metadata.get('Reviewed', ''),
                'similarity_score': 1.0,  # Query protein has 100% similarity to itself
                'protein_sequence': query_sequence,
                'is_query': True
            }
            rows.append(query_row)
            
            # Add top k matches, sorted by similarity score (highest first)
            sorted_matches = sorted(
                top_k_similarity_scores.items(),
                key=lambda x: x[1],
                reverse=True
            )
            
            for uniprot_id, similarity_score in sorted_matches:
                metadata = top_k_metadata.get(uniprot_id, {'Entry': uniprot_id})
                sequence = top_k_sequences.get(uniprot_id, '')
                
                match_row = {
                    'uniprot_id': uniprot_id,
                    'protein_name': metadata.get('Protein names', ''),
                    'organism': metadata.get('Organism', ''),
                    'ec_number': metadata.get('EC number', ''),
                    'protein_families': metadata.get('Protein families', ''),
                    'reviewed': metadata.get('Reviewed', ''),
                    'similarity_score': similarity_score,
                    'protein_sequence': sequence,
                    'is_query': False
                }
                rows.append(match_row)
            
            # Create DataFrame and save to TSV
            df = pd.DataFrame(rows)
            tsv_path = os.path.join(output_dir, f"network_analysis_{safe_protein_id}_{protein_index}.tsv")
            df.to_csv(tsv_path, index=False, sep='\t')
            logger.info(f"Top k matches TSV saved to {tsv_path}")
            
            return tsv_path
            
        except Exception as e:
            logger.error(f"Failed to create top k matches TSV: {e}", exc_info=True)
            return None


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
