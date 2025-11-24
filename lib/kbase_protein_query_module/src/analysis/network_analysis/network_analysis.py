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

# Set up path for script execution - add lib directory to path
# File is at: src/analysis/network_analysis/network_analysis.py
_script_file = os.path.abspath(__file__)
_current_dir = os.path.dirname(_script_file)  # analysis/network_analysis
_src_dir = os.path.dirname(os.path.dirname(_current_dir))  # src
_module_dir = os.path.dirname(_src_dir) # kbase_protein_query_module
_lib_dir = os.path.dirname(_module_dir) # lib

if _lib_dir not in sys.path:
    sys.path.insert(0, _lib_dir)

# Also add src to path for direct imports if needed (though full package imports are preferred)
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
    HAS_GRAPH_DEPS = True
except Exception:
    HAS_GRAPH_DEPS = False


# --- Utility module imports ---
# When running as script, use absolute imports; when imported as module, use relative
if __name__ == "__main__":
    # Script execution - use full package imports
    from kbase_protein_query_module.src.util.storage.storage import ProteinStorage
    from kbase_protein_query_module.src.util.embeddings.generator import ProteinEmbeddingGenerator
    from kbase_protein_query_module.src.util.uniprot.api import fetch_metadata, fetch_protein_sequence
else:
    # Module import - use relative imports
    try:
        from ...util.storage.storage import ProteinStorage
        from ...util.embeddings.generator import ProteinEmbeddingGenerator
        from ...util.uniprot.api import fetch_metadata, fetch_protein_sequence
    except ImportError:
        # Fallback to absolute if relative fails
        from kbase_protein_query_module.src.util.storage.storage import ProteinStorage
        from kbase_protein_query_module.src.util.embeddings.generator import ProteinEmbeddingGenerator
        from kbase_protein_query_module.src.util.uniprot.api import fetch_metadata, fetch_protein_sequence


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
        
        # Initialize Storage
        # Note: ProteinStorage now uses KBUtilLib/PLM API and doesn't require local embeddings file
        plm_api_url = self.config.get('plm_api_url', "https://kbase.us/services/llm_homology_api")
        self.storage = ProteinStorage(plm_api_url=plm_api_url, config=self.config)
        
        self.fetch_metadata = fetch_metadata
        self.fetch_protein_sequence = fetch_protein_sequence
        
        logger.info("NetworkAnalysis initialized")
    
    def run_network_analysis(self, input_data: Dict[str, Any]):
        """Run network analysis for a list of proteins."""
        start_time = time.time()
        
        if input_data is None:
            raise ValueError("Input data is required for network analysis")
        
        try:
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
            
            # Check if we have proteins dict (from input processing)
            proteins = input_data.get('proteins', {})
            
            # Handle single protein sequence input (for testing or simple calls)
            if not proteins and 'protein_sequence' in input_data:
                proteins = {
                    'query_protein': {
                        'sequence': input_data['protein_sequence'],
                        'source': 'user_input',
                        'original_id': 'Query'
                    }
                }

            if not isinstance(proteins, dict):
                raise ValueError(f"Expected proteins to be a dict, got {type(proteins)}")
            
            if proteins and len(proteins) > 0:
                # Process multiple proteins - run analysis for each protein separately
                logger.info(f"Processing {len(proteins)} proteins for network analysis")
                all_saved_files = []
                all_network_stats = []
                protein_results = []
                
                protein_items = list(proteins.items())
                for idx, (protein_id, protein_data) in enumerate(protein_items):
                    try:
                        original_id = protein_data.get('original_id', protein_id)
                        logger.info(f"Processing protein {idx + 1}/{len(proteins)}: {original_id} ({protein_id})")
                        
                        # Convert to format expected by _process_single_protein
                        protein = {
                            'protein_id': protein_id,
                            'sequence': protein_data.get('sequence', ''),
                            'source': protein_data.get('source', 'unknown'),
                            'original_id': original_id
                        }
                        
                        result = self._process_single_protein(protein, output_dir, idx)
                        if result and result.get('success'):
                            all_saved_files.extend(result.get('output_files', []))
                            all_network_stats.append(result.get('network_stats', {}))
                            protein_results.append({
                                'protein_id': protein_id,
                                'original_id': original_id,
                                'success': True,
                                'output_files': result.get('output_files', [])
                            })
                        else:
                            logger.warning(f"Failed to process protein {original_id}: {result.get('error_message', 'Unknown error') if result else 'No result'}")
                            protein_results.append({
                                'protein_id': protein_id,
                                'original_id': original_id,
                                'success': False,
                                'error_message': result.get('error_message', 'Unknown error') if result else 'No result'
                            })
                    except Exception as e:
                        logger.error(f"Error processing protein {protein_id}: {e}", exc_info=True)
                        protein_results.append({
                            'protein_id': protein_id,
                            'original_id': protein_data.get('original_id', protein_id),
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
                raise ValueError("No proteins provided for network analysis")
            
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
            
            # Generate embedding for query protein sequence (local fallback)
            local_query_embedding = self.embedding_generator.generate_embedding(protein_sequence, pooling="mean")
            query_embedding = local_query_embedding
            
            # Perform similarity search using ProteinStorage
            # Construct proteins dict for query
            proteins_query = {
                protein_id: {
                    'id': protein_id,
                    'sequence': protein_sequence
                }
            }
            
            # Query PLM API
            # Note: return_embeddings=True ensures we get embeddings for hits AND query if available
            search_results = self.storage.query_similar_proteins(
                proteins=proteins_query,
                max_hits=self.k_neighbors,
                similarity_threshold=self.similarity_threshold,
                return_embeddings=True
            )
            
            if protein_id not in search_results:
                raise ValueError(f"No search results returned for protein {protein_id}")
            
            result_entry = search_results[protein_id]
            
            # Use server-provided query embedding if available (avoids dimension mismatch)
            if 'query_embedding' in result_entry:
                query_embedding = result_entry['query_embedding']
                logger.info(f"Using server-provided query embedding (shape: {query_embedding.shape})")
            elif query_embedding is None or query_embedding.size == 0:
                 raise ValueError(f"Failed to generate embedding for protein {protein_id}")
            
            hits = result_entry.get('hits', [])
            
            if not hits:
                logger.warning(f"No similar proteins found for {protein_id}")
                # We can still proceed with just the query protein if desired, or raise error
                # For network analysis, a single node is boring but valid
            
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
                
            # Process hits
            top_k_similar_proteins_metadata = {}
            top_k_sequences = {}
            top_k_similarity_scores = {}
            top_k_embeddings = {}
            
            # Collect UniProt IDs for batch metadata fetch
            uniprot_ids = [hit.get('id') or hit.get('uniprot_id') for hit in hits if hit.get('id') or hit.get('uniprot_id')]
            
            # Fetch metadata in batch
            if uniprot_ids:
                try:
                    metadata_list = self.fetch_metadata(uniprot_ids)
                    for meta in metadata_list:
                        uid = meta.get('Entry')
                        if uid:
                            top_k_similar_proteins_metadata[uid] = meta
                except Exception as e:
                    logger.warning(f"Batch metadata fetch failed: {e}")

            for hit in hits:
                uniprot_id = hit.get('id') or hit.get('uniprot_id')
                if not uniprot_id:
                    continue
                    
                score = hit.get('score') or hit.get('plm_score', 0.0)
                top_k_similarity_scores[uniprot_id] = float(score)
                
                # Ensure metadata exists
                if uniprot_id not in top_k_similar_proteins_metadata:
                    top_k_similar_proteins_metadata[uniprot_id] = {'Entry': uniprot_id}
                
                # Fetch sequence (one by one for now, or could use batch if available)
                # KBUtilLib has get_uniprot_sequences but we are using util.uniprot.api here
                # We'll stick to util.uniprot.api for consistency within this module for now
                try:
                    seq = self.fetch_protein_sequence(uniprot_id)
                    top_k_sequences[uniprot_id] = seq if seq else ''
                except Exception as e:
                    logger.warning(f"Could not fetch sequence for {uniprot_id}: {e}")
                    top_k_sequences[uniprot_id] = ''
                
                # Get embedding
                # Try to get from hit first
                embedding = hit.get('embedding')
                if embedding is None:
                    # Try getting from storage cache
                    embedding = self.storage.get_embedding(uniprot_id)
                
                if embedding is not None:
                    top_k_embeddings[uniprot_id] = embedding
            
            # Combine embeddings: query + top k
            query_protein_id_key = query_protein_metadata.get('Entry', protein_id)
            protein_ids = [query_protein_id_key]
            embeddings_list = [query_embedding]
            all_metadata = {query_protein_id_key: query_protein_metadata}
            
            for uniprot_id, embedding in top_k_embeddings.items():
                protein_ids.append(uniprot_id)
                embeddings_list.append(embedding)
                all_metadata[uniprot_id] = top_k_similar_proteins_metadata.get(uniprot_id, {'Entry': uniprot_id})
            
            if len(embeddings_list) == 0:
                raise ValueError(f"No embeddings available for protein {protein_id}")
            
            # Stack embeddings
            # Ensure all are numpy arrays of same shape
            try:
                # Check shapes
                shapes = [e.shape for e in embeddings_list]
                if len(set(shapes)) > 1:
                    logger.warning(f"Embedding shapes mismatch: {shapes}. Truncating or skipping.")
                    # Simple fix: filter out mismatches (assuming query is correct)
                    target_shape = query_embedding.shape
                    valid_indices = [i for i, e in enumerate(embeddings_list) if e.shape == target_shape]
                    embeddings_list = [embeddings_list[i] for i in valid_indices]
                    protein_ids = [protein_ids[i] for i in valid_indices]
                
                embeddings = np.vstack([np.asarray(v, dtype=np.float32) for v in embeddings_list])
            except Exception as e:
                raise ValueError(f"Failed to stack embeddings: {e}")
            
            # Create visualization and save results
            saved_files = []
            network_stats = {}
            
            if HAS_GRAPH_DEPS:
                global NetworkVisualizer
                if NetworkVisualizer is None:
                    try:
                        from .network_visualizer import NetworkVisualizer
                    except ImportError:
                        # Fallback for script execution where relative import might fail
                        # Try absolute import assuming src is in path
                        try:
                            from analysis.network_analysis.network_visualizer import NetworkVisualizer
                        except ImportError:
                            # Try direct import if in same directory
                            from network_visualizer import NetworkVisualizer
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
                except Exception as e:
                    logger.warning(f"Could not copy/rename HTML file: {e}")
                    if os.path.exists(html_path):
                        saved_files.append(html_path)
            else:
                saved_files.append(html_path)
        
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
        # Test sequence (Human Insulin)
        test_sequence = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        
        input_data = {
            'input_type': 'protein_sequence',
            'protein_sequence': test_sequence,
            'output_dir': output_dir
        }
        
        print(f"Starting NetworkAnalysis self-test with output dir: {output_dir}")
        analysis = NetworkAnalysis(config={})
        result = analysis.run_network_analysis(input_data)
        
        if not isinstance(result, dict) or result.get('success') is not True:
            # Check if it failed due to missing dependencies (e.g. KBPLMUtils) which is expected in some envs
            error = result.get('error_message', '')
            if "KBPLMUtils not available" in error:
                print(f"ANALYSIS_SKIP: {error}")
            else:
                raise RuntimeError(f"Network analysis failed: {result}")
        else:
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
