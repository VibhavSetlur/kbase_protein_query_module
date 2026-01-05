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

import sys
import os
import logging
import time
from typing import List, Dict, Tuple, Optional, Any, Union

import numpy as np
import pandas as pd

# --- Third-party package imports (optional) ---
try:
    import networkx as nx
    HAS_GRAPH_DEPS = True
except ImportError:
    HAS_GRAPH_DEPS = False

# --- Utility module imports ---
# Handle both script execution and module import
if __name__ == "__main__":
    # Script execution - use full package imports
    # Add lib directory to path if not present
    _script_file = os.path.abspath(__file__)
    _current_dir = os.path.dirname(_script_file)
    _src_dir = os.path.dirname(os.path.dirname(_current_dir))
    _module_dir = os.path.dirname(_src_dir)
    _lib_dir = os.path.dirname(_module_dir)
    
    if _lib_dir not in sys.path:
        sys.path.insert(0, _lib_dir)

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
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the NetworkAnalysis class.
        
        Args:
            config: Configuration dictionary.
        """
        self.config = config or {}
        self.analysis_name = "network_analysis"
        
        # Initialize network parameters 
        self.k_neighbors = self.config.get('k_neighbors', 10)
        self.similarity_threshold = self.config.get('similarity_threshold', 0.995)
        
        # Initialize utilities
        self.embedding_generator = ProteinEmbeddingGenerator()
        
        # Initialize Storage
        plm_api_url = self.config.get('plm_api_url', "https://kbase.us/services/llm_homology_api")
        self.storage = ProteinStorage(plm_api_url=plm_api_url, config=self.config)
        
        self.fetch_metadata = fetch_metadata
        self.fetch_protein_sequence = fetch_protein_sequence
        
        logger.info("NetworkAnalysis initialized")
    
    def run_network_analysis(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run network analysis for a list of proteins.
        
        Args:
            input_data: Dictionary containing input proteins and parameters.
            
        Returns:
            Dictionary containing analysis results.
        """
        start_time = time.time()
        
        if not input_data:
            raise ValueError("Input data is required for network analysis")
        
        try:
            output_dir = self._prepare_output_dir(input_data)
            proteins = self._get_proteins_from_input(input_data)
            
            logger.info(f"Processing {len(proteins)} proteins for network analysis")
            
            results = self._process_all_proteins(proteins, output_dir)
            
            execution_time = time.time() - start_time
            return {
                'success': True,
                'output_files': results['all_saved_files'],
                'processing_time': execution_time,
                'network_stats': results['all_network_stats'],
                'proteins_processed': len(proteins),
                'protein_results': results['protein_results']
            }
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Network analysis failed: {str(e)}", exc_info=True)
            return {
                'success': False,
                'output_files': [],
                'processing_time': execution_time,
                'error_message': str(e)
            }

    def _prepare_output_dir(self, input_data: Dict[str, Any]) -> str:
        """Prepare and verify output directory."""
        output_dir = input_data.get('output_dir', 'outputs')
        try:
            os.makedirs(output_dir, exist_ok=True)
            if not os.path.isdir(output_dir):
                raise ValueError(f"Output directory is not a directory: {output_dir}")
            if not os.access(output_dir, os.W_OK):
                raise ValueError(f"Output directory is not writable: {output_dir}")
            logger.info(f"Network analysis output directory: {output_dir}")
            return output_dir
        except Exception as e:
            logger.error(f"Failed to create or verify output directory {output_dir}: {e}")
            raise

    def _get_proteins_from_input(self, input_data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Extract proteins dictionary from input data."""
        proteins = input_data.get('proteins', {})
        
        # Handle single protein sequence input
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
        
        if not proteins:
            raise ValueError("No proteins provided for network analysis")
            
        return proteins

    def _process_all_proteins(self, proteins: Dict[str, Dict[str, Any]], output_dir: str) -> Dict[str, Any]:
        """Process all proteins and collect results."""
        all_saved_files = []
        all_network_stats = []
        protein_results = []
        
        for idx, (protein_id, protein_data) in enumerate(proteins.items()):
            try:
                original_id = protein_data.get('original_id', protein_id)
                logger.info(f"Processing protein {idx + 1}/{len(proteins)}: {original_id}")
                
                protein = {
                    'protein_id': protein_id,
                    'sequence': protein_data.get('sequence', ''),
                    'source': protein_data.get('source', 'unknown'),
                    'original_id': original_id
                }
                
                result = self._process_single_protein(protein, output_dir, idx)
                
                if result.get('success'):
                    all_saved_files.extend(result.get('output_files', []))
                    all_network_stats.append(result.get('network_stats', {}))
                    protein_results.append({
                        'protein_id': protein_id,
                        'original_id': original_id,
                        'success': True,
                        'output_files': result.get('output_files', [])
                    })
                else:
                    error_msg = result.get('error_message', 'Unknown error')
                    logger.warning(f"Failed to process protein {original_id}: {error_msg}")
                    protein_results.append({
                        'protein_id': protein_id,
                        'original_id': original_id,
                        'success': False,
                        'error_message': error_msg
                    })
            except Exception as e:
                logger.error(f"Error processing protein {protein_id}: {e}", exc_info=True)
                protein_results.append({
                    'protein_id': protein_id,
                    'original_id': protein_data.get('original_id', protein_id),
                    'success': False,
                    'error_message': str(e)
                })
                
        return {
            'all_saved_files': all_saved_files,
            'all_network_stats': all_network_stats,
            'protein_results': protein_results
        }

    def _process_single_protein(self, protein: Dict[str, Any], output_dir: str, protein_index: int) -> Dict[str, Any]:
        """Process a single protein for network analysis."""
        try:
            protein_id = protein.get('protein_id', f'PROTEIN_{protein_index}')
            protein_sequence = protein.get('sequence', '')
            
            if not protein_sequence:
                raise ValueError(f"No sequence provided for protein {protein_id}")
            
            # 1. Perform Similarity Search
            search_results = self._perform_similarity_search(protein_id, protein_sequence)
            result_entry = search_results[protein_id]
            
            # 2. Get Query Embedding
            query_embedding = self._get_query_embedding(protein_sequence, result_entry)
            
            # 3. Process Hits and Metadata
            hits = result_entry.get('hits', [])
            hits_data = self._process_hits_and_metadata(hits, protein_id, protein.get('source', 'unknown'))
            
            # 4. Prepare Embeddings
            embeddings, protein_ids, all_metadata = self._prepare_embeddings(
                query_embedding, 
                hits_data['top_k_embeddings'], 
                hits_data['query_metadata'], 
                hits_data['top_k_metadata'],
                protein_id
            )
            
            # 5. Visualize and Save
            return self._visualize_and_save(
                embeddings, protein_ids, all_metadata, query_embedding,
                hits_data, protein_sequence, output_dir, protein_index, protein_id
            )
            
        except Exception as e:
            logger.error(f"Error processing single protein: {e}", exc_info=True)
            return {
                'success': False,
                'output_files': [],
                'error_message': str(e),
                'protein_id': protein.get('protein_id', 'unknown')
            }

    def _perform_similarity_search(self, protein_id: str, sequence: str) -> Dict[str, Any]:
        """Perform similarity search using ProteinStorage."""
        proteins_query = {
            protein_id: {
                'id': protein_id,
                'sequence': sequence
            }
        }
        
        results = self.storage.query_similar_proteins(
            proteins=proteins_query,
            max_hits=self.k_neighbors,
            similarity_threshold=self.similarity_threshold,
            return_embeddings=True
        )
        
        # Debug logging for missing query ID investigation
        if protein_id in results:
            hits = results[protein_id].get('hits', [])
            logger.info(f"Search returned {len(hits)} hits for {protein_id}")
            hit_ids = [h.get('id') or h.get('uniprot_id') for h in hits]
            logger.info(f"Hit IDs: {hit_ids}")
            
            # Check if query ID is in hits (it might be filtered out if it's a perfect match and we exclude self, or if the ID format differs)
            # If the user provides a sequence that matches a UniProt ID, we expect that ID to be in the hits.
        else:
            logger.warning(f"No results entry for {protein_id}")
        
        if protein_id not in results:
            raise ValueError(f"No search results returned for protein {protein_id}")
            
        return results

    def _get_query_embedding(self, sequence: str, result_entry: Dict[str, Any]) -> np.ndarray:
        """Get query embedding from results or generate locally."""
        if 'query_embedding' in result_entry:
            logger.info("Using server-provided query embedding")
            return result_entry['query_embedding']
            
        logger.info("Generating local query embedding")
        return self.embedding_generator.generate_embedding(sequence, pooling="mean")

    def _process_hits_and_metadata(self, hits: List[Dict[str, Any]], protein_id: str, source: str) -> Dict[str, Any]:
        """Process hits, fetch metadata and sequences."""
        # Get query metadata
        query_metadata = self._get_query_metadata(protein_id, source)
        
        top_k_metadata = {}
        top_k_sequences = {}
        top_k_scores = {}
        top_k_embeddings = {}
        
        # Collect IDs
        uniprot_ids = [h.get('id') or h.get('uniprot_id') for h in hits if h.get('id') or h.get('uniprot_id')]
        
        # Batch fetch metadata
        if uniprot_ids:
            try:
                metadata_list = self.fetch_metadata(uniprot_ids)
                for meta in metadata_list:
                    if uid := meta.get('Entry'):
                        # Ensure we keep all relevant fields
                        top_k_metadata[uid] = {
                            'Entry': uid,
                            'Protein names': meta.get('Protein names', 'N/A'),
                            'Organism': meta.get('Organism', 'N/A'),
                            'EC number': meta.get('EC number', 'N/A'),
                            'Protein families': meta.get('Protein families', 'N/A'),
                            'Reviewed': meta.get('Reviewed', 'N/A'),
                            'Function [CC]': meta.get('Function [CC]', 'N/A'),
                            'Gene Names': meta.get('Gene Names', 'N/A'),
                            'Length': meta.get('Length', 'N/A')
                        }
            except Exception as e:
                logger.warning(f"Batch metadata fetch failed: {e}")

        for hit in hits:
            uniprot_id = hit.get('id') or hit.get('uniprot_id')
            if not uniprot_id:
                continue
                
            top_k_scores[uniprot_id] = float(hit.get('score') or hit.get('plm_score', 0.0))
            
            if uniprot_id not in top_k_metadata:
                # Try individual fetch if missed in batch
                try:
                    meta_list = self.fetch_metadata([uniprot_id])
                    if meta_list and (meta := meta_list[0]):
                         top_k_metadata[uniprot_id] = {
                            'Entry': uniprot_id,
                            'Protein names': meta.get('Protein names', 'N/A'),
                            'Organism': meta.get('Organism', 'N/A'),
                            'EC number': meta.get('EC number', 'N/A'),
                            'Protein families': meta.get('Protein families', 'N/A'),
                            'Reviewed': meta.get('Reviewed', 'N/A'),
                            'Function [CC]': meta.get('Function [CC]', 'N/A'),
                            'Gene Names': meta.get('Gene Names', 'N/A'),
                            'Length': meta.get('Length', 'N/A')
                        }
                    else:
                        top_k_metadata[uniprot_id] = {'Entry': uniprot_id}
                except Exception as e:
                    logger.warning(f"Individual metadata fetch failed for {uniprot_id}: {e}")
                    top_k_metadata[uniprot_id] = {'Entry': uniprot_id}
            
            # Fetch sequence
            try:
                seq = self.fetch_protein_sequence(uniprot_id)
                top_k_sequences[uniprot_id] = seq if seq else ''
            except Exception as e:
                logger.warning(f"Could not fetch sequence for {uniprot_id}: {e}")
                top_k_sequences[uniprot_id] = ''
            
            # Get embedding
            embedding = hit.get('embedding')
            if embedding is None:
                embedding = self.storage.get_embedding(uniprot_id)
            
            if embedding is not None:
                top_k_embeddings[uniprot_id] = embedding
                
        return {
            'query_metadata': query_metadata,
            'top_k_metadata': top_k_metadata,
            'top_k_sequences': top_k_sequences,
            'top_k_scores': top_k_scores,
            'top_k_embeddings': top_k_embeddings
        }

    def _get_query_metadata(self, protein_id: str, source: str) -> Dict[str, Any]:
        """Get metadata for query protein."""
        if source == 'uniprot' and protein_id:
            try:
                metadata_list = self.fetch_metadata([protein_id])
                if metadata_list:
                    return metadata_list[0]
            except Exception as e:
                logger.warning(f"Could not fetch metadata for query protein {protein_id}: {e}")
        
        return {
            'Entry': protein_id,
            'Protein names': f'Query Protein {protein_id}',
            'Organism': 'N/A',
            'EC number': 'N/A',
            'Protein families': 'N/A',
            'Reviewed': 'N/A'
        }

    def _prepare_embeddings(self, query_embedding: np.ndarray, top_k_embeddings: Dict[str, np.ndarray],
                          query_metadata: Dict[str, Any], top_k_metadata: Dict[str, Any],
                          protein_id: str) -> Tuple[np.ndarray, List[str], Dict[str, Any]]:
        """Stack embeddings and prepare metadata map."""
        query_id_key = query_metadata.get('Entry', protein_id)
        protein_ids = [query_id_key]
        embeddings_list = [query_embedding]
        all_metadata = {query_id_key: query_metadata}
        
        for uniprot_id, embedding in top_k_embeddings.items():
            protein_ids.append(uniprot_id)
            embeddings_list.append(embedding)
            all_metadata[uniprot_id] = top_k_metadata.get(uniprot_id, {'Entry': uniprot_id})
            
        if not embeddings_list:
             raise ValueError(f"No embeddings available for protein {protein_id}")

        # Stack embeddings
        try:
            # Filter by shape matching query embedding
            target_shape = query_embedding.shape
            valid_indices = [i for i, e in enumerate(embeddings_list) if e.shape == target_shape]
            
            if len(valid_indices) < len(embeddings_list):
                logger.warning(f"Filtered {len(embeddings_list) - len(valid_indices)} embeddings due to shape mismatch")
                
            filtered_embeddings = [embeddings_list[i] for i in valid_indices]
            filtered_ids = [protein_ids[i] for i in valid_indices]
            
            embeddings = np.vstack([np.asarray(v, dtype=np.float32) for v in filtered_embeddings])
            return embeddings, filtered_ids, all_metadata
            
        except Exception as e:
            raise ValueError(f"Failed to stack embeddings: {e}")

    def _visualize_and_save(self, embeddings: np.ndarray, protein_ids: List[str], 
                          all_metadata: Dict[str, Any], query_embedding: np.ndarray,
                          hits_data: Dict[str, Any], query_sequence: str,
                          output_dir: str, protein_index: int, protein_id: str) -> Dict[str, Any]:
        """Create visualization and save all results."""
        saved_files = []
        network_stats = {}
        query_id_key = protein_ids[0]
        
        # 1. Visualization
        if HAS_GRAPH_DEPS:
            try:
                visualizer = self._get_visualizer()
                viz_results = visualizer.create_interactive_visualization(
                    embeddings=embeddings,
                    protein_ids=protein_ids,
                    metadata=all_metadata,
                    query_embedding=query_embedding,
                    query_protein_id=query_id_key,
                    query_sequence=query_sequence, # Added to display in sidebar
                    output_dir=output_dir
                )
                
                viz_files = self.save_results(viz_results, query_id_key, output_dir, protein_index)
                saved_files.extend(viz_files)
                network_stats = viz_results.get('network_statistics', {})
            except Exception as e:
                logger.error(f"Visualization failed: {e}")
        else:
            logger.warning("Graph dependencies not available, skipping visualization")
            
        # 2. Top K TSV
        tsv_path = self._create_top_k_matches_tsv(
            query_id_key, query_sequence,
            hits_data['query_metadata'],
            hits_data['top_k_metadata'],
            hits_data['top_k_sequences'],
            hits_data['top_k_scores'],
            output_dir, protein_index
        )
        if tsv_path:
            saved_files.append(tsv_path)
            
        return {
            'success': True,
            'output_files': saved_files,
            'network_stats': network_stats,
            'protein_id': protein_id
        }

    def _get_visualizer(self):
        """Get NetworkVisualizer instance."""
        global NetworkVisualizer
        if NetworkVisualizer is None:
            try:
                from .network_visualizer import NetworkVisualizer
            except ImportError:
                try:
                    from analysis.network_analysis.network_visualizer import NetworkVisualizer
                except ImportError:
                    from network_visualizer import NetworkVisualizer
                    
        return NetworkVisualizer({
            'k_neighbors': self.k_neighbors,
            'similarity_threshold': self.similarity_threshold
        })

    def save_results(self, visualization_results: Dict[str, Any], 
                    query_protein_id: str, output_dir: str, protein_index: int = 0) -> List[str]:
        """Save analysis results and return file paths."""
        saved_files = []
        safe_protein_id = query_protein_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
        
        # Save HTML
        if html_path := visualization_results.get('html_path'):
            if os.path.exists(html_path):
                new_html_path = os.path.join(output_dir, f"network_visualization_{safe_protein_id}_{protein_index}.html")
                if os.path.abspath(html_path) != os.path.abspath(new_html_path):
                    try:
                        import shutil
                        shutil.copy2(html_path, new_html_path)
                        saved_files.append(new_html_path)
                    except Exception as e:
                        logger.warning(f"Could not copy HTML file: {e}")
                        saved_files.append(html_path)
                else:
                    saved_files.append(html_path)
        
        # Save Network Stats & Edges
        if network_graph := visualization_results.get('network_graph'):
            self._save_network_stats(network_graph, query_protein_id, output_dir, safe_protein_id, protein_index, saved_files)
            self._save_network_edges(network_graph, query_protein_id, output_dir, safe_protein_id, protein_index, saved_files)
            
        return saved_files

    def _save_network_stats(self, graph, query_id, output_dir, safe_id, idx, saved_files):
        """Save network statistics to TSV."""
        try:
            stats = [{
                'protein_id': node,
                'degree': graph.degree(node),
                'clustering_coefficient': nx.clustering(graph, node),
                'is_query_protein': (node == query_id)
            } for node in graph.nodes()]
            
            path = os.path.join(output_dir, f"network_statistics_{safe_id}_{idx}.tsv")
            pd.DataFrame(stats).to_csv(path, index=False, sep='\t')
            saved_files.append(path)
        except Exception as e:
            logger.warning(f"Failed to generate statistics TSV: {e}")

    def _save_network_edges(self, graph, query_id, output_dir, safe_id, idx, saved_files):
        """Save network edges to TSV."""
        try:
            edges = [{
                'protein_1': u,
                'protein_2': v,
                'similarity_weight': data.get('weight', 0),
                'edge_type': 'query_connection' if query_id in (u, v) else 'protein_connection'
            } for u, v, data in graph.edges(data=True)]
            
            path = os.path.join(output_dir, f"network_edges_{safe_id}_{idx}.tsv")
            pd.DataFrame(edges).to_csv(path, index=False, sep='\t')
            saved_files.append(path)
        except Exception as e:
            logger.warning(f"Failed to generate edges TSV: {e}")
    
    def _create_top_k_matches_tsv(self, query_protein_id: str, query_sequence: str,
                                  query_metadata: Dict[str, Any],
                                  top_k_metadata: Dict[str, Dict[str, Any]],
                                  top_k_sequences: Dict[str, str],
                                  top_k_similarity_scores: Dict[str, float],
                                  output_dir: str, protein_index: int) -> Optional[str]:
        """Create a TSV file with top k matches and UniProt data for a protein."""
        try:
            safe_protein_id = query_protein_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
            
            rows = [{
                'uniprot_id': query_protein_id,
                'protein_name': query_metadata.get('Protein names', ''),
                'organism': query_metadata.get('Organism', ''),
                'ec_number': query_metadata.get('EC number', ''),
                'protein_families': query_metadata.get('Protein families', ''),
                'reviewed': query_metadata.get('Reviewed', ''),
                'gene_names': query_metadata.get('Gene Names', ''),
                'length': query_metadata.get('Length', ''),
                'similarity_score': 1.0,
                'protein_sequence': query_sequence,
                'is_query': True
            }]
            
            sorted_matches = sorted(top_k_similarity_scores.items(), key=lambda x: x[1], reverse=True)
            
            for uniprot_id, score in sorted_matches:
                # Skip if this hit IS the query protein (prevent duplicate rows)
                if uniprot_id == query_protein_id or uniprot_id == safe_protein_id:
                    continue

                meta = top_k_metadata.get(uniprot_id, {'Entry': uniprot_id})
                rows.append({
                    'uniprot_id': uniprot_id,
                    'protein_name': meta.get('Protein names', ''),
                    'organism': meta.get('Organism', ''),
                    'ec_number': meta.get('EC number', ''),
                    'protein_families': meta.get('Protein families', ''),
                    'reviewed': meta.get('Reviewed', ''),
                    'gene_names': meta.get('Gene Names', ''),
                    'length': meta.get('Length', ''),
                    'similarity_score': score,
                    'protein_sequence': top_k_sequences.get(uniprot_id, ''),
                    'is_query': False
                })
            
            path = os.path.join(output_dir, f"network_analysis_{safe_protein_id}_{protein_index}.tsv")
            pd.DataFrame(rows).to_csv(path, index=False, sep='\t')
            return path
            
        except Exception as e:
            logger.error(f"Failed to create top k matches TSV: {e}", exc_info=True)
            return None


def main():
    """Self-test for network analysis."""
    import shutil
    ok = True
    # Create readable timestamp
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    base_output_dir = os.path.join(os.getcwd(), 'test_local', 'output')
    output_dir = os.path.join(base_output_dir, f"run_{timestamp}")
    
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
        # User requested to keep files for review
        pass
    return 0 if ok else 1


if __name__ == "__main__":
    main()
