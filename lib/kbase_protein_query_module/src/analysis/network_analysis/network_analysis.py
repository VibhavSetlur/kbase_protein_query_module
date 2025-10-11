"""
Network Analysis Stage

This stage provides comprehensive network analysis for protein similarity search results,
including interactive network visualizations, CSV data outputs, and statistics, following 
the BaseStage/StageResult interface used throughout the pipeline.

Features:
- Interactive Plotly network visualization with Kamada-Kawai layout
- CSV outputs for network statistics and protein similarity data
- Network properties and clustering analysis
- Top similar proteins analysis with metadata
- Robust edge selection for connectivity
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union, Any
import logging
from tqdm import tqdm
import os
import json
import time

# Optional imports with fallbacks
try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    nx = None

try:
    from sklearn.metrics.pairwise import cosine_similarity
    from sklearn.cluster import AgglomerativeClustering
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    cosine_similarity = None
    AgglomerativeClustering = None

# Simple result class for compatibility
class StageResult:
    def __init__(self, success: bool, data: Any = None, error_message: str = None, 
                 metadata: Dict[str, Any] = None, artifacts: List[Dict[str, Any]] = None):
        self.success = success
        self.data = data
        self.error_message = error_message
        self.metadata = metadata or {}
        self.artifacts = artifacts or []

logger = logging.getLogger(__name__)

# Visualization imports
try:
    import plotly.graph_objects as go
    import plotly.colors
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    go = None
    logger.warning("Plotly not available. Install with: pip install plotly")


# ==============================================================================
# HELPER FUNCTIONS FOR NETWORK VISUALIZATION
# ==============================================================================

def wrap_text(text, max_chars=80):
    """Wraps text for clean tooltips."""
    if not isinstance(text, str) or pd.isna(text) or text.strip() == '':
        return 'N/A'
    words = text.split()
    lines, current_line, current_length = [], [], 0
    for word in words:
        if current_length + len(word) + 1 <= max_chars:
            current_line.append(word)
            current_length += len(word) + 1
        else:
            lines.append(' '.join(current_line))
            current_line, current_length = [word], len(word)
    if current_line:
        lines.append(' '.join(current_line))
    return '<br>'.join(lines)

def compute_cosine_similarity_matrix(embeddings):
    """Compute cosine similarity matrix for embeddings (shape [N, D])."""
    # Use numpy for basic cosine similarity computation
    # This doesn't require sklearn, but we'll use it if available for consistency
    normed = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
    sim_matrix = np.dot(normed, normed.T)
    return sim_matrix

def build_robust_network_edges(sim_matrix, ids, k_neighbors=5, similarity_threshold=0.1):
    """
    Build a network with robust edge selection to ensure connectivity.
    
    Args:
        sim_matrix: Cosine similarity matrix
        ids: List of protein IDs
        k_neighbors: Number of top neighbors to consider per node
        similarity_threshold: Minimum similarity to include an edge (lowered for connectivity)
    
    Returns:
        NetworkX graph with selected edges
    """
    if not NETWORKX_AVAILABLE:
        raise ImportError("NetworkX is required for network analysis but not available")
    
    N = sim_matrix.shape[0]
    G = nx.Graph() if NETWORKX_AVAILABLE else None
    
    # Add all nodes
    for i in range(N):
        G.add_node(ids[i])
    
    # For each node, find its top k neighbors above threshold
    for i in range(N):
        # Get similarities to all other nodes
        similarities = sim_matrix[i]
        
        # Find indices of nodes with similarity above threshold (excluding self)
        high_sim_indices = []
        for j in range(N):
            if i != j and similarities[j] >= similarity_threshold:
                high_sim_indices.append(j)
        
        # Sort by similarity and take top k
        if high_sim_indices:
            high_sim_indices.sort(key=lambda j: similarities[j], reverse=True)
            top_k_indices = high_sim_indices[:k_neighbors]
            
            # Add edges to top k neighbors
            for j in top_k_indices:
                weight = similarities[j]
                G.add_edge(ids[i], ids[j], weight=weight)
        else:
            # If no neighbors above threshold, connect to top k regardless of threshold
            # This ensures connectivity
            all_indices = [j for j in range(N) if j != i]
            all_indices.sort(key=lambda j: similarities[j], reverse=True)
            top_k_indices = all_indices[:k_neighbors]
            
            for j in top_k_indices:
                weight = similarities[j]
                G.add_edge(ids[i], ids[j], weight=weight)
    
    return G

def create_kamada_kawai_layout(G, seed=42):
    """Create a Kamada-Kawai layout for the graph using similarity weights."""
    if not NETWORKX_AVAILABLE:
        raise ImportError("NetworkX is required for network layout but not available")
    
    # Kamada-Kawai layout is excellent for showing edge relationships clearly
    # It positions nodes based on their graph-theoretic distances
    # Perfect for revealing natural groupings without forcing clusters
    
    # Create a copy of the graph with similarity-based distances
    G_weighted = G.copy()
    
    # Convert similarity to distance: distance = 1 - similarity
    # This way, similar proteins (high similarity) will be closer together
    for u, v, data in G_weighted.edges(data=True):
        similarity = data.get('weight', 0)
        # Convert similarity to distance: distance = 1 - similarity
        # Add small epsilon to avoid zero distance
        distance = max(0.1, 1.0 - similarity)
        G_weighted[u][v]['weight'] = distance
    
    # Use Kamada-Kawai layout with the distance-based weights
    # This layout is excellent for showing edge relationships and natural groupings
    pos = nx.kamada_kawai_layout(G_weighted, weight='weight')
    return pos


class NetworkAnalysis:
    """Network analysis for protein similarity search results."""
    
    def __init__(self, config: Dict[str, Any] = None):
        """Initialize the NetworkAnalysis class."""
        self.config = config or {}
        
        # Check dependencies
        if not NETWORKX_AVAILABLE:
            raise ImportError("NetworkX is required for network analysis but not available")
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn is required for network analysis but not available")
    
    def analyze(self, proteins: List[Any], **kwargs) -> Dict[str, Any]:
        """Main analysis method called by AnalysisManager."""
        try:
            # Convert proteins to the expected format
            input_data = {
                'proteins': proteins,
                'similarity_results': kwargs.get('similarity_results', []),
                'analysis_config': kwargs.get('analysis_config', {})
            }
            
            # Run the analysis
            result = self.run(input_data)
            
            # Convert StageResult to dictionary format expected by AnalysisManager
            return {
                'success': result.success,
                'data': result.data,
                'error_message': result.error_message,
                'metadata': result.metadata,
                'artifacts': result.artifacts
            }
        except Exception as e:
            return {
                'success': False,
                'data': None,
                'error_message': str(e),
                'metadata': {},
                'artifacts': []
            }

    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.k_neighbors = self.config.get('k_neighbors', 8)
        
        # Initialize util components
        self.family_assignment = None
        self.similarity_search = None
        self.embedding_generator = None
        self._initialize_utils()
        self.similarity_threshold = self.config.get('similarity_threshold', 0.1)
        self.mutual_knn = self.config.get('mutual_knn', True)
        self.min_network_size = self.config.get('min_network_size', 5)
        self.max_network_size = self.config.get('max_network_size', 100)
    
    def _initialize_utils(self):
        """Initialize utility components."""
        try:
            from ...util.family_assignment import FamilyAssignment
            from ...util.similarity_search import SimilaritySearch
            from ...util.embeddings import ProteinEmbeddingGenerator
            
            self.family_assignment = FamilyAssignment(self.config)
            self.similarity_search = SimilaritySearch(self.config)
            self.embedding_generator = ProteinEmbeddingGenerator()
            
            logger.info("Network analysis utilities initialized")
            
        except Exception as e:
            logger.warning(f"Could not initialize all utilities: {e}")

    def get_stage_name(self) -> str:
        return "network_analysis"

    def get_required_inputs(self) -> List[str]:
        return ['similarity_results', 'embeddings', 'protein_ids', 'metadata_df', 'query_embedding', 'query_protein_id']

    def get_optional_inputs(self) -> List[str]:
        return ['network_config']

    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        required = self.get_required_inputs()
        for key in required:
            if key not in input_data:
                return False
        return True

    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'network_analysis': {
                'type': 'object',
                'properties': {
                    'query_protein_id': {'type': 'string'},
                    'similarity_table': {'type': 'object'},
                    'top_similar_proteins': {'type': 'array'},
                    'network_visualization': {'type': 'object'},
                    'network_properties': {'type': 'object'},
                    'network_statistics': {'type': 'object'},
                    'clustering_analysis': {'type': 'object'}
                }
            }
        }

    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        try:
            start_time = time.time()
            
            # Check for multi-protein input
            if 'multi_protein_data' in input_data:
                return self._run_multi_protein_analysis(input_data, start_time)
            
            # Single protein analysis
            similar_proteins = input_data['similarity_results']
            embeddings = input_data['embeddings']
            protein_ids = input_data['protein_ids']
            metadata_df = input_data['metadata_df']
            query_embedding = input_data['query_embedding']
            query_protein_id = input_data['query_protein_id']

            # Generate the interactive network visualization and CSV outputs
            results = self.create_interactive_network_visualization(
                embeddings=embeddings,
                protein_ids=protein_ids,
                metadata_df=metadata_df,
                query_embedding=query_embedding,
                query_protein_id=query_protein_id,
                output_dir=input_data.get('output_dir', 'test/outputs')
            )

            execution_time = time.time() - start_time

            return StageResult(
                success=True,
                output_data={'network_analysis': results},
                metadata={
                    'k_neighbors': self.k_neighbors,
                    'similarity_threshold': self.similarity_threshold,
                    'mutual_knn': self.mutual_knn,
                    'execution_time': execution_time
                },
                execution_time=execution_time
            )
        except Exception as e:
            logger.error(f"Network analysis failed: {str(e)}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=0.0,
                error_message=str(e)
            )

    def create_interactive_network_visualization(self,
                                               embeddings: np.ndarray,
                                               protein_ids: list,
                                               metadata_df: pd.DataFrame,
                                               query_embedding: np.ndarray = None,
                                               query_protein_id: str = None,
                                               output_dir: str = 'test/outputs',
                                               id_column: str = 'Entry',
                                               multi_query_proteins: List[str] = None) -> Dict[str, Any]:
        """
        Generates interactive network visualization and saves CSV data outputs.
        Supports both single and multi-protein query analysis.
        
        Args:
            embeddings: Protein embeddings array (N x D)
            protein_ids: List of protein IDs corresponding to embeddings
            metadata_df: DataFrame containing protein metadata
            query_embedding: Query protein embedding vector (1 x D) - for single protein
            query_protein_id: Query protein ID - for single protein
            output_dir: Directory to save outputs
            id_column: Column name for protein IDs in metadata
            multi_query_proteins: List of query protein IDs for multi-protein analysis
            
        Returns:
            Dictionary containing visualization figure and file paths
        """
        logger.info("Creating interactive network visualization...")
        
        # Clean protein_ids to remove None/NaN values
        original_count = len(protein_ids)
        valid_indices = []
        valid_protein_ids = []
        
        for i, protein_id in enumerate(protein_ids):
            if protein_id is not None and str(protein_id).strip() != '' and str(protein_id).lower() != 'nan':
                valid_indices.append(i)
                valid_protein_ids.append(protein_id)
            else:
                logger.warning(f"Skipping invalid protein ID at index {i}: {protein_id}")
        
        if len(valid_protein_ids) != original_count:
            logger.info(f"Filtered {original_count - len(valid_protein_ids)} invalid protein IDs")
            logger.info(f"Using {len(valid_protein_ids)} valid protein IDs")
            
            # Filter embeddings to match valid protein IDs
            embeddings = embeddings[valid_indices]
            protein_ids = valid_protein_ids
        
        # Merge query embedding if provided
        if query_embedding is not None:
            if query_embedding.ndim == 1:
                query_embedding = query_embedding.reshape(1, -1)
            
            # Handle multi-protein queries
            if multi_query_proteins is not None:
                # Multi-protein mode: all specified proteins are queries
                query_proteins = multi_query_proteins
                logger.info(f"Multi-protein mode: {len(query_proteins)} query proteins specified")
                # Use first query protein as primary for file naming
                if query_protein_id is None:
                    query_protein_id = f"multi_query_{len(query_proteins)}_proteins"
            else:
                # Single protein mode
                if query_protein_id is None:
                    query_protein_id = "QUERY_PROTEIN"
                query_proteins = [query_protein_id]
                
                # Check if query protein ID is already in the protein_ids list
                if query_protein_id in protein_ids:
                    logger.info(f"Query protein '{query_protein_id}' is already in the dataset.")
                    embeddings = embeddings  # Keep embeddings as is
                else:
                    # Query protein is new, add it to embeddings and protein_ids
                    logger.info(f"Adding new query protein '{query_protein_id}' to the dataset.")
                    embeddings = np.vstack([embeddings, query_embedding])
                    protein_ids.append(query_protein_id)
        else:
            # No query embedding provided, use the last protein as query
            query_protein_id = protein_ids[-1]
            query_proteins = [query_protein_id]
            logger.info(f"No query embedding provided, using last protein as query: '{query_protein_id}'")

        # Create mapping from protein IDs to metadata
        metadata_dict = metadata_df.set_index(id_column).to_dict('index')
        
        # Build Network
        logger.info("Building protein network...")
        
        # Compute similarity matrix
        sim_matrix = compute_cosine_similarity_matrix(embeddings)
        
        # Build network with robust edge selection
        G = build_robust_network_edges(sim_matrix, protein_ids, self.k_neighbors, self.similarity_threshold)
        
        # Create Layout and Prepare Data
        logger.info("Creating network layout...")
        
        # Create Kamada-Kawai layout
        pos = create_kamada_kawai_layout(G)
        
        # Save CSV outputs
        csv_files = self._save_network_csv_data(G, sim_matrix, protein_ids, metadata_dict, query_protein_id, output_dir, query_proteins)
        
        # Create the interactive visualization
        fig = self._create_plotly_network_figure(G, pos, protein_ids, metadata_dict, query_protein_id, sim_matrix, query_proteins)
        
        # Save the raw network visualization as standalone HTML
        html_filename = f"network_visualization_{query_protein_id}_{int(time.time())}.html"
        html_path = os.path.join(output_dir, html_filename)
        html_path = save_plot_html(fig, html_path)
        
        logger.info(f"Network visualization report saved to: {html_path}")
        logger.info(f"CSV data files saved: {csv_files}")
        
        return {
            'visualization_figure': fig,
            'html_path': html_path,
            'csv_files': csv_files,
            'network_properties': {
                'num_nodes': len(G.nodes()),
                'num_edges': len(G.edges()),
                'density': nx.density(G),
                'connected_components': nx.number_connected_components(G)
            },
            'query_protein_id': query_protein_id,
            'k_neighbors': self.k_neighbors,
            'similarity_threshold': self.similarity_threshold,
            'visualization_type': 'plotly'
        }

    def _save_network_csv_data(self, G: nx.Graph, sim_matrix: np.ndarray, 
                             protein_ids: List[str], metadata_dict: Dict, 
                             query_protein_id: str, output_dir: str,
                             query_proteins: List[str] = None) -> Dict[str, str]:
        """Save network analysis data as CSV files."""
        csv_files = {}
        
        # 1. Network statistics CSV
        network_stats = []
        if query_proteins is None:
            query_proteins = [query_protein_id]
        
        for node in G.nodes():
            node_idx = protein_ids.index(node)
            
            # Calculate max similarity to any query protein
            max_similarity = 0
            connected_to_queries = []
            for qp in query_proteins:
                if qp in protein_ids:
                    q_idx = protein_ids.index(qp)
                    sim = sim_matrix[q_idx, node_idx] if q_idx < len(sim_matrix) and node_idx < len(sim_matrix) else 0
                    max_similarity = max(max_similarity, sim)
                    if G.has_edge(node, qp):
                        connected_to_queries.append(qp)
            
            # Get metadata for this node
            metadata = metadata_dict.get(node, {})
            
            network_stats.append({
                'protein_id': node,
                'protein_name': metadata.get('Protein names', 'N/A'),
                'organism': metadata.get('Organism', 'N/A'),
                'ec_number': metadata.get('EC number', 'N/A'),
                'family': metadata.get('Protein families', 'N/A'),
                'reviewed': metadata.get('Reviewed', 'N/A'),
                'degree': G.degree(node),
                'clustering_coefficient': nx.clustering(G, node),
                'max_similarity_to_queries': max_similarity,
                'is_query_protein': (node in query_proteins),
                'connected_to_queries': len(connected_to_queries),
                'connected_query_proteins': ';'.join(connected_to_queries) if connected_to_queries else 'None'
            })
        
        network_stats_df = pd.DataFrame(network_stats)
        network_stats_path = os.path.join(output_dir, f"network_statistics_{query_protein_id}_{int(time.time())}.csv")
        network_stats_df.to_csv(network_stats_path, index=False)
        csv_files['network_statistics'] = network_stats_path
        
        # 2. Top similar proteins CSV  
        # For multi-protein, use the first query protein as reference
        primary_query = query_proteins[0] if query_proteins else query_protein_id
        if primary_query in protein_ids:
            query_idx = protein_ids.index(primary_query)
            query_similarities = [(protein_ids[i], sim_matrix[query_idx, i]) for i in range(len(protein_ids)) if i != query_idx]
            query_similarities.sort(key=lambda x: x[1], reverse=True)
        else:
            query_similarities = []
        
        top_similar = []
        for rank, (protein_id, similarity) in enumerate(query_similarities[:50], 1):  # Top 50
            metadata = metadata_dict.get(protein_id, {})
            top_similar.append({
                'rank': rank,
                'protein_id': protein_id,
                'similarity_score': similarity,
                'protein_name': metadata.get('Protein names', 'N/A'),
                'organism': metadata.get('Organism', 'N/A'),
                'ec_number': metadata.get('EC number', 'N/A'),
                'family': metadata.get('Protein families', 'N/A'),
                'reviewed': metadata.get('Reviewed', 'N/A'),
                'function': metadata.get('Function [CC]', 'N/A')
            })
        
        top_similar_df = pd.DataFrame(top_similar)
        top_similar_path = os.path.join(output_dir, f"top_similar_proteins_{query_protein_id}_{int(time.time())}.csv")
        top_similar_df.to_csv(top_similar_path, index=False)
        csv_files['top_similar_proteins'] = top_similar_path
        
        # 3. Network edges CSV
        edges_data = []
        for u, v, data in G.edges(data=True):
            edges_data.append({
                'protein_1': u,
                'protein_2': v,
                'similarity_weight': data.get('weight', 0),
                'edge_type': 'query_connection' if (u == query_protein_id or v == query_protein_id) else 'protein_connection'
            })
        
        edges_df = pd.DataFrame(edges_data)
        edges_path = os.path.join(output_dir, f"network_edges_{query_protein_id}_{int(time.time())}.csv")
        edges_df.to_csv(edges_path, index=False)
        csv_files['network_edges'] = edges_path
        
        return csv_files

    def _create_plotly_network_figure(self, G: nx.Graph, pos: Dict, protein_ids: List[str], 
                                    metadata_dict: Dict, query_protein_id: str, 
                                    sim_matrix: np.ndarray, query_proteins: List[str] = None) -> go.Figure:
        """Create the Plotly interactive network figure."""
        # Prepare data for visualization
        node_x, node_y = [], []
        node_customdata = []
        node_colors, node_symbols, node_sizes = [], [], []
        edge_x, edge_y = [], []
        
        # Get connected components for metadata
        connected_components = list(nx.connected_components(G))
        node_to_component = {}
        for i, component in enumerate(connected_components):
            for node in component:
                node_to_component[node] = i
        
        # Add edges to edge lists
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        # Setup query proteins
        if query_proteins is None:
            query_proteins = [query_protein_id]
        
        # Define colors for multiple query proteins
        query_colors = ['gold', 'orange', 'purple', 'green', 'cyan', 'magenta']
        
        for node in G.nodes():
            node_x.append(pos[node][0])
            node_y.append(pos[node][1])
            
            # Check if this is a query protein
            is_query = (node in query_proteins)
            query_index = query_proteins.index(node) if is_query else -1
            
            # Get metadata for this node
            if is_query:
                # Always use custom metadata for the query protein
                protein_metadata = {
                    'Protein names': 'Query Protein',
                    'Organism': 'N/A',
                    'EC number': 'N/A',
                    'Function [CC]': 'Query protein - no metadata available',
                    'Protein families': 'N/A',
                    'Reviewed': 'N/A'
                }
            else:
                # Get metadata from the dictionary, with fallback for missing entries
                protein_metadata = metadata_dict.get(node, {})
                if not protein_metadata:
                    # Fallback metadata for proteins not in the metadata DataFrame
                    protein_metadata = {
                        'Protein names': f'Protein {node}',
                        'Organism': 'N/A',
                        'EC number': 'N/A',
                        'Function [CC]': 'No metadata available',
                        'Protein families': 'N/A',
                        'Reviewed': 'N/A'
                    }
            
            # Get component info
            component_id = node_to_component.get(node, -1)
            component_size = len(connected_components[component_id]) if component_id >= 0 else 1
            
            # Get node degree
            node_degree = G.degree(node)
            
            # Get maximum similarity to any query protein from similarity matrix
            node_idx = protein_ids.index(node)
            similarity_to_query = 0
            for qp in query_proteins:
                if qp in protein_ids:
                    q_idx = protein_ids.index(qp)
                    if q_idx < len(sim_matrix) and node_idx < len(sim_matrix):
                        sim = sim_matrix[q_idx, node_idx]
                        similarity_to_query = max(similarity_to_query, sim)
            
            # Store all necessary info in customdata for hover
            protein_name = str(protein_metadata.get('Protein names', 'N/A')) if protein_metadata.get('Protein names') else 'N/A'
            organism = str(protein_metadata.get('Organism', 'N/A')) if protein_metadata.get('Organism') else 'N/A'
            ec_number = str(protein_metadata.get('EC number', 'N/A')) if protein_metadata.get('EC number') else 'N/A'
            function_text = wrap_text(protein_metadata.get('Function [CC]', '')) if protein_metadata.get('Function [CC]') else 'N/A'
            family = str(protein_metadata.get('Protein families', 'N/A')) if protein_metadata.get('Protein families') else 'N/A'
            reviewed = str(protein_metadata.get('Reviewed', 'N/A')) if protein_metadata.get('Reviewed') else 'N/A'
            
            node_customdata.append([
                str(node),                                    # 0: ID
                protein_name,                                 # 1: Protein Name
                organism,                                     # 2: Organism
                ec_number,                                    # 3: EC Number
                function_text,                                # 4: Function
                family,                                       # 5: Family
                reviewed,                                     # 6: Reviewed
                component_id,                                 # 7: Component ID
                component_size,                               # 8: Component Size
                node_degree,                                  # 9: Node Degree
                is_query,                                     # 10: Is Query Protein
                similarity_to_query                           # 11: Similarity to Query
            ])
            
            # Assign colors and symbols based on query protein status
            if is_query:
                # Use different colors for different query proteins
                color = query_colors[query_index % len(query_colors)]
                node_colors.append(color)
                node_symbols.append('star')  # Star symbol for query proteins
                node_sizes.append(15)  # Larger size for query proteins
            else:
                # Check if this protein is connected to any query protein
                connected_to_any_query = any(G.has_edge(node, qp) for qp in query_proteins if qp in protein_ids)
                
                if connected_to_any_query:
                    node_colors.append('red')  # Red for proteins connected to any query
                    node_symbols.append('circle')
                    node_sizes.append(8)  # Slightly larger for connected proteins
                else:
                    node_colors.append('#CCCCCC')  # Gray for proteins not connected to any query
                    node_symbols.append('circle')
                    node_sizes.append(6)  # Regular size for other proteins

        # Create Figure and Traces
        fig = go.Figure()

        # Edge Trace
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y,
            mode='lines',
            line=dict(width=0.5, color='#888'),
            opacity=0.6,
            hoverinfo='none',
            showlegend=False
        ))

        # Node Trace
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            customdata=node_customdata,
            marker=dict(
                size=node_sizes,
                color=node_colors,
                symbol=node_symbols,
                line=dict(width=2, color='black'),
                opacity=0.9
            ),
            hovertemplate=(
                "<b>Protein ID:</b> %{customdata[0]}<br>"
                "<b>Name:</b> %{customdata[1]}<br>"
                "<b>Organism:</b> %{customdata[2]}<br>"
                "<b>EC Number:</b> %{customdata[3]}<br>"
                "<b>Family:</b> %{customdata[5]}<br>"
                "<b>Reviewed:</b> %{customdata[6]}<br>"
                "<b>Similarity to Query:</b> %{customdata[11]:.3f}<br>"
                "<b>Component:</b> %{customdata[7]} (Size: %{customdata[8]})<br>"
                "<b>Degree:</b> %{customdata[9]}<br>"
                "<b>Type:</b> %{marker.symbol}<br>"
                "<hr><b>Function:</b><br>%{customdata[4]}"
                "<extra></extra>"
            ),
            hoverinfo='all',
            showlegend=False
        ))

        # Add legend for different node types
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=12, symbol='star', color='red'),
            name='Query Protein',
            showlegend=True
        ))
        
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=8, symbol='circle', color='red'),
            name='Proteins Connected to Query',
            showlegend=True
        ))
        
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=6, symbol='circle', color='#CCCCCC'),
            name='Other Proteins',
            showlegend=True
        ))

        # Update layout
        fig.update_layout(
            title=f'<b>Interactive Protein Network (k={self.k_neighbors}, threshold={self.similarity_threshold})</b>',
            height=800,
            plot_bgcolor='white',
            hovermode='closest',
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            margin=dict(l=10, r=10, t=40, b=10),
            legend=dict(
                x=0.02,
                y=0.98,
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='black',
                borderwidth=1
            ),
            modebar=dict(
                orientation='v',
                bgcolor='rgba(255,255,255,0.8)',
                color='black'
            )
        )
        
        # Add instruction annotation for interactive click functionality
        fig.add_annotation(
            text="Click on a node to highlight its edges. Double-click to reset.",
            xref="paper", yref="paper",
            x=0.5, y=-0.05,
            showarrow=False,
            font=dict(size=12, color="gray"),
            align="center"
        )

        return fig

    def analyze_similarity_search_results(self,
                                        similar_proteins: List[Dict],
                                        embeddings: np.ndarray,
                                        protein_ids: List[str],
                                        metadata_df: pd.DataFrame,
                                        query_embedding: np.ndarray,
                                        query_protein_id: str) -> Dict[str, Any]:
        logger.info(f"Analyzing similarity search results for {len(similar_proteins)} proteins")
        results = {
            'query_protein_id': query_protein_id,
            'similarity_table': self._create_similarity_table(similar_proteins, metadata_df),
            'top_similar_proteins': self._get_top_similar_proteins(similar_proteins, 10),
            'network_visualization': None,
            'network_properties': {},
            'network_statistics': None,
            'clustering_analysis': {}
        }

        G = self._build_network_from_similar_proteins(
            similar_proteins, embeddings, protein_ids,
            query_embedding, query_protein_id
        )

        if G is not None:
            results['network_properties'] = self._analyze_network_properties(G)
            results['network_statistics'] = self._get_network_statistics(G)
            results['clustering_analysis'] = self._perform_clustering_analysis(
                similar_proteins, embeddings, protein_ids
            )
            if PLOTLY_AVAILABLE:
                results['network_visualization'] = self._create_network_visualization(
                    G, similar_proteins, metadata_df, query_protein_id
                )

        return results

    def _create_similarity_table(self, similar_proteins: List[Dict],
                                 metadata_df: pd.DataFrame) -> pd.DataFrame:
        table_data = []
        for i, protein in enumerate(similar_proteins):
            protein_id = protein['protein_id']
            similarity_score = protein.get('similarity_score', 0.0)
            if protein_id in metadata_df.index:
                metadata = metadata_df.loc[protein_id]
                protein_name = metadata.get('Protein names', 'N/A')
                organism = metadata.get('Organism', 'N/A')
                ec_number = metadata.get('EC number', 'N/A')
                function_text = metadata.get('Function [CC]', 'N/A')
                family = metadata.get('Protein families', 'N/A')
                reviewed = metadata.get('Reviewed', 'N/A')
            else:
                protein_name = f'Protein {protein_id}'
                organism = 'N/A'
                ec_number = 'N/A'
                function_text = 'No metadata available'
                family = 'N/A'
                reviewed = 'N/A'
            table_data.append({
                'Rank': i + 1,
                'Protein ID': protein_id,
                'Similarity Score': similarity_score,
                'Protein Name': protein_name,
                'Organism': organism,
                'EC Number': ec_number,
                'Family': family,
                'Reviewed': reviewed,
                'Function': function_text
            })
        return pd.DataFrame(table_data)

    def _get_top_similar_proteins(self, similar_proteins: List[Dict],
                                   top_n: int = 10) -> List[Dict]:
        sorted_proteins = sorted(similar_proteins,
                                 key=lambda x: x.get('similarity_score', 0),
                                 reverse=True)
        return sorted_proteins[:top_n]

    def _build_network_from_similar_proteins(self,
                                             similar_proteins: List[Dict],
                                             embeddings: np.ndarray,
                                             protein_ids: List[str],
                                             query_embedding: np.ndarray,
                                             query_protein_id: str) -> Optional[nx.Graph]:
        try:
            logger.info(f"Building network from {len(similar_proteins)} similar proteins...")
            similar_ids = [p['protein_id'] for p in similar_proteins]
            similar_indices = []
            for pid in similar_ids:
                if pid in protein_ids:
                    similar_indices.append(protein_ids.index(pid))
            if not similar_indices:
                logger.warning("No similar proteins found in the full protein list")
                similar_indices = list(range(min(50, len(protein_ids))))
            similar_embeddings = embeddings[similar_indices]
            similar_protein_ids = [protein_ids[i] for i in similar_indices]
            if query_protein_id not in similar_protein_ids:
                similar_embeddings = np.vstack([similar_embeddings, query_embedding])
                similar_protein_ids.append(query_protein_id)
            return self._build_mutual_knn_network(similar_embeddings, similar_protein_ids)
        except Exception as e:
            logger.error(f"Failed to build network: {e}")
            return None

    def _build_mutual_knn_network(self, embeddings: np.ndarray,
                                  protein_ids: List[str]) -> nx.Graph:
        logger.info("Building mutual k-NN network...")
        embeddings_norm = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        similarity_matrix = cosine_similarity(embeddings_norm)
        G = nx.Graph()
        for protein_id in protein_ids:
            G.add_node(protein_id)
        n_nodes = len(protein_ids)
        for i in range(n_nodes):
            similarities = similarity_matrix[i]
            neighbor_indices = np.argsort(similarities)[::-1][1:self.k_neighbors+1]
            for j in neighbor_indices:
                if similarities[j] >= self.similarity_threshold:
                    if self.mutual_knn:
                        j_similarities = similarity_matrix[j]
                        j_neighbors = np.argsort(j_similarities)[::-1][1:self.k_neighbors+1]
                        if i in j_neighbors and j_similarities[i] >= self.similarity_threshold:
                            G.add_edge(protein_ids[i], protein_ids[j],
                                       weight=similarities[j])
                    else:
                        G.add_edge(protein_ids[i], protein_ids[j],
                                   weight=similarities[j])
        if len(G.nodes()) < self.min_network_size:
            logger.warning(f"Network too small ({len(G.nodes())} nodes), adding more edges")
            self._expand_network(G, similarity_matrix, protein_ids)
        if len(G.nodes()) > self.max_network_size:
            logger.warning(f"Network too large ({len(G.nodes())} nodes), trimming")
            self._trim_network(G)
        logger.info(f"Built network with {len(G.nodes())} nodes and {len(G.edges())} edges")
        return G

    def _expand_network(self, G: nx.Graph, similarity_matrix: np.ndarray,
                        protein_ids: List[str]):
        threshold = self.similarity_threshold
        while len(G.nodes()) < self.min_network_size and threshold > 0.01:
            threshold *= 0.9
            for i in range(len(protein_ids)):
                for j in range(i+1, len(protein_ids)):
                    similarity = similarity_matrix[i, j]
                    if similarity >= threshold:
                        G.add_edge(protein_ids[i], protein_ids[j], weight=similarity)

    def _trim_network(self, G: nx.Graph):
        while len(G.nodes()) > self.max_network_size:
            degrees = dict(G.degree())
            min_degree_node = min(degrees, key=degrees.get)
            G.remove_node(min_degree_node)

    def _analyze_network_properties(self, G: nx.Graph) -> Dict:
        properties = {
            'num_nodes': len(G.nodes()),
            'num_edges': len(G.edges()),
            'density': nx.density(G),
            'average_degree': np.mean([d for n, d in G.degree()]),
            'connected_components': [list(component) for component in nx.connected_components(G)],
            'num_components': nx.number_connected_components(G),
            'largest_component_size': len(max(nx.connected_components(G), key=len)),
            'average_clustering': nx.average_clustering(G),
            'average_shortest_path': None
        }
        largest_cc = max(nx.connected_components(G), key=len)
        if len(largest_cc) > 1:
            largest_cc_graph = G.subgraph(largest_cc)
            try:
                properties['average_shortest_path'] = nx.average_shortest_path_length(largest_cc_graph)
            except nx.NetworkXError:
                properties['average_shortest_path'] = None
        return properties

    def _get_network_statistics(self, G: nx.Graph) -> pd.DataFrame:
        stats = []
        for node in G.nodes():
            node_stats = {
                'protein_id': node,
                'degree': G.degree(node),
                'clustering_coefficient': nx.clustering(G, node),
                'betweenness_centrality': nx.betweenness_centrality(G)[node],
                'closeness_centrality': nx.closeness_centrality(G)[node],
                'eigenvector_centrality': nx.eigenvector_centrality_numpy(G)[node]
            }
            stats.append(node_stats)
        return pd.DataFrame(stats)

    def _perform_clustering_analysis(self, similar_proteins: List[Dict],
                                     embeddings: np.ndarray,
                                     protein_ids: List[str]) -> Dict:
        try:
            similar_ids = [p['protein_id'] for p in similar_proteins]
            similar_indices = []
            for pid in similar_ids:
                if pid in protein_ids:
                    similar_indices.append(protein_ids.index(pid))
            if len(similar_indices) < 2:
                return {'clusters': [], 'cluster_assignments': {}}
            similar_embeddings = embeddings[similar_indices]
            similar_protein_ids = [protein_ids[i] for i in similar_indices]
            n_clusters = min(5, len(similar_indices))
            clustering = AgglomerativeClustering(n_clusters=n_clusters)
            cluster_labels = clustering.fit_predict(similar_embeddings)
            clusters = {}
            cluster_assignments = {}
            for i, label in enumerate(cluster_labels):
                protein_id = similar_protein_ids[i]
                if label not in clusters:
                    clusters[label] = []
                clusters[label].append(protein_id)
                cluster_assignments[protein_id] = label
            return {
                'clusters': clusters,
                'cluster_assignments': cluster_assignments,
                'n_clusters': n_clusters
            }
        except Exception as e:
            logger.error(f"Failed to perform clustering analysis: {e}")
            return {'clusters': [], 'cluster_assignments': {}}

    def _create_network_visualization(self, G: nx.Graph,
                                      similar_proteins: List[Dict],
                                      metadata_df: pd.DataFrame,
                                      query_protein_id: str) -> Optional[object]:
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly not available for network visualization")
            return None
        try:
            pos = nx.kamada_kawai_layout(G)
            node_x, node_y = [], []
            node_customdata = []
            node_colors, node_symbols, node_sizes = [], [], []
            edge_x, edge_y = [], []
            connected_components = list(nx.connected_components(G))
            node_to_component = {}
            for i, component in enumerate(connected_components):
                for node in component:
                    node_to_component[node] = i
            for edge in G.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
            for node in G.nodes():
                node_x.append(pos[node][0])
                node_y.append(pos[node][1])
                is_query = (node == query_protein_id)
                if is_query:
                    protein_metadata = {
                        'Protein names': 'Query Protein',
                        'Organism': 'N/A',
                        'EC number': 'N/A',
                        'Function [CC]': 'Query protein - no metadata available',
                        'Protein families': 'N/A',
                        'Reviewed': 'N/A'
                    }
                else:
                    protein_metadata = metadata_df.loc[node].to_dict() if node in metadata_df.index else {
                        'Protein names': f'Protein {node}',
                        'Organism': 'N/A',
                        'EC number': 'N/A',
                        'Function [CC]': 'No metadata available',
                        'Protein families': 'N/A',
                        'Reviewed': 'N/A'
                    }
                component_id = node_to_component.get(node, -1)
                component_size = len(connected_components[component_id]) if component_id >= 0 else 1
                node_degree = G.degree(node)
                protein_name = str(protein_metadata.get('Protein names', 'N/A'))
                organism = str(protein_metadata.get('Organism', 'N/A'))
                ec_number = str(protein_metadata.get('EC number', 'N/A'))
                function_text = str(protein_metadata.get('Function [CC]', 'N/A'))
                family = str(protein_metadata.get('Protein families', 'N/A'))
                reviewed = str(protein_metadata.get('Reviewed', 'N/A'))
                node_customdata.append([
                    str(node), protein_name, organism, ec_number, function_text,
                    family, reviewed, component_id, component_size, node_degree, is_query
                ])
                if is_query:
                    node_colors.append('red')
                    node_symbols.append('star')
                    node_sizes.append(12)
                else:
                    node_colors.append('#CCCCCC')
                    node_symbols.append('circle')
                    node_sizes.append(6)
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=edge_x, y=edge_y,
                mode='lines',
                line=dict(width=0.5, color='#888'),
                opacity=0.6,
                hoverinfo='none',
                showlegend=False
            ))
            fig.add_trace(go.Scatter(
                x=node_x, y=node_y,
                mode='markers',
                customdata=node_customdata,
                marker=dict(
                    size=node_sizes,
                    color=node_colors,
                    symbol=node_symbols,
                    line=dict(width=2, color='black'),
                    opacity=0.9
                ),
                hovertemplate=(
                    "<b>Protein ID:</b> %{customdata[0]}<br>"
                    "<b>Name:</b> %{customdata[1]}<br>"
                    "<b>Organism:</b> %{customdata[2]}<br>"
                    "<b>EC Number:</b> %{customdata[3]}<br>"
                    "<b>Family:</b> %{customdata[5]}<br>"
                    "<b>Reviewed:</b> %{customdata[6]}<br>"
                    "<b>Component:</b> %{customdata[7]} (Size: %{customdata[8]})<br>"
                    "<b>Degree:</b> %{customdata[9]}<br>"
                    "<b>Type:</b> %{marker.symbol}<br>"
                    "<hr><b>Function:</b><br>%{customdata[4]}"
                    "<extra></extra>"
                ),
                hoverinfo='all',
                showlegend=False
            ))
            fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                                     marker=dict(size=12, symbol='star', color='red'),
                                     name='Query Protein', showlegend=True))
            fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                                     marker=dict(size=6, symbol='circle', color='#CCCCCC'),
                                     name='Similar Proteins', showlegend=True))
            fig.update_layout(
                title=f'<b>Protein Similarity Network (k={self.k_neighbors}, threshold={self.similarity_threshold})</b>',
                height=600,
                plot_bgcolor='white',
                hovermode='closest',
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                margin=dict(l=10, r=10, t=40, b=10),
                legend=dict(x=0.02, y=0.98, bgcolor='rgba(255,255,255,0.8)', bordercolor='black', borderwidth=1)
            )
            return fig
        except Exception as e:
            logger.error(f"Failed to create network visualization: {e}")
            return None
    
    def _run_multi_protein_analysis(self, input_data: Dict[str, Any], start_time: float) -> StageResult:
        """Run network analysis for multiple query proteins."""
        multi_protein_data = input_data['multi_protein_data']
        output_dir = input_data.get('output_dir', 'test/outputs')
        
        logger.info(f"Running multi-protein network analysis for {len(multi_protein_data)} proteins")
        
        all_results = {}
        
        # Process each query protein individually  
        for query_protein_id, protein_data in multi_protein_data.items():
            logger.info(f"Processing network for query protein: {query_protein_id}")
            
            individual_results = self.create_interactive_network_visualization(
                embeddings=protein_data['embeddings'],
                protein_ids=protein_data['protein_ids'],
                metadata_df=protein_data['metadata_df'],
                query_embedding=protein_data['query_embedding'],
                query_protein_id=query_protein_id,
                output_dir=os.path.join(output_dir, f"individual_networks")
            )
            
            all_results[query_protein_id] = individual_results
        
        # Generate multi-protein CSV data
        multi_csv_files = self._generate_multi_protein_csv_data(all_results, output_dir)
        
        execution_time = time.time() - start_time
        
        return StageResult(
            success=True,
            output_data={
                'multi_protein_network_analysis': {
                    'individual_networks': all_results,
                    'csv_data': multi_csv_files,
                    'total_query_proteins': len(multi_protein_data)
                }
            },
            metadata={
                'analysis_type': 'multi_protein_network',
                'query_proteins': list(multi_protein_data.keys()),
                'execution_time': execution_time
            }
        )
    
    def _generate_multi_protein_csv_data(self, all_results: Dict[str, Any], output_dir: str) -> Dict[str, str]:
        """Generate comprehensive CSV data for multi-protein analysis."""
        csv_files = {}
        os.makedirs(output_dir, exist_ok=True)
        
        # Multi-protein summary CSV
        summary_data = []
        for query_protein_id, results in all_results.items():
            network_props = results.get('network_properties', {})
            summary_data.append({
                'query_protein_id': query_protein_id,
                'network_nodes': network_props.get('num_nodes', 0),
                'network_edges': network_props.get('num_edges', 0),
                'network_density': network_props.get('density', 0),
                'connected_components': network_props.get('connected_components', 0),
                'visualization_file': os.path.basename(results.get('html_path', ''))
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_path = os.path.join(output_dir, f"multi_protein_summary_{int(time.time())}.csv")
        summary_df.to_csv(summary_path, index=False)
        csv_files['multi_protein_summary'] = summary_path
        
        return csv_files


