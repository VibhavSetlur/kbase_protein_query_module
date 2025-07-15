"""
Protein Network Builder Module

This module handles the construction and visualization of localized protein networks using
mutual k-nearest neighbors and other network construction methods.
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import List, Dict, Tuple, Optional, Union
import logging
from tqdm import tqdm
import os
from sklearn.metrics.pairwise import cosine_similarity

logger = logging.getLogger(__name__)

# Visualization imports
try:
    import plotly.graph_objects as go
    import plotly.colors
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly not available. Install with: pip install plotly")




# ==============================================================================
# 1. HELPER FUNCTIONS (for layout and text wrapping)
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
    N = sim_matrix.shape[0]
    G = nx.Graph()
    
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


# ==============================================================================
# 2. MAIN VISUALIZATION FUNCTION
# ==============================================================================

def visualize_interactive_protein_network(
    embeddings: np.ndarray,
    protein_ids: list,
    metadata_df: pd.DataFrame,
    k_neighbors: int = 8,
    similarity_threshold: float = 0.1,
    id_column: str = 'Entry',
    query_embedding: np.ndarray = None,
    query_protein_id: str = None,
    output_file: Optional[str] = None
):
    """
    Generates a highly interactive visualization of protein network based on embeddings.

    Args:
        embeddings: Protein embeddings array (N x D)
        protein_ids: List of protein IDs corresponding to embeddings
        metadata_df: DataFrame containing protein metadata
        k_neighbors: Number of top neighbors to connect per node
        similarity_threshold: Minimum similarity to include an edge
        id_column: Column name for protein IDs in metadata
        query_embedding: Query protein embedding vector (1 x D). If provided, will be merged with embeddings.
        query_protein_id: Query protein ID
        output_file: Optional path to save the visualization as HTML

    Features:
        - Shows all proteins in the dataset
        - Rich hover-tooltips with detailed protein information
        - Robust edge selection to avoid dense connections
        - Interactive zoom, pan, and hover
        - Color-coded nodes based on connectivity
    """
    
    if not PLOTLY_AVAILABLE:
        logger.error("Plotly not available. Cannot create interactive visualization.")
        return None, None
    
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
        
        # If query_protein_id is not provided, use a default name
        if query_protein_id is None:
            query_protein_id = "QUERY_PROTEIN"
        
        # Check if query protein ID is already in the protein_ids list
        if query_protein_id in protein_ids:
            # Query protein is already in the dataset, don't add it again
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
        logger.info(f"No query embedding provided, using last protein as query: '{query_protein_id}'")
    
    # Warn if query protein ID is in metadata but not in protein_ids (shouldn't happen now)
    # Check if the id_column exists in metadata_df.columns, otherwise use index
    if id_column is not None and id_column in metadata_df.columns:
        if query_protein_id in set(metadata_df[id_column]) and query_protein_id not in protein_ids:
            logger.warning(f"Query protein ID '{query_protein_id}' is present in metadata but not in protein_ids.")
        # Create mapping from protein IDs to metadata using the specified column
        metadata_dict = metadata_df.set_index(id_column).to_dict('index')
    else:
        # If id_column doesn't exist or is None, assume protein IDs are in the index
        if query_protein_id in metadata_df.index and query_protein_id not in protein_ids:
            logger.warning(f"Query protein ID '{query_protein_id}' is present in metadata but not in protein_ids.")
        # Create mapping from protein IDs to metadata using the index
        metadata_dict = metadata_df.to_dict('index')
    
    # --- 3. Build Network ---
    logger.info("Building protein network...")
    
    # Compute similarity matrix
    sim_matrix = compute_cosine_similarity_matrix(embeddings)
    
    # Build network with robust edge selection
    G = build_robust_network_edges(sim_matrix, protein_ids, k_neighbors, similarity_threshold)
    
    # --- 4. Create Layout and Prepare Data ---
    logger.info("Creating network layout...")
    
    # Create Kamada-Kawai layout
    pos = create_kamada_kawai_layout(G)

    # If output_file is specified, ensure directory exists (if any)
    if output_file:
        output_dir = os.path.dirname(output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
    
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
    
    # Add nodes to node lists
    for node in G.nodes():
        node_x.append(pos[node][0])
        node_y.append(pos[node][1])
        
        # Check if this is the query protein
        is_query = (node == query_protein_id)
        
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
        
        # Get similarity to query protein from similarity matrix
        query_idx = protein_ids.index(query_protein_id)
        node_idx = protein_ids.index(node)
        similarity_to_query = sim_matrix[query_idx, node_idx]
        
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
        
        # Assign colors and symbols based on connection to query protein
        if is_query:
            node_colors.append('red')  # Red for query protein
            node_symbols.append('star')  # Star symbol for query protein
            node_sizes.append(12)  # Larger size for query protein
        else:
            # Check if this protein is connected to the query protein
            is_connected_to_query = G.has_edge(node, query_protein_id)
            
            if is_connected_to_query:
                node_colors.append('red')  # Bright red for proteins connected to query
                node_symbols.append('circle')
                node_sizes.append(8)  # Slightly larger for connected proteins
            else:
                node_colors.append('#CCCCCC')  # Gray for proteins not connected to query
                node_symbols.append('circle')
                node_sizes.append(6)  # Regular size for other proteins

    # --- 5. Create Figure and Traces ---
    logger.info("Creating interactive visualization...")
    
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
        title=f'<b>Interactive Protein Network (k={k_neighbors}, threshold={similarity_threshold})</b>',
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

    # Save to file if specified
    if output_file:
        output_dir = os.path.dirname(output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        fig.write_html(output_file)
        logger.info(f"Interactive visualization saved to: {output_file}")

    # --- 6. Display the Figure ---
    logger.info("Interactive protein network created.")
    logger.info("Features:")
    logger.info("- Query protein highlighted as red star (larger size)")
    logger.info("- Proteins connected to query shown in red")
    logger.info("- Other proteins shown in gray")
    logger.info("- Kamada-Kawai layout for optimal edge visibility")
    logger.info("- Natural groupings revealed without forced clustering")
    logger.info("- Clear edge relationships and protein connections")
    logger.info("- Click on nodes to highlight their edges")
    logger.info("- Double-click to reset edge highlighting")
    logger.info("- Hover over nodes to see detailed protein information")
    logger.info("- Edges show protein similarities")
    logger.info("- Interactive zoom, pan, and hover")
    logger.info("- Legend shows node types")
    
    return fig, G


class DynamicNetworkBuilder:
    """
    Builds dynamic, localized protein networks from similarity search results.
    
    This class constructs targeted networks around query proteins using
    various network construction methods including mutual k-nearest neighbors.
    """
    
    def __init__(self, k_neighbors: int = 8, 
                 similarity_threshold: float = 0.1,
                 mutual_knn: bool = True,
                 min_network_size: int = 5,
                 max_network_size: int = 100):
        """
        Initialize the network builder.
        
        Args:
            k_neighbors: Number of neighbors for k-NN construction
            similarity_threshold: Minimum similarity to include an edge
            mutual_knn: Whether to use mutual k-nearest neighbors
            min_network_size: Minimum number of nodes in network
            max_network_size: Maximum number of nodes in network
        """
        self.k_neighbors = k_neighbors
        self.similarity_threshold = similarity_threshold
        self.mutual_knn = mutual_knn
        self.min_network_size = min_network_size
        self.max_network_size = max_network_size
    
    def build_mutual_knn_network(self, embeddings: np.ndarray,
                               protein_ids: List[str],
                               query_embedding: Optional[np.ndarray] = None,
                               query_protein_id: Optional[str] = None) -> nx.Graph:
        """
        Build a mutual k-nearest neighbors network.
        
        Args:
            embeddings: Protein embeddings array
            protein_ids: List of protein IDs
            query_embedding: Optional query protein embedding
            query_protein_id: Optional query protein ID
            
        Returns:
            NetworkX graph
        """
        logger.info("Building mutual k-NN network...")
        
        # Add query protein if provided
        if query_embedding is not None:
            embeddings = np.vstack([embeddings, query_embedding])
            protein_ids = protein_ids + [query_protein_id or "QUERY_PROTEIN"]
        
        # Normalize embeddings
        embeddings_norm = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        
        # Compute similarity matrix
        similarity_matrix = cosine_similarity(embeddings_norm)
        
        # Build mutual k-NN graph
        G = nx.Graph()
        
        # Add all nodes
        for protein_id in protein_ids:
            G.add_node(protein_id)
        
        # Find k-nearest neighbors for each node
        n_nodes = len(protein_ids)
        for i in range(n_nodes):
            # Get similarities to all other nodes
            similarities = similarity_matrix[i]
            
            # Find top k neighbors (excluding self)
            neighbor_indices = np.argsort(similarities)[::-1][1:self.k_neighbors+1]
            
            for j in neighbor_indices:
                if similarities[j] >= self.similarity_threshold:
                    # Check if this is a mutual edge
                    if self.mutual_knn:
                        # Check if j also has i as one of its top k neighbors
                        j_similarities = similarity_matrix[j]
                        j_neighbors = np.argsort(j_similarities)[::-1][1:self.k_neighbors+1]
                        
                        if i in j_neighbors and j_similarities[i] >= self.similarity_threshold:
                            G.add_edge(protein_ids[i], protein_ids[j], 
                                     weight=similarities[j])
                    else:
                        # Add directed edge
                        G.add_edge(protein_ids[i], protein_ids[j], 
                                 weight=similarities[j])
        
        # Ensure network size constraints
        if len(G.nodes()) < self.min_network_size:
            logger.warning(f"Network too small ({len(G.nodes())} nodes), adding more edges")
            self._expand_network(G, similarity_matrix, protein_ids)
        
        if len(G.nodes()) > self.max_network_size:
            logger.warning(f"Network too large ({len(G.nodes())} nodes), trimming")
            self._trim_network(G, query_protein_id)
        
        logger.info(f"Built network with {len(G.nodes())} nodes and {len(G.edges())} edges")
        return G
    
    def build_threshold_network(self, embeddings: np.ndarray,
                              protein_ids: List[str],
                              query_embedding: Optional[np.ndarray] = None,
                              query_protein_id: Optional[str] = None) -> nx.Graph:
        """
        Build a network using similarity threshold.
        
        Args:
            embeddings: Protein embeddings array
            protein_ids: List of protein IDs
            query_embedding: Optional query protein embedding
            query_protein_id: Optional query protein ID
            
        Returns:
            NetworkX graph
        """
        logger.info("Building threshold-based network...")
        
        # Add query protein if provided
        if query_embedding is not None:
            embeddings = np.vstack([embeddings, query_embedding])
            protein_ids = protein_ids + [query_protein_id or "QUERY_PROTEIN"]
        
        # Normalize embeddings
        embeddings_norm = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        
        # Compute similarity matrix
        similarity_matrix = cosine_similarity(embeddings_norm)
        
        # Build graph
        G = nx.Graph()
        
        # Add all nodes
        for protein_id in protein_ids:
            G.add_node(protein_id)
        
        # Add edges above threshold
        n_nodes = len(protein_ids)
        for i in range(n_nodes):
            for j in range(i+1, n_nodes):
                similarity = similarity_matrix[i, j]
                if similarity >= self.similarity_threshold:
                    G.add_edge(protein_ids[i], protein_ids[j], weight=similarity)
        
        logger.info(f"Built network with {len(G.nodes())} nodes and {len(G.edges())} edges")
        return G
    
    def build_hybrid_network(self, embeddings: np.ndarray,
                           protein_ids: List[str],
                           query_embedding: Optional[np.ndarray] = None,
                           query_protein_id: Optional[str] = None) -> nx.Graph:
        """
        Build a hybrid network combining k-NN and threshold methods.
        
        Args:
            embeddings: Protein embeddings array
            protein_ids: List of protein IDs
            query_embedding: Optional query protein embedding
            query_protein_id: Optional query protein ID
            
        Returns:
            NetworkX graph
        """
        logger.info("Building hybrid network...")
        
        # Add query protein if provided
        if query_embedding is not None:
            embeddings = np.vstack([embeddings, query_embedding])
            protein_ids = protein_ids + [query_protein_id or "QUERY_PROTEIN"]
        
        # Normalize embeddings
        embeddings_norm = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        
        # Compute similarity matrix
        similarity_matrix = cosine_similarity(embeddings_norm)
        
        # Build graph
        G = nx.Graph()
        
        # Add all nodes
        for protein_id in protein_ids:
            G.add_node(protein_id)
        
        # Add edges using hybrid approach
        n_nodes = len(protein_ids)
        for i in range(n_nodes):
            # Get similarities to all other nodes
            similarities = similarity_matrix[i]
            
            # Find top k neighbors
            neighbor_indices = np.argsort(similarities)[::-1][1:self.k_neighbors+1]
            
            for j in neighbor_indices:
                similarity = similarities[j]
                
                # Add edge if it meets threshold OR is in top k
                if similarity >= self.similarity_threshold:
                    G.add_edge(protein_ids[i], protein_ids[j], weight=similarity)
        
        logger.info(f"Built hybrid network with {len(G.nodes())} nodes and {len(G.edges())} edges")
        return G
    
    def _expand_network(self, G: nx.Graph, similarity_matrix: np.ndarray, 
                       protein_ids: List[str]):
        """Expand network by adding more edges if too small."""
        # Add edges with lower threshold until minimum size is reached
        threshold = self.similarity_threshold
        while len(G.nodes()) < self.min_network_size and threshold > 0.01:
            threshold *= 0.9
            
            for i in range(len(protein_ids)):
                for j in range(i+1, len(protein_ids)):
                    similarity = similarity_matrix[i, j]
                    if similarity >= threshold:
                        G.add_edge(protein_ids[i], protein_ids[j], weight=similarity)
    
    def _trim_network(self, G: nx.Graph, query_protein_id: Optional[str] = None):
        """Trim network to maximum size while preserving query connectivity."""
        if query_protein_id is None:
            # Remove nodes with lowest degree
            while len(G.nodes()) > self.max_network_size:
                degrees = dict(G.degree())
                min_degree_node = min(degrees, key=degrees.get)
                G.remove_node(min_degree_node)
        else:
            # Keep query protein and its neighbors, remove others by degree
            query_neighbors = set(G.neighbors(query_protein_id))
            query_neighbors.add(query_protein_id)
            
            # Remove nodes outside query neighborhood
            nodes_to_remove = [n for n in G.nodes() if n not in query_neighbors]
            
            # Sort by degree (remove lowest degree first)
            nodes_to_remove.sort(key=lambda n: G.degree(n))
            
            # Remove nodes until we reach max size
            for node in nodes_to_remove:
                if len(G.nodes()) <= self.max_network_size:
                    break
                G.remove_node(node)
    
    def build_network_from_similar_proteins(self, 
                                          similar_proteins: List[Dict],
                                          embeddings: np.ndarray,
                                          protein_ids: List[str],
                                          query_embedding: np.ndarray,
                                          query_protein_id: str,
                                          method: str = "mutual_knn") -> nx.Graph:
        """
        Build network from similarity search results.
        
        Args:
            similar_proteins: List of similar proteins from search
            embeddings: Full embeddings array
            protein_ids: Full list of protein IDs
            query_embedding: Query protein embedding
            query_protein_id: Query protein ID
            method: Network construction method
            
        Returns:
            NetworkX graph
        """
        logger.info(f"Building network from {len(similar_proteins)} similar proteins...")
        
        # Extract embeddings and IDs for similar proteins
        similar_ids = [p['protein_id'] for p in similar_proteins]
        similar_indices = [protein_ids.index(pid) for pid in similar_ids if pid in protein_ids]
        similar_embeddings = embeddings[similar_indices]
        
        # Build network using specified method
        if method == "mutual_knn":
            return self.build_mutual_knn_network(
                similar_embeddings, similar_ids, query_embedding, query_protein_id
            )
        elif method == "threshold":
            return self.build_threshold_network(
                similar_embeddings, similar_ids, query_embedding, query_protein_id
            )
        elif method == "hybrid":
            return self.build_hybrid_network(
                similar_embeddings, similar_ids, query_embedding, query_protein_id
            )
        else:
            raise ValueError(f"Unknown network construction method: {method}")
    
    def create_interactive_visualization(self, 
                                       embeddings: np.ndarray,
                                       protein_ids: List[str],
                                       metadata_df: pd.DataFrame,
                                       query_embedding: Optional[np.ndarray] = None,
                                       query_protein_id: Optional[str] = None,
                                       output_file: Optional[str] = None) -> Tuple[Optional[object], Optional[nx.Graph]]:
        """
        Create an interactive visualization of the protein network.
        
        Args:
            embeddings: Protein embeddings array
            protein_ids: List of protein IDs
            metadata_df: DataFrame containing protein metadata
            query_embedding: Optional query protein embedding
            query_protein_id: Optional query protein ID
            output_file: Optional path to save the visualization as HTML
            
        Returns:
            Tuple of (Plotly figure, NetworkX graph)
        """
        return visualize_interactive_protein_network(
            embeddings=embeddings,
            protein_ids=protein_ids,
            metadata_df=metadata_df,
            k_neighbors=self.k_neighbors,
            similarity_threshold=self.similarity_threshold,
            id_column=None,  # Use index instead of column name
            query_embedding=query_embedding,
            query_protein_id=query_protein_id,
            output_file=output_file
        )
    
    def analyze_network_properties(self, G: nx.Graph) -> Dict:
        """
        Analyze network properties.
        
        Args:
            G: NetworkX graph
            
        Returns:
            Dictionary of network properties
        """
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
        
        # Calculate average shortest path for largest component
        largest_cc = max(nx.connected_components(G), key=len)
        if len(largest_cc) > 1:
            largest_cc_graph = G.subgraph(largest_cc)
            try:
                properties['average_shortest_path'] = nx.average_shortest_path_length(largest_cc_graph)
            except nx.NetworkXError:
                properties['average_shortest_path'] = None
        
        return properties
    
    def get_network_statistics(self, G: nx.Graph) -> pd.DataFrame:
        """
        Get detailed network statistics.
        
        Args:
            G: NetworkX graph
            
        Returns:
            DataFrame with node-level statistics
        """
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


def create_localized_network(query_embedding: np.ndarray,
                           query_protein_id: str,
                           similar_proteins: List[Dict],
                           embeddings: np.ndarray,
                           protein_ids: List[str],
                           method: str = "mutual_knn",
                           **kwargs) -> Tuple[nx.Graph, Dict]:
    """
    Create a localized network around a query protein.
    
    Args:
        query_embedding: Query protein embedding
        query_protein_id: Query protein ID
        similar_proteins: List of similar proteins from search
        embeddings: Full embeddings array
        protein_ids: Full list of protein IDs
        method: Network construction method
        **kwargs: Additional parameters for network builder
        
    Returns:
        Tuple of (NetworkX graph, network properties)
    """
    # Create network builder
    builder = DynamicNetworkBuilder(**kwargs)
    
    # Build network
    G = builder.build_network_from_similar_proteins(
        similar_proteins, embeddings, protein_ids, 
        query_embedding, query_protein_id, method
    )
    
    # Analyze network properties
    properties = builder.analyze_network_properties(G)
    
    return G, properties


def save_network(G: nx.Graph, output_file: str, properties: Optional[Dict] = None):
    """
    Save network to file.
    
    Args:
        G: NetworkX graph
        output_file: Path to output file
        properties: Optional network properties
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save network as GraphML
    nx.write_graphml(G, output_file)
    
    # Save properties if provided
    if properties is not None:
        properties_file = output_file.replace('.graphml', '_properties.json')
        import json
        with open(properties_file, 'w') as f:
            json.dump(properties, f, indent=2)
    
    logger.info(f"Saved network to {output_file}")


def load_network(input_file: str) -> Tuple[nx.Graph, Optional[Dict]]:
    """
    Load network from file.
    
    Args:
        input_file: Path to input file
        
    Returns:
        Tuple of (NetworkX graph, network properties)
    """
    # Load network
    G = nx.read_graphml(input_file)
    
    # Load properties if available
    properties = None
    properties_file = input_file.replace('.graphml', '_properties.json')
    if os.path.exists(properties_file):
        import json
        with open(properties_file, 'r') as f:
            properties = json.load(f)
    
    return G, properties 