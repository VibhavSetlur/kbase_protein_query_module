"""
Network Visualization Module

This module provides interactive network visualization capabilities for protein similarity networks.
It creates Plotly-based interactive visualizations with proper layout algorithms and styling.
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import List, Dict, Tuple, Optional, Union, Any
import logging
import os
import time

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

# NetworkX import
try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    nx = None
    logger.error("NetworkX is required for network visualization but not available")


class NetworkVisualizer:
    """
    Creates interactive network visualizations for protein similarity networks.
    
    This class handles the creation of Plotly-based interactive visualizations
    with proper layout algorithms, styling, and metadata integration.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        Initialize the Network Visualizer.
        
        Args:
            config: Configuration dictionary for visualization parameters
        """
        self.config = config or {}
        self.k_neighbors = self.config.get('k_neighbors', 8)
        self.similarity_threshold = self.config.get('similarity_threshold', 0.1)
        
        if not PLOTLY_AVAILABLE:
            raise ImportError("Plotly is required for network visualization but not available")
        if not NETWORKX_AVAILABLE:
            raise ImportError("NetworkX is required for network visualization but not available")
    
    def create_interactive_visualization(self,
                                       embeddings: np.ndarray,
                                       protein_ids: List[str],
                                       metadata_df: pd.DataFrame,
                                       query_embedding: Optional[np.ndarray] = None,
                                       query_protein_id: Optional[str] = None,
                                       output_dir: str = 'outputs',
                                       id_column: str = 'Entry') -> Dict[str, Any]:
        """
        Create an interactive network visualization.
        
        Args:
            embeddings: Protein embeddings array (N x D)
            protein_ids: List of protein IDs corresponding to embeddings
            metadata_df: DataFrame containing protein metadata
            query_embedding: Optional query protein embedding vector
            query_protein_id: Optional query protein ID
            output_dir: Directory to save the visualization
            id_column: Column name for protein IDs in metadata
            
        Returns:
            Dictionary containing visualization figure and metadata
        """
        logger.info("Creating interactive network visualization...")
        
        # Clean and validate protein IDs
        embeddings, protein_ids = self._clean_protein_ids(embeddings, protein_ids)
        
        # Merge query embedding if provided
        if query_embedding is not None:
            embeddings, protein_ids, query_protein_id = self._merge_query_protein(
                embeddings, protein_ids, query_embedding, query_protein_id
            )
        else:
            query_protein_id = protein_ids[-1] if protein_ids else "unknown"
        
        # Create metadata mapping
        metadata_dict = self._create_metadata_mapping(metadata_df, id_column)
        
        # Build network
        G = self._build_network(embeddings, protein_ids)
        
        # Create layout
        pos = self._create_layout(G)
        
        # Create visualization
        fig = self._create_plotly_figure(G, pos, protein_ids, metadata_dict, query_protein_id)
        
        # Save visualization
        html_path = self._save_visualization(fig, query_protein_id, output_dir)
        
        # Generate network statistics
        network_stats = self._generate_network_statistics(G, protein_ids, metadata_dict, query_protein_id)
        
        return {
            'visualization_figure': fig,
            'html_path': html_path,
            'network_graph': G,
            'network_properties': {
                'num_nodes': len(G.nodes()),
                'num_edges': len(G.edges()),
                'density': nx.density(G),
                'connected_components': nx.number_connected_components(G)
            },
            'network_statistics': network_stats,
            'query_protein_id': query_protein_id,
            'k_neighbors': self.k_neighbors,
            'similarity_threshold': self.similarity_threshold
        }
    
    def _clean_protein_ids(self, embeddings: np.ndarray, protein_ids: List[str]) -> Tuple[np.ndarray, List[str]]:
        """Clean protein IDs to remove None/NaN values."""
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
            embeddings = embeddings[valid_indices]
            protein_ids = valid_protein_ids
        
        return embeddings, protein_ids
    
    def _merge_query_protein(self, embeddings: np.ndarray, protein_ids: List[str], 
                           query_embedding: np.ndarray, query_protein_id: Optional[str]) -> Tuple[np.ndarray, List[str], str]:
        """Merge query protein into the dataset."""
        if query_embedding.ndim == 1:
            query_embedding = query_embedding.reshape(1, -1)
        
        if query_protein_id is None:
            query_protein_id = "QUERY_PROTEIN"
        
        if query_protein_id in protein_ids:
            logger.info(f"Query protein '{query_protein_id}' is already in the dataset.")
            return embeddings, protein_ids, query_protein_id
        else:
            logger.info(f"Adding new query protein '{query_protein_id}' to the dataset.")
            embeddings = np.vstack([embeddings, query_embedding])
            protein_ids.append(query_protein_id)
            return embeddings, protein_ids, query_protein_id
    
    def _create_metadata_mapping(self, metadata_df: pd.DataFrame, id_column: str) -> Dict[str, Dict]:
        """Create mapping from protein IDs to metadata."""
        if id_column and id_column in metadata_df.columns:
            return metadata_df.set_index(id_column).to_dict('index')
        else:
            return metadata_df.to_dict('index')
    
    def _build_network(self, embeddings: np.ndarray, protein_ids: List[str]) -> nx.Graph:
        """Build a network from embeddings using similarity."""
        # Compute similarity matrix
        sim_matrix = self._compute_similarity_matrix(embeddings)
        
        # Build network with robust edge selection
        G = nx.Graph()
        
        # Add all nodes
        for protein_id in protein_ids:
            G.add_node(protein_id)
        
        # Add edges based on similarity
        N = len(protein_ids)
        for i in range(N):
            similarities = sim_matrix[i]
            
            # Find top k neighbors above threshold
            high_sim_indices = []
            for j in range(N):
                if i != j and similarities[j] >= self.similarity_threshold:
                    high_sim_indices.append(j)
            
            if high_sim_indices:
                high_sim_indices.sort(key=lambda j: similarities[j], reverse=True)
                top_k_indices = high_sim_indices[:self.k_neighbors]
                
                for j in top_k_indices:
                    weight = similarities[j]
                    G.add_edge(protein_ids[i], protein_ids[j], weight=weight)
            else:
                # Ensure connectivity by connecting to top k regardless of threshold
                all_indices = [j for j in range(N) if j != i]
                all_indices.sort(key=lambda j: similarities[j], reverse=True)
                top_k_indices = all_indices[:self.k_neighbors]
                
                for j in top_k_indices:
                    weight = similarities[j]
                    G.add_edge(protein_ids[i], protein_ids[j], weight=weight)
        
        return G
    
    def _compute_similarity_matrix(self, embeddings: np.ndarray) -> np.ndarray:
        """Compute cosine similarity matrix for embeddings."""
        normed = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        sim_matrix = np.dot(normed, normed.T)
        return sim_matrix
    
    def _create_layout(self, G: nx.Graph) -> Dict[str, Tuple[float, float]]:
        """Create Kamada-Kawai layout for the graph."""
        # Create weighted copy for layout
        G_weighted = G.copy()
        
        # Convert similarity to distance
        for u, v, data in G_weighted.edges(data=True):
            similarity = data.get('weight', 0)
            distance = max(0.1, 1.0 - similarity)
            G_weighted[u][v]['weight'] = distance
        
        # Use Kamada-Kawai layout
        pos = nx.kamada_kawai_layout(G_weighted, weight='weight')
        return pos
    
    def _create_plotly_figure(self, G: nx.Graph, pos: Dict, protein_ids: List[str], 
                            metadata_dict: Dict, query_protein_id: str) -> go.Figure:
        """Create the Plotly interactive network figure."""
        # Prepare data for visualization
        node_x, node_y = [], []
        node_customdata = []
        node_colors, node_symbols, node_sizes = [], [], []
        edge_x, edge_y = [], []
        
        # Get connected components
        connected_components = list(nx.connected_components(G))
        node_to_component = {}
        for i, component in enumerate(connected_components):
            for node in component:
                node_to_component[node] = i
        
        # Add edges
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        # Add nodes
        for node in G.nodes():
            node_x.append(pos[node][0])
            node_y.append(pos[node][1])
            
            is_query = (node == query_protein_id)
            
            # Get metadata
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
                protein_metadata = metadata_dict.get(node, {})
                if not protein_metadata:
                    protein_metadata = {
                        'Protein names': f'Protein {node}',
                        'Organism': 'N/A',
                        'EC number': 'N/A',
                        'Function [CC]': 'No metadata available',
                        'Protein families': 'N/A',
                        'Reviewed': 'N/A'
                    }
            
            # Prepare custom data for hover
            component_id = node_to_component.get(node, -1)
            component_size = len(connected_components[component_id]) if component_id >= 0 else 1
            node_degree = G.degree(node)
            
            node_customdata.append([
                str(node),                                    # 0: ID
                str(protein_metadata.get('Protein names', 'N/A')),  # 1: Name
                str(protein_metadata.get('Organism', 'N/A')),        # 2: Organism
                str(protein_metadata.get('EC number', 'N/A')),       # 3: EC Number
                self._wrap_text(protein_metadata.get('Function [CC]', 'N/A')),  # 4: Function
                str(protein_metadata.get('Protein families', 'N/A')), # 5: Family
                str(protein_metadata.get('Reviewed', 'N/A')),        # 6: Reviewed
                component_id,                                 # 7: Component ID
                component_size,                               # 8: Component Size
                node_degree,                                  # 9: Node Degree
                is_query,                                     # 10: Is Query Protein
                1.0 if is_query else 0.0                     # 11: Similarity to Query (placeholder)
            ])
            
            # Assign colors and symbols
            if is_query:
                node_colors.append('red')
                node_symbols.append('star')
                node_sizes.append(15)
            else:
                is_connected_to_query = G.has_edge(node, query_protein_id)
                if is_connected_to_query:
                    node_colors.append('red')
                    node_symbols.append('circle')
                    node_sizes.append(8)
                else:
                    node_colors.append('#CCCCCC')
                    node_symbols.append('circle')
                    node_sizes.append(6)
        
        # Create figure
        fig = go.Figure()
        
        # Add edge trace
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y,
            mode='lines',
            line=dict(width=0.5, color='#888'),
            opacity=0.6,
            hoverinfo='none',
            showlegend=False
        ))
        
        # Add node trace
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
        
        # Add legend
        self._add_legend(fig)
        
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
            )
        )
        
        return fig
    
    def _add_legend(self, fig: go.Figure):
        """Add legend traces to the figure."""
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
    
    def _save_visualization(self, fig: go.Figure, query_protein_id: str, output_dir: str) -> str:
        """Save the visualization as HTML file."""
        os.makedirs(output_dir, exist_ok=True)
        html_filename = f"network_visualization_{query_protein_id}_{int(time.time())}.html"
        html_path = os.path.join(output_dir, html_filename)
        fig.write_html(html_path)
        logger.info(f"Network visualization saved to: {html_path}")
        return html_path
    
    def _generate_network_statistics(self, G: nx.Graph, protein_ids: List[str], 
                                   metadata_dict: Dict, query_protein_id: str) -> Dict[str, Any]:
        """Generate network statistics."""
        stats = {
            'total_nodes': len(G.nodes()),
            'total_edges': len(G.edges()),
            'density': nx.density(G),
            'connected_components': nx.number_connected_components(G),
            'average_degree': np.mean([d for n, d in G.degree()]) if G.nodes() else 0,
            'query_protein_degree': G.degree(query_protein_id) if query_protein_id in G.nodes() else 0
        }
        
        return stats
    
    def _wrap_text(self, text: str, max_chars: int = 80) -> str:
        """Wrap text for clean tooltips."""
        if not isinstance(text, str) or text.strip() == '':
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
