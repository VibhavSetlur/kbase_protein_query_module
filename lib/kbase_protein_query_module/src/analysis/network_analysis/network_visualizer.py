"""
Network Visualization Module

Creates interactive Plotly-based network visualizations for protein similarity networks.
Simple, clean, and functional.
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import List, Dict, Tuple, Optional, Union, Any
import logging
import os
import time

logger = logging.getLogger(__name__)

# Required dependencies
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    logger.error("Plotly is required. Install with: pip install plotly")

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    nx = None
    logger.error("NetworkX is required. Install with: pip install networkx")


class NetworkVisualizer:
    """
    Creates interactive network visualizations for protein similarity networks.
    
    Simple, clean, and functional visualization with interactive hover metadata.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """Initialize the Network Visualizer."""
        if not HAS_PLOTLY:
            raise ImportError("Plotly is required for network visualization")
        if not HAS_NETWORKX:
            raise ImportError("NetworkX is required for network visualization")
        
        self.config = config or {}
        self.k_neighbors = self.config.get('k_neighbors', 10)
        self.similarity_threshold = self.config.get('similarity_threshold', 0.1)
    
    def create_interactive_visualization(self,
                                       embeddings: np.ndarray,
                                       protein_ids: List[str],
                                       metadata: Union[pd.DataFrame, Dict[str, Dict]],
                                       query_embedding: Optional[np.ndarray] = None,
                                       query_protein_id: Optional[str] = None,
                                       output_dir: str = 'outputs',
                                       id_column: str = 'Entry') -> Dict[str, Any]:
        """
        Create an interactive network visualization.
        
        Args:
            embeddings: Protein embeddings array (N x D)
            protein_ids: List of protein IDs
            metadata: DataFrame or dict with protein metadata
            query_embedding: Query protein embedding (optional)
            query_protein_id: Query protein ID (optional)
            output_dir: Directory to save visualization
            id_column: Column name for protein IDs in metadata
        
        Returns:
            Dictionary with visualization results and network graph
        """
        logger.info("Creating interactive network visualization...")
        
        # Prepare data
        embeddings, protein_ids = self._validate_and_clean_ids(embeddings, protein_ids)
        embeddings, protein_ids, query_protein_id = self._prepare_query_protein(
            embeddings, protein_ids, query_embedding, query_protein_id
        )
        metadata_dict = self._prepare_metadata(metadata, id_column)
        
        # Build network
        network_graph = self._build_network(embeddings, protein_ids)
        node_positions = self._compute_layout(network_graph)
        
        # Compute similarity to query for hover
        similarity_scores = self._compute_similarity_to_query(
            embeddings, protein_ids, query_embedding
        )
        
        # Create visualization
        fig = self._create_plotly_figure(
            network_graph, node_positions, protein_ids, metadata_dict,
            query_protein_id, similarity_scores
        )
        
        # Save visualization
        html_path = self._save_html(fig, query_protein_id, output_dir)
        
        # Generate statistics
        network_stats = self._compute_network_statistics(network_graph, query_protein_id)
        
        return {
            'visualization_figure': fig,
            'html_path': html_path,
            'network_graph': network_graph,
            'network_properties': {
                'num_nodes': len(network_graph.nodes()),
                'num_edges': len(network_graph.edges()),
                'density': nx.density(network_graph),
                'connected_components': nx.number_connected_components(network_graph)
            },
            'network_statistics': network_stats,
            'query_protein_id': query_protein_id,
            'k_neighbors': self.k_neighbors,
            'similarity_threshold': self.similarity_threshold
        }
    
    def _validate_and_clean_ids(self, embeddings: np.ndarray, protein_ids: List[str]) -> Tuple[np.ndarray, List[str]]:
        """Remove invalid protein IDs and corresponding embeddings."""
        valid_indices = []
        valid_ids = []
        
        for i, protein_id in enumerate(protein_ids):
            if protein_id and str(protein_id).strip() and str(protein_id).lower() != 'nan':
                valid_indices.append(i)
                valid_ids.append(protein_id)
            else:
                logger.warning(f"Removing invalid protein ID at index {i}: {protein_id}")
        
        if len(valid_ids) < len(protein_ids):
            logger.info(f"Filtered {len(protein_ids) - len(valid_ids)} invalid protein IDs")
            embeddings = embeddings[valid_indices]
            protein_ids = valid_ids
        
        return embeddings, protein_ids
    
    def _prepare_query_protein(self, embeddings: np.ndarray, protein_ids: List[str],
                              query_embedding: Optional[np.ndarray],
                              query_protein_id: Optional[str]) -> Tuple[np.ndarray, List[str], str]:
        """Add query protein to dataset if not already present."""
        if query_embedding is None:
            query_protein_id = protein_ids[0] if protein_ids else "QUERY_PROTEIN"
            return embeddings, protein_ids, query_protein_id
        
        # Ensure query embedding is 2D
        if query_embedding.ndim == 1:
            query_embedding = query_embedding.reshape(1, -1)
        
        if query_protein_id is None:
            query_protein_id = "QUERY_PROTEIN"
        
        # Add query protein if not already in dataset
        if query_protein_id not in protein_ids:
            embeddings = np.vstack([embeddings, query_embedding])
            protein_ids.append(query_protein_id)
            logger.info(f"Added query protein '{query_protein_id}' to network")
        else:
            logger.info(f"Query protein '{query_protein_id}' already in dataset")
        
        return embeddings, protein_ids, query_protein_id
    
    def _prepare_metadata(self, metadata: Union[pd.DataFrame, Dict[str, Dict]], id_column: str) -> Dict[str, Dict]:
        """Convert metadata to dictionary mapping protein IDs to metadata."""
        if isinstance(metadata, dict):
            return metadata
        elif isinstance(metadata, pd.DataFrame):
            if id_column and id_column in metadata.columns:
                return metadata.set_index(id_column).to_dict('index')
            else:
                return metadata.to_dict('index')
        else:
            return {}
    
    def _build_network(self, embeddings: np.ndarray, protein_ids: List[str]) -> nx.Graph:
        """Build network graph from protein embeddings using cosine similarity."""
        # Compute similarity matrix
        normalized_embeddings = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        similarity_matrix = np.dot(normalized_embeddings, normalized_embeddings.T)
        
        # Create graph
        graph = nx.Graph()
        graph.add_nodes_from(protein_ids)
        
        # Add edges based on similarity
        num_proteins = len(protein_ids)
        for i in range(num_proteins):
            similarities = similarity_matrix[i]
            
            # Find neighbors above threshold
            candidate_indices = [
                j for j in range(num_proteins)
                if i != j and similarities[j] >= self.similarity_threshold
            ]
            
            if candidate_indices:
                # Sort by similarity and take top k
                candidate_indices.sort(key=lambda j: similarities[j], reverse=True)
                top_neighbors = candidate_indices[:self.k_neighbors]
            else:
                # If no neighbors above threshold, connect to top k anyway for connectivity
                all_indices = [j for j in range(num_proteins) if j != i]
                all_indices.sort(key=lambda j: similarities[j], reverse=True)
                top_neighbors = all_indices[:self.k_neighbors]
            
            # Add edges
            for j in top_neighbors:
                similarity = float(similarities[j])
                graph.add_edge(protein_ids[i], protein_ids[j], weight=similarity)
        
        logger.info(f"Built network with {len(graph.nodes())} nodes and {len(graph.edges())} edges")
        return graph
    
    def _compute_layout(self, graph: nx.Graph) -> Dict[str, Tuple[float, float]]:
        """Compute node positions using Kamada-Kawai layout algorithm."""
        # Create weighted copy (convert similarity to distance)
        weighted_graph = graph.copy()
        for u, v, data in weighted_graph.edges(data=True):
            similarity = data.get('weight', 0)
            distance = max(0.1, 1.0 - similarity)
            weighted_graph[u][v]['weight'] = distance
        
        # Compute layout
        positions = nx.kamada_kawai_layout(weighted_graph, weight='weight')
        return positions
    
    def _compute_similarity_to_query(self, embeddings: np.ndarray, protein_ids: List[str],
                                    query_embedding: Optional[np.ndarray]) -> Dict[str, float]:
        """Compute cosine similarity of each protein to the query protein."""
        if query_embedding is None:
            return {}
        
        try:
            # Normalize query embedding
            query_vec = query_embedding.astype(np.float32)
            if query_vec.ndim > 1:
                query_vec = query_vec.reshape(-1)
            query_vec = query_vec / (np.linalg.norm(query_vec) + 1e-8)
            
            # Normalize all embeddings
            normalized_embeddings = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
            
            # Compute similarities
            similarities = np.dot(normalized_embeddings, query_vec)
            
            # Create dictionary
            similarity_dict = {}
            for protein_id, similarity in zip(protein_ids, similarities):
                similarity_dict[protein_id] = float(similarity)
            
            return similarity_dict
        except Exception as e:
            logger.warning(f"Failed to compute similarity to query: {e}")
            return {}
    
    def _create_plotly_figure(self, graph: nx.Graph, positions: Dict[str, Tuple[float, float]],
                             protein_ids: List[str], metadata_dict: Dict[str, Dict],
                             query_protein_id: str, similarity_scores: Dict[str, float]) -> go.Figure:
        """Create interactive Plotly figure with enhanced hover information."""
        # Prepare edge traces
        edge_x, edge_y = [], []
        for u, v in graph.edges():
            x0, y0 = positions[u]
            x1, y1 = positions[v]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        # Prepare node data
        node_x, node_y = [], []
        node_hover_data = []
        node_colors, node_symbols, node_sizes = [], [], []
        
        # Get connected components for display
        components = list(nx.connected_components(graph))
        node_to_component = {}
        for comp_id, component in enumerate(components):
            for node in component:
                node_to_component[node] = comp_id
        
        # Process each node
        for node in graph.nodes():
            node_x.append(positions[node][0])
            node_y.append(positions[node][1])
            
            is_query = (node == query_protein_id)
            is_connected_to_query = graph.has_edge(node, query_protein_id) if query_protein_id in graph.nodes() else False
            
            # Get metadata
            node_metadata = metadata_dict.get(node, {})
            if not node_metadata:
                node_metadata = {
                    'Protein names': f'Protein {node}',
                    'Organism': 'N/A',
                    'EC number': 'N/A',
                    'Protein families': 'N/A',
                    'Reviewed': 'N/A',
                    'Function [CC]': 'No metadata available'
                }
            elif is_query:
                # Override for query protein
                node_metadata['Protein names'] = 'Query Protein'
                node_metadata['Function [CC]'] = 'Query protein - user input'
            
            # Compute node properties
            component_id = node_to_component.get(node, -1)
            component_size = len(components[component_id]) if component_id >= 0 else 1
            node_degree = graph.degree(node)
            clustering_coef = nx.clustering(graph, node)
            similarity_to_query = similarity_scores.get(node, 0.0)
            if is_query:
                similarity_to_query = 1.0
            
            # Prepare hover data as list for Plotly (rich metadata for interactive display)
            hover_data = [
                str(node),  # 0: protein_id
                str(node_metadata.get('Protein names', 'N/A')),  # 1: protein_name
                str(node_metadata.get('Organism', 'N/A')),  # 2: organism
                str(node_metadata.get('EC number', 'N/A')),  # 3: ec_number
                str(node_metadata.get('Protein families', 'N/A')),  # 4: protein_family
                str(node_metadata.get('Reviewed', 'N/A')),  # 5: reviewed
                self._format_text(node_metadata.get('Function [CC]', 'N/A')),  # 6: function
                f"{similarity_to_query:.3f}",  # 7: similarity_to_query
                str(node_degree),  # 8: node_degree
                f"{clustering_coef:.3f}",  # 9: clustering_coefficient
                str(component_id),  # 10: component_id
                str(component_size),  # 11: component_size
                'Query Protein' if is_query else ('Connected to Query' if is_connected_to_query else 'Other Protein')  # 12: node_type
            ]
            node_hover_data.append(hover_data)
            
            # Set node appearance
            if is_query:
                node_colors.append('red')
                node_symbols.append('star')
                node_sizes.append(20)
            elif is_connected_to_query:
                node_colors.append('#FF6B6B')
                node_symbols.append('circle')
                node_sizes.append(12)
            else:
                node_colors.append('#4ECDC4')
                node_symbols.append('circle')
                node_sizes.append(8)
        
        # Create figure
        fig = go.Figure()
        
        # Add edges
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y,
            mode='lines',
            line=dict(width=1, color='#CCCCCC'),
            opacity=0.4,
            hoverinfo='none',
            showlegend=False
        ))
        
        # Add nodes with enhanced hover
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            customdata=node_hover_data,
            marker=dict(
                size=node_sizes,
                color=node_colors,
                symbol=node_symbols,
                line=dict(width=2, color='black'),
                opacity=0.9
            ),
            hovertemplate=self._create_hover_template(),
            hoverinfo='all',
            showlegend=False
        ))
        
        # Add legend
        self._add_legend(fig)
        
        # Update layout
        fig.update_layout(
            title=dict(
                text=f'<b>Interactive Protein Similarity Network</b><br><sub>k={self.k_neighbors}, threshold={self.similarity_threshold:.2f}</sub>',
                x=0.5,
                xanchor='center'
            ),
            height=900,
            plot_bgcolor='white',
            hovermode='closest',
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            margin=dict(l=20, r=20, t=80, b=20),
            legend=dict(
                x=0.02,
                y=0.98,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='black',
                borderwidth=1,
                font=dict(size=12)
            )
        )
        
        return fig
    
    def _create_hover_template(self) -> str:
        """Create rich hover template for interactive metadata display."""
        return (
            "<b>Protein ID:</b> %{customdata[0]}<br>"
            "<b>Name:</b> %{customdata[1]}<br>"
            "<b>Organism:</b> %{customdata[2]}<br>"
            "<b>EC Number:</b> %{customdata[3]}<br>"
            "<b>Family:</b> %{customdata[4]}<br>"
            "<b>Reviewed:</b> %{customdata[5]}<br>"
            "<hr>"
            "<b>Network Properties:</b><br>"
            "  Similarity to Query: %{customdata[7]}<br>"
            "  Node Degree: %{customdata[8]}<br>"
            "  Clustering Coefficient: %{customdata[9]}<br>"
            "  Component: %{customdata[10]} (Size: %{customdata[11]})<br>"
            "  Node Type: %{customdata[12]}<br>"
            "<hr>"
            "<b>Function:</b><br>%{customdata[6]}<br>"
            "<extra></extra>"
        )
    
    def _format_text(self, text: str, max_length: int = 100) -> str:
        """Format text for display in hover tooltips."""
        if not text or not isinstance(text, str) or text.strip() == '':
            return 'N/A'
        
        # Clean up text
        text = text.strip()
        
        # Wrap long text
        if len(text) > max_length:
            words = text.split()
            lines = []
            current_line = []
            current_length = 0
            
            for word in words:
                word_length = len(word) + 1
                if current_length + word_length <= max_length:
                    current_line.append(word)
                    current_length += word_length
                else:
                    if current_line:
                        lines.append(' '.join(current_line))
                    current_line = [word]
                    current_length = word_length
            
            if current_line:
                lines.append(' '.join(current_line))
            
            text = '<br>'.join(lines)
        
        return text
    
    def _add_legend(self, fig: go.Figure):
        """Add legend to the figure."""
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=15, symbol='star', color='red'),
            name='Query Protein',
            showlegend=True
        ))
        
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=12, symbol='circle', color='#FF6B6B'),
            name='Connected to Query',
            showlegend=True
        ))
        
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=8, symbol='circle', color='#4ECDC4'),
            name='Other Proteins',
            showlegend=True
        ))
    
    def _save_html(self, fig: go.Figure, query_protein_id: str, output_dir: str) -> str:
        """Save visualization as HTML file."""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = int(time.time())
        html_filename = f"network_visualization_{query_protein_id}_{timestamp}.html"
        html_path = os.path.join(output_dir, html_filename)
        
        # Save with embedded plotly.js for robust offline viewing
        # Using 'include' embeds the full library for offline use
        fig.write_html(html_path, include_plotlyjs='include', config={'displayModeBar': True})
        logger.info(f"Network visualization saved to: {html_path}")
        
        return html_path
    
    def _compute_network_statistics(self, graph: nx.Graph, query_protein_id: str) -> Dict[str, Any]:
        """Compute network statistics."""
        if not graph.nodes():
            return {}
        
        degrees = [d for _, d in graph.degree()]
        
        stats = {
            'total_nodes': len(graph.nodes()),
            'total_edges': len(graph.edges()),
            'density': nx.density(graph),
            'connected_components': nx.number_connected_components(graph),
            'average_degree': np.mean(degrees) if degrees else 0,
            'max_degree': max(degrees) if degrees else 0,
            'min_degree': min(degrees) if degrees else 0,
            'query_protein_degree': graph.degree(query_protein_id) if query_protein_id in graph.nodes() else 0
        }
        
        return stats
