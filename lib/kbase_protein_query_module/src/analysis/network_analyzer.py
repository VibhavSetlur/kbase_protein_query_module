"""
Network Analyzer Module

This module provides comprehensive network analysis for protein similarity search results,
including similarity tables, network visualization, and network statistics.

Features:
- Similarity search results in table format
- Interactive network visualization
- Network statistics and properties
- Top similar proteins analysis
- Network clustering and community detection
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import List, Dict, Tuple, Optional, Union, Any
import logging
from tqdm import tqdm
import os
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import AgglomerativeClustering
import json

logger = logging.getLogger(__name__)

# Visualization imports
try:
    import plotly.graph_objects as go
    import plotly.colors
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly not available. Install with: pip install plotly")


class NetworkAnalyzer:
    """
    Comprehensive network analyzer for protein similarity search results.
    
    This class provides functionality to:
    - Generate similarity tables from search results
    - Create interactive network visualizations
    - Analyze network properties and statistics
    - Identify top similar proteins
    - Perform network clustering analysis
    """
    
    def __init__(self, k_neighbors: int = 8, 
                 similarity_threshold: float = 0.1,
                 mutual_knn: bool = True,
                 min_network_size: int = 5,
                 max_network_size: int = 100):
        """
        Initialize the network analyzer.
        
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
    
    def analyze_similarity_search_results(self, 
                                        similar_proteins: List[Dict],
                                        embeddings: np.ndarray,
                                        protein_ids: List[str],
                                        metadata_df: pd.DataFrame,
                                        query_embedding: np.ndarray,
                                        query_protein_id: str) -> Dict[str, Any]:
        """
        Comprehensive analysis of similarity search results.
        
        Args:
            similar_proteins: List of similar proteins from search
            embeddings: Full embeddings array
            protein_ids: Full list of protein IDs
            metadata_df: DataFrame containing protein metadata
            query_embedding: Query protein embedding
            query_protein_id: Query protein ID
            
        Returns:
            Dictionary containing all analysis results
        """
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
        
        # Build network
        G = self._build_network_from_similar_proteins(
            similar_proteins, embeddings, protein_ids, 
            query_embedding, query_protein_id
        )
        
        if G is not None:
            # Analyze network properties
            results['network_properties'] = self._analyze_network_properties(G)
            
            # Get network statistics
            results['network_statistics'] = self._get_network_statistics(G)
            
            # Perform clustering analysis
            results['clustering_analysis'] = self._perform_clustering_analysis(
                similar_proteins, embeddings, protein_ids
            )
            
            # Create network visualization
            if PLOTLY_AVAILABLE:
                results['network_visualization'] = self._create_network_visualization(
                    G, similar_proteins, metadata_df, query_protein_id
                )
        
        return results
    
    def _create_similarity_table(self, similar_proteins: List[Dict], 
                               metadata_df: pd.DataFrame) -> pd.DataFrame:
        """
        Create a comprehensive similarity table from search results.
        
        Args:
            similar_proteins: List of similar proteins from search
            metadata_df: DataFrame containing protein metadata
            
        Returns:
            DataFrame with similarity information and metadata
        """
        table_data = []
        
        for i, protein in enumerate(similar_proteins):
            protein_id = protein['protein_id']
            similarity_score = protein.get('similarity_score', 0.0)
            
            # Get metadata for this protein
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
        """
        Get the top N most similar proteins.
        
        Args:
            similar_proteins: List of similar proteins from search
            top_n: Number of top proteins to return
            
        Returns:
            List of top similar proteins
        """
        # Sort by similarity score (descending)
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
        """
        Build a network from similarity search results.
        
        Args:
            similar_proteins: List of similar proteins from search
            embeddings: Full embeddings array
            protein_ids: Full list of protein IDs
            query_embedding: Query protein embedding
            query_protein_id: Query protein ID
            
        Returns:
            NetworkX graph or None if building fails
        """
        try:
            logger.info(f"Building network from {len(similar_proteins)} similar proteins...")
            
            # Extract embeddings and IDs for similar proteins
            similar_ids = [p['protein_id'] for p in similar_proteins]
            
            # Find indices of similar proteins in the full protein list
            similar_indices = []
            for pid in similar_ids:
                if pid in protein_ids:
                    similar_indices.append(protein_ids.index(pid))
            
            if not similar_indices:
                logger.warning("No similar proteins found in the full protein list")
                # Use a subset of the full embeddings
                similar_indices = list(range(min(50, len(protein_ids))))
            
            similar_embeddings = embeddings[similar_indices]
            similar_protein_ids = [protein_ids[i] for i in similar_indices]
            
            # Add query protein
            if query_protein_id not in similar_protein_ids:
                similar_embeddings = np.vstack([similar_embeddings, query_embedding])
                similar_protein_ids.append(query_protein_id)
            
            # Build network using mutual k-NN
            return self._build_mutual_knn_network(similar_embeddings, similar_protein_ids)
            
        except Exception as e:
            logger.error(f"Failed to build network: {e}")
            return None
    
    def _build_mutual_knn_network(self, embeddings: np.ndarray,
                                protein_ids: List[str]) -> nx.Graph:
        """
        Build a mutual k-nearest neighbors network.
        
        Args:
            embeddings: Protein embeddings array
            protein_ids: List of protein IDs
            
        Returns:
            NetworkX graph
        """
        logger.info("Building mutual k-NN network...")
        
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
            self._trim_network(G)
        
        logger.info(f"Built network with {len(G.nodes())} nodes and {len(G.edges())} edges")
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
    
    def _trim_network(self, G: nx.Graph):
        """Trim network to maximum size while preserving connectivity."""
        # Remove nodes with lowest degree
        while len(G.nodes()) > self.max_network_size:
            degrees = dict(G.degree())
            min_degree_node = min(degrees, key=degrees.get)
            G.remove_node(min_degree_node)
    
    def _analyze_network_properties(self, G: nx.Graph) -> Dict:
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
    
    def _get_network_statistics(self, G: nx.Graph) -> pd.DataFrame:
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
    
    def _perform_clustering_analysis(self, similar_proteins: List[Dict],
                                   embeddings: np.ndarray,
                                   protein_ids: List[str]) -> Dict:
        """
        Perform clustering analysis on similar proteins.
        
        Args:
            similar_proteins: List of similar proteins from search
            embeddings: Full embeddings array
            protein_ids: Full list of protein IDs
            
        Returns:
            Dictionary containing clustering results
        """
        try:
            # Extract embeddings for similar proteins
            similar_ids = [p['protein_id'] for p in similar_proteins]
            similar_indices = []
            for pid in similar_ids:
                if pid in protein_ids:
                    similar_indices.append(protein_ids.index(pid))
            
            if len(similar_indices) < 2:
                return {'clusters': [], 'cluster_assignments': {}}
            
            similar_embeddings = embeddings[similar_indices]
            similar_protein_ids = [protein_ids[i] for i in similar_indices]
            
            # Perform hierarchical clustering
            n_clusters = min(5, len(similar_indices))  # Max 5 clusters
            clustering = AgglomerativeClustering(n_clusters=n_clusters)
            cluster_labels = clustering.fit_predict(similar_embeddings)
            
            # Organize results
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
        """
        Create an interactive network visualization.
        
        Args:
            G: NetworkX graph
            similar_proteins: List of similar proteins from search
            metadata_df: DataFrame containing protein metadata
            query_protein_id: Query protein ID
            
        Returns:
            Plotly figure object or None if visualization fails
        """
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly not available for network visualization")
            return None
        
        try:
            # Create Kamada-Kawai layout
            pos = nx.kamada_kawai_layout(G)
            
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
                
                # Get component info
                component_id = node_to_component.get(node, -1)
                component_size = len(connected_components[component_id]) if component_id >= 0 else 1
                
                # Get node degree
                node_degree = G.degree(node)
                
                # Store metadata for hover
                protein_name = str(protein_metadata.get('Protein names', 'N/A'))
                organism = str(protein_metadata.get('Organism', 'N/A'))
                ec_number = str(protein_metadata.get('EC number', 'N/A'))
                function_text = protein_metadata.get('Function [CC]', 'N/A')
                family = str(protein_metadata.get('Protein families', 'N/A'))
                reviewed = str(protein_metadata.get('Reviewed', 'N/A'))
                
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
                    is_query                                      # 10: Is Query Protein
                ])
                
                # Assign colors and symbols
                if is_query:
                    node_colors.append('red')
                    node_symbols.append('star')
                    node_sizes.append(12)
                else:
                    node_colors.append('#CCCCCC')
                    node_symbols.append('circle')
                    node_sizes.append(6)
            
            # Create figure
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
                marker=dict(size=6, symbol='circle', color='#CCCCCC'),
                name='Similar Proteins',
                showlegend=True
            ))
            
            # Update layout
            fig.update_layout(
                title=f'<b>Protein Similarity Network (k={self.k_neighbors}, threshold={self.similarity_threshold})</b>',
                height=600,
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
            
        except Exception as e:
            logger.error(f"Failed to create network visualization: {e}")
            return None
    
    def save_analysis_results(self, results: Dict[str, Any], output_dir: str):
        """
        Save analysis results to files.
        
        Args:
            results: Analysis results dictionary
            output_dir: Directory to save results
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Save similarity table
        if 'similarity_table' in results:
            similarity_file = os.path.join(output_dir, 'similarity_table.csv')
            results['similarity_table'].to_csv(similarity_file, index=False)
            logger.info(f"Saved similarity table to {similarity_file}")
        
        # Save network statistics
        if 'network_statistics' in results and results['network_statistics'] is not None:
            stats_file = os.path.join(output_dir, 'network_statistics.csv')
            results['network_statistics'].to_csv(stats_file, index=False)
            logger.info(f"Saved network statistics to {stats_file}")
        
        # Save network properties
        if 'network_properties' in results:
            props_file = os.path.join(output_dir, 'network_properties.json')
            with open(props_file, 'w') as f:
                json.dump(results['network_properties'], f, indent=2)
            logger.info(f"Saved network properties to {props_file}")
        
        # Save clustering analysis
        if 'clustering_analysis' in results:
            cluster_file = os.path.join(output_dir, 'clustering_analysis.json')
            with open(cluster_file, 'w') as f:
                json.dump(results['clustering_analysis'], f, indent=2)
            logger.info(f"Saved clustering analysis to {cluster_file}")
        
        # Save network visualization
        if 'network_visualization' in results and results['network_visualization'] is not None:
            viz_file = os.path.join(output_dir, 'network_visualization.html')
            results['network_visualization'].write_html(viz_file)
            logger.info(f"Saved network visualization to {viz_file}")
        
        # Save top similar proteins
        if 'top_similar_proteins' in results:
            top_file = os.path.join(output_dir, 'top_similar_proteins.json')
            with open(top_file, 'w') as f:
                json.dump(results['top_similar_proteins'], f, indent=2)
            logger.info(f"Saved top similar proteins to {top_file}")


def create_network_analyzer(**kwargs) -> NetworkAnalyzer:
    """
    Factory function to create a network analyzer.
    
    Args:
        **kwargs: Parameters for NetworkAnalyzer
        
    Returns:
        NetworkAnalyzer instance
    """
    return NetworkAnalyzer(**kwargs)
