"""
Network Analysis Stage

This stage provides comprehensive network analysis for protein similarity search results,
including tables, visualizations, and statistics, following the BaseStage/StageResult
interface used throughout the pipeline.

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

from ..base_stage import BaseStage, StageResult

logger = logging.getLogger(__name__)

# Visualization imports
try:
    import plotly.graph_objects as go
    import plotly.colors
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly not available. Install with: pip install plotly")


class NetworkAnalysisStage(BaseStage):
    """Pipeline stage for network analysis of similarity search results."""

    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config or {})
        self.k_neighbors = self.config.get('k_neighbors', 8)
        self.similarity_threshold = self.config.get('similarity_threshold', 0.1)
        self.mutual_knn = self.config.get('mutual_knn', True)
        self.min_network_size = self.config.get('min_network_size', 5)
        self.max_network_size = self.config.get('max_network_size', 100)

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
            similar_proteins = input_data['similarity_results']
            embeddings = input_data['embeddings']
            protein_ids = input_data['protein_ids']
            metadata_df = input_data['metadata_df']
            query_embedding = input_data['query_embedding']
            query_protein_id = input_data['query_protein_id']

            results = self.analyze_similarity_search_results(
                similar_proteins=similar_proteins,
                embeddings=embeddings,
                protein_ids=protein_ids,
                metadata_df=metadata_df,
                query_embedding=query_embedding,
                query_protein_id=query_protein_id
            )

            return StageResult(
                success=True,
                output_data={'network_analysis': results},
                metadata={
                    'k_neighbors': self.k_neighbors,
                    'similarity_threshold': self.similarity_threshold,
                    'mutual_knn': self.mutual_knn
                },
                execution_time=0.0
            )
        except Exception as e:
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=0.0,
                error_message=str(e)
            )

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


