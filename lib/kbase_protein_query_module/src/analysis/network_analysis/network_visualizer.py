"""
Network Visualization Module

Creates professional, interactive dashboards for protein similarity networks.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union, Any
import logging
import os
import time
import json

logger = logging.getLogger(__name__)

# Required dependencies
try:
    import plotly.graph_objects as go
    from plotly.utils import PlotlyJSONEncoder
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
    Creates research-level interactive dashboards for protein similarity networks.
    
    Generates a standalone HTML file containing:
    - Interactive 2D Network Graph (Plotly)
    - Sortable/Searchable Data Table
    - Detailed Metadata Sidebar
    - Statistical Analysis
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the Network Visualizer.
        
        Args:
            config: Configuration dictionary.
        """
        if not HAS_PLOTLY:
            raise ImportError("Plotly is required for network visualization")
        if not HAS_NETWORKX:
            raise ImportError("NetworkX is required for network visualization")
        
        self.config = config or {}
        self.k_neighbors = self.config.get('k_neighbors', 10)
        self.similarity_threshold = self.config.get('similarity_threshold', 0.995)
    
    def create_interactive_visualization(self,
                                       embeddings: np.ndarray,
                                       protein_ids: List[str],
                                       metadata: Union[pd.DataFrame, Dict[str, Dict]],
                                       query_embedding: Optional[np.ndarray] = None,
                                       query_protein_id: Optional[str] = None,
                                       query_sequence: Optional[str] = None,
                                       output_dir: str = 'outputs',
                                       id_column: str = 'Entry') -> Dict[str, Any]:
        """
        Create the full interactive dashboard.
        
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
        logger.info("Creating interactive network dashboard...")
        
        # Prepare data
        embeddings, protein_ids = self._validate_and_clean_ids(embeddings, protein_ids)
        embeddings, protein_ids, query_protein_id = self._prepare_query_protein(
            embeddings, protein_ids, query_embedding, query_protein_id
        )
        metadata_dict = self._prepare_metadata(metadata, id_column)
        
        # Build network
        network_graph = self._build_network(embeddings, protein_ids)
        node_positions = self._compute_layout(network_graph)
        
        # Compute similarity to query
        similarity_scores = self._compute_similarity_to_query(
            embeddings, protein_ids, query_embedding
        )
        
        # Generate Network Statistics
        network_stats = self._compute_network_statistics(network_graph, query_protein_id)
        
        # Prepare Dashboard Data Payload
        dashboard_data = self._prepare_dashboard_data(
            network_graph, node_positions, protein_ids, metadata_dict,
            query_protein_id, query_sequence, similarity_scores, network_stats
        )
        
        # Create Plotly Figure (for fallback/embedding check, not primary display)
        # Note: We are now generating custom HTML that re-creates the plot using PlotlyJS
        # but we can still create the python figure object if needed.
        # For this implementation, we will generate the figure as JSON and embed it.
        fig = self._create_plotly_figure_object(
            network_graph, node_positions, protein_ids, metadata_dict,
            query_protein_id, similarity_scores
        )
        
        # Save Dashboard
        html_path = self._save_dashboard_html(dashboard_data, fig, query_protein_id, output_dir)
        
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
    
    # ... Helper methods reused from previous implementation ...
    
    def _validate_and_clean_ids(self, embeddings: np.ndarray, protein_ids: List[str]) -> Tuple[np.ndarray, List[str]]:
        """Remove invalid protein IDs."""
        valid_indices = []
        valid_ids = []
        for i, protein_id in enumerate(protein_ids):
            if protein_id and str(protein_id).strip() and str(protein_id).lower() != 'nan':
                valid_indices.append(i)
                valid_ids.append(protein_id)
        
        if len(valid_ids) < len(protein_ids):
            logger.info(f"Filtered invalid IDs. Original: {len(protein_ids)}, Valid: {len(valid_ids)}")
            embeddings = embeddings[valid_indices]
            protein_ids = valid_ids
        return embeddings, protein_ids

    def _prepare_query_protein(self, embeddings, protein_ids, query_embedding, query_protein_id):
        """Add query protein if missing."""
        if query_embedding is None:
            query_protein_id = protein_ids[0] if protein_ids else "QUERY_PROTEIN"
            return embeddings, protein_ids, query_protein_id
        
        if query_embedding.ndim == 1:
            query_embedding = query_embedding.reshape(1, -1)
        
        query_protein_id = query_protein_id or "QUERY_PROTEIN"
        if query_protein_id not in protein_ids:
            embeddings = np.vstack([embeddings, query_embedding])
            protein_ids.append(query_protein_id)
        return embeddings, protein_ids, query_protein_id

    def _prepare_metadata(self, metadata, id_column):
        """Convert metadata to dict."""
        if isinstance(metadata, dict): return metadata
        if isinstance(metadata, pd.DataFrame):
            if id_column in metadata.columns: return metadata.set_index(id_column).to_dict('index')
            return metadata.to_dict('index')
        return {}

    def _build_network(self, embeddings, protein_ids):
        """
        Build similarity network using dynamic k-NN with query priority.
        Optimized for visual clarity (avoiding hairballs).
        """
        normalized_embeddings = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        similarity_matrix = np.dot(normalized_embeddings, normalized_embeddings.T)
        
        graph = nx.Graph()
        graph.add_nodes_from(protein_ids)
        
        num_proteins = len(protein_ids)
        if num_proteins < 2:
            return graph
            
        # Dynamic Thresholding Calculation
        upper_tri = np.triu(similarity_matrix, k=1)
        non_zero_sims = upper_tri[upper_tri > 0]
        
        dynamic_threshold = self.similarity_threshold
        if len(non_zero_sims) > 0:
            mean_sim = np.mean(non_zero_sims)
            std_sim = np.std(non_zero_sims)
            # Stricter Statistical Threshold: Mean + 1.0 * StdDev for cleaner graph
            stat_threshold = mean_sim + 1.0 * std_sim
            # Ensure floor is high enough (0.8) to avoid noise
            dynamic_threshold = min(self.similarity_threshold, max(0.8, stat_threshold))
        
        logger.info(f"Building network with dynamic threshold: {dynamic_threshold:.4f}")

        for i in range(num_proteins):
            sim_row = similarity_matrix[i]
            neighbor_indices = np.argsort(sim_row)[::-1]
            edges_added = 0
            
            for j in neighbor_indices:
                if i == j: continue
                
                sim = float(sim_row[j])
                
                # Logic:
                # 1. Query Priority: Keep top K for query to show context.
                # 2. Backbone: Keep top 2 for everyone else (Minimum Spanning-ish).
                # 3. Strong Edges: Keep if very high similarity.
                
                is_strong = sim >= dynamic_threshold
                is_query_edge = (i == 0 or j == 0) # Assuming 0 is query usually, or checking ID later
                
                # Determine max edges for this node type
                # Query gets full K (context)
                # Others get strict limit (1) to create tree-like branches unless sim is very high
                max_edges = self.k_neighbors if is_query_edge else 1 
                
                if (is_strong or edges_added < max_edges):
                    if not graph.has_edge(protein_ids[i], protein_ids[j]):
                        graph.add_edge(protein_ids[i], protein_ids[j], weight=sim)
                    edges_added += 1
                
                # Hard stop to prevent starbursts on common domains
                if edges_added >= self.k_neighbors and not is_query_edge:
                    break
                    
        return graph

    def _compute_layout(self, graph):
        """
        Compute layout optimized for biological networks.
        Uses Fruchterman-Reingold (Spring) layout with adjusted k for separation.
        """
        if not graph.nodes():
            return {}
            
        # Spring layout often separates clusters better than Kamada-Kawai for dense graphs
        # k parameter controls node spacing. Higher = more spread.
        # k = 1/sqrt(n) is default. We scale it up.
        k_val = 2.0 / np.sqrt(len(graph.nodes()) or 1)
        
        pos = nx.spring_layout(
            graph, 
            weight='weight', 
            k=k_val, 
            iterations=100, 
            seed=42
        )
        return pos

    def _compute_similarity_to_query(self, embeddings, protein_ids, query_embedding):
        """Compute cosine similarity to query."""
        if query_embedding is None: return {}
        q_vec = query_embedding.flatten()
        q_vec = q_vec / (np.linalg.norm(q_vec) + 1e-8)
        norm_emb = embeddings / (np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-8)
        sims = np.dot(norm_emb, q_vec)
        return {pid: float(sim) for pid, sim in zip(protein_ids, sims)}

    def _compute_network_statistics(self, graph, query_protein_id):
        """Compute statistical metrics."""
        if not graph.nodes(): return {}
        degrees = [d for _, d in graph.degree()]
        return {
            'total_nodes': len(graph.nodes()),
            'total_edges': len(graph.edges()),
            'density': nx.density(graph),
            'connected_components': nx.number_connected_components(graph),
            'avg_degree': np.mean(degrees) if degrees else 0,
            'max_degree': max(degrees) if degrees else 0,
            'query_degree': graph.degree(query_protein_id) if query_protein_id in graph.nodes() else 0
        }

    def _prepare_dashboard_data(self, graph, positions, protein_ids, metadata_dict, 
                              query_protein_id, query_sequence, similarity_scores, stats):
        """Serialize all data for the frontend."""
        nodes = []
        for pid in protein_ids:
            meta = metadata_dict.get(pid, {})
            is_query = (pid == query_protein_id)
            
            # Inject Sequence for Query
            sequence = meta.get('Sequence', '') 
            if is_query and query_sequence:
                sequence = query_sequence
                
            # Infer length if missing
            length = meta.get('Length', 'N/A')
            if (length == 'N/A' or not length) and sequence:
                length = str(len(sequence))
                meta['Length'] = length
            
            nodes.append({
                'id': pid,
                'x': positions[pid][0],
                'y': positions[pid][1],
                'is_query': is_query,
                'sequence': sequence,
                'similarity': similarity_scores.get(pid, 1.0 if is_query else 0.0),
                'degree': graph.degree(pid),
                'cluster_coeff': nx.clustering(graph, pid),
                'component': 0, # Simplified for now
                'meta': {
                    'name': meta.get('Protein names', 'N/A'),
                    'organism': meta.get('Organism', 'N/A'),
                    'length': meta.get('Length', 'N/A'),
                    'gene': meta.get('Gene Names', 'N/A'),
                    'ec': meta.get('EC number', 'N/A'),
                    'reviewed': meta.get('Reviewed', 'N/A'),
                    'function': meta.get('Function [CC]', 'N/A'),
                    'families': meta.get('Protein families', 'N/A')
                }
            })
            
        edges = []
        for u, v, data in graph.edges(data=True):
            edges.append({
                'source': u,
                'target': v,
                'weight': data.get('weight', 0.0)
            })
            
        return {
            'nodes': nodes,
            'edges': edges,
            'stats': stats,
            'query_id': query_protein_id,
            'generated_at': time.strftime("%Y-%m-%d %H:%M:%S")
        }

    def _create_plotly_figure_object(self, graph, positions, protein_ids, metadata_dict, query_protein_id, similarity_scores):
        """Create the Plotly Figure object (for JSON serialization)."""
        edge_x, edge_y = [], []
        for u, v in graph.edges():
            x0, y0 = positions[u]
            x1, y1 = positions[v]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            
        node_x, node_y, node_color, node_size, node_text = [], [], [], [], []
        
        for node in graph.nodes():
            node_x.append(positions[node][0])
            node_y.append(positions[node][1])
            is_query = (node == query_protein_id)
            node_color.append('#ef4444' if is_query else '#22d3ee') # Red or Cyan
            node_size.append(20 if is_query else 10)
            meta = metadata_dict.get(node, {})
            node_text.append(f"{node}<br>{meta.get('Protein names', 'N/A')}")
            
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y, mode='lines', 
            line=dict(width=1, color='#cbd5e1'), hoverinfo='none'
        ))
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y, mode='markers',
            marker=dict(size=node_size, color=node_color, line=dict(width=1, color='#0f172a')),
            text=node_text, hoverinfo='text'
        ))
        fig.update_layout(
            showlegend=False,
            plot_bgcolor='white',
            margin=dict(l=0, r=0, t=0, b=0),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
        )
        return fig

    def _save_dashboard_html(self, data: Dict, fig: go.Figure, query_id: str, output_dir: str) -> str:
        """Generate and save the HTML dashboard."""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = int(time.time())
        filename = f"network_analysis_{query_id}_{timestamp}.html"
        path = os.path.join(output_dir, filename)
        
        # Serialize data
        json_data = json.dumps(data, cls=PlotlyJSONEncoder)
        
        # HTML Template
        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Protein Network Analysis | {query_id}</title>
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <!-- Font -->
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
    <!-- Tailwind CSS (via CDN for portability) -->
    <script src="https://cdn.tailwindcss.com"></script>
    <style>
        body {{ font-family: 'Inter', sans-serif; background-color: #f8fafc; }}
        .tab-active {{ border-bottom: 2px solid #2563eb; color: #2563eb; font-weight: 600; }}
        .tab-inactive {{ color: #64748b; cursor: pointer; }}
        .tab-inactive:hover {{ color: #334155; }}
        .custom-scrollbar::-webkit-scrollbar {{ width: 6px; }}
        .custom-scrollbar::-webkit-scrollbar-track {{ background: #f1f5f9; }}
        .custom-scrollbar::-webkit-scrollbar-thumb {{ background-color: #cbd5e1; border-radius: 3px; }}
        tr:hover td {{ background-color: #f1f5f9; }}
    </style>
</head>
<body class="h-screen flex flex-col overflow-hidden text-slate-800">

    <!-- Header -->
    <header class="bg-white border-b border-slate-200 px-6 py-3 flex justify-between items-center shadow-sm z-10">
        <div>
            <h1 class="text-xl font-bold text-slate-900 flex items-center gap-2">
                <svg class="w-6 h-6 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19.428 15.428a2 2 0 00-1.022-.547l-2.384-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z"></path></svg>
                Protein Query Analysis
            </h1>
            <p class="text-sm text-slate-500">Query ID: <span class="font-mono bg-slate-100 px-1 rounded">{query_id}</span> • {data['stats']['total_nodes']} Proteins • {data['stats']['total_edges']} Interactions</p>
        </div>
        <div class="flex gap-4">
            <button onclick="switchTab('network')" id="tab-btn-network" class="tab-active px-4 py-2">Network View</button>
            <button onclick="switchTab('table')" id="tab-btn-table" class="tab-inactive px-4 py-2">Data Table</button>
        </div>
    </header>

    <!-- Content Area -->
    <main class="flex-1 relative overflow-hidden">
        
        <!-- Network Tab -->
        <div id="tab-network" class="absolute inset-0 flex">
            <!-- Graph Area -->
            <div class="flex-1 relative bg-white" id="plotly-graph"></div>
            
            <!-- Sidebar -->
            <div class="w-96 bg-slate-50 border-l border-slate-200 flex flex-col shadow-lg z-20">
                <div class="p-4 border-b border-slate-200 bg-white">
                    <h2 class="font-semibold text-slate-900 mb-2">Selected Protein</h2>
                    <div id="node-details" class="text-sm">
                        <div class="p-3 bg-blue-50 text-blue-700 rounded-md text-center">
                            Click a node to view details
                        </div>
                    </div>
                </div>
                
                <div class="flex-1 overflow-y-auto custom-scrollbar p-4 space-y-6">
                    <!-- Controls -->
                    <div>
                        <h3 class="text-xs font-bold text-slate-400 uppercase tracking-wider mb-3">Controls</h3>
                        <div class="space-y-2">
                             <label class="flex items-center gap-2 cursor-pointer">
                                <input type="checkbox" id="show-labels" class="rounded text-blue-600" onchange="updateGraph()">
                                <span class="text-sm text-slate-600">Show Node Labels</span>
                            </label>
                            <label class="flex items-center gap-2 cursor-pointer">
                                <input type="checkbox" id="color-edges" class="rounded text-blue-600" onchange="updateGraph()">
                                <span class="text-sm text-slate-600">Color Edges by Similarity</span>
                            </label>
                        </div>
                    </div>
                
                    <!-- Stats -->
                    <div>
                        <h3 class="text-xs font-bold text-slate-400 uppercase tracking-wider mb-3">Network Statistics</h3>
                        <div class="grid grid-cols-2 gap-3 text-sm">
                            <div class="bg-white p-3 rounded border border-slate-200">
                                <span class="block text-slate-500 text-xs">Density</span>
                                <span class="font-mono font-medium">{data['stats']['density']:.3f}</span>
                            </div>
                            <div class="bg-white p-3 rounded border border-slate-200">
                                <span class="block text-slate-500 text-xs">Avg Degree</span>
                                <span class="font-mono font-medium">{data['stats']['avg_degree']:.1f}</span>
                            </div>
                             <div class="bg-white p-3 rounded border border-slate-200">
                                <span class="block text-slate-500 text-xs">Clustering Coeff</span>
                                <span class="font-mono font-medium" id="global-clustering">-</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Table Tab -->
        <div id="tab-table" class="absolute inset-0 bg-white overflow-auto hidden p-6">
            <div class="max-w-7xl mx-auto">
                <div class="mb-4 flex justify-between items-center">
                    <h2 class="text-lg font-bold">Protein Metadata Table</h2>
                    <input type="text" id="table-search" placeholder="Search..." class="border border-slate-300 rounded px-3 py-1.5 text-sm w-64 focus:outline-none focus:ring-2 focus:ring-blue-500">
                </div>
                <div class="border border-slate-200 rounded-lg overflow-hidden">
                    <table class="w-full text-sm text-left">
                        <thead class="bg-slate-50 text-slate-500 font-medium border-b border-slate-200">
                            <tr>
                                <th class="px-4 py-3 cursor-pointer hover:bg-slate-100" onclick="sortTable('id')">ID</th>
                                <th class="px-4 py-3 cursor-pointer hover:bg-slate-100" onclick="sortTable('name')">Name</th>
                                <th class="px-4 py-3 cursor-pointer hover:bg-slate-100" onclick="sortTable('organism')">Organism</th>
                                <th class="px-4 py-3 cursor-pointer hover:bg-slate-100" onclick="sortTable('similarity')">Similarity to Query</th>
                                <th class="px-4 py-3">Length</th>
                                <th class="px-4 py-3">Gene</th>
                                <th class="px-4 py-3">Actions</th>
                            </tr>
                        </thead>
                        <tbody id="table-body" class="divide-y divide-slate-100">
                            <!-- Populated by JS -->
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    </main>

    <script>
        // Data Payload
        const DATA = {json_data};
        
        // --- State ---
        let currentSort = {{ field: 'similarity', dir: 'desc' }};
        
        // --- Initialization ---
        function init() {{
            renderGraph();
            renderTable(DATA.nodes);
            document.getElementById('global-clustering').innerText = calculateAvgClustering().toFixed(3);
            
            // Search Listener
            document.getElementById('table-search').addEventListener('input', (e) => {{
                const term = e.target.value.toLowerCase();
                const filtered = DATA.nodes.filter(n => 
                    n.id.toLowerCase().includes(term) || 
                    n.meta.name.toLowerCase().includes(term) ||
                    n.meta.organism.toLowerCase().includes(term)
                );
                renderTable(filtered);
            }});
        }}
        
        // --- Tab Switching ---
        function switchTab(tab) {{
            ['network', 'table'].forEach(t => {{
                document.getElementById(`tab-${{t}}`).classList.toggle('hidden', t !== tab);
                const btn = document.getElementById(`tab-btn-${{t}}`);
                if (t === tab) {{
                    btn.className = 'tab-active px-4 py-2';
                }} else {{
                    btn.className = 'tab-inactive px-4 py-2';
                }}
            }});
            if(tab === 'network') Plotly.relayout('plotly-graph', {{autosize: true}});
        }}
        
        // --- Graph Rendering ---
        function renderGraph() {{
            const showLabels = document.getElementById('show-labels').checked;
            
            const colorEdges = document.getElementById('color-edges').checked;
            
            const traces = [];
            
            // Edge Traces - Pseudo-Gradient using Dynamic Bins
            if (colorEdges) {{
                // 1. Calculate Min/Max Weight dynamically
                let minW = 1.0, maxW = 0.0;
                DATA.edges.forEach(e => {{
                     if (e.weight < minW) minW = e.weight;
                     if (e.weight > maxW) maxW = e.weight;
                }});
                
                // Avoid divide by zero if all weights are same
                if (Math.abs(maxW - minW) < 0.001) {{
                     minW = 0; maxW = 1; 
                }}

                // 2. Create Bins (e.g., 10 steps)
                const numBins = 10;
                const bins = [];
                for(let i=0; i<numBins; i++) {{
                     const t = i / (numBins - 1); // 0.0 to 1.0
                     // Interpolate color: Blue (#60a5fa) to Red (#ef4444)
                     // Simple approach: mix RGB
                     // Blue: [96, 165, 250], Red: [239, 68, 68]
                     const r = Math.round(96 + (239 - 96) * t);
                     const g = Math.round(165 + (68 - 165) * t);
                     const b = Math.round(250 + (68 - 250) * t);
                     const color = `rgb(${{r}},${{g}},${{b}})`;
                     
                     // Range for this bin
                     const wStart = minW + t * (maxW - minW);
                     // Make bins cover the range slightly overlapping or continuous
                     // We actually just assign edges to nearest bin or bucket them
                     // Let's bucket them
                     const bucketMin = minW + (i / numBins) * (maxW - minW);
                     const bucketMax = minW + ((i + 1) / numBins) * (maxW - minW);
                     
                     bins.push({{
                         min: bucketMin, 
                         max: i === numBins - 1 ? 1.1 : bucketMax, // Ensure last bin catches max
                         color: color,
                         width: 1.5 + t, // Thicker as it gets stronger
                         opacity: 0.5 + (0.5 * t) // More opaque as strong
                     }});
                }}
                
                bins.forEach(bin => {{
                    const binX = [], binY = [];
                    DATA.edges.forEach(edge => {{
                        // Check if edge falls in this bin
                        if (edge.weight >= bin.min && edge.weight < bin.max) {{
                            const u = DATA.nodes.find(n => n.id === edge.source);
                            const v = DATA.nodes.find(n => n.id === edge.target);
                            if(u && v) {{
                                binX.push(u.x, v.x, null);
                                binY.push(u.y, v.y, null);
                            }}
                        }}
                    }});
                    if(binX.length > 0) {{
                        traces.push({{
                            x: binX, y: binY, mode: 'lines',
                            line: {{ width: bin.width, color: bin.color }},
                            hoverinfo: 'none',
                            type: 'scatter',
                            opacity: bin.opacity
                        }});
                    }}
                }});
            }} else {{
                // Default Single Trace
                const edgeX = [], edgeY = [];
                DATA.edges.forEach(edge => {{
                    const u = DATA.nodes.find(n => n.id === edge.source);
                    const v = DATA.nodes.find(n => n.id === edge.target);
                    if(u && v) {{
                        edgeX.push(u.x, v.x, null);
                        edgeY.push(u.y, v.y, null);
                    }}
                }});
                traces.push({{
                    x: edgeX, y: edgeY, mode: 'lines',
                    line: {{ width: 1, color: '#e2e8f0' }},
                    hoverinfo: 'none',
                    type: 'scatter'
                }});
            }}
            
            const nodeTrace = {{
                x: DATA.nodes.map(n => n.x),
                y: DATA.nodes.map(n => n.y),
                text: DATA.nodes.map(n => n.meta.name),
                mode: showLabels ? 'markers+text' : 'markers',
                textposition: 'top center',
                hoverinfo: 'text',
                marker: {{
                    size: DATA.nodes.map(n => n.is_query ? 20 : 10 + (n.similarity * 5)),
                    color: DATA.nodes.map(n => n.is_query ? '#ef4444' : '#0ea5e9'),
                    line: {{ width: 1, color: 'white' }}
                }},
                type: 'scatter'
            }};
            traces.push(nodeTrace);
            
            const layout = {{
                showlegend: false,
                hovermode: 'closest',
                margin: {{l:0, r:0, t:0, b:0}},
                xaxis: {{visible: false}},
                yaxis: {{visible: false}}
            }};
            
            Plotly.newPlot('plotly-graph', traces, layout, {{responsive: true}});
            
            document.getElementById('plotly-graph').on('plotly_click', function(data){{
                const pt = data.points[0];
                const nodeIdx = pt.pointIndex;
                const node = DATA.nodes[nodeIdx];
                showNodeDetails(node);
            }});
        }}
        
        function updateGraph() {{
            renderGraph();
        }}
        
        // --- Sidebar Details ---
        function isValidUniprotId(id) {{
            if (!id || typeof id !== 'string') return false;
            // Reject 'query_protein' or similar placeholders
            if (id.toLowerCase().includes('query')) return false;
            // Basic UniProt ID regex (alphanumeric 6-12 chars)
            return /^[A-Z0-9]{{6,12}}$/i.test(id);
        }}

        function showNodeDetails(node) {{
            const html = `
                <div class="space-y-4">
                    <div>
                        <div class="flex items-center gap-2 mb-1">
                            <span class="text-xs font-bold px-2 py-0.5 rounded ${{node.is_query ? 'bg-red-100 text-red-700' : 'bg-blue-100 text-blue-700'}}">
                                ${{node.is_query ? 'QUERY INPUT' : 'HIT'}}
                            </span>
                            <span class="text-xs text-slate-400 font-mono">${{node.id}}</span>
                        </div>
                        <h3 class="font-bold text-slate-800 leading-snug">${{node.meta.name}}</h3>
                        <p class="text-sm text-slate-600 italic mt-1">${{node.meta.organism}}</p>
                    </div>
                    
                    <div class="grid grid-cols-2 gap-x-4 gap-y-2 text-sm">
                        <div class="text-slate-500">Similarity</div>
                        <div class="font-mono text-slate-900">${{ (node.similarity*100).toFixed(1) }}%</div>
                        
                        <div class="text-slate-500">Length</div>
                        <div class="text-slate-900">${{node.meta.length}} aa</div>
                        
                        ${{node.sequence ? `
                        <div class="col-span-2 mt-2">
                             <div class="text-xs font-bold text-slate-500 uppercase mb-1">Sequence Preview</div>
                             <div class="bg-slate-100 p-2 rounded font-mono text-xs break-all border border-slate-200" style="max-height: 80px; overflow-y: auto;">
                                 ${{node.sequence}}
                             </div>
                        </div>
                        ` : ''}}
                        
                        <div class="text-slate-500">Gene</div>
                        <div class="text-slate-900 truncate" title="${{node.meta.gene}}">${{node.meta.gene || '-'}}</div>
                        
                        ${{ isValidUniprotId(node.id) ? `
                        <div class="text-slate-500">UniProt</div>
                        <a href="https://www.uniprot.org/uniprotkb/${{node.id}}" target="_blank" class="text-blue-600 hover:underline flex items-center gap-1">
                            View Entry 
                            <svg class="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path></svg>
                        </a>
                        ` : `
                        <div class="text-slate-500">Source</div>
                        <div class="text-slate-700 italic">Sequence Input</div>
                        ` }}
                    </div>
                    
                    <div class="bg-slate-50 p-3 rounded text-sm border border-slate-100">
                        <strong class="block text-slate-500 text-xs uppercase mb-1">Function</strong>
                        <p class="text-slate-700 line-clamp-6" title="${{node.meta.function}}">${{node.meta.function || 'No data available'}}</p>
                    </div>
                </div>
            `;
            document.getElementById('node-details').innerHTML = html;
        }}
        
        // --- Table Rendering ---
        function renderTable(nodes) {{
            const tbody = document.getElementById('table-body');
            
            // Sort
            nodes.sort((a, b) => {{
                let valA = a[currentSort.field] || a.meta[currentSort.field];
                let valB = b[currentSort.field] || b.meta[currentSort.field];
                
                if (typeof valA === 'string') valA = valA.toLowerCase();
                if (typeof valB === 'string') valB = valB.toLowerCase();
                
                if (valA < valB) return currentSort.dir === 'asc' ? -1 : 1;
                if (valA > valB) return currentSort.dir === 'asc' ? 1 : -1;
                return 0;
            }});
            
            tbody.innerHTML = nodes.map(n => `
                <tr class="group transition-colors">
                    <td class="px-4 py-3 font-mono font-medium text-slate-600 group-hover:text-blue-600">${{n.id}}</td>
                    <td class="px-4 py-3 font-medium text-slate-900">${{n.meta.name}}</td>
                    <td class="px-4 py-3 text-slate-600 italic">${{n.meta.organism}}</td>
                    <td class="px-4 py-3">
                        <div class="flex items-center gap-2">
                            <div class="w-16 h-1.5 bg-slate-100 rounded-full overflow-hidden">
                                <div class="h-full bg-blue-500" style="width: ${{ (n.similarity*100) }}%"></div>
                            </div>
                            <span class="font-mono text-xs">${{ (n.similarity).toFixed(3) }}</span>
                        </div>
                    </td>
                    <td class="px-4 py-3 text-slate-600">${{n.meta.length}}</td>
                    <td class="px-4 py-3 text-slate-600 text-xs truncate max-w-[100px]" title="${{n.meta.gene}}">${{n.meta.gene || '-'}}</td>
                    <td class="px-4 py-3">
                        <button onclick="viewInNetwork('${{n.id}}')" class="text-xs bg-slate-100 hover:bg-slate-200 text-slate-700 px-2 py-1 rounded">View</button>
                    </td>
                </tr>
            `).join('');
        }}
        
        function sortTable(field) {{
            if (currentSort.field === field) {{
                currentSort.dir = currentSort.dir === 'asc' ? 'desc' : 'asc';
            }} else {{
                currentSort.field = field;
                currentSort.dir = 'desc';
            }}
            renderTable(DATA.nodes);
        }}
        
        function viewInNetwork(id) {{
            switchTab('network');
            const node = DATA.nodes.find(n => n.id === id);
            if(node) showNodeDetails(node);
        }}
        
        function calculateAvgClustering() {{
            const sum = DATA.nodes.reduce((acc, n) => acc + (n.cluster_coeff || 0), 0);
            return sum / DATA.nodes.length;
        }}
        
        // Start
        init();
    </script>
</body>
</html>
"""
        with open(path, 'w', encoding='utf-8') as f:
            f.write(html)
            
        logger.info(f"Dashboard saved to {path}")
        return path

