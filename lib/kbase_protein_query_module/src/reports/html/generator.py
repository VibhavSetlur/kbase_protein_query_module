"""
Comprehensive HTML Report Generator for KBase Protein Query Module

This module generates comprehensive HTML reports with:
- Multi-layer design (dashboard + detail layers)
- Interactive tables and visualizations
- KBase UI compliance
- Plotly/D3.js integration
- Single and multi-protein analysis integration
"""

import os
import time
import json
import logging
from typing import Dict, List, Any, Optional, Tuple
import numpy as np
import pandas as pd
from pathlib import Path

from ...analysis.sequence_analyzer import ProteinSequenceAnalyzer
from ...processing.networks.builder import DynamicNetworkBuilder

logger = logging.getLogger(__name__)

class HTMLReportGenerator:
    """
    Generates comprehensive HTML reports with multi-layer design.
    
    Features:
    - Dashboard layer with aggregated visualizations
    - Detail layer for individual protein analysis
    - Multi-protein analysis integration
    - Interactive tables and charts
    - KBase UI compliance
    """
    
    def __init__(self, output_dir: str = None):
        """Initialize the HTML report generator."""
        self.output_dir = output_dir or "html_reports"
        self.sequence_analyzer = ProteinSequenceAnalyzer()
        self.network_builder = DynamicNetworkBuilder()
        
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
    
    def generate_comprehensive_report(self, 
                                   pipeline_results: Dict[str, Any],
                                   protein_id: str = None,
                                   sequence: str = None) -> Dict[str, Any]:
        """
        Generate a comprehensive HTML report with multi-layer design.
        
        Args:
            pipeline_results: Results from the protein analysis pipeline
            protein_id: Protein identifier
            sequence: Protein sequence string
            
        Returns:
            Dictionary containing report information and file paths
        """
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # Extract pipeline data
        pipeline_data = pipeline_results.get('pipeline_data', {})
        stage_results = pipeline_results.get('results', {})
        
        # Generate dashboard data
        dashboard_data = self._generate_dashboard_data(pipeline_data, stage_results)
        
        # Generate HTML content
        html_content = self._generate_html_content(
            pipeline_results, dashboard_data, protein_id, timestamp
        )
        
        # Save HTML file
        report_filename = f"protein_analysis_report_{int(time.time())}.html"
        report_path = os.path.join(self.output_dir, report_filename)
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # Generate network visualizations
        network_viz_path = None
        if 'network_analysis' in stage_results:
            try:
                network_viz_path = self._generate_network_visualizations(
                    stage_results['network_analysis'], protein_id
                )
            except Exception as e:
                logger.warning(f"Could not generate network visualizations: {e}")
        
        return {
            'html_path': report_path,
            'network_viz_path': network_viz_path,
            'dashboard_data': dashboard_data,
            'timestamp': timestamp
        }
    
    def _generate_dashboard_data(self, pipeline_data: Dict[str, Any], 
                               stage_results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate data for the dashboard layer."""
        dashboard_data = {
            'summary_stats': {},
            'protein_table': [],
            'family_distribution': {},
            'confidence_distribution': {},
            'sequence_length_distribution': {},
            'similarity_distribution': {},
            'network_stats': {},
            'multi_protein_stats': {}
        }
        
        # Extract protein records
        protein_records = pipeline_data.get('protein_records', [])
        if protein_records:
            # Generate protein table
            for record in protein_records:
                protein_info = {
                    'protein_id': record.protein_id,
                    'sequence_length': len(record.sequence),
                    'source': record.source,
                    'family_id': 'Unknown',
                    'confidence': 0.0,
                    'similar_proteins_count': 0,
                    'network_size': 0
                }
                
                # Add family assignment if available
                family_assignments = pipeline_data.get('family_assignments', {})
                if record.protein_id in family_assignments:
                    assignment = family_assignments[record.protein_id]
                    if assignment.get('status') == 'success':
                        protein_info['family_id'] = assignment.get('family_id', 'Unknown')
                        protein_info['confidence'] = assignment.get('confidence', 0.0)
                
                # Add similarity search results if available
                similarity_results = pipeline_data.get('similarity_results', {})
                if record.protein_id in similarity_results:
                    sim_result = similarity_results[record.protein_id]
                    if sim_result.get('status') == 'success':
                        protein_info['similar_proteins_count'] = sim_result.get('total_found', 0)
                
                # Add network analysis results if available
                networks = pipeline_data.get('networks', {})
                if record.protein_id in networks:
                    network_data = networks[record.protein_id]
                    if network_data.get('status') == 'success':
                        network_props = network_data.get('network_properties', {})
                        protein_info['network_size'] = network_props.get('node_count', 0)
                
                dashboard_data['protein_table'].append(protein_info)
            
            # Calculate summary statistics
            dashboard_data['summary_stats'] = {
                'total_proteins': len(protein_records),
                'avg_sequence_length': np.mean([len(r.sequence) for r in protein_records]),
                'total_families': len(set(info['family_id'] for info in dashboard_data['protein_table'] if info['family_id'] != 'Unknown')),
                'avg_confidence': np.mean([info['confidence'] for info in dashboard_data['protein_table']]),
                'avg_similar_proteins': np.mean([info['similar_proteins_count'] for info in dashboard_data['protein_table']]),
                'avg_network_size': np.mean([info['network_size'] for info in dashboard_data['protein_table']])
            }
            
            # Generate distributions
            dashboard_data['family_distribution'] = self._calculate_family_distribution(dashboard_data['protein_table'])
            dashboard_data['confidence_distribution'] = self._calculate_confidence_distribution(dashboard_data['protein_table'])
            dashboard_data['sequence_length_distribution'] = self._calculate_sequence_length_distribution(protein_records)
            dashboard_data['similarity_distribution'] = self._calculate_similarity_distribution(similarity_results)
        
        # Add network statistics
        if 'network_analysis' in stage_results:
            network_stats = stage_results['network_analysis'].get('output_data', {}).get('network_stats', {})
            dashboard_data['network_stats'] = {
                'total_networks': network_stats.get('total_networks', 0),
                'successful_networks': network_stats.get('successful_networks', 0),
                'avg_node_count': network_stats.get('network_properties', {}).get('node_counts', {}).get('mean', 0),
                'avg_edge_count': network_stats.get('network_properties', {}).get('edge_counts', {}).get('mean', 0),
                'avg_density': network_stats.get('network_properties', {}).get('densities', {}).get('mean', 0)
            }
        
        # Add multi-protein analysis statistics
        if 'multi_protein_analysis' in stage_results:
            multi_protein_data = stage_results['multi_protein_analysis'].get('output_data', {}).get('multi_protein_analysis', {})
            dashboard_data['multi_protein_stats'] = {
                'alignment_length': multi_protein_data.get('metadata', {}).get('alignment_length', 0),
                'conservation_score': multi_protein_data.get('metadata', {}).get('conservation_score', 0),
                'num_proteins': multi_protein_data.get('metadata', {}).get('num_proteins', 0)
            }
        
        return dashboard_data
    
    def _calculate_family_distribution(self, protein_table: List[Dict]) -> Dict[str, int]:
        """Calculate family distribution."""
        family_counts = {}
        for protein in protein_table:
            family_id = protein['family_id']
            family_counts[family_id] = family_counts.get(family_id, 0) + 1
        return family_counts
    
    def _calculate_confidence_distribution(self, protein_table: List[Dict]) -> Dict[str, int]:
        """Calculate confidence score distribution."""
        confidence_ranges = {
            '0.0-0.2': 0, '0.2-0.4': 0, '0.4-0.6': 0, '0.6-0.8': 0, '0.8-1.0': 0
        }
        
        for protein in protein_table:
            confidence = protein['confidence']
            if confidence < 0.2:
                confidence_ranges['0.0-0.2'] += 1
            elif confidence < 0.4:
                confidence_ranges['0.2-0.4'] += 1
            elif confidence < 0.6:
                confidence_ranges['0.4-0.6'] += 1
            elif confidence < 0.8:
                confidence_ranges['0.6-0.8'] += 1
            else:
                confidence_ranges['0.8-1.0'] += 1
        
        return confidence_ranges
    
    def _calculate_sequence_length_distribution(self, protein_records: List) -> Dict[str, int]:
        """Calculate sequence length distribution."""
        length_ranges = {
            '0-100': 0, '100-200': 0, '200-500': 0, '500-1000': 0, '1000+': 0
        }
        
        for record in protein_records:
            length = len(record.sequence)
            if length <= 100:
                length_ranges['0-100'] += 1
            elif length <= 200:
                length_ranges['100-200'] += 1
            elif length <= 500:
                length_ranges['200-500'] += 1
            elif length <= 1000:
                length_ranges['500-1000'] += 1
            else:
                length_ranges['1000+'] += 1
        
        return length_ranges
    
    def _calculate_similarity_distribution(self, similarity_results: Dict) -> Dict[str, int]:
        """Calculate similarity score distribution."""
        similarity_ranges = {
            '0.0-0.2': 0, '0.2-0.4': 0, '0.4-0.6': 0, '0.6-0.8': 0, '0.8-1.0': 0
        }
        
        all_similarities = []
        for result in similarity_results.values():
            if result.get('status') == 'success':
                for protein in result.get('similar_proteins', []):
                    all_similarities.append(protein.get('similarity', 0))
        
        for similarity in all_similarities:
            if similarity < 0.2:
                similarity_ranges['0.0-0.2'] += 1
            elif similarity < 0.4:
                similarity_ranges['0.2-0.4'] += 1
            elif similarity < 0.6:
                similarity_ranges['0.4-0.6'] += 1
            elif similarity < 0.8:
                similarity_ranges['0.6-0.8'] += 1
            else:
                similarity_ranges['0.8-1.0'] += 1
        
        return similarity_ranges
    
    def _generate_html_content(self, pipeline_results: Dict[str, Any],
                             dashboard_data: Dict[str, Any],
                             protein_id: str, timestamp: str) -> str:
        """Generate the complete HTML content."""
        
        # Generate dashboard layer
        dashboard_html = self._generate_dashboard_layer(dashboard_data)
        
        # Generate detail layer templates
        detail_templates = self._generate_detail_layer_templates(pipeline_results)
        
        # Generate multi-protein analysis section
        multi_protein_section = self._generate_multi_protein_section(pipeline_results)
        
        # Generate chart scripts
        chart_scripts = self._generate_chart_scripts(dashboard_data)
        
        # Complete HTML structure
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Comprehensive Protein Query Analysis Report</title>
    
    <!-- KBase UI Styles -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.7.2/font/bootstrap-icons.css" rel="stylesheet">
    
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    
    <!-- DataTables -->
    <link href="https://cdn.datatables.net/1.11.5/css/dataTables.bootstrap5.min.css" rel="stylesheet">
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.11.5/js/dataTables.bootstrap5.min.js"></script>
    
    <!-- D3.js for advanced visualizations -->
    <script src="https://d3js.org/d3.v7.min.js"></script>
    
    <style>
        .dashboard-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 15px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }}
        
        .stat-card {{
            background: white;
            border-radius: 10px;
            padding: 15px;
            margin-bottom: 15px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
            border-left: 4px solid #667eea;
        }}
        
        .chart-container {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        
        .protein-table {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        
        .multi-protein-section {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        
        .detail-modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.5);
        }}
        
        .detail-content {{
            background-color: white;
            margin: 5% auto;
            padding: 20px;
            border-radius: 10px;
            width: 90%;
            max-width: 1200px;
            max-height: 80vh;
            overflow-y: auto;
        }}
        
        .close {{
            color: #aaa;
            float: right;
            font-size: 28px;
            font-weight: bold;
            cursor: pointer;
        }}
        
        .close:hover {{
            color: black;
        }}
        
        .tab-content {{
            padding: 20px 0;
        }}
        
        .nav-tabs .nav-link {{
            color: #667eea;
        }}
        
        .nav-tabs .nav-link.active {{
            background-color: #667eea;
            color: white;
            border-color: #667eea;
        }}
        
        .bioinformatics-links {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }}
        
        .bioinformatics-link {{
            background: #f8f9fa;
            padding: 10px;
            border-radius: 5px;
            text-decoration: none;
            color: #333;
            border: 1px solid #dee2e6;
            transition: all 0.3s ease;
        }}
        
        .bioinformatics-link:hover {{
            background: #667eea;
            color: white;
            text-decoration: none;
        }}
        
        .conservation-bar {{
            height: 20px;
            background: linear-gradient(90deg, #ff6b6b 0%, #feca57 50%, #48dbfb 100%);
            border-radius: 10px;
            margin: 5px 0;
        }}
        
        .tree-container {{
            width: 100%;
            height: 400px;
            border: 1px solid #dee2e6;
            border-radius: 5px;
            overflow: auto;
        }}
    </style>
</head>
<body>
    <div class="container-fluid">
        <!-- Header -->
        <div class="row">
            <div class="col-12">
                <div class="dashboard-card">
                    <h1><i class="bi bi-graph-up"></i> Comprehensive Protein Query Analysis Report</h1>
                    <p class="mb-0">Generated on {timestamp}</p>
                    <p class="mb-0">Analysis ID: {protein_id or 'Bulk Analysis'}</p>
                </div>
            </div>
        </div>
        
        <!-- Main Navigation -->
        <div class="row">
            <div class="col-12">
                <ul class="nav nav-tabs" id="mainTabs" role="tablist">
                    <li class="nav-item" role="presentation">
                        <button class="nav-link active" id="dashboard-tab" data-bs-toggle="tab" data-bs-target="#dashboard" type="button" role="tab">
                            <i class="bi bi-speedometer2"></i> Dashboard
                        </button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="multi-protein-tab" data-bs-toggle="tab" data-bs-target="#multi-protein" type="button" role="tab">
                            <i class="bi bi-diagram-3"></i> Multi-Protein Analysis
                        </button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="bioinformatics-tab" data-bs-toggle="tab" data-bs-target="#bioinformatics" type="button" role="tab">
                            <i class="bi bi-link-45deg"></i> Bioinformatics Links
                        </button>
                    </li>
                </ul>
            </div>
        </div>
        
        <!-- Tab Content -->
        <div class="tab-content" id="mainTabContent">
            <!-- Dashboard Tab -->
            <div class="tab-pane fade show active" id="dashboard" role="tabpanel">
                {dashboard_html}
            </div>
            
            <!-- Multi-Protein Analysis Tab -->
            <div class="tab-pane fade" id="multi-protein" role="tabpanel">
                {multi_protein_section}
            </div>
            
            <!-- Bioinformatics Links Tab -->
            <div class="tab-pane fade" id="bioinformatics" role="tabpanel">
                {self._generate_bioinformatics_links_section(pipeline_results)}
            </div>
        </div>
        
        <!-- Detail Layer Modals -->
        {detail_templates}
    </div>
    
    <!-- Scripts -->
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    {chart_scripts}
    
    <script>
        // Initialize DataTable
        $(document).ready(function() {{
            $('#proteinTable').DataTable({{
                pageLength: 25,
                order: [[1, 'desc']],
                responsive: true,
                dom: 'Bfrtip',
                buttons: ['copy', 'csv', 'excel']
            }});
        }});
        
        // Detail layer functionality
        function showProteinDetail(proteinId) {{
            const modal = document.getElementById('detailModal');
            const content = document.getElementById('detailContent');
            
            // Load protein-specific content
            content.innerHTML = document.getElementById('detailTemplate_' + proteinId).innerHTML;
            
            modal.style.display = 'block';
        }}
        
        function closeDetailModal() {{
            document.getElementById('detailModal').style.display = 'none';
        }}
        
        // Close modal when clicking outside
        window.onclick = function(event) {{
            const modal = document.getElementById('detailModal');
            if (event.target == modal) {{
                modal.style.display = 'none';
            }}
        }}
        
        // Initialize phylogenetic tree visualization
        function initializePhylogeneticTree(treeData) {{
            if (!treeData || !treeData.tree) return;
            
            const container = document.getElementById('phylogeneticTree');
            if (!container) return;
            
            // Create tree visualization using D3.js
            const width = container.clientWidth;
            const height = 400;
            
            const tree = d3.tree().size([height, width - 160]);
            const root = d3.hierarchy(treeData.tree);
            
            const svg = d3.select(container)
                .append('svg')
                .attr('width', width)
                .attr('height', height);
            
            const g = svg.append('g')
                .attr('transform', 'translate(80,0)');
            
            const link = g.selectAll('.link')
                .data(tree(root).links())
                .enter().append('path')
                .attr('class', 'link')
                .attr('d', d3.linkHorizontal()
                    .x(d => d.y)
                    .y(d => d.x));
            
            const node = g.selectAll('.node')
                .data(root.descendants())
                .enter().append('g')
                .attr('class', 'node')
                .attr('transform', d => `translate(${{d.y}},${{d.x}})`);
            
            node.append('circle')
                .attr('r', 3);
            
            node.append('text')
                .attr('dy', '.31em')
                .attr('x', d => d.children ? -8 : 8)
                .style('text-anchor', d => d.children ? 'end' : 'start')
                .text(d => d.data.name);
        }}
    </script>
</body>
</html>
        """
        
        return html_content
    
    def _generate_dashboard_layer(self, dashboard_data: Dict[str, Any]) -> str:
        """Generate the dashboard layer HTML."""
        
        # Summary statistics cards
        stats_html = ""
        if dashboard_data['summary_stats']:
            stats = dashboard_data['summary_stats']
            multi_protein_stats = dashboard_data.get('multi_protein_stats', {})
            
            stats_html = f"""
            <div class="row mb-4">
                <div class="col-md-2">
                    <div class="stat-card">
                        <h4>{stats.get('total_proteins', 0)}</h4>
                        <p class="text-muted mb-0">Total Proteins</p>
                    </div>
                </div>
                <div class="col-md-2">
                    <div class="stat-card">
                        <h4>{stats.get('total_families', 0)}</h4>
                        <p class="text-muted mb-0">Protein Families</p>
                    </div>
                </div>
                <div class="col-md-2">
                    <div class="stat-card">
                        <h4>{stats.get('avg_sequence_length', 0):.0f}</h4>
                        <p class="text-muted mb-0">Avg Length</p>
                    </div>
                </div>
                <div class="col-md-2">
                    <div class="stat-card">
                        <h4>{stats.get('avg_confidence', 0):.2f}</h4>
                        <p class="text-muted mb-0">Avg Confidence</p>
                    </div>
                </div>
                <div class="col-md-2">
                    <div class="stat-card">
                        <h4>{multi_protein_stats.get('conservation_score', 0):.2f}</h4>
                        <p class="text-muted mb-0">Conservation</p>
                    </div>
                </div>
                <div class="col-md-2">
                    <div class="stat-card">
                        <h4>{multi_protein_stats.get('alignment_length', 0)}</h4>
                        <p class="text-muted mb-0">Alignment Length</p>
                    </div>
                </div>
            </div>
            """
        
        # Protein table
        table_html = ""
        if dashboard_data['protein_table']:
            table_rows = ""
            for protein in dashboard_data['protein_table']:
                table_rows += f"""
                <tr onclick="showProteinDetail('{protein['protein_id']}')" style="cursor: pointer;">
                    <td>{protein['protein_id']}</td>
                    <td>{protein['sequence_length']}</td>
                    <td>{protein['family_id']}</td>
                    <td>{protein['confidence']:.3f}</td>
                    <td>{protein['similar_proteins_count']}</td>
                    <td>{protein['network_size']}</td>
                    <td>{protein['source']}</td>
                </tr>
                """
            
            table_html = f"""
            <div class="protein-table">
                <h3><i class="bi bi-table"></i> Protein Analysis Summary</h3>
                <p>Click on any row to view detailed analysis for that protein.</p>
                <table id="proteinTable" class="table table-striped table-hover">
                    <thead>
                        <tr>
                            <th>Protein ID</th>
                            <th>Length</th>
                            <th>Family</th>
                            <th>Confidence</th>
                            <th>Similar Proteins</th>
                            <th>Network Size</th>
                            <th>Source</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </div>
            """
        
        # Chart containers
        charts_html = """
        <div class="row">
            <div class="col-md-6">
                <div class="chart-container">
                    <h4><i class="bi bi-pie-chart"></i> Family Distribution</h4>
                    <div id="familyChart"></div>
                </div>
            </div>
            <div class="col-md-6">
                <div class="chart-container">
                    <h4><i class="bi bi-bar-chart"></i> Confidence Distribution</h4>
                    <div id="confidenceChart"></div>
                </div>
            </div>
        </div>
        <div class="row">
            <div class="col-md-6">
                <div class="chart-container">
                    <h4><i class="bi bi-bar-chart"></i> Sequence Length Distribution</h4>
                    <div id="lengthChart"></div>
                </div>
            </div>
            <div class="col-md-6">
                <div class="chart-container">
                    <h4><i class="bi bi-bar-chart"></i> Similarity Distribution</h4>
                    <div id="similarityChart"></div>
                </div>
            </div>
        </div>
        """
        
        return stats_html + table_html + charts_html
    
    def _generate_detail_layer_templates(self, pipeline_results: Dict[str, Any]) -> str:
        """Generate detail layer templates for each protein."""
        pipeline_data = pipeline_results.get('pipeline_data', {})
        protein_records = pipeline_data.get('protein_records', [])
        
        templates = ""
        for record in protein_records:
            protein_id = record.protein_id
            template = self._generate_single_protein_detail_template(protein_id, record, pipeline_data)
            templates += f"""
            <div id="detailTemplate_{protein_id}" style="display: none;">
                {template}
            </div>
            """
        
        # Add modal container
        modal_html = """
        <div id="detailModal" class="detail-modal">
            <div class="detail-content">
                <span class="close" onclick="closeDetailModal()">&times;</span>
                <div id="detailContent"></div>
            </div>
        </div>
        """
        
        return templates + modal_html
    
    def _generate_single_protein_detail_template(self, protein_id: str, record, 
                                               pipeline_data: Dict[str, Any]) -> str:
        """Generate detail template for a single protein."""
        
        # Get analysis results for this protein
        sequence_analysis = pipeline_data.get('sequence_analyses', {}).get(protein_id, {})
        family_assignment = pipeline_data.get('family_assignments', {}).get(protein_id, {})
        similarity_result = pipeline_data.get('similarity_results', {}).get(protein_id, {})
        network_data = pipeline_data.get('networks', {}).get(protein_id, {})
        
        # Generate detail content
        detail_content = f"""
        <div class="container-fluid">
            <h2><i class="bi bi-protein"></i> Detailed Analysis: {protein_id}</h2>
            
            <ul class="nav nav-tabs" id="detailTabs" role="tablist">
                <li class="nav-item" role="presentation">
                    <button class="nav-link active" id="overview-tab" data-bs-toggle="tab" data-bs-target="#overview" type="button" role="tab">Overview</button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link" id="sequence-tab" data-bs-toggle="tab" data-bs-target="#sequence" type="button" role="tab">Sequence Analysis</button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link" id="similarity-tab" data-bs-toggle="tab" data-bs-target="#similarity" type="button" role="tab">Similarity Search</button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link" id="network-tab" data-bs-toggle="tab" data-bs-target="#network" type="button" role="tab">Network Analysis</button>
                </li>
            </ul>
            
            <div class="tab-content" id="detailTabContent">
                <!-- Overview Tab -->
                <div class="tab-pane fade show active" id="overview" role="tabpanel">
                    <div class="row">
                        <div class="col-md-6">
                            <h4>Basic Information</h4>
                            <table class="table">
                                <tr><td>Protein ID:</td><td>{protein_id}</td></tr>
                                <tr><td>Sequence Length:</td><td>{len(record.sequence)}</td></tr>
                                <tr><td>Source:</td><td>{record.source}</td></tr>
                                <tr><td>Family ID:</td><td>{family_assignment.get('family_id', 'Unknown')}</td></tr>
                                <tr><td>Confidence:</td><td>{family_assignment.get('confidence', 0):.3f}</td></tr>
                            </table>
                        </div>
                        <div class="col-md-6">
                            <h4>Sequence Preview</h4>
                            <div style="background: #f8f9fa; padding: 10px; border-radius: 5px; font-family: monospace; font-size: 12px;">
                                {record.sequence[:100]}{'...' if len(record.sequence) > 100 else ''}
                            </div>
                        </div>
                    </div>
                </div>
                
                <!-- Sequence Analysis Tab -->
                <div class="tab-pane fade" id="sequence" role="tabpanel">
                    {self._generate_sequence_analysis_content(sequence_analysis)}
                </div>
                
                <!-- Similarity Search Tab -->
                <div class="tab-pane fade" id="similarity" role="tabpanel">
                    {self._generate_similarity_content(similarity_result)}
                </div>
                
                <!-- Network Analysis Tab -->
                <div class="tab-pane fade" id="network" role="tabpanel">
                    {self._generate_network_content(network_data)}
                </div>
            </div>
        </div>
        """
        
        return detail_content
    
    def _generate_sequence_analysis_content(self, sequence_analysis: Dict[str, Any]) -> str:
        """Generate sequence analysis content."""
        if not sequence_analysis or sequence_analysis.get('status') == 'error':
            return "<p>Sequence analysis not available.</p>"
        
        content = "<div class='row'>"
        
        # Amino acid composition
        if 'amino_acid_composition' in sequence_analysis:
            aa_comp = sequence_analysis['amino_acid_composition']
            content += """
            <div class="col-md-6">
                <h4>Amino Acid Composition</h4>
                <div id="aaCompositionChart"></div>
            </div>
            """
        
        # Physicochemical properties
        if 'physicochemical_properties' in sequence_analysis:
            props = sequence_analysis['physicochemical_properties']
            content += f"""
            <div class="col-md-6">
                <h4>Physicochemical Properties</h4>
                <table class="table">
                    <tr><td>Molecular Weight:</td><td>{props.get('molecular_weight', 0):.2f}</td></tr>
                    <tr><td>Isoelectric Point:</td><td>{props.get('isoelectric_point', 0):.2f}</td></tr>
                    <tr><td>Instability Index:</td><td>{props.get('instability_index', 0):.2f}</td></tr>
                    <tr><td>Aliphatic Index:</td><td>{props.get('aliphatic_index', 0):.2f}</td></tr>
                    <tr><td>Extinction Coefficient:</td><td>{props.get('extinction_coefficient', 0):.2f}</td></tr>
                </table>
            </div>
            """
        
        content += "</div>"
        return content
    
    def _generate_similarity_content(self, similarity_result: Dict[str, Any]) -> str:
        """Generate similarity search content."""
        if not similarity_result or similarity_result.get('status') != 'success':
            return "<p>Similarity search results not available.</p>"
        
        similar_proteins = similarity_result.get('similar_proteins', [])
        
        if not similar_proteins:
            return "<p>No similar proteins found.</p>"
        
        table_rows = ""
        for protein in similar_proteins[:20]:  # Show top 20
            table_rows += f"""
            <tr>
                <td>{protein['protein_id']}</td>
                <td>{protein['similarity']:.3f}</td>
                <td>{protein['family_id']}</td>
                <td>{protein['rank']}</td>
            </tr>
            """
        
        return f"""
        <h4>Top Similar Proteins</h4>
        <p>Found {len(similar_proteins)} similar proteins</p>
        <table class="table table-striped">
            <thead>
                <tr>
                    <th>Protein ID</th>
                    <th>Similarity</th>
                    <th>Family</th>
                    <th>Rank</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>
        """
    
    def _generate_network_content(self, network_data: Dict[str, Any]) -> str:
        """Generate network analysis content."""
        if not network_data or network_data.get('status') != 'success':
            return "<p>Network analysis results not available.</p>"
        
        network_props = network_data.get('network_properties', {})
        
        return f"""
        <h4>Network Properties</h4>
        <div class="row">
            <div class="col-md-6">
                <table class="table">
                    <tr><td>Node Count:</td><td>{network_props.get('node_count', 0)}</td></tr>
                    <tr><td>Edge Count:</td><td>{network_props.get('edge_count', 0)}</td></tr>
                    <tr><td>Density:</td><td>{network_props.get('density', 0):.3f}</td></tr>
                    <tr><td>Diameter:</td><td>{network_props.get('diameter', 0)}</td></tr>
                    <tr><td>Clustering Coefficient:</td><td>{network_props.get('clustering_coefficient', 0):.3f}</td></tr>
                </table>
            </div>
            <div class="col-md-6">
                <div id="networkChart"></div>
            </div>
        </div>
        """
    
    def _generate_chart_scripts(self, dashboard_data: Dict[str, Any]) -> str:
        """Generate JavaScript for interactive charts."""
        
        # Family distribution chart
        family_chart = ""
        if dashboard_data['family_distribution']:
            families = list(dashboard_data['family_distribution'].keys())
            counts = list(dashboard_data['family_distribution'].values())
            
            family_chart = f"""
            var familyData = [{{
                values: {counts},
                labels: {json.dumps(families)},
                type: 'pie',
                marker: {{
                    colors: ['#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0', '#9966FF', '#FF9F40']
                }}
            }}];
            
            var familyLayout = {{
                title: 'Protein Family Distribution',
                height: 400
            }};
            
            Plotly.newPlot('familyChart', familyData, familyLayout);
            """
        
        # Confidence distribution chart
        confidence_chart = ""
        if dashboard_data['confidence_distribution']:
            ranges = list(dashboard_data['confidence_distribution'].keys())
            counts = list(dashboard_data['confidence_distribution'].values())
            
            confidence_chart = f"""
            var confidenceData = [{{
                x: {json.dumps(ranges)},
                y: {counts},
                type: 'bar',
                marker: {{
                    color: '#36A2EB'
                }}
            }}];
            
            var confidenceLayout = {{
                title: 'Confidence Score Distribution',
                xaxis: {{ title: 'Confidence Range' }},
                yaxis: {{ title: 'Number of Proteins' }},
                height: 400
            }};
            
            Plotly.newPlot('confidenceChart', confidenceData, confidenceLayout);
            """
        
        # Sequence length distribution chart
        length_chart = ""
        if dashboard_data['sequence_length_distribution']:
            ranges = list(dashboard_data['sequence_length_distribution'].keys())
            counts = list(dashboard_data['sequence_length_distribution'].values())
            
            length_chart = f"""
            var lengthData = [{{
                x: {json.dumps(ranges)},
                y: {counts},
                type: 'bar',
                marker: {{
                    color: '#FFCE56'
                }}
            }}];
            
            var lengthLayout = {{
                title: 'Sequence Length Distribution',
                xaxis: {{ title: 'Length Range' }},
                yaxis: {{ title: 'Number of Proteins' }},
                height: 400
            }};
            
            Plotly.newPlot('lengthChart', lengthData, lengthLayout);
            """
        
        # Similarity distribution chart
        similarity_chart = ""
        if dashboard_data['similarity_distribution']:
            ranges = list(dashboard_data['similarity_distribution'].keys())
            counts = list(dashboard_data['similarity_distribution'].values())
            
            similarity_chart = f"""
            var similarityData = [{{
                x: {json.dumps(ranges)},
                y: {counts},
                type: 'bar',
                marker: {{
                    color: '#4BC0C0'
                }}
            }}];
            
            var similarityLayout = {{
                title: 'Similarity Score Distribution',
                xaxis: {{ title: 'Similarity Range' }},
                yaxis: {{ title: 'Number of Similar Proteins' }},
                height: 400
            }};
            
            Plotly.newPlot('similarityChart', similarityData, similarityLayout);
            """
        
        return f"""
        <script>
            // Chart initialization
            document.addEventListener('DOMContentLoaded', function() {{
                {family_chart}
                {confidence_chart}
                {length_chart}
                {similarity_chart}
            }});
        </script>
        """
    
    def _generate_multi_protein_section(self, pipeline_results: Dict[str, Any]) -> str:
        """Generate the multi-protein analysis section."""
        multi_protein_data = pipeline_results.get('results', {}).get('multi_protein_analysis', {}).get('output_data', {})
        
        if not multi_protein_data:
            return "<p>Multi-protein analysis results not available.</p>"
        
        content = f"""
        <div class="multi-protein-section">
            <h3><i class="bi bi-diagram-3"></i> Multi-Protein Analysis</h3>
            <p>This section provides an overview of the analysis across multiple proteins.</p>
            
            <div class="row">
                <div class="col-md-6">
                    <h4>Alignment Length</h4>
                    <p>{multi_protein_data.get('metadata', {}).get('alignment_length', 0)}</p>
                </div>
                <div class="col-md-6">
                    <h4>Conservation Score</h4>
                    <p>{multi_protein_data.get('metadata', {}).get('conservation_score', 0):.2f}</p>
                </div>
                <div class="col-md-6">
                    <h4>Number of Proteins Analyzed</h4>
                    <p>{multi_protein_data.get('metadata', {}).get('num_proteins', 0)}</p>
                </div>
            </div>
            
            <h4>Protein Alignment Visualization</h4>
            <div id="alignmentChart"></div>
            
            <h4>Phylogenetic Tree</h4>
            <div id="phylogeneticTree" class="tree-container"></div>
        </div>
        """
        
        return content
    
    def _generate_bioinformatics_links_section(self, pipeline_results: Dict[str, Any]) -> str:
        """Generate the bioinformatics links section."""
        links = [
            ("KBase Protein Atlas", "https://www.kbase.us/atlas/"),
            ("UniProt Knowledgebase", "https://www.uniprot.org/"),
            ("NCBI Protein Data Bank", "https://www.ncbi.nlm.nih.gov/protein"),
            ("Pfam", "https://pfam.xfam.org/"),
            ("InterPro", "https://www.ebi.ac.uk/interpro/"),
            ("SWISS-PROT", "https://www.uniprot.org/sprot/"),
            ("TrEMBL", "https://www.uniprot.org/trembl/")
        ]
        
        links_html = ""
        for name, url in links:
            links_html += f"""
            <a href="{url}" target="_blank" class="bioinformatics-link">
                <i class="bi bi-link-45deg"></i> {name}
            </a>
            """
        
        return f"""
        <div class="multi-protein-section">
            <h3><i class="bi bi-link-45deg"></i> Bioinformatics Links</h3>
            <p>Useful resources for protein analysis and data retrieval.</p>
            <div class="bioinformatics-links">
                {links_html}
            </div>
        </div>
        """
    
    def _generate_network_visualizations(self, network_analysis_result: Dict[str, Any], 
                                       protein_id: str) -> Optional[str]:
        """Generate network visualization files."""
        try:
            # This would generate separate network visualization files
            # For now, return None as network visualizations are embedded in the main report
            return None
        except Exception as e:
            logger.warning(f"Could not generate network visualizations: {e}")
            return None
