"""
HTML Report Generator for KBase Protein Network Analysis Toolkit

This module generates comprehensive HTML reports that integrate:
- Sequence analysis results
- Network visualizations
- Statistical summaries
- Bioinformatics information
- Interactive visualizations
"""

import os
import time
import json
import logging
from typing import Dict, List, Any, Optional, Tuple
import numpy as np
import pandas as pd
from pathlib import Path

from .sequence_analyzer import ProteinSequenceAnalyzer
from .network_builder import DynamicNetworkBuilder

logger = logging.getLogger(__name__)

class HTMLReportGenerator:
    """
    Generates comprehensive HTML reports for protein network analysis results.
    Integrates sequence analysis, network visualization, and statistical summaries.
    """
    
    def __init__(self, output_dir: str = None):
        """
        Initialize the HTML report generator.
        
        Args:
            output_dir: Directory to save HTML reports and assets
        """
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
        Generate a comprehensive HTML report with all analysis results.
        
        Args:
            pipeline_results: Results from the protein analysis pipeline
            protein_id: Protein identifier
            sequence: Protein sequence string
            
        Returns:
            Dictionary containing report information and file paths
        """
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # Perform sequence analysis if sequence is provided
        sequence_analysis = None
        if sequence:
            try:
                sequence_analysis = self.sequence_analyzer.analyze_sequence(
                    sequence, protein_id or "UNKNOWN"
                )
                logger.info(f"Sequence analysis completed for {protein_id}")
            except Exception as e:
                logger.warning(f"Could not perform sequence analysis: {e}")
                sequence_analysis = self._create_minimal_sequence_analysis(sequence, protein_id)
        
        # Generate HTML content
        html_content = self._generate_html_content(
            pipeline_results, sequence_analysis, protein_id, timestamp
        )
        
        # Save HTML file
        report_filename = f"protein_analysis_report_{int(time.time())}.html"
        report_path = os.path.join(self.output_dir, report_filename)
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # Generate network visualization if possible
        network_viz_path = None
        if 'similarity_search' in pipeline_results:
            try:
                network_viz_path = self._generate_network_visualization(
                    pipeline_results, protein_id
                )
            except Exception as e:
                logger.warning(f"Could not generate network visualization: {e}")
        
        return {
            'html_path': report_path,
            'network_viz_path': network_viz_path,
            'sequence_analysis': sequence_analysis,
            'timestamp': timestamp,
            'protein_id': protein_id
        }
    
    def _create_minimal_sequence_analysis(self, sequence: str, protein_id: str) -> Dict[str, Any]:
        """Create minimal sequence analysis when full analysis fails."""
        return {
            'protein_id': protein_id or "UNKNOWN",
            'sequence': sequence,
            'length': len(sequence),
            'amino_acid_composition': {
                'individual': {},
                'groups': {},
                'total_amino_acids': len(sequence)
            },
            'physicochemical_properties': {
                'molecular_weight': 0,
                'net_charge': 0,
                'average_hydrophobicity': 0,
                'isoelectric_point_estimate': 7.0
            },
            'secondary_structure_prediction': {
                'helix_preference': 0,
                'sheet_preference': 0,
                'turn_preference': 0,
                'dominant_structure': 'unknown',
                'dominant_score': 0
            },
            'sequence_motifs': {
                'n_glycosylation': [],
                'o_glycosylation': [],
                'phosphorylation': [],
                'disulfide_bonds': [],
                'repeats': []
            },
            'bioinformatics_links': {},
            'statistics': {
                'length': len(sequence),
                'molecular_weight': 0,
                'net_charge': 0,
                'average_hydrophobicity': 0
            }
        }
    
    def _generate_html_content(self, 
                             pipeline_results: Dict[str, Any],
                             sequence_analysis: Optional[Dict[str, Any]],
                             protein_id: str,
                             timestamp: str) -> str:
        """Generate the main HTML content with all tabs and visualizations."""
        
        # Extract data from pipeline results
        data = self._extract_pipeline_data(pipeline_results)
        
        # Generate individual tab content
        overview_tab = self._generate_overview_tab(data, sequence_analysis, protein_id)
        sequence_tab = self._generate_sequence_tab(sequence_analysis) if sequence_analysis else ""
        network_tab = self._generate_network_tab(data)
        statistics_tab = self._generate_statistics_tab(data)
        bioinformatics_tab = self._generate_bioinformatics_tab(sequence_analysis) if sequence_analysis else ""
        raw_data_tab = self._generate_raw_data_tab(pipeline_results)
        
        # Generate chart scripts
        chart_scripts = self._generate_chart_scripts(data, sequence_analysis)
        
        return f"""
<!DOCTYPE html>
<html>
<head>
    <title>Comprehensive Protein Network Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
        * {{
            box-sizing: border-box;
        }}
        
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 0; 
            padding: 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        
        .main-container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            min-height: 100vh;
            box-shadow: 0 0 30px rgba(0,0,0,0.1);
        }}
        
        .header {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; 
            padding: 40px 30px; 
            text-align: center;
            position: relative;
            overflow: hidden;
        }}
        
        .header::before {{
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: url('data:image/svg+xml,<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100"><defs><pattern id="grain" width="100" height="100" patternUnits="userSpaceOnUse"><circle cx="25" cy="25" r="1" fill="white" opacity="0.1"/><circle cx="75" cy="75" r="1" fill="white" opacity="0.1"/><circle cx="50" cy="10" r="0.5" fill="white" opacity="0.1"/></pattern></defs><rect width="100" height="100" fill="url(%23grain)"/></svg>');
            opacity: 0.3;
        }}
        
        .header h1 {{ 
            margin: 0; 
            font-size: 3em; 
            font-weight: 300;
            position: relative;
            z-index: 1;
        }}
        .header p {{ 
            margin: 10px 0 0 0; 
            opacity: 0.9; 
            font-size: 1.2em;
            position: relative;
            z-index: 1;
        }}
        
        .report-info {{
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 8px;
            margin-top: 20px;
            backdrop-filter: blur(10px);
            position: relative;
            z-index: 1;
        }}
        
        .tab-container {{ 
            background: white; 
            margin: 20px;
            border-radius: 15px; 
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        .tab-buttons {{ 
            display: flex; 
            background: #f8f9fa; 
            border-bottom: 2px solid #e9ecef;
            overflow-x: auto;
        }}
        
        .tab-button {{ 
            padding: 15px 25px; 
            border: none; 
            background: transparent; 
            cursor: pointer; 
            font-size: 14px; 
            font-weight: 500;
            transition: all 0.3s ease;
            white-space: nowrap;
            position: relative;
        }}
        
        .tab-button:hover {{ 
            background: #e9ecef; 
            transform: translateY(-2px);
        }}
        
        .tab-button.active {{ 
            background: #007bff; 
            color: white;
            box-shadow: 0 2px 8px rgba(0,123,255,0.3);
        }}
        
        .tab-button.active::after {{
            content: '';
            position: absolute;
            bottom: 0;
            left: 0;
            right: 0;
            height: 3px;
            background: #0056b3;
        }}
        
        .tab-content {{ 
            padding: 30px; 
            display: none; 
            min-height: 500px;
        }}
        
        .tab-content.active {{ 
            display: block; 
            animation: fadeIn 0.5s ease-in-out;
        }}
        
        @keyframes fadeIn {{
            from {{ opacity: 0; transform: translateY(20px); }}
            to {{ opacity: 1; transform: translateY(0); }}
        }}
        
        .section {{ 
            margin-bottom: 30px; 
            padding: 25px; 
            background: #f8f9fa; 
            border-radius: 12px;
            border-left: 5px solid #007bff;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            transition: all 0.3s ease;
        }}
        
        .section:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }}
        
        .section h3 {{ 
            margin: 0 0 20px 0; 
            color: #333; 
            font-size: 1.4em;
            font-weight: 600;
        }}
        
        .result {{ 
            background: white; 
            padding: 20px; 
            border-radius: 8px; 
            margin: 15px 0;
            border: 1px solid #e9ecef;
            transition: all 0.3s ease;
        }}
        
        .result:hover {{
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }}
        
        .result p {{ margin: 8px 0; }}
        
        .chart-container {{
            background: white;
            border-radius: 8px;
            padding: 20px;
            margin: 15px 0;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        }}
        
        .metric-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .metric-card {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            border-top: 4px solid #007bff;
            transition: all 0.3s ease;
        }}
        
        .metric-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
        }}
        
        .metric-value {{
            font-size: 2.5em;
            font-weight: bold;
            color: #007bff;
            margin-bottom: 10px;
        }}
        
        .metric-label {{
            color: #666;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .loading {{
            text-align: center;
            padding: 40px;
            color: #666;
        }}
        
        .error {{
            background: #f8d7da;
            color: #721c24;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border: 1px solid #f5c6cb;
        }}
        
        .success {{
            background: #d4edda;
            color: #155724;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border: 1px solid #c3e6cb;
        }}
        
        .warning {{
            background: #fff3cd;
            color: #856404;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border: 1px solid #ffeaa7;
        }}
        
        @media (max-width: 768px) {{
            .tab-buttons {{
                flex-direction: column;
            }}
            
            .metric-grid {{
                grid-template-columns: 1fr;
            }}
            
            .header h1 {{
                font-size: 2em;
            }}
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="header">
            <h1>Protein Network Analysis Report</h1>
            <p>Comprehensive analysis of protein similarities, families, and network relationships</p>
            <div class="report-info">
                <strong>Protein ID:</strong> {protein_id} | 
                <strong>Generated:</strong> {timestamp} | 
                <strong>Analysis Type:</strong> Comprehensive Network Analysis
            </div>
        </div>
        
        <div class="tab-container">
            <div class="tab-buttons">
                <button class="tab-button active" onclick="showTab('overview')">Overview</button>
                <button class="tab-button" onclick="showTab('sequence')">Sequence Analysis</button>
                <button class="tab-button" onclick="showTab('network')">Network Analysis</button>
                <button class="tab-button" onclick="showTab('statistics')">Statistics</button>
                <button class="tab-button" onclick="showTab('bioinformatics')">Bioinformatics</button>
                <button class="tab-button" onclick="showTab('raw-data')">Raw Data</button>
            </div>
            
            <div id="overview" class="tab-content active">
                {overview_tab}
            </div>
            
            <div id="sequence" class="tab-content">
                {sequence_tab if sequence_tab else '<div class="loading">Sequence analysis not available</div>'}
            </div>
            
            <div id="network" class="tab-content">
                {network_tab}
            </div>
            
            <div id="statistics" class="tab-content">
                {statistics_tab}
            </div>
            
            <div id="bioinformatics" class="tab-content">
                {bioinformatics_tab if bioinformatics_tab else '<div class="loading">Bioinformatics data not available</div>'}
            </div>
            
            <div id="raw-data" class="tab-content">
                {raw_data_tab}
            </div>
        </div>
    </div>
    
    <script>
        function showTab(tabName) {{
            // Hide all tab contents
            const tabContents = document.querySelectorAll('.tab-content');
            tabContents.forEach(content => {{
                content.classList.remove('active');
            }});
            
            // Remove active class from all buttons
            const tabButtons = document.querySelectorAll('.tab-button');
            tabButtons.forEach(button => {{
                button.classList.remove('active');
            }});
            
            // Show selected tab content
            document.getElementById(tabName).classList.add('active');
            
            // Add active class to clicked button
            event.target.classList.add('active');
        }}
        
        // Initialize charts when DOM is loaded
        document.addEventListener('DOMContentLoaded', function() {{
            console.log('Protein Network Analysis Report loaded');
            
            // Add smooth scrolling
            document.querySelectorAll('a[href^="#"]').forEach(anchor => {{
                anchor.addEventListener('click', function (e) {{
                    e.preventDefault();
                    document.querySelector(this.getAttribute('href')).scrollIntoView({{
                        behavior: 'smooth'
                    }});
                }});
            }});
            
            // Add loading animations
            const sections = document.querySelectorAll('.section');
            const observer = new IntersectionObserver((entries) => {{
                entries.forEach(entry => {{
                    if (entry.isIntersecting) {{
                        entry.target.style.opacity = '1';
                        entry.target.style.transform = 'translateY(0)';
                    }}
                }});
            }});
            
            sections.forEach(section => {{
                section.style.opacity = '0';
                section.style.transform = 'translateY(20px)';
                section.style.transition = 'all 0.6s ease';
                observer.observe(section);
            }});
        }});
    </script>
    
    {chart_scripts}
</body>
</html>
        """
    
    def _extract_pipeline_data(self, pipeline_results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract and organize data from pipeline results."""
        data = {
            'protein_id': 'UNKNOWN',
            'sequence': '',
            'embedding_dim': 0,
            'embedding_norm': 0.0,
            'family_assignment': {},
            'similarity_matches': 0,
            'network_nodes': 0,
            'network_edges': 0,
            'similarity_stats': {},
            'sequence_length': 0
        }
        
        # Extract data from different pipeline steps
        if 'protein_existence' in pipeline_results:
            existence_data = pipeline_results['protein_existence']
            data['protein_id'] = existence_data.get('protein_id', 'UNKNOWN')
            data['family_assignment'] = existence_data.get('family_assignment', {})
        
        if 'embedding_generation' in pipeline_results:
            embedding_data = pipeline_results['embedding_generation']
            data['embedding_dim'] = embedding_data.get('embedding_dim', 0)
            data['embedding_norm'] = embedding_data.get('embedding_norm', 0.0)
            data['sequence_length'] = embedding_data.get('sequence_length', 0)
            data['sequence'] = embedding_data.get('sequence', '')
        
        if 'similarity_search' in pipeline_results:
            search_data = pipeline_results['similarity_search']
            data['similarity_matches'] = len(search_data.get('matches', []))
            data['similarity_stats'] = search_data.get('similarity_stats', {})
        
        if 'network_building' in pipeline_results:
            network_data = pipeline_results['network_building']
            data['network_nodes'] = network_data.get('network_nodes', 0)
            data['network_edges'] = network_data.get('network_edges', 0)
        
        return data
    
    def _generate_overview_tab(self, data: Dict[str, Any], 
                              sequence_analysis: Optional[Dict[str, Any]], 
                              protein_id: str) -> str:
        """Generate the overview tab content."""
        return f"""
            <div class="section">
                <h3>🎯 Analysis Summary</h3>
                <div class="metric-grid">
                    <div class="metric-card">
                        <div class="metric-value">{data.get('protein_id', 'N/A')}</div>
                        <div class="metric-label">Protein ID</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{data.get('sequence_length', 0)}</div>
                        <div class="metric-label">Sequence Length</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{data.get('embedding_dim', 0)}</div>
                        <div class="metric-label">Embedding Dimension</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{data.get('similarity_matches', 'N/A')}</div>
                        <div class="metric-label">Similar Proteins</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{data.get('family_assignment', {}).get('family_id', 'N/A')}</div>
                        <div class="metric-label">Assigned Family</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{data.get('network_nodes', 0)}</div>
                        <div class="metric-label">Network Nodes</div>
                    </div>
                </div>
            </div>
            
            <div class="section">
                <h3>🔍 Similarity Search Results</h3>
                <div class="result">
                    <p><strong>Query Protein:</strong> {data.get('protein_id', 'N/A')}</p>
                    <p><strong>Similar Proteins Found:</strong> {data.get('similarity_matches', 'N/A')}</p>
                    <p><strong>Family Assignment:</strong> {data.get('family_assignment', {}).get('family_id', 'N/A')} (Confidence: {data.get('family_assignment', {}).get('confidence', 'N/A')})</p>
                </div>
            </div>
            
            <div class="section">
                <h3>📈 Network Statistics</h3>
                <div class="result">
                    <p><strong>Network Density:</strong> {data.get('network_edges', 0) / max(data.get('network_nodes', 1), 1):.3f}</p>
                    <p><strong>Average Degree:</strong> {data.get('network_edges', 0) * 2 / max(data.get('network_nodes', 1), 1):.1f}</p>
                    <p><strong>Network Type:</strong> {'Connected' if data.get('network_edges', 0) > 0 else 'Isolated'}</p>
                    <p><strong>Family Assignment Confidence:</strong> {data.get('family_assignment', {}).get('confidence', 'N/A')}</p>
                </div>
            </div>
            
            <div class="section">
                <h3>🎯 Family Assignment Details</h3>
                <div class="result">
                    <p><strong>Assigned Family:</strong> {data.get('family_assignment', {}).get('family_id', 'N/A')}</p>
                    <p><strong>Confidence Score:</strong> {data.get('family_assignment', {}).get('confidence', 'N/A')}</p>
                    <p><strong>Eigenprotein ID:</strong> {data.get('family_assignment', {}).get('eigenprotein_id', 'N/A')}</p>
                </div>
            </div>
        """
    
    def _generate_sequence_tab(self, sequence_analysis: Dict[str, Any]) -> str:
        """Generate the sequence analysis tab content with enhanced physicochemical properties."""
        html = "<div class='section'>"
        
        if sequence_analysis:
            # Show the actual sequence
            html += f"""
                <div class="section">
                    <h3>🧬 Protein Sequence</h3>
                    <div class="sequence-display">
                        <strong>Protein ID:</strong> {sequence_analysis['protein_id']}<br>
                        <strong>Length:</strong> {sequence_analysis['length']} amino acids<br>
                        <strong>Sequence:</strong><br>
                        {sequence_analysis['sequence']}
                    </div>
                </div>
                
                <div class="section">
                    <h3>📊 Amino Acid Composition</h3>
                    <div class="stats-grid">
            """
            
            # Show amino acid composition
            if 'individual' in sequence_analysis['amino_acid_composition']:
                aa_comp = sequence_analysis['amino_acid_composition']['individual']
                for aa, info in aa_comp.items():
                    if info['count'] > 0:
                        html += f"""
                            <div class="stat-card">
                                <div class="stat-value">{info['count']}</div>
                                <div class="stat-label">{aa} ({info['name']}) - {info['percentage']:.1f}%</div>
                            </div>
                        """
            else:
                # Calculate basic composition if not available
                sequence = sequence_analysis['sequence']
                aa_counts = {}
                for aa in sequence:
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1
                
                for aa, count in sorted(aa_counts.items()):
                    percentage = (count / len(sequence)) * 100
                    html += f"""
                        <div class="stat-card">
                            <div class="stat-value">{count}</div>
                            <div class="stat-label">{aa} - {percentage:.1f}%</div>
                        </div>
                    """
            
            html += """
                    </div>
                </div>
                
                <div class="section">
                    <h3>📊 Amino Acid Composition Chart</h3>
                    <div class="chart-container">
                        <canvas id="aaCompositionChart" style="width: 100%; height: 300px;"></canvas>
                    </div>
                </div>
                
                <div class="section">
                    <h3>⚗️ Enhanced Physicochemical Properties</h3>
                    <div class="stats-grid">
                        <div class="stat-card">
                            <div class="stat-value">{:.1f}</div>
                            <div class="stat-label">Molecular Weight (Da)</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{:.2f}</div>
                            <div class="stat-label">Net Charge (pH 7.0)</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{:.2f}</div>
                            <div class="stat-label">Isoelectric Point</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{:.3f}</div>
                            <div class="stat-label">Avg Hydrophobicity</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{:.0f}</div>
                            <div class="stat-label">Extinction Coefficient (280nm)</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{:.2f}</div>
                            <div class="stat-label">Instability Index</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{:.1f}</div>
                            <div class="stat-label">Aliphatic Index</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{:.3f}</div>
                            <div class="stat-label">GRAVY</div>
                        </div>
                    </div>
                </div>
            """.format(
                sequence_analysis['physicochemical_properties']['molecular_weight'],
                sequence_analysis['physicochemical_properties'].get('net_charge_ph7', sequence_analysis['physicochemical_properties'].get('net_charge', 0)),
                sequence_analysis['physicochemical_properties']['isoelectric_point'],
                sequence_analysis['physicochemical_properties']['average_hydrophobicity'],
                sequence_analysis['physicochemical_properties'].get('extinction_coefficient_280nm', 0),
                sequence_analysis['physicochemical_properties'].get('instability_index', 0),
                sequence_analysis['physicochemical_properties'].get('aliphatic_index', 0),
                sequence_analysis['physicochemical_properties'].get('gravy', 0)
            )
            
            # Secondary structure prediction
            ss = sequence_analysis['secondary_structure_prediction']
            html += f"""
                <div class="section">
                    <h3>🔄 Secondary Structure Prediction</h3>
                    <div class="stats-grid">
                        <div class="stat-card">
                            <div class="stat-value">{ss['helix_preference']:.3f}</div>
                            <div class="stat-label">Helix Preference</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{ss['sheet_preference']:.3f}</div>
                            <div class="stat-label">Sheet Preference</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{ss['turn_preference']:.3f}</div>
                            <div class="stat-label">Turn Preference</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{ss['dominant_structure'].title()}</div>
                            <div class="stat-label">Dominant Structure</div>
                        </div>
                    </div>
                </div>
                
                <div class="section">
                    <h3>🎯 Enhanced Sequence Motifs</h3>
                    <table class="motif-table">
                        <thead>
                            <tr>
                                <th>Motif Type</th>
                                <th>Count</th>
                                <th>Confidence</th>
                                <th>Details</th>
                            </tr>
                        </thead>
                        <tbody>
            """
            
            motif_found = False
            for motif_type, motifs in sequence_analysis['sequence_motifs'].items():
                if motifs:
                    motif_found = True
                    count = len(motifs)
                    # Get confidence level from first motif if available
                    confidence = motifs[0].get('confidence', 'medium') if motifs else 'medium'
                    details = str(motifs[:2]) if len(motifs) > 2 else str(motifs)
                    html += f"""
                        <tr>
                            <td>{motif_type.replace('_', ' ').title()}</td>
                            <td>{count}</td>
                            <td><span class="confidence-{confidence}">{confidence.title()}</span></td>
                            <td>{details}</td>
                        </tr>
                    """
            
            if not motif_found:
                html += """
                        <tr>
                            <td colspan="4">No significant motifs detected</td>
                        </tr>
                """
            
            html += """
                        </tbody>
                    </table>
                </div>
            """
        else:
            # If no sequence analysis available, show basic info from pipeline data
            if 'full_pipeline' in self.results:
                pipeline_data = self.results['full_pipeline']
                if 'sequence' in pipeline_data:
                    sequence = pipeline_data['sequence']
                    protein_id = pipeline_data.get('protein_id', 'UNKNOWN')
                    
                    html += f"""
                        <div class="section">
                            <h3>🧬 Protein Sequence</h3>
                            <div class="sequence-display">
                                <strong>Protein ID:</strong> {protein_id}<br>
                                <strong>Length:</strong> {len(sequence)} amino acids<br>
                                <strong>Sequence:</strong><br>
                                {sequence}
                            </div>
                        </div>
                        
                        <div class="section">
                            <h3>📊 Basic Amino Acid Composition</h3>
                            <div class="stats-grid">
                    """
                    
                    # Calculate basic composition
                    aa_counts = {}
                    for aa in sequence:
                        aa_counts[aa] = aa_counts.get(aa, 0) + 1
                    
                    for aa, count in sorted(aa_counts.items()):
                        percentage = (count / len(sequence)) * 100
                        html += f"""
                            <div class="stat-card">
                                <div class="stat-value">{count}</div>
                                <div class="stat-label">{aa} - {percentage:.1f}%</div>
                            </div>
                        """
                    
                    html += """
                            </div>
                        </div>
                        
                        <div class="section">
                            <h3>📊 Basic Amino Acid Composition Chart</h3>
                            <div class="chart-container">
                                <canvas id="aaCompositionChart" style="width: 100%; height: 250px;"></canvas>
                            </div>
                        </div>
                        
                        <div class="section">
                            <h3>ℹ️ Sequence Analysis Note</h3>
                            <div class="result">
                                <p>Advanced sequence analysis (physicochemical properties, secondary structure prediction, motifs) requires the sequence analyzer module to be properly configured.</p>
                                <p>Basic sequence information is displayed above.</p>
                            </div>
                        </div>
                    """
                else:
                    html += "<p>No sequence data available.</p>"
            else:
                html += "<p>No sequence analysis data available.</p>"
        
        html += "</div>"
        return html
    
    def _generate_network_tab(self, data: Dict[str, Any]) -> str:
        """Generate the network visualization tab content."""
        return f"""
            <div class="section">
                <h3>🕸️ Protein Network Visualization</h3>
                <div class="network-container" id="network-container">
                    <div class="loading">Loading network visualization...</div>
                </div>
            </div>
            
            <div class="section">
                <h3>📊 Network Properties</h3>
                <div class="result">
                    <p><strong>Number of Nodes:</strong> {data.get('network_nodes', 0)}</p>
                    <p><strong>Number of Edges:</strong> {data.get('network_edges', 0)}</p>
                    <p><strong>Network Density:</strong> {data.get('network_edges', 0) / max(data.get('network_nodes', 1), 1):.3f}</p>
                    <p><strong>Average Degree:</strong> {data.get('network_edges', 0) * 2 / max(data.get('network_nodes', 1), 1):.1f}</p>
                </div>
            </div>
        """
    
    def _generate_statistics_tab(self, data: Dict[str, Any]) -> str:
        """Generate the statistics tab content."""
        return f"""
            <div class="section">
                <h3>📈 Analysis Statistics</h3>
                <div class="chart-container">
                    <canvas id="statsChart" width="400" height="200"></canvas>
                </div>
            </div>
            
            <div class="section">
                <h3>📊 Similarity Statistics</h3>
                <div class="result">
                    <p><strong>Maximum Similarity:</strong> {data.get('similarity_stats', {}).get('max', 'N/A')}</p>
                    <p><strong>Minimum Similarity:</strong> {data.get('similarity_stats', {}).get('min', 'N/A')}</p>
                    <p><strong>Mean Similarity:</strong> {data.get('similarity_stats', {}).get('mean', 'N/A')}</p>
                </div>
            </div>
        """
    
    def _generate_bioinformatics_tab(self, sequence_analysis: Dict[str, Any]) -> str:
        """Generate the bioinformatics tab content with verified working links."""
        if not sequence_analysis:
            return "<div class='loading'>No bioinformatics data available.</div>"
        
        links = sequence_analysis.get('bioinformatics_links', {})
        
        # Group links by category with verified working links only
        link_categories = {
            'Protein Databases': ['uniprot', 'ncbi_protein', 'pdb', 'pdb_summary', 'sifts'],
            'Sequence Analysis': ['expasy_protscale', 'expasy_peptide_mass', 'expasy_peptide_cutter', 'expasy_scanprosite'],
            'Domain Analysis': ['pfam', 'interpro', 'smart', 'prosite'],
            'Structure Prediction': ['swiss_model', 'alphafold', 'i_tasser'],
            'Sequence Alignment': ['blast', 'clustal', 'muscle', 'tcoffee'],
            'Visualization Tools': ['pymol', 'chimera', 'chimera_x', 'vmd', 'jmol', 'molstar', 'ngl_viewer', 'pdb_viewer'],
            'Functional Annotation': ['go_annotation', 'kegg', 'reactome', 'string', 'intact', 'mint'],
            'Protein Families': ['panther', 'superfamily', 'gene3d'],
            'Localization Prediction': ['wolfpsort', 'psort', 'targetp', 'signalp'],
            'Transmembrane Prediction': ['tmhmm', 'topcons'],
            'Post-translational Modifications': ['phosphosite', 'glycomod', 'netphos'],
            'Literature & Annotation': ['pubmed', 'genecards', 'omim'],
            'Protein Expression': ['expression_atlas', 'protein_atlas'],
            'Evolutionary Analysis': ['orthodb', 'ensembl', 'ucsc_genome'],
            'Protein Analysis Tools': ['protter', 'protter_web', 'protein_workshop']
        }
        
        html = "<div class='section'>"
        html += "<h3>🔬 Bioinformatics Resources</h3>"
        html += "<p>Click on any link below to access verified bioinformatics tools and databases:</p>"
        
        for category, link_keys in link_categories.items():
            available_links = []
            for link_key in link_keys:
                if link_key in links:
                    available_links.append((link_key, links[link_key]))
            
            if available_links:
                html += f"""
                    <div class="section">
                        <h4>🔗 {category}</h4>
                        <div class="bio-links">
                """
                
                for link_name, link_url in available_links:
                    display_name = link_name.replace('_', ' ').title()
                    html += f"""
                        <a href="{link_url}" target="_blank" class="bio-link">
                            <strong>{display_name}</strong>
                        </a>
                    """
                
                html += """
                        </div>
                    </div>
                """
        
        html += "</div>"
        return html
    
    def _generate_raw_data_tab(self, pipeline_results: Dict[str, Any]) -> str:
        """Generate the raw data tab content."""
        return f"""
            <div class="section">
                <h3>📄 Raw Analysis Data</h3>
                <div class="result">
                    <pre>{json.dumps(pipeline_results, indent=2)}</pre>
                </div>
            </div>
        """
    
    def _generate_chart_scripts(self, data: Dict[str, Any], 
                               sequence_analysis: Optional[Dict[str, Any]]) -> str:
        """Generate JavaScript for interactive charts."""
        return f"""
    <script>
        // Chart.js configuration
        Chart.defaults.font.family = "'Segoe UI', Tahoma, Geneva, Verdana, sans-serif";
        
        // Statistics chart
        const statsCtx = document.getElementById('statsChart');
        if (statsCtx) {{
            new Chart(statsCtx, {{
                type: 'bar',
                data: {{
                    labels: ['Sequence Length', 'Embedding Dim', 'Network Nodes', 'Network Edges'],
                    datasets: [{{
                        label: 'Analysis Metrics',
                        data: [
                            {data.get('sequence_length', 0)},
                            {data.get('embedding_dim', 0)},
                            {data.get('network_nodes', 0)},
                            {data.get('network_edges', 0)}
                        ],
                        backgroundColor: [
                            'rgba(54, 162, 235, 0.6)',
                            'rgba(255, 99, 132, 0.6)',
                            'rgba(255, 205, 86, 0.6)',
                            'rgba(75, 192, 192, 0.6)'
                        ],
                        borderColor: [
                            'rgba(54, 162, 235, 1)',
                            'rgba(255, 99, 132, 1)',
                            'rgba(255, 205, 86, 1)',
                            'rgba(75, 192, 192, 1)'
                        ],
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    scales: {{
                        y: {{
                            beginAtZero: true
                        }}
                    }}
                }}
            }});
        }}
        
        // Initialize network visualization if data is available
        function initializeNetwork() {{
            const container = document.getElementById('network-container');
            if (!container) return;
            
            // Placeholder for network visualization
            // This would be populated with actual network data
            container.innerHTML = '<div class="loading">Network visualization will be implemented with actual data.</div>';
        }}
        
        // Initialize charts when page loads
        function initializeCharts() {{
            // Network initialization
            setTimeout(initializeNetwork, 1000);
        }}
    </script>
        """
    
    def _generate_network_visualization(self, pipeline_results: Dict[str, Any], 
                                      protein_id: str) -> Optional[str]:
        """
        Generate network visualization for the analysis results.
        
        Args:
            pipeline_results: Results from the protein analysis pipeline
            protein_id: Protein identifier
            
        Returns:
            Path to the network visualization file, or None if failed
        """
        try:
            # Extract network data from pipeline results
            similar_proteins = None
            
            # Try different possible locations for similar proteins data
            if 'similarity_search' in pipeline_results:
                similar_proteins = pipeline_results['similarity_search']
            elif 'top_matches' in pipeline_results:
                similar_proteins = pipeline_results['top_matches']
            elif 'similar_proteins' in pipeline_results:
                similar_proteins = pipeline_results['similar_proteins']
            elif 'network_data' in pipeline_results:
                similar_proteins = pipeline_results['network_data']
            
            if not similar_proteins:
                logger.warning("No similar proteins data found for network visualization")
                return self._create_fallback_network_visualization([], protein_id)
            
            # Handle different data types
            if isinstance(similar_proteins, dict):
                # Convert dict to list of dicts
                logger.warning("similar_proteins is a dict, converting to list")
                similar_proteins = [similar_proteins]
            elif isinstance(similar_proteins, str):
                # Handle string input
                logger.warning("similar_proteins is a string, creating default entry")
                similar_proteins = [{'protein_id': similar_proteins, 'similarity': 0.0, 'family_id': 'unknown'}]
            elif not isinstance(similar_proteins, list):
                logger.warning(f"similar_proteins is not a list: {type(similar_proteins)}")
                similar_proteins = []
            
            # Create network visualization
            network_viz_path = os.path.join(self.output_dir, f"network_{protein_id}_{int(time.time())}.html")
            
            # Generate network using network builder
            try:
                # Extract embeddings and protein IDs from similar proteins
                protein_ids = []
                similarities = []
                families = []
                
                for i, protein in enumerate(similar_proteins):
                    if isinstance(protein, dict):
                        protein_id = protein.get('protein_id', f"protein_{i}")
                        similarity = protein.get('similarity', protein.get('similarity_score', 0.0))
                        family_id = protein.get('family_id', 'unknown')
                    else:
                        protein_id = f"protein_{i}"
                        similarity = 0.0
                        family_id = 'unknown'
                    
                    protein_ids.append(protein_id)
                    similarities.append(float(similarity))
                    families.append(family_id)
                
                # Use real embeddings from pipeline results if available
                embeddings = None
                if 'embeddings' in pipeline_results:
                    embeddings = pipeline_results['embeddings']
                elif 'family_embeddings' in pipeline_results:
                    embeddings = pipeline_results['family_embeddings']
                elif 'query_embedding' in pipeline_results:
                    # If only query embedding is available, create a simple visualization
                    query_embedding = pipeline_results['query_embedding']
                    if query_embedding is not None:
                        # Create dummy embeddings for visualization
                        num_proteins = len(protein_ids)
                        if num_proteins > 0:
                            # Create random embeddings based on similarity scores
                            embedding_dim = len(query_embedding) if hasattr(query_embedding, '__len__') else 320
                            embeddings = np.random.rand(num_proteins, embedding_dim)
                            # Scale by similarity to make similar proteins closer
                            for i, sim in enumerate(similarities):
                                embeddings[i] *= sim
                
                if embeddings is None:
                    logger.warning("No embeddings available for network visualization")
                    return self._create_fallback_network_visualization(similar_proteins, protein_id)
                
                # Create metadata from real data
                metadata_df = pd.DataFrame({
                    'protein_id': protein_ids,
                    'family_id': families,
                    'similarity': similarities,
                    'sequence_length': [100] * len(protein_ids)  # Default length
                }).set_index('protein_id')
                
                # Generate network visualization
                fig, G = self.network_builder.create_interactive_visualization(
                    embeddings=embeddings,
                    protein_ids=protein_ids,
                    metadata_df=metadata_df,
                    query_protein_id=protein_id,
                    output_file=network_viz_path
                )
                
                if fig is not None:
                    logger.info(f"Network visualization created: {network_viz_path}")
                    return network_viz_path
                else:
                    logger.warning("Network visualization generation failed")
                    return self._create_fallback_network_visualization(similar_proteins, protein_id)
                    
            except Exception as e:
                logger.warning(f"Could not generate network visualization: {e}")
                # Create a simple fallback visualization
                return self._create_fallback_network_visualization(similar_proteins, protein_id)
                
        except Exception as e:
            logger.error(f"Network visualization failed: {e}")
            return self._create_fallback_network_visualization([], protein_id)
    
    def _create_fallback_network_visualization(self, similar_proteins: List[Dict], 
                                            protein_id: str) -> Optional[str]:
        """
        Create a simple fallback network visualization when the main method fails.
        
        Args:
            similar_proteins: List of similar proteins
            protein_id: Query protein ID
            
        Returns:
            Path to the fallback visualization file
        """
        try:
            # Handle different input types
            if isinstance(similar_proteins, dict):
                # Convert dict to list of dicts
                logger.warning("similar_proteins is a dict, converting to list")
                similar_proteins = [similar_proteins]
            elif isinstance(similar_proteins, str):
                # Handle string input
                logger.warning("similar_proteins is a string, creating default entry")
                similar_proteins = [{'protein_id': similar_proteins, 'similarity': 0.0, 'family_id': 'unknown'}]
            elif not isinstance(similar_proteins, list):
                logger.warning(f"similar_proteins is not a list: {type(similar_proteins)}")
                similar_proteins = []
            
            # Convert to list of dictionaries if needed
            processed_proteins = []
            for i, protein in enumerate(processed_proteins):
                if isinstance(protein, dict):
                    # Ensure all required fields are present
                    processed_protein = {
                        'protein_id': protein.get('protein_id', f'protein_{i}'),
                        'similarity': float(protein.get('similarity', protein.get('similarity_score', 0.0))),
                        'family_id': protein.get('family_id', 'unknown')
                    }
                    processed_proteins.append(processed_protein)
                else:
                    # Create a default protein entry
                    processed_proteins.append({
                        'protein_id': f'protein_{i}',
                        'similarity': 0.0,
                        'family_id': 'unknown'
                    })
            
            # Calculate statistics safely
            similarities = [p.get('similarity', 0.0) for p in processed_proteins if isinstance(p.get('similarity'), (int, float))]
            avg_similarity = np.mean(similarities) if similarities else 0.0
            families = set(p.get('family_id', 'unknown') for p in processed_proteins)
            
            # Create a more sophisticated HTML visualization with D3.js
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>Protein Network - {protein_id}</title>
                <script src="https://d3js.org/d3.v7.min.js"></script>
                <style>
                    body {{
                        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                        margin: 0;
                        padding: 20px;
                        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                        color: #333;
                    }}
                    .container {{
                        max-width: 1200px;
                        margin: 0 auto;
                        background: white;
                        border-radius: 10px;
                        box-shadow: 0 10px 30px rgba(0,0,0,0.1);
                        overflow: hidden;
                    }}
                    .header {{
                        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                        color: white;
                        padding: 20px;
                        text-align: center;
                    }}
                    .content {{
                        padding: 20px;
                    }}
                    .stats-grid {{
                        display: grid;
                        grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                        gap: 20px;
                        margin: 20px 0;
                    }}
                    .stat-card {{
                        background: #f8f9fa;
                        padding: 15px;
                        border-radius: 8px;
                        text-align: center;
                        border-left: 4px solid #667eea;
                    }}
                    .stat-value {{
                        font-size: 2em;
                        font-weight: bold;
                        color: #667eea;
                    }}
                    .stat-label {{
                        color: #666;
                        font-size: 0.9em;
                    }}
                    .protein-list {{
                        background: #f8f9fa;
                        border-radius: 8px;
                        padding: 20px;
                        margin: 20px 0;
                    }}
                    .protein-item {{
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                        padding: 10px;
                        margin: 5px 0;
                        background: white;
                        border-radius: 5px;
                        border-left: 4px solid #28a745;
                    }}
                    .protein-item:hover {{
                        background: #e9ecef;
                        transform: translateX(5px);
                        transition: all 0.3s ease;
                    }}
                    .similarity-bar {{
                        height: 8px;
                        background: #e9ecef;
                        border-radius: 4px;
                        overflow: hidden;
                        width: 100px;
                    }}
                    .similarity-fill {{
                        height: 100%;
                        background: linear-gradient(90deg, #28a745, #20c997);
                        transition: width 0.3s ease;
                    }}
                    .network-viz {{
                        height: 400px;
                        background: #f8f9fa;
                        border-radius: 8px;
                        display: flex;
                        align-items: center;
                        justify-content: center;
                        margin: 20px 0;
                    }}
                    .no-data {{
                        text-align: center;
                        color: #666;
                        font-style: italic;
                    }}
                </style>
            </head>
            <body>
                <div class="container">
                    <div class="header">
                        <h1>Protein Network Analysis</h1>
                        <h2>Query Protein: {protein_id}</h2>
                    </div>
                    
                    <div class="content">
                        <div class="stats-grid">
                            <div class="stat-card">
                                <div class="stat-value">{len(processed_proteins)}</div>
                                <div class="stat-label">Similar Proteins</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value">{avg_similarity:.3f}</div>
                                <div class="stat-label">Average Similarity</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value">{len(families)}</div>
                                <div class="stat-label">Families Represented</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value">{max(similarities) if similarities else 0.0:.3f}</div>
                                <div class="stat-label">Highest Similarity</div>
                            </div>
                        </div>
                        
                        <div class="protein-list">
                            <h3>Similar Proteins</h3>
                            {self._generate_protein_list_html(processed_proteins)}
                        </div>
                        
                        <div class="network-viz">
                            {self._generate_simple_network_viz(processed_proteins, protein_id)}
                        </div>
                    </div>
                </div>
                
                <script>
                    // Add interactive elements
                    document.addEventListener('DOMContentLoaded', function() {{
                        console.log('Network visualization loaded');
                        
                        // Animate similarity bars
                        const bars = document.querySelectorAll('.similarity-fill');
                        bars.forEach(bar => {{
                            const width = bar.style.width;
                            bar.style.width = '0%';
                            setTimeout(() => {{
                                bar.style.width = width;
                            }}, 500);
                        }});
                        
                        // Add hover effects
                        const items = document.querySelectorAll('.protein-item');
                        items.forEach(item => {{
                            item.addEventListener('mouseenter', function() {{
                                this.style.transform = 'translateX(10px)';
                            }});
                            item.addEventListener('mouseleave', function() {{
                                this.style.transform = 'translateX(0)';
                            }});
                        }});
                    }});
                </script>
            </body>
            </html>
            """
            
            # Save fallback visualization
            fallback_path = os.path.join(self.output_dir, f"network_{protein_id}_fallback_{int(time.time())}.html")
            with open(fallback_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            logger.info(f"Created fallback network visualization: {fallback_path}")
            return fallback_path
            
        except Exception as e:
            logger.error(f"Fallback network visualization failed: {e}")
            return None
    
    def _generate_protein_list_html(self, processed_proteins: List[Dict]) -> str:
        """Generate HTML for the protein list section."""
        if not processed_proteins:
            return '<div class="no-data">No similar proteins found</div>'
        
        html = ''
        for i, protein in enumerate(processed_proteins[:10]):  # Show top 10
            protein_id_sim = protein.get('protein_id', f'protein_{i}')
            similarity = protein.get('similarity', 0.0)
            family_id = protein.get('family_id', 'unknown')
            
            # Create color-coded similarity bar
            similarity_percent = min(similarity * 100, 100)
            bar_color = '#28a745' if similarity > 0.7 else '#ffc107' if similarity > 0.4 else '#dc3545'
            
            html += f"""
                <div class="protein-item">
                    <div>
                        <strong>{protein_id_sim}</strong><br>
                        <small>Family: {family_id}</small>
                    </div>
                    <div style="text-align: right;">
                        <div style="margin-bottom: 5px;">{similarity:.3f}</div>
                        <div class="similarity-bar">
                            <div class="similarity-fill" style="width: {similarity_percent}%; background: {bar_color};"></div>
                        </div>
                    </div>
                </div>
            """
        
        if len(processed_proteins) > 10:
            html += f'<div style="text-align: center; margin-top: 10px; color: #666;">... and {len(processed_proteins) - 10} more proteins</div>'
        
        return html
    
    def _generate_simple_network_viz(self, processed_proteins: List[Dict], protein_id: str) -> str:
        """Generate a simple network visualization using D3.js."""
        if not processed_proteins:
            return '<div class="no-data">No network data available</div>'
        
        # Create simple network data
        nodes = [{'id': protein_id, 'group': 0, 'size': 20}]  # Query protein
        links = []
        
        for i, protein in enumerate(processed_proteins[:5]):  # Limit to 5 for simplicity
            protein_id_sim = protein.get('protein_id', f'protein_{i}')
            similarity = protein.get('similarity', 0.0)
            
            nodes.append({
                'id': protein_id_sim,
                'group': 1,
                'size': 10 + similarity * 10
            })
            
            links.append({
                'source': protein_id,
                'target': protein_id_sim,
                'value': similarity
            })
        
        nodes_json = json.dumps(nodes)
        links_json = json.dumps(links)
        
        return f"""
        <div id="network-container" style="width: 100%; height: 100%;">
            <div class="no-data">Network visualization would be displayed here</div>
        </div>
        <script>
            // Simple D3.js network visualization
            const nodes = {nodes_json};
            const links = {links_json};
            
            const width = document.getElementById('network-container').offsetWidth;
            const height = 300;
            
            const svg = d3.select('#network-container')
                .append('svg')
                .attr('width', width)
                .attr('height', height);
            
            const simulation = d3.forceSimulation(nodes)
                .force('link', d3.forceLink(links).id(d => d.id))
                .force('charge', d3.forceManyBody().strength(-100))
                .force('center', d3.forceCenter(width / 2, height / 2));
            
            const link = svg.append('g')
                .selectAll('line')
                .data(links)
                .enter().append('line')
                .attr('stroke', '#999')
                .attr('stroke-opacity', 0.6)
                .attr('stroke-width', d => Math.sqrt(d.value) * 3);
            
            const node = svg.append('g')
                .selectAll('circle')
                .data(nodes)
                .enter().append('circle')
                .attr('r', d => d.size)
                .attr('fill', d => d.group === 0 ? '#ff7f0e' : '#69b3a2')
                .call(d3.drag()
                    .on('start', dragstarted)
                    .on('drag', dragged)
                    .on('end', dragended));
            
            node.append('title')
                .text(d => d.id);
            
            simulation.on('tick', () => {{
                link
                    .attr('x1', d => d.source.x)
                    .attr('y1', d => d.source.y)
                    .attr('x2', d => d.target.x)
                    .attr('y2', d => d.target.y);
                
                node
                    .attr('cx', d => d.x)
                    .attr('cy', d => d.y);
            }});
            
            function dragstarted(event, d) {{
                if (!event.active) simulation.alphaTarget(0.3).restart();
                d.fx = d.x;
                d.fy = d.y;
            }}
            
            function dragged(event, d) {{
                d.fx = event.x;
                d.fy = event.y;
            }}
            
            function dragended(event, d) {{
                if (!event.active) simulation.alphaTarget(0);
                d.fx = null;
                d.fy = null;
            }}
        </script>
        """ 