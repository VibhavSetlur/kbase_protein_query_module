"""
HTML Report Generator

This module creates comprehensive HTML reports for protein analysis pipeline outputs.
It integrates with the visualization generator to create rich, interactive reports.
"""

import os
import time
import json
import logging
from typing import Dict, List, Any, Optional, Union
from pathlib import Path

from .visualization_generator import VisualizationGenerator

logger = logging.getLogger(__name__)

class HTMLReportGenerator:
    """
    Comprehensive HTML report generator for protein analysis pipeline.
    
    Features:
    - Single and multi-protein analysis reports
    - Network visualization integration
    - Sequence analysis summaries
    - Embedding generation results
    - Family assignment information
    - Similarity search results
    - Professional responsive design
    """
    
    def __init__(self, output_dir: str = None):
        """Initialize the HTML report generator."""
        if output_dir:
            self.output_dir = output_dir
        else:
            # Use test/outputs for testing, or scratch space for KBase
            if os.path.exists('test/outputs'):
                self.output_dir = 'test/outputs'
            else:
                self.output_dir = os.environ.get('SCRATCH_DIR', '/kb/module/work/tmp')
        
        os.makedirs(self.output_dir, exist_ok=True)
        self.visualization_generator = VisualizationGenerator(output_dir)
        
    def generate_comprehensive_report(self, pipeline_results: Dict[str, Any], 
                                   protein_id: str, sequence: str = None) -> Dict[str, Any]:
        """
        Generate a comprehensive HTML report for protein analysis results.
        
        Args:
            pipeline_results: Dictionary containing all pipeline stage results
            protein_id: The protein identifier
            sequence: Optional protein sequence
            
        Returns:
            Dictionary with report paths and metadata
        """
        try:
            # Generate timestamp for unique filenames
            timestamp = int(time.time())
            
            # Create report filename
            report_filename = f"protein_analysis_report_{timestamp}.html"
            report_path = os.path.join(self.output_dir, report_filename)
            
            # Generate the HTML content
            html_content = self._generate_html_content(pipeline_results, protein_id, sequence)
            
            # Write the HTML file
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            # Generate network visualization if available
            network_viz_path = None
            if 'network_building' in pipeline_results:
                network_viz_path = self._generate_network_visualization(
                    pipeline_results['network_building'], protein_id, timestamp
                )
            
            # Generate dashboard data
            dashboard_data = self._extract_dashboard_data(pipeline_results)
            
            result = {
                'html_path': report_path,
                'network_viz_path': network_viz_path,
                'dashboard_data': dashboard_data,
                'timestamp': timestamp,
                'protein_id': protein_id
            }
            
            logger.info(f"Generated comprehensive report: {report_path}")
            return result
            
        except Exception as e:
            logger.error(f"Error generating comprehensive report: {e}")
            raise
    
    def _generate_html_content(self, pipeline_results: Dict[str, Any], 
                             protein_id: str, sequence: str = None) -> str:
        """Generate the main HTML content for the report."""
        
        # Extract key information
        protein_exists = pipeline_results.get('protein_existence', {}).get('exists', False)
        family_id = pipeline_results.get('family_assignment', {}).get('family_id', 'Unknown')
        embedding_info = pipeline_results.get('embedding_generation', {})
        similarity_info = pipeline_results.get('similarity_search', {})
        network_info = pipeline_results.get('network_building', {})
        
        # Generate HTML content
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Protein Analysis Report - {protein_id}</title>
    <style>
        body {{
            font-family: 'Times New Roman', serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .header {{
            text-align: center;
            border-bottom: 3px solid #2c3e50;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        .section {{
            margin-bottom: 30px;
            padding: 20px;
            border: 1px solid #ddd;
            border-radius: 5px;
            background-color: #fafafa;
        }}
        .section h2 {{
            color: #2c3e50;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
        }}
        .status {{
            padding: 10px;
            border-radius: 5px;
            font-weight: bold;
            text-align: center;
        }}
        .status.exists {{
            background-color: #d4edda;
            color: #155724;
            border: 1px solid #c3e6cb;
        }}
        .status.not-exists {{
            background-color: #f8d7da;
            color: #721c24;
            border: 1px solid #f5c6cb;
        }}
        .data-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        .data-card {{
            background: white;
            padding: 15px;
            border-radius: 5px;
            border: 1px solid #ddd;
        }}
        .metric {{
            display: flex;
            justify-content: space-between;
            margin-bottom: 10px;
            padding: 8px;
            background-color: #f8f9fa;
            border-radius: 3px;
        }}
        .metric-label {{
            font-weight: bold;
            color: #495057;
        }}
        .metric-value {{
            color: #6c757d;
        }}
        .sequence {{
            font-family: 'Courier New', monospace;
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            border: 1px solid #dee2e6;
            word-break: break-all;
            font-size: 14px;
        }}
        .network-info {{
            background-color: #e3f2fd;
            padding: 15px;
            border-radius: 5px;
            border: 1px solid #bbdefb;
        }}
        .similarity-list {{
            list-style: none;
            padding: 0;
        }}
        .similarity-item {{
            display: flex;
            justify-content: space-between;
            padding: 8px;
            margin-bottom: 5px;
            background-color: #f8f9fa;
            border-radius: 3px;
            border-left: 4px solid #3498db;
        }}
        .footer {{
            text-align: center;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #6c757d;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Comprehensive Protein Network Analysis Report</h1>
            <p>Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="section">
            <h2>Protein Information</h2>
            <div class="status {'exists' if protein_exists else 'not-exists'}">
                {f'Protein {protein_id} found in database' if protein_exists else f'Protein {protein_id} not found in database'}
            </div>
            <div class="metric">
                <span class="metric-label">Protein ID:</span>
                <span class="metric-value">{protein_id}</span>
            </div>
            <div class="metric">
                <span class="metric-label">Family ID:</span>
                <span class="metric-value">{family_id}</span>
            </div>
            {f'<div class="metric"><span class="metric-label">Sequence Length:</span><span class="metric-value">{len(sequence)}</span></div>' if sequence else ''}
        </div>
        
        {f'<div class="section"><h2>Protein Sequence</h2><div class="sequence">{sequence}</div></div>' if sequence else ''}
        
        {self._generate_embedding_section(embedding_info) if embedding_info else ''}
        
        {self._generate_family_section(pipeline_results.get('family_assignment', {})) if pipeline_results.get('family_assignment') else ''}
        
        {self._generate_similarity_section(similarity_info) if similarity_info else ''}
        
        {self._generate_network_section(network_info) if network_info else ''}
        
        <div class="footer">
            <p>Report generated by KBase Protein Query Module</p>
            <p>Analysis completed successfully</p>
        </div>
    </div>
</body>
</html>
        """
        
        return html_content
    
    def _generate_embedding_section(self, embedding_info: Dict[str, Any]) -> str:
        """Generate HTML section for embedding information."""
        if not embedding_info:
            return ""
            
        return f"""
        <div class="section">
            <h2>Embedding Generation</h2>
            <div class="data-grid">
                <div class="data-card">
                    <div class="metric">
                        <span class="metric-label">Model:</span>
                        <span class="metric-value">{embedding_info.get('model_name', 'Unknown')}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">Embedding Dimension:</span>
                        <span class="metric-value">{embedding_info.get('embedding_dim', 'Unknown')}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">Sequence Length:</span>
                        <span class="metric-value">{embedding_info.get('sequence_length', 'Unknown')}</span>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_family_section(self, family_info: Dict[str, Any]) -> str:
        """Generate HTML section for family assignment information."""
        if not family_info:
            return ""
            
        return f"""
        <div class="section">
            <h2>Family Assignment</h2>
            <div class="data-grid">
                <div class="data-card">
                    <div class="metric">
                        <span class="metric-label">Family ID:</span>
                        <span class="metric-value">{family_info.get('family_id', 'Unknown')}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">Confidence:</span>
                        <span class="metric-value">{family_info.get('confidence', 'Unknown')}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">Eigenprotein ID:</span>
                        <span class="metric-value">{family_info.get('eigenprotein_id', 'Unknown')}</span>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_similarity_section(self, similarity_info: Dict[str, Any]) -> str:
        """Generate HTML section for similarity search results."""
        if not similarity_info:
            return ""
            
        matches = similarity_info.get('matches', [])
        stats = similarity_info.get('similarity_stats', {})
        
        matches_html = ""
        for match in matches:
            matches_html += f"""
                <li class="similarity-item">
                    <span>{match.get('protein_id', 'Unknown')}</span>
                    <span>{match.get('similarity', 'Unknown'):.3f}</span>
                </li>
            """
        
        return f"""
        <div class="section">
            <h2>Similarity Search Results</h2>
            <div class="data-grid">
                <div class="data-card">
                    <h3>Top Matches</h3>
                    <ul class="similarity-list">
                        {matches_html}
                    </ul>
                </div>
                <div class="data-card">
                    <h3>Similarity Statistics</h3>
                    <div class="metric">
                        <span class="metric-label">Maximum:</span>
                        <span class="metric-value">{stats.get('max', 'Unknown'):.3f}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">Minimum:</span>
                        <span class="metric-value">{stats.get('min', 'Unknown'):.3f}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">Mean:</span>
                        <span class="metric-value">{stats.get('mean', 'Unknown'):.3f}</span>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_network_section(self, network_info: Dict[str, Any]) -> str:
        """Generate HTML section for network information."""
        if not network_info:
            return ""
            
        return f"""
        <div class="section">
            <h2>Network Analysis</h2>
            <div class="network-info">
                <div class="data-grid">
                    <div class="data-card">
                        <div class="metric">
                            <span class="metric-label">Network Nodes:</span>
                            <span class="metric-value">{network_info.get('network_nodes', 'Unknown')}</span>
                        </div>
                        <div class="metric">
                            <span class="metric-label">Network Edges:</span>
                            <span class="metric-value">{network_info.get('network_edges', 'Unknown')}</span>
                        </div>
                        <div class="metric">
                            <span class="metric-label">Network Density:</span>
                            <span class="metric-value">{network_info.get('network_density', 'Unknown'):.3f}</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_network_visualization(self, network_info: Dict[str, Any], 
                                      protein_id: str, timestamp: int) -> Optional[str]:
        """Generate network visualization using the visualization generator."""
        try:
            # This would integrate with the visualization generator
            # For now, return None as the main HTML report includes network info
            return None
        except Exception as e:
            logger.warning(f"Could not generate network visualization: {e}")
            return None
    
    def _extract_dashboard_data(self, pipeline_results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract data for dashboard components."""
        dashboard_data = {
            'protein_count': 1,  # Single protein analysis
            'analysis_stages': list(pipeline_results.keys()),
            'timestamp': time.time(),
            'status': 'completed'
        }
        
        # Add specific metrics
        if 'embedding_generation' in pipeline_results:
            dashboard_data['embedding_dim'] = pipeline_results['embedding_generation'].get('embedding_dim')
            
        if 'similarity_search' in pipeline_results:
            dashboard_data['similarity_matches'] = len(pipeline_results['similarity_search'].get('matches', []))
            
        if 'network_building' in pipeline_results:
            dashboard_data['network_nodes'] = pipeline_results['network_building'].get('network_nodes')
            dashboard_data['network_edges'] = pipeline_results['network_building'].get('network_edges')
            
        return dashboard_data
    
    def generate_multi_protein_report(self, pipeline_results: Dict[str, Any], 
                                    protein_ids: List[str]) -> Dict[str, Any]:
        """
        Generate a comprehensive HTML report for multi-protein analysis.
        
        Args:
            pipeline_results: Dictionary containing all pipeline stage results
            protein_ids: List of protein identifiers
            
        Returns:
            Dictionary with report paths and metadata
        """
        try:
            timestamp = int(time.time())
            report_filename = f"multi_protein_report_{timestamp}.html"
            report_path = os.path.join(self.output_dir, report_filename)
            
            # Generate multi-protein HTML content
            html_content = self._generate_multi_protein_html_content(pipeline_results, protein_ids)
            
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            result = {
                'html_path': report_path,
                'timestamp': timestamp,
                'protein_count': len(protein_ids),
                'protein_ids': protein_ids
            }
            
            logger.info(f"Generated multi-protein report: {report_path}")
            return result
            
        except Exception as e:
            logger.error(f"Error generating multi-protein report: {e}")
            raise
    
    def _generate_multi_protein_html_content(self, pipeline_results: Dict[str, Any], 
                                           protein_ids: List[str]) -> str:
        """Generate HTML content for multi-protein analysis."""
        
        # This would be similar to single protein but adapted for multiple proteins
        # For now, return a basic multi-protein report
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Multi-Protein Analysis Report</title>
    <style>
        body {{
            font-family: 'Times New Roman', serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .header {{
            text-align: center;
            border-bottom: 3px solid #2c3e50;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        .section {{
            margin-bottom: 30px;
            padding: 20px;
            border: 1px solid #ddd;
            border-radius: 5px;
            background-color: #fafafa;
        }}
        .protein-list {{
            list-style: none;
            padding: 0;
        }}
        .protein-item {{
            padding: 15px;
            margin-bottom: 10px;
            background-color: white;
            border-radius: 5px;
            border: 1px solid #ddd;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Multi-Protein Analysis Report</h1>
            <p>Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="section">
            <h2>Analysis Summary</h2>
            <p>Analyzed {len(protein_ids)} proteins</p>
            <ul class="protein-list">
                {''.join([f'<li class="protein-item">{pid}</li>' for pid in protein_ids])}
            </ul>
        </div>
        
        <div class="section">
            <h2>Pipeline Results</h2>
            <p>Analysis stages completed: {', '.join(pipeline_results.keys())}</p>
        </div>
    </div>
</body>
</html>
        """
        
        return html_content
