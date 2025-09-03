"""
Visualization Generator

This module creates comprehensive HTML reports for visualization outputs from analysis stages.
It supports multiple visualization libraries and formats, converting them to interactive HTML reports.
"""

import os
import time
import json
import logging
import base64
from io import BytesIO
from typing import Dict, List, Any, Optional, Union
from pathlib import Path

logger = logging.getLogger(__name__)

# Try to import visualization libraries with graceful fallbacks
try:
    import plotly.graph_objects as go
    import plotly.io as pio
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly not available. Install with: pip install plotly")

try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("Matplotlib not available. Install with: pip install matplotlib")

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    logger.warning("Seaborn not available. Install with: pip install seaborn")

try:
    import altair as alt
    ALTAIR_AVAILABLE = True
except ImportError:
    ALTAIR_AVAILABLE = False
    logger.warning("Altair not available. Install with: pip install altair")

try:
    import bokeh.plotting as bk
    from bokeh.embed import file_html
    from bokeh.resources import CDN
    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False
    logger.warning("Bokeh not available. Install with: pip install bokeh")

class VisualizationGenerator:
    """
    Universal visualization generator that converts various visualization formats to HTML reports.
    
    Supported Formats:
    - Plotly: Interactive web-based visualizations
    - Matplotlib: Static plots with optional interactive features
    - Seaborn: Statistical visualizations
    - Altair: Declarative statistical visualizations  
    - Bokeh: Interactive web visualizations
    - NetworkX: Network graphs and trees
    - Custom: User-defined visualization formats
    
    Features:
    - Universal format detection and conversion
    - Automatic HTML embedding with proper CDN links
    - Data file integration (CSV, H5, JSON)
    - Rich metadata and analysis summaries
    - Professional responsive design
    - Future-ready for MSA trees, phylogenetic analysis, etc.
    """
    
    def __init__(self, output_dir: str = None):
        """Initialize the visualization generator."""
        if output_dir:
            self.output_dir = output_dir
        else:
            # Use test/outputs for testing, or scratch space for KBase
            if os.path.exists('test/outputs'):
                self.output_dir = 'test/outputs'
            else:
                self.output_dir = os.environ.get('SCRATCH_DIR', '/kb/module/work/tmp')
        
        os.makedirs(self.output_dir, exist_ok=True)
        self.report_metadata = {
            'generator_version': '2.0.0',
            'creation_time': time.strftime('%Y-%m-%d %H:%M:%S'),
            'reports_created': [],
            'supported_formats': {
                'plotly': PLOTLY_AVAILABLE,
                'matplotlib': MATPLOTLIB_AVAILABLE,
                'seaborn': SEABORN_AVAILABLE,
                'altair': ALTAIR_AVAILABLE,
                'bokeh': BOKEH_AVAILABLE
            }
        }
    
    # ==========================================
    # UNIVERSAL VISUALIZATION CONVERSION METHODS
    # ==========================================
    
    def detect_visualization_type(self, visualization_obj: Any) -> str:
        """
        Automatically detect the type of visualization object.
        
        Args:
            visualization_obj: Any visualization object from supported libraries
            
        Returns:
            String indicating the detected format
        """
        if hasattr(visualization_obj, 'to_html') and hasattr(visualization_obj, 'data'):
            # Plotly figure
            return 'plotly'
        elif hasattr(visualization_obj, 'savefig') or str(type(visualization_obj)) == "<class 'matplotlib.figure.Figure'>":
            # Matplotlib figure
            return 'matplotlib'
        elif hasattr(visualization_obj, 'to_json') and hasattr(visualization_obj, 'mark'):
            # Altair chart
            return 'altair'
        elif hasattr(visualization_obj, 'output_backend'):
            # Bokeh plot
            return 'bokeh'
        elif isinstance(visualization_obj, str) and visualization_obj.startswith(('<svg', '<?xml')):
            # SVG string
            return 'svg'
        elif isinstance(visualization_obj, str) and '<html>' in visualization_obj.lower():
            # Pre-rendered HTML
            return 'html'
        else:
            return 'unknown'
    
    def convert_visualization_to_html(self, visualization_obj: Any, 
                                    visualization_type: str = None) -> str:
        """
        Convert any supported visualization to HTML string.
        
        Args:
            visualization_obj: Visualization object from any supported library
            visualization_type: Optional explicit type specification
            
        Returns:
            HTML string containing the embedded visualization
        """
        if visualization_type is None:
            visualization_type = self.detect_visualization_type(visualization_obj)
        
        logger.info(f"Converting {visualization_type} visualization to HTML")
        
        if visualization_type == 'plotly' and PLOTLY_AVAILABLE:
            return self._convert_plotly_to_html(visualization_obj)
        elif visualization_type == 'matplotlib' and MATPLOTLIB_AVAILABLE:
            return self._convert_matplotlib_to_html(visualization_obj)
        elif visualization_type == 'altair' and ALTAIR_AVAILABLE:
            return self._convert_altair_to_html(visualization_obj)
        elif visualization_type == 'bokeh' and BOKEH_AVAILABLE:
            return self._convert_bokeh_to_html(visualization_obj)
        elif visualization_type == 'svg':
            return self._convert_svg_to_html(visualization_obj)
        elif visualization_type == 'html':
            return str(visualization_obj)
        else:
            logger.warning(f"Unsupported visualization type: {visualization_type}")
            return f"<p>Unsupported visualization type: {visualization_type}</p>"
    
    def _convert_plotly_to_html(self, fig) -> str:
        """Convert Plotly figure to HTML."""
        try:
            return fig.to_html(include_plotlyjs='cdn', div_id="plotly-visualization")
        except Exception as e:
            logger.error(f"Failed to convert Plotly figure: {e}")
            return "<p>Error converting Plotly visualization</p>"
    
    def _convert_matplotlib_to_html(self, fig) -> str:
        """Convert Matplotlib figure to HTML with base64 embedding."""
        try:
            buffer = BytesIO()
            fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
            buffer.seek(0)
            img_base64 = base64.b64encode(buffer.getvalue()).decode()
            plt.close(fig)  # Clean up
            
            return f'''
            <div class="matplotlib-container">
                <img src="data:image/png;base64,{img_base64}" 
                     alt="Matplotlib visualization" 
                     style="max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);">
            </div>
            '''
        except Exception as e:
            logger.error(f"Failed to convert Matplotlib figure: {e}")
            return "<p>Error converting Matplotlib visualization</p>"
    
    def _convert_altair_to_html(self, chart) -> str:
        """Convert Altair chart to HTML."""
        try:
            return chart.to_html()
        except Exception as e:
            logger.error(f"Failed to convert Altair chart: {e}")
            return "<p>Error converting Altair visualization</p>"
    
    def _convert_bokeh_to_html(self, plot) -> str:
        """Convert Bokeh plot to HTML."""
        try:
            return file_html(plot, CDN, "Bokeh Visualization")
        except Exception as e:
            logger.error(f"Failed to convert Bokeh plot: {e}")
            return "<p>Error converting Bokeh visualization</p>"
    
    def _convert_svg_to_html(self, svg_string: str) -> str:
        """Convert SVG string to HTML."""
        return f'''
        <div class="svg-container" style="text-align: center; padding: 20px;">
            {svg_string}
        </div>
        '''
    
    def create_universal_visualization_report(self, 
                                            visualization_obj: Any,
                                            title: str = "Visualization Report",
                                            description: str = "Analysis Visualization",
                                            data_files: Dict[str, str] = None,
                                            metadata: Dict[str, Any] = None) -> str:
        """
        Create a comprehensive HTML report for any visualization type.
        
        Args:
            visualization_obj: Any supported visualization object
            title: Report title
            description: Description of the analysis
            data_files: Dictionary of data files {type: path}
            metadata: Additional metadata to include
            
        Returns:
            Path to generated HTML report
        """
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        viz_type = self.detect_visualization_type(visualization_obj)
        
        # Convert visualization to HTML
        viz_html = self.convert_visualization_to_html(visualization_obj)
        
        # Prepare data files section
        data_files_html = ""
        if data_files:
            data_files_html = '<h3 class="section-title"><i class="fas fa-database"></i> Data Files</h3><div class="file-grid">'
            for file_type, file_path in data_files.items():
                filename = os.path.basename(file_path)
                icon = "fas fa-table" if file_type.endswith(('.csv', '.tsv')) else "fas fa-file-alt"
                data_files_html += f'''
                    <a href="{filename}" download class="file-link">
                        <i class="{icon}"></i>
                        {file_type.replace("_", " ").title()}
                    </a>
                '''
            data_files_html += '</div>'
        
        # Prepare metadata section
        metadata_html = ""
        if metadata:
            metadata_html = f'''
            <div class="metadata-section">
                <h3><i class="fas fa-info-circle"></i> Analysis Information</h3>
                <div class="metadata-grid">
                    {''.join([f'<div class="metadata-item"><span class="label">{k}:</span> <span class="value">{v}</span></div>' for k, v in metadata.items()])}
                </div>
            </div>
            '''
        
        # Create comprehensive HTML report
        html_content = f'''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0; padding: 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .main-container {{
            max-width: 1600px; margin: 0 auto; padding: 20px;
        }}
        .report-card {{
            background: white; border-radius: 12px;
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
            overflow: hidden; margin-bottom: 20px;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; padding: 30px; text-align: center;
        }}
        .header h1 {{ margin: 0; font-size: 2.5em; font-weight: 300; }}
        .header .subtitle {{ font-size: 1.1em; opacity: 0.9; margin-top: 10px; }}
        .content {{ padding: 30px; }}
        .visualization-container {{
            background: #f8f9fa; border-radius: 12px; padding: 30px;
            margin: 30px 0; box-shadow: inset 0 2px 10px rgba(0,0,0,0.05);
        }}
        .section-title {{
            font-size: 1.5em; font-weight: 600; margin-bottom: 20px;
            color: #2c3e50; display: flex; align-items: center; gap: 10px;
        }}
        .metadata-section {{
            background: #f0f4f8; border-radius: 12px; padding: 25px; margin: 30px 0;
        }}
        .metadata-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 15px; margin-top: 20px;
        }}
        .metadata-item {{
            background: white; padding: 15px; border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .label {{ font-weight: 600; color: #2c3e50; display: block; margin-bottom: 5px; }}
        .value {{ color: #34495e; }}
        .file-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px; margin-top: 20px;
        }}
        .file-link {{
            display: flex; align-items: center; gap: 10px; padding: 15px 20px;
            background: #1976d2; color: white; text-decoration: none;
            border-radius: 8px; font-weight: 500; transition: all 0.3s ease;
        }}
        .file-link:hover {{ background: #1565c0; transform: translateY(-2px); color: white; }}
        .viz-info {{
            background: #e3f2fd; border-radius: 12px; padding: 20px; margin: 20px 0;
            border-left: 4px solid #1976d2;
        }}
        /* Responsive design */
        @media (max-width: 768px) {{
            .main-container {{ padding: 10px; }}
            .content {{ padding: 15px; }}
            .header h1 {{ font-size: 2em; }}
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="report-card">
            <div class="header">
                <h1><i class="fas fa-chart-line"></i> {title}</h1>
                <div class="subtitle">{description} | Generated on {timestamp}</div>
            </div>
            
            <div class="content">
                <div class="viz-info">
                    <h4><i class="fas fa-info-circle"></i> Visualization Information</h4>
                    <p><strong>Format:</strong> {viz_type.title()}</p>
                    <p><strong>Generated:</strong> {timestamp}</p>
                    <p><strong>Interactive:</strong> {'Yes' if viz_type in ['plotly', 'bokeh', 'altair'] else 'No'}</p>
                </div>
                
                {metadata_html}
                
                <div class="visualization-container">
                    <h3 class="section-title"><i class="fas fa-chart-area"></i> Visualization</h3>
                    {viz_html}
                </div>
                
                {data_files_html}
                
                <div class="instructions" style="background: white; border-radius: 12px; padding: 30px; margin-top: 30px;">
                    <h3 class="section-title"><i class="fas fa-question-circle"></i> How to Use This Report</h3>
                    <ul style="color: #555; line-height: 1.6;">
                        <li><strong>Interactive Features:</strong> {'Use mouse to zoom, pan, and explore data points' if viz_type in ['plotly', 'bokeh', 'altair'] else 'This is a static visualization'}</li>
                        <li><strong>Data Access:</strong> Download raw data files using the links above</li>
                        <li><strong>Export:</strong> Right-click on visualizations to save images</li>
                        <li><strong>Sharing:</strong> This HTML file is self-contained and can be shared</li>
                    </ul>
                </div>
            </div>
        </div>
    </div>
</body>
</html>
        '''
        
        # Save HTML file
        safe_title = "".join(c for c in title if c.isalnum() or c in (' ', '-', '_')).rstrip()
        html_filename = f"{safe_title.replace(' ', '_').lower()}_{int(time.time())}.html"
        html_path = os.path.join(self.output_dir, html_filename)
        
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # Track this report
        self.report_metadata['reports_created'].append({
            'type': 'universal_visualization',
            'format': viz_type,
            'filename': html_filename,
            'title': title,
            'timestamp': timestamp,
            'data_files': len(data_files) if data_files else 0,
            'has_metadata': bool(metadata)
        })
        
        logger.info(f"Universal visualization report saved to: {html_path}")
        return html_path
    
    # ==========================================
    # SPECIALIZED VISUALIZATION METHODS
    # ==========================================
    
    def create_msa_tree_report(self, 
                             tree_visualization: Any,
                             alignment_data: Dict[str, Any] = None,
                             sequence_data: Dict[str, str] = None,
                             title: str = "MSA Phylogenetic Tree",
                             protein_id: str = None) -> str:
        """
        Create an HTML report for Multiple Sequence Alignment (MSA) tree visualizations.
        Future-ready for phylogenetic analysis integration.
        
        Args:
            tree_visualization: Tree visualization object (Plotly, Matplotlib, or custom)
            alignment_data: MSA alignment information and statistics
            sequence_data: Dictionary of sequence_id -> sequence
            title: Report title
            protein_id: Query protein identifier
            
        Returns:
            Path to generated HTML report
        """
        logger.info("Creating MSA tree visualization report")
        
        # Detect and convert tree visualization
        viz_type = self.detect_visualization_type(tree_visualization)
        tree_html = self.convert_visualization_to_html(tree_visualization)
        
        # Prepare MSA-specific metadata
        msa_metadata = {
            'Analysis Type': 'Multiple Sequence Alignment & Phylogeny',
            'Query Protein': protein_id or 'Unknown',
            'Tree Format': viz_type.title(),
            'Sequences Analyzed': len(sequence_data) if sequence_data else 'Unknown'
        }
        
        if alignment_data:
            msa_metadata.update({
                'Alignment Length': alignment_data.get('alignment_length', 'Unknown'),
                'Conservation Score': alignment_data.get('conservation_score', 'Unknown'),
                'Gap Percentage': f"{alignment_data.get('gap_percentage', 0):.1f}%" if 'gap_percentage' in alignment_data else 'Unknown'
            })
        
        # Prepare sequence data section
        sequence_html = ""
        if sequence_data:
            sequence_html = f'''
            <div class="sequences-section">
                <h3 class="section-title"><i class="fas fa-dna"></i> Aligned Sequences</h3>
                <div class="sequence-grid">
                    <div class="sequence-count">Total sequences: {len(sequence_data)}</div>
                    <div class="sequence-preview">
                        <h4>Sequence Preview (first 5):</h4>
                        <div class="sequence-list">
                            {''.join([f'<div class="sequence-item"><strong>{seq_id[:20]}...</strong><br><code>{seq[:80]}...</code></div>' for i, (seq_id, seq) in enumerate(sequence_data.items()) if i < 5])}
                        </div>
                    </div>
                </div>
            </div>
            '''
        
        # Use universal visualization report with MSA-specific enhancements
        return self.create_universal_visualization_report(
            visualization_obj=tree_visualization,
            title=title,
            description=f"Phylogenetic analysis and MSA tree for protein family analysis",
            metadata=msa_metadata
        )
    
    def create_statistical_plot_report(self,
                                     plot_object: Any,
                                     analysis_type: str = "Statistical Analysis",
                                     data_summary: Dict[str, Any] = None,
                                     protein_id: str = None) -> str:
        """
        Create an HTML report for statistical plots (Seaborn, Matplotlib stats, etc.).
        
        Args:
            plot_object: Statistical plot object
            analysis_type: Type of statistical analysis
            data_summary: Summary statistics and data information
            protein_id: Associated protein identifier
            
        Returns:
            Path to generated HTML report
        """
        logger.info(f"Creating {analysis_type} statistical plot report")
        
        # Prepare statistical metadata
        stats_metadata = {
            'Analysis Type': analysis_type,
            'Protein ID': protein_id or 'Multiple proteins',
            'Plot Type': self.detect_visualization_type(plot_object).title()
        }
        
        if data_summary:
            stats_metadata.update(data_summary)
        
        return self.create_universal_visualization_report(
            visualization_obj=plot_object,
            title=f"{analysis_type} Report",
            description=f"Statistical analysis and visualization for protein data",
            metadata=stats_metadata
        )
    
    def create_network_report(self, network_result: Dict[str, Any], 
                            protein_id: str = None, 
                            additional_metadata: Dict[str, Any] = None) -> str:
        """
        Create a comprehensive HTML report for network visualization.
        
        Args:
            network_result: Result from network analysis stage
            protein_id: Protein identifier for naming
            additional_metadata: Additional metadata to include in report
            
        Returns:
            Path to generated HTML file
        """
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # Get the Plotly figure
        fig = network_result.get('visualization_figure')
        if not fig:
            logger.warning("No visualization figure found in network result")
            return None
        
        # Get CSV file links
        csv_files = network_result.get('csv_files', {})
        
        # Generate CSV links HTML with enhanced styling
        csv_links_html = ""
        if csv_files:
            for file_type, file_path in csv_files.items():
                filename = os.path.basename(file_path)
                icon = "fas fa-table" if file_type.endswith("csv") else "fas fa-file-alt"
                csv_links_html += f'''
                    <a href="{filename}" download class="file-link">
                        <i class="{icon}"></i>
                        {file_type.replace("_", " ").title()}
                    </a>
                '''
        else:
            csv_links_html = '<p class="no-files">No data files available</p>'
        
        # Get network properties
        props = network_result.get('network_properties', {})
        
        # Add additional metadata if provided
        metadata_section = ""
        if additional_metadata:
            metadata_section = f"""
            <div class="metadata-section">
                <h3>📋 Analysis Metadata</h3>
                <div class="metadata-grid">
                    {''.join([f'<div class="metadata-item"><span class="label">{k}:</span> <span class="value">{v}</span></div>' for k, v in additional_metadata.items()])}
                </div>
            </div>
            """

        # Create comprehensive HTML content
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Protein Network Analysis Report - {protein_id or 'Analysis'}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css">
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
            max-width: 1600px;
            margin: 0 auto;
            padding: 20px;
        }}
        .report-card {{
            background: white;
            border-radius: 12px;
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
            overflow: hidden;
            margin-bottom: 20px;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
            font-weight: 300;
        }}
        .header .subtitle {{
            font-size: 1.1em;
            opacity: 0.9;
            margin-top: 10px;
        }}
        .content {{
            padding: 30px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            color: white;
            padding: 25px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            transition: transform 0.3s ease;
        }}
        .stat-card:hover {{
            transform: translateY(-5px);
        }}
        .stat-value {{
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        .stat-label {{
            font-size: 0.9em;
            opacity: 0.9;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        .visualization-container {{
            background: #f8f9fa;
            border-radius: 12px;
            padding: 30px;
            margin: 30px 0;
            box-shadow: inset 0 2px 10px rgba(0,0,0,0.05);
        }}
        .section-title {{
            font-size: 1.5em;
            font-weight: 600;
            margin-bottom: 20px;
            color: #2c3e50;
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .data-files {{
            background: #e3f2fd;
            border-radius: 12px;
            padding: 25px;
            margin: 30px 0;
        }}
        .file-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 20px;
        }}
        .file-link {{
            display: flex;
            align-items: center;
            gap: 10px;
            padding: 15px 20px;
            background: #1976d2;
            color: white;
            text-decoration: none;
            border-radius: 8px;
            font-weight: 500;
            transition: all 0.3s ease;
        }}
        .file-link:hover {{
            background: #1565c0;
            transform: translateY(-2px);
            color: white;
        }}
        .metadata-section {{
            background: #f0f4f8;
            border-radius: 12px;
            padding: 25px;
            margin: 30px 0;
        }}
        .metadata-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 15px;
            margin-top: 20px;
        }}
        .metadata-item {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .label {{
            font-weight: 600;
            color: #2c3e50;
            display: block;
            margin-bottom: 5px;
        }}
        .value {{
            color: #34495e;
        }}
        .instructions {{
            background: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 12px;
            padding: 20px;
            margin: 30px 0;
        }}
        .instructions h4 {{
            margin-top: 0;
            color: #856404;
        }}
        .instructions ul {{
            margin-bottom: 0;
            color: #856404;
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="report-card">
            <div class="header">
                <h1><i class="fas fa-project-diagram"></i> Protein Network Analysis Report</h1>
                <div class="subtitle">Generated on {timestamp} | Query: {protein_id or 'Multiple proteins'}</div>
            </div>
            
            <div class="content">
                <div class="stats-grid">
                    <div class="stat-card">
                        <div class="stat-value">{props.get('num_nodes', 0)}</div>
                        <div class="stat-label">Network Nodes</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">{props.get('num_edges', 0)}</div>
                        <div class="stat-label">Network Edges</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">{props.get('density', 0):.3f}</div>
                        <div class="stat-label">Network Density</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">{props.get('connected_components', 0)}</div>
                        <div class="stat-label">Connected Components</div>
                    </div>
                </div>
                
                {metadata_section}
                
                <div class="visualization-container">
                    <h3 class="section-title"><i class="fas fa-chart-network"></i> Interactive Network Visualization</h3>
                    <div id="networkPlot"></div>
                </div>
                
                <div class="data-files">
                    <h3 class="section-title"><i class="fas fa-database"></i> Data Files</h3>
                    <div class="file-grid">
                        {csv_links_html}
                    </div>
                </div>
                
                <div class="instructions">
                    <h4><i class="fas fa-info-circle"></i> How to Use This Report</h4>
                    <ul>
                        <li><strong>Interactive Visualization:</strong> Hover over nodes to see detailed protein information</li>
                        <li><strong>Navigation:</strong> Use mouse wheel to zoom, click and drag to pan</li>
                        <li><strong>Data Access:</strong> Click on data file links above to download CSV files</li>
                        <li><strong>Network Analysis:</strong> Red nodes show proteins connected to your query</li>
                        <li><strong>Toolbar:</strong> Use the visualization toolbar for additional controls</li>
                    </ul>
                </div>
            </div>
        </div>
    </div>
    
    <script>
        // Embed the Plotly figure with enhanced configuration
        var figureData = {fig.to_json()};
        var config = {{
            displayModeBar: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            responsive: true
        }};
        Plotly.newPlot('networkPlot', figureData.data, figureData.layout, config);
        
        // Add resize handler for responsiveness
        window.addEventListener('resize', function() {{
            Plotly.Plots.resize('networkPlot');
        }});
    </script>
</body>
</html>
        """
        
        # Save HTML file
        html_filename = f"network_report_{protein_id or 'analysis'}_{int(time.time())}.html"
        html_path = os.path.join(self.output_dir, html_filename)
        
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # Track this report
        self.report_metadata['reports_created'].append({
            'type': 'network_analysis',
            'filename': html_filename,
            'protein_id': protein_id,
            'timestamp': timestamp,
            'csv_files': len(csv_files),
            'has_metadata': bool(additional_metadata)
        })
        
        logger.info(f"Network analysis report saved to: {html_path}")
        return html_path
    
    def create_pipeline_summary_report(self, pipeline_results: Dict[str, Any], 
                                      pipeline_metadata: Dict[str, Any] = None) -> str:
        """
        Create a comprehensive HTML summary of all pipeline results and data files.
        
        Args:
            pipeline_results: Complete pipeline results
            pipeline_metadata: Additional pipeline metadata
            
        Returns:
            Path to generated HTML summary file
        """
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # Collect all output files systematically
        all_files = []
        analysis_summaries = {}
        
        # Network analysis files
        if 'network_analysis' in pipeline_results:
            network_data = pipeline_results['network_analysis']
            csv_files = network_data.get('csv_files', {})
            
            analysis_summaries['network_analysis'] = {
                'title': 'Network Analysis',
                'description': 'Protein similarity network analysis with interactive visualization',
                'files_count': len(csv_files),
                'has_visualization': bool(network_data.get('html_path'))
            }
            
            for file_type, file_path in csv_files.items():
                all_files.append({
                    'type': 'CSV',
                    'category': 'Network Analysis',
                    'description': file_type.replace('_', ' ').title(),
                    'filename': os.path.basename(file_path),
                    'size': os.path.getsize(file_path) if os.path.exists(file_path) else 0
                })
            
            html_path = network_data.get('html_path')
            if html_path and os.path.exists(html_path):
                all_files.append({
                    'type': 'HTML',
                    'category': 'Network Visualization',
                    'description': 'Interactive Network Plot',
                    'filename': os.path.basename(html_path),
                    'size': os.path.getsize(html_path)
                })
        
        # Sequence analysis files
        if 'sequence_analysis' in pipeline_results:
            seq_data = pipeline_results['sequence_analysis']
            csv_files = seq_data.get('csv_files', {})
            
            analysis_summaries['sequence_analysis'] = {
                'title': 'Sequence Analysis',
                'description': 'Comprehensive protein sequence analysis and properties',
                'files_count': len(csv_files),
                'proteins_analyzed': seq_data.get('analysis_stats', {}).get('total_analyses', 0)
            }
            
            for file_type, file_path in csv_files.items():
                all_files.append({
                    'type': 'CSV',
                    'category': 'Sequence Analysis',
                    'description': file_type.replace('_', ' ').title(),
                    'filename': os.path.basename(file_path),
                    'size': os.path.getsize(file_path) if os.path.exists(file_path) else 0
                })
        
        # Create comprehensive summary HTML
        analysis_cards_html = ""
        for analysis_type, summary in analysis_summaries.items():
            analysis_cards_html += f"""
            <div class="analysis-card">
                <h3>{summary['title']}</h3>
                <p>{summary['description']}</p>
                <div class="card-stats">
                    <span class="stat">📁 {summary['files_count']} files</span>
                    {f'<span class="stat">🔬 {summary.get("proteins_analyzed", 0)} proteins</span>' if 'proteins_analyzed' in summary else ''}
                    {f'<span class="stat">📊 Interactive viz</span>' if summary.get('has_visualization') else ''}
                </div>
            </div>
            """
        
        # Create file table HTML
        file_table_html = ""
        if all_files:
            file_table_html = """
            <div class="files-table">
                <table>
                    <thead>
                        <tr>
                            <th>File Type</th>
                            <th>Category</th>
                            <th>Description</th>
                            <th>Filename</th>
                            <th>Size</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for file_info in sorted(all_files, key=lambda x: x['category']):
                size_str = f"{file_info['size']:,} bytes" if file_info['size'] > 0 else "Unknown"
                file_table_html += f"""
                        <tr>
                            <td><span class="file-type-badge {file_info['type'].lower()}">{file_info['type']}</span></td>
                            <td>{file_info['category']}</td>
                            <td>{file_info['description']}</td>
                            <td><a href="{file_info['filename']}" download class="file-link">{file_info['filename']}</a></td>
                            <td>{size_str}</td>
                        </tr>
                """
            
            file_table_html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Pipeline metadata section
        metadata_html = ""
        if pipeline_metadata:
            metadata_html = f"""
            <div class="metadata-section">
                <h3><i class="fas fa-info-circle"></i> Pipeline Information</h3>
                <div class="metadata-grid">
                    {''.join([f'<div class="metadata-item"><span class="label">{k}:</span> <span class="value">{v}</span></div>' for k, v in pipeline_metadata.items()])}
                </div>
            </div>
            """
        
        # Create comprehensive HTML
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Protein Analysis Pipeline - Complete Report</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css">
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
            padding: 20px;
        }}
        .report-header {{
            background: white;
            border-radius: 12px;
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
            padding: 40px;
            text-align: center;
            margin-bottom: 30px;
        }}
        .report-header h1 {{
            margin: 0;
            font-size: 2.5em;
            color: #2c3e50;
            margin-bottom: 10px;
        }}
        .report-header .subtitle {{
            color: #7f8c8d;
            font-size: 1.1em;
        }}
        .analysis-overview {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .analysis-card {{
            background: white;
            border-radius: 12px;
            padding: 25px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            transition: transform 0.3s ease;
        }}
        .analysis-card:hover {{
            transform: translateY(-5px);
        }}
        .analysis-card h3 {{
            margin: 0 0 10px 0;
            color: #2c3e50;
            font-size: 1.3em;
        }}
        .analysis-card p {{
            color: #7f8c8d;
            margin-bottom: 15px;
        }}
        .card-stats {{
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
        }}
        .stat {{
            background: #e3f2fd;
            color: #1976d2;
            padding: 5px 12px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 500;
        }}
        .files-section {{
            background: white;
            border-radius: 12px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            padding: 30px;
            margin-bottom: 30px;
        }}
        .section-title {{
            font-size: 1.5em;
            font-weight: 600;
            margin-bottom: 20px;
            color: #2c3e50;
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .files-table table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }}
        .files-table th,
        .files-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }}
        .files-table th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #2c3e50;
        }}
        .file-type-badge {{
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.8em;
            font-weight: 600;
            text-transform: uppercase;
        }}
        .file-type-badge.csv {{
            background: #e8f5e8;
            color: #2e7d32;
        }}
        .file-type-badge.html {{
            background: #fff3e0;
            color: #ef6c00;
        }}
        .file-type-badge.json {{
            background: #e3f2fd;
            color: #1976d2;
        }}
        .file-link {{
            color: #1976d2;
            text-decoration: none;
            font-weight: 500;
        }}
        .file-link:hover {{
            text-decoration: underline;
        }}
        .metadata-section {{
            background: #f8f9fa;
            border-radius: 12px;
            padding: 25px;
            margin-bottom: 30px;
        }}
        .metadata-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin-top: 20px;
        }}
        .metadata-item {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .label {{
            font-weight: 600;
            color: #2c3e50;
            display: block;
            margin-bottom: 5px;
        }}
        .value {{
            color: #34495e;
        }}
        .instructions {{
            background: white;
            border-radius: 12px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            padding: 30px;
        }}
        .instructions ul {{
            color: #555;
            line-height: 1.6;
        }}
        .instructions li {{
            margin-bottom: 10px;
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="report-header">
            <h1><i class="fas fa-dna"></i> Protein Analysis Pipeline Report</h1>
            <div class="subtitle">Complete analysis results generated on {timestamp}</div>
        </div>
        
        <div class="analysis-overview">
            {analysis_cards_html}
        </div>
        
        {metadata_html}
        
        <div class="files-section">
            <h2 class="section-title"><i class="fas fa-folder-open"></i> Generated Files</h2>
            {file_table_html}
        </div>
        
        <div class="instructions">
            <h2 class="section-title"><i class="fas fa-question-circle"></i> How to Use These Results</h2>
            <ul>
                <li><strong>CSV Files:</strong> Open with Excel, R, Python, or any data analysis tool</li>
                <li><strong>HTML Reports:</strong> Open directly in your web browser for interactive visualizations</li>
                <li><strong>Network Visualizations:</strong> Explore protein similarity networks with zoom, pan, and hover features</li>
                <li><strong>Sequence Analysis:</strong> Access detailed protein properties, motifs, and physicochemical characteristics</li>
                <li><strong>Data Integration:</strong> All files use consistent protein IDs for easy cross-referencing</li>
            </ul>
        </div>
    </div>
</body>
</html>
        """
        
        # Save summary HTML
        summary_path = os.path.join(self.output_dir, f"pipeline_report_{int(time.time())}.html")
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # Track this report
        self.report_metadata['reports_created'].append({
            'type': 'pipeline_summary',
            'filename': os.path.basename(summary_path),
            'timestamp': timestamp,
            'files_included': len(all_files),
            'analyses_included': len(analysis_summaries)
        })
        
        logger.info(f"Comprehensive pipeline report saved to: {summary_path}")
        return summary_path
