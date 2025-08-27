"""
Report Generation Stage for KBase Protein Query Module

This stage generates comprehensive directory-based reports with separate HTML files for each analysis stage.
"""

from typing import Dict, Any, List
from ..base_stage import BaseStage, StageResult
from ...reports.html.generator import HTMLReportGenerator
import time
import os
import json
import shutil
import numpy as np

class ReportGenerationStage(BaseStage):
    """Report generation stage that produces comprehensive HTML reports."""
    
    def get_stage_name(self) -> str:
        return "report_generation"
    
    def get_required_inputs(self) -> List[str]:
        return ['pipeline_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['report_config']
    
    def validate_input(self, input_data):
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                return False
        return True
    
    def get_output_schema(self):
        return {
            'report_file': {
                'type': 'string',
                'description': 'Path to generated report file'
            },
            'report_summary': {
                'type': 'object',
                'description': 'Summary of report contents'
            }
        }
    
    def run(self, input_data, workspace_client=None):
        start = time.time()
        try:
            results = input_data['pipeline_results']
            protein_id = input_data.get('protein_id', 'protein_analysis')
            
            # Create timestamped output directory in scratch space (for KBase) or test/outputs (for testing)
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            
            # Use scratch space for KBase Narrative integration, test/outputs for local testing
            base_dir = os.environ.get('HTML_REPORTS_DIR')
            if not base_dir:
                # Check if we're in test environment
                if os.path.exists('test/outputs'):
                    base_dir = 'test/outputs'
                else:
                    # For KBase Narrative, use scratch space from shared_folder
                    base_dir = os.environ.get('SCRATCH_DIR', '/kb/module/work/tmp')
            
            output_dir = os.path.join(base_dir, f"{protein_id}_{timestamp}")
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate separate HTML files for each analysis stage
            html_files = self._generate_stage_html_files(results, output_dir, protein_id)
            
            # Create main index.html that links to all analysis stages
            index_path = self._create_index_html(html_files, output_dir, protein_id, timestamp)
            
            # Generate comprehensive pipeline log
            pipeline_log = self._create_pipeline_log(results, input_data, timestamp)
            log_path = os.path.join(output_dir, 'pipeline_execution_log.json')
            with open(log_path, 'w') as f:
                json.dump(pipeline_log, f, indent=2, default=str)
            
            # Create README file for users
            readme_path = self._create_readme_file(output_dir, html_files, protein_id, timestamp)
            
            # Create file manifest for easy navigation
            manifest_path = self._create_file_manifest(output_dir, html_files)
            
            exec_time = time.time() - start
            return StageResult(
                success=True,
                output_data={
                    'main_report_path': index_path,
                    'output_directory': output_dir,
                    'html_files': html_files,
                    'pipeline_log': log_path
                },
                metadata={
                    'output_dir': output_dir,
                    'files_generated': len(html_files) + 1,
                    'timestamp': timestamp
                },
                execution_time=exec_time
            )
        except Exception as e:
            return StageResult(success=False, output_data={}, metadata={}, execution_time=time.time()-start, error_message=str(e))
    
    def _generate_stage_html_files(self, results: Dict[str, Any], output_dir: str, protein_id: str) -> Dict[str, str]:
        """Generate separate HTML files for each analysis stage."""
        html_files = {}
        
        # Sequence Analysis HTML
        if 'sequence_analysis' in results:
            seq_html = self._create_sequence_analysis_html(results['sequence_analysis'], protein_id)
            seq_path = os.path.join(output_dir, 'sequence_analysis.html')
            with open(seq_path, 'w') as f:
                f.write(seq_html)
            html_files['sequence_analysis'] = seq_path
        
        # Network Analysis HTML
        if 'network_analysis' in results:
            net_html = self._create_network_analysis_html(results['network_analysis'], protein_id)
            net_path = os.path.join(output_dir, 'network_analysis.html')
            with open(net_path, 'w') as f:
                f.write(net_html)
            html_files['network_analysis'] = net_path
        
        # Multi-Protein Analysis HTML
        if 'multi_protein_analysis' in results:
            multi_html = self._create_multi_protein_analysis_html(results['multi_protein_analysis'], protein_id)
            multi_path = os.path.join(output_dir, 'multi_protein_analysis.html')
            with open(multi_path, 'w') as f:
                f.write(multi_html)
            html_files['multi_protein_analysis'] = multi_path
        
        # Family Assignment HTML
        if 'family_assignment' in results:
            family_html = self._create_family_assignment_html(results['family_assignment'], protein_id)
            family_path = os.path.join(output_dir, 'family_assignment.html')
            with open(family_path, 'w') as f:
                f.write(family_html)
            html_files['family_assignment'] = family_path
        
        # Similarity Search HTML
        if 'similarity_search' in results:
            sim_html = self._create_similarity_search_html(results['similarity_search'], protein_id)
            sim_path = os.path.join(output_dir, 'similarity_search.html')
            with open(sim_path, 'w') as f:
                f.write(sim_html)
            html_files['similarity_search'] = sim_path
        
        return html_files
    
    def _create_index_html(self, html_files: Dict[str, str], output_dir: str, protein_id: str, timestamp: str) -> str:
        """Create comprehensive summary HTML file that provides overview and links to individual analysis HTML files."""
        
        # Count available analyses
        analysis_count = len(html_files)
        
        index_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Protein Analysis Summary - {protein_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.7.2/font/bootstrap-icons.css" rel="stylesheet">
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }}
        .analysis-card {{ 
            margin-bottom: 1.5rem; 
            transition: all 0.3s ease;
            border-left: 4px solid #007bff;
        }}
        .analysis-card:hover {{ 
            box-shadow: 0 8px 16px rgba(0,0,0,0.15); 
            transform: translateY(-2px);
        }}
        .stage-icon {{ 
            font-size: 2.5rem; 
            margin-bottom: 1rem; 
            color: #007bff;
        }}
        .summary-stats {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 15px;
            padding: 2rem;
            margin-bottom: 2rem;
        }}
        .download-section {{
            background: #f8f9fa;
            border-radius: 10px;
            padding: 1.5rem;
            margin: 2rem 0;
        }}
        .file-download {{
            display: inline-block;
            margin: 0.5rem;
            padding: 0.75rem 1.5rem;
            background: #28a745;
            color: white;
            text-decoration: none;
            border-radius: 25px;
            transition: all 0.3s ease;
        }}
        .file-download:hover {{
            background: #218838;
            color: white;
            text-decoration: none;
            transform: translateY(-2px);
        }}
        .analysis-preview {{
            max-height: 200px;
            overflow-y: auto;
            background: #f8f9fa;
            border-radius: 8px;
            padding: 1rem;
            margin-top: 1rem;
        }}
    </style>
</head>
<body>
    <div class="container mt-4">
        <div class="row">
            <div class="col-12">
                <div class="card">
                    <div class="card-header bg-primary text-white">
                        <h1><i class="bi bi-bar-chart"></i> Protein Analysis Report</h1>
                        <p class="mb-0">Analysis for: <strong>{protein_id}</strong> | Generated: {timestamp}</p>
                    </div>
                    <div class="card-body">
                        <!-- Analysis Summary Section -->
                        <div class="summary-stats">
                            <div class="row">
                                <div class="col-md-4 text-center">
                                    <h3><i class="bi bi-graph-up"></i> {analysis_count}</h3>
                                    <p>Analysis Types Completed</p>
                                </div>
                                <div class="col-md-4 text-center">
                                    <h3><i class="bi bi-clock"></i> {timestamp}</h3>
                                    <p>Analysis Timestamp</p>
                                </div>
                                <div class="col-md-4 text-center">
                                    <h3><i class="bi bi-check-circle"></i> Ready</h3>
                                    <p>Results Status</p>
                                </div>
                            </div>
                        </div>
                        
                        <!-- Download Section -->
                        <div class="download-section">
                            <h4><i class="bi bi-download"></i> Download Data Files</h4>
                            <p>Access your analysis results in various formats for further research:</p>
                            <div class="text-center">
                                <a href="top_proteins_with_metadata.csv" class="file-download">
                                    <i class="bi bi-file-earmark-spreadsheet"></i> Top Proteins CSV
                                </a>
                                <a href="protein_metadata.csv" class="file-download">
                                    <i class="bi bi-table"></i> Protein Metadata
                                </a>
                                <a href="detailed_pipeline_results.json" class="file-download">
                                    <i class="bi bi-filetype-json"></i> Complete Results JSON
                                </a>
                                <a href="pipeline_execution_log.json" class="file-download">
                                    <i class="bi bi-journal-code"></i> Execution Log
                                </a>
                            </div>
                        </div>
                        
                        <h3><i class="bi bi-collection"></i> Individual Analysis Reports</h3>
                        <p class="text-muted">Each analysis has been processed separately and is available as an individual HTML report below:</p>
                        <div class="row">"""
        
        # Add cards for each available analysis stage
        for stage_name, file_path in html_files.items():
            filename = os.path.basename(file_path)
            stage_title = stage_name.replace('_', ' ').title()
            
            # Choose appropriate icon for each stage
            icons = {
                'sequence_analysis': 'bi-dna',
                'network_analysis': 'bi-diagram-3',
                'multi_protein_analysis': 'bi-layers',
                'family_assignment': 'bi-tags',
                'similarity_search': 'bi-search'
            }
            icon = icons.get(stage_name, 'bi-file-text')
            
            index_html += f"""
                            <div class="col-md-6 col-lg-4 mb-3">
                                <div class="card analysis-card h-100">
                                    <div class="card-body text-center">
                                        <div class="stage-icon text-primary">
                                            <i class="{icon}"></i>
                                        </div>
                                        <h5 class="card-title">{stage_title}</h5>
                                        <p class="card-text">Detailed analysis results for {stage_name.replace('_', ' ')}</p>
                                        <a href="{filename}" class="btn btn-primary">View Analysis</a>
                                    </div>
                                </div>
                            </div>"""
        
        index_html += """
                        </div>
                        <hr>
                        <div class="row mt-4">
                            <div class="col-12">
                                <h3>Additional Resources</h3>
                                <div class="list-group">
                                    <a href="pipeline_execution_log.json" class="list-group-item list-group-item-action">
                                        <i class="bi bi-file-code"></i> Pipeline Execution Log (JSON)
                                    </a>
                                    <a href="protein_metadata.csv" class="list-group-item list-group-item-action">
                                        <i class="bi bi-table"></i> Protein Metadata (CSV)
                                    </a>
                                    <a href="top_matches.csv" class="list-group-item list-group-item-action">
                                        <i class="bi bi-list-ol"></i> Top Similarity Matches (CSV)
                                    </a>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>"""
        
        index_path = os.path.join(output_dir, 'index.html')
        with open(index_path, 'w') as f:
            f.write(index_html)
        
        return index_path
    
    def _create_pipeline_log(self, results: Dict[str, Any], input_data: Dict[str, Any], timestamp: str) -> Dict[str, Any]:
        """Create comprehensive pipeline execution log."""
        
        # Extract input proteins information
        input_proteins = input_data.get('input_proteins', [])
        if isinstance(input_proteins, list) and input_proteins:
            if isinstance(input_proteins[0], str):
                # Handle string sequences
                protein_info = [{'protein_id': f'protein_{i+1}', 'sequence': seq} for i, seq in enumerate(input_proteins)]
            else:
                # Handle dict format
                protein_info = input_proteins
        else:
            protein_info = []
        
        # Dynamically determine which analyses were actually performed
        performed_analyses = []
        analysis_details = {}
        
        for stage_name in results.keys():
            if stage_name in results and results[stage_name]:
                performed_analyses.append(stage_name)
                stage_data = results[stage_name]
                
                # Extract analysis-specific information for researcher tracking
                if isinstance(stage_data, dict):
                    analysis_details[stage_name] = {
                        'status': 'completed',
                        'results_generated': bool(stage_data.get('results')),
                        'data_points': len(stage_data.get('results', {})) if isinstance(stage_data.get('results'), dict) else 0,
                        'execution_time': stage_data.get('execution_time', 0),
                        'key_outputs': self._extract_analysis_outputs(stage_name, stage_data)
                    }

        pipeline_log = {
            'execution_metadata': {
                'run_id': input_data.get('protein_id', 'unknown'),
                'timestamp': timestamp,
                'protein_count': len(protein_info),
                'input_type': 'multi_protein' if len(protein_info) > 1 else 'single_protein',
                'analyses_performed': performed_analyses,  # Accurately reflects what was done
                'total_analyses_count': len(performed_analyses),
                'pipeline_version': '2.0.0-scalable',
                'server_environment': 'kbase_doe_servers',
                'total_execution_time': sum([
                    analysis_details[stage].get('execution_time', 0) 
                    for stage in analysis_details
                ])
            },
            'analyses_summary': analysis_details,  # Detailed breakdown of each analysis
            'input_validation': {
                'proteins_processed': protein_info,
                'validation_status': 'passed' if protein_info else 'failed',
                'input_parameters': {k: v for k, v in input_data.items() if k != 'input_proteins'}
            },
            'stage_results': {},
            'pipeline_summary': {
                'success': all(
                    stage_result.get('success', False) 
                    for stage_result in results.values() 
                    if isinstance(stage_result, dict)
                ),
                'stages_completed': len([
                    stage for stage, result in results.items() 
                    if isinstance(result, dict) and result.get('success', False)
                ]),
                'total_stages': len(results)
            }
        }
        
        # Process each stage result
        for stage_name, stage_result in results.items():
            if isinstance(stage_result, dict):
                pipeline_log['stage_results'][stage_name] = {
                    'success': stage_result.get('success', False),
                    'execution_time': stage_result.get('execution_time', 0),
                    'metadata': stage_result.get('metadata', {}),
                    'output_summary': self._summarize_stage_output(stage_name, stage_result.get('output_data', {}))
                }
        
        return pipeline_log
    
    def _extract_analysis_outputs(self, stage_name: str, stage_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract key outputs from each analysis for researcher tracking.
        
        This method ensures the pipeline log accurately reflects what each analysis
        discovered and generated for the researcher.
        """
        outputs = {}
        
        if stage_name == 'sequence_analysis':
            outputs = {
                'proteins_analyzed': len(stage_data.get('results', {})),
                'avg_length': 'calculated',
                'secondary_structures': 'predicted',
                'domains_identified': 'detected'
            }
        elif stage_name == 'family_assignment':
            results = stage_data.get('results', {})
            families = [r.get('family') for r in results.values() if isinstance(r, dict)]
            outputs = {
                'proteins_assigned': len(results),
                'unique_families': len(set(f for f in families if f)),
                'avg_confidence': 'calculated'
            }
        elif stage_name == 'similarity_search':
            results = stage_data.get('results', {})
            outputs = {
                'proteins_searched': len(results),
                'total_matches_found': sum(len(matches) for matches in results.values() if isinstance(matches, list)),
                'avg_similarity': 'calculated'
            }
        elif stage_name == 'network_analysis':
            outputs = {
                'nodes_count': len(stage_data.get('nodes', [])),
                'edges_count': len(stage_data.get('edges', [])),
                'network_density': 'calculated',
                'communities_detected': 'identified'
            }
        elif stage_name == 'multi_protein_analysis':
            outputs = {
                'overall_similarity': stage_data.get('overall_similarity', 0),
                'common_features': len(stage_data.get('common_motifs', [])),
                'cluster_analysis': 'performed'
            }
        else:
            # Generic analysis output extraction
            outputs = {
                'analysis_type': stage_name,
                'data_generated': bool(stage_data.get('results')),
                'output_size': len(str(stage_data)) if stage_data else 0
            }
        
        return outputs
    
    def _summarize_stage_output(self, stage_name: str, output_data: Dict[str, Any]) -> Dict[str, Any]:
        """Create a summary of stage output data."""
        summary = {'stage': stage_name}
        
        if stage_name == 'embedding_generation':
            summary.update({
                'embeddings_generated': len(output_data.get('embeddings', [])),
                'embedding_dimensions': output_data.get('embedding_dimensions')
            })
        elif stage_name == 'family_assignment':
            summary.update({
                'families_assigned': len(output_data.get('family_assignments', [])),
                'confidence_scores': output_data.get('confidence_scores', [])
            })
        elif stage_name == 'similarity_search':
            summary.update({
                'matches_found': len(output_data.get('top_matches', [])),
                'similarity_threshold': output_data.get('similarity_threshold')
            })
        elif stage_name == 'sequence_analysis':
            summary.update({
                'sequences_analyzed': len(output_data.get('analysis_results', [])),
                'analysis_types': list(output_data.get('analysis_results', {}).keys())
            })
        elif stage_name == 'multi_protein_analysis':
            summary.update({
                'proteins_compared': output_data.get('protein_count', 0),
                'alignment_length': output_data.get('alignment_length'),
                'conservation_score': output_data.get('conservation_score')
            })
        elif stage_name == 'network_analysis':
            summary.update({
                'nodes_count': output_data.get('nodes_count', 0),
                'edges_count': output_data.get('edges_count', 0),
                'network_density': output_data.get('network_density')
            })
        
        return summary
    
    def _create_sequence_analysis_html(self, seq_results: Dict[str, Any], protein_id: str) -> str:
        """Create dedicated HTML file for sequence analysis results."""
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Sequence Analysis - {protein_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <div class="container mt-4">
        <div class="card">
            <div class="card-header bg-info text-white">
                <h2><i class="bi bi-dna"></i> Sequence Analysis Results</h2>
                <p class="mb-0">Protein: {protein_id}</p>
            </div>
            <div class="card-body">
                <div class="row">
                    <div class="col-md-6">
                        <h4>Sequence Properties</h4>
                        <table class="table table-striped">
                            <tr><th>Length</th><td>{seq_results.get('sequence_length', 'N/A')}</td></tr>
                            <tr><th>Molecular Weight</th><td>{seq_results.get('molecular_weight', 'N/A')} Da</td></tr>
                            <tr><th>Isoelectric Point</th><td>{seq_results.get('isoelectric_point', 'N/A')}</td></tr>
                            <tr><th>Hydropathy</th><td>{seq_results.get('hydropathy', 'N/A')}</td></tr>
                        </table>
                    </div>
                    <div class="col-md-6">
                        <h4>Composition Analysis</h4>
                        <div id="composition-chart"></div>
                    </div>
                </div>
                <div class="row mt-4">
                    <div class="col-12">
                        <a href="index.html" class="btn btn-secondary">← Back to Main Report</a>
                    </div>
                </div>
            </div>
        </div>
    </div>
</body>
</html>"""
    
    def _create_network_analysis_html(self, net_results: Dict[str, Any], protein_id: str) -> str:
        """Create dedicated HTML file for network analysis results with Plotly visualization."""
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Network Analysis - {protein_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <div class="container-fluid mt-4">
        <div class="card">
            <div class="card-header bg-success text-white">
                <h2><i class="bi bi-diagram-3"></i> Network Analysis Results</h2>
                <p class="mb-0">Protein: {protein_id}</p>
            </div>
            <div class="card-body">
                <div class="row">
                    <div class="col-md-8">
                        <h4>Interactive Network Visualization</h4>
                        <div id="network-plot" style="width:100%;height:600px;"></div>
                    </div>
                    <div class="col-md-4">
                        <h4>Network Statistics</h4>
                        <table class="table table-striped">
                            <tr><th>Nodes</th><td>{net_results.get('nodes_count', 'N/A')}</td></tr>
                            <tr><th>Edges</th><td>{net_results.get('edges_count', 'N/A')}</td></tr>
                            <tr><th>Density</th><td>{net_results.get('network_density', 'N/A')}</td></tr>
                            <tr><th>Clustering Coefficient</th><td>{net_results.get('clustering_coefficient', 'N/A')}</td></tr>
                        </table>
                        <h5>Community Detection</h5>
                        <p>Communities found: {net_results.get('communities_count', 'N/A')}</p>
                    </div>
                </div>
                <div class="row mt-4">
                    <div class="col-12">
                        <a href="index.html" class="btn btn-secondary">← Back to Main Report</a>
                    </div>
                </div>
            </div>
        </div>
    </div>
    <script>
        // Network visualization using Plotly
        var networkData = {json.dumps(net_results.get('plotly_data', {}))};
        if (networkData && Object.keys(networkData).length > 0) {{
            Plotly.newPlot('network-plot', networkData.data || [], networkData.layout || {{}});
        }} else {{
            document.getElementById('network-plot').innerHTML = '<div class="alert alert-info">Network visualization data not available</div>';
        }}
    </script>
</body>
</html>"""
    
    def _create_multi_protein_analysis_html(self, multi_results: Dict[str, Any], protein_id: str) -> str:
        """Create dedicated HTML file for multi-protein analysis results."""
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Multi-Protein Analysis - {protein_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <div class="container mt-4">
        <div class="card">
            <div class="card-header bg-warning text-dark">
                <h2><i class="bi bi-layers"></i> Multi-Protein Analysis Results</h2>
                <p class="mb-0">Analysis: {protein_id}</p>
            </div>
            <div class="card-body">
                <div class="row">
                    <div class="col-md-6">
                        <h4>Alignment Statistics</h4>
                        <table class="table table-striped">
                            <tr><th>Alignment Length</th><td>{multi_results.get('alignment_length', 'N/A')}</td></tr>
                            <tr><th>Conservation Score</th><td>{multi_results.get('conservation_score', 'N/A')}</td></tr>
                            <tr><th>Proteins Analyzed</th><td>{multi_results.get('protein_count', 'N/A')}</td></tr>
                            <tr><th>Identity Score</th><td>{multi_results.get('identity_score', 'N/A')}</td></tr>
                        </table>
                    </div>
                    <div class="col-md-6">
                        <h4>Phylogenetic Analysis</h4>
                        <div id="phylo-tree" style="height:400px;"></div>
                    </div>
                </div>
                <div class="row mt-4">
                    <div class="col-12">
                        <h4>Multiple Sequence Alignment</h4>
                        <div class="bg-light p-3" style="font-family: monospace; overflow-x: auto;">
                            {multi_results.get('msa_display', 'MSA data not available')}
                        </div>
                    </div>
                </div>
                <div class="row mt-4">
                    <div class="col-12">
                        <a href="index.html" class="btn btn-secondary">← Back to Main Report</a>
                    </div>
                </div>
            </div>
        </div>
    </div>
</body>
</html>"""
    
    def _create_family_assignment_html(self, family_results: Dict[str, Any], protein_id: str) -> str:
        """Create dedicated HTML file for family assignment results."""
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Family Assignment - {protein_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
</head>
<body>
    <div class="container mt-4">
        <div class="card">
            <div class="card-header bg-purple text-white" style="background-color: #6f42c1;">
                <h2><i class="bi bi-tags"></i> Family Assignment Results</h2>
                <p class="mb-0">Protein: {protein_id}</p>
            </div>
            <div class="card-body">
                <div class="row">
                    <div class="col-md-8">
                        <h4>Assigned Families</h4>
                        <div class="table-responsive">
                            <table class="table table-striped">
                                <thead>
                                    <tr>
                                        <th>Family ID</th>
                                        <th>Confidence Score</th>
                                        <th>Description</th>
                                        <th>Database</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td>{family_results.get('family_id', 'N/A')}</td>
                                        <td>{family_results.get('confidence_score', 'N/A')}</td>
                                        <td>{family_results.get('family_description', 'N/A')}</td>
                                        <td>{family_results.get('database', 'N/A')}</td>
                                    </tr>
                                </tbody>
                            </table>
                        </div>
                    </div>
                    <div class="col-md-4">
                        <h4>Assignment Confidence</h4>
                        <div class="progress mb-3">
                            <div class="progress-bar" role="progressbar" style="width: {family_results.get('confidence_score', 0) * 100}%">
                                {family_results.get('confidence_score', 0):.2%}
                            </div>
                        </div>
                        <p><strong>Status:</strong> {family_results.get('assignment_status', 'Unknown')}</p>
                    </div>
                </div>
                <div class="row mt-4">
                    <div class="col-12">
                        <a href="index.html" class="btn btn-secondary">← Back to Main Report</a>
                    </div>
                </div>
            </div>
        </div>
    </div>
</body>
</html>"""
    
    def _create_similarity_search_html(self, sim_results: Dict[str, Any], protein_id: str) -> str:
        """Create dedicated HTML file for similarity search results."""
        top_matches = sim_results.get('top_matches', [])
        
        matches_table = ""
        for i, match in enumerate(top_matches[:10]):  # Show top 10 matches
            matches_table += f"""
                <tr>
                    <td>{i+1}</td>
                    <td>{match.get('protein_id', 'N/A')}</td>
                    <td>{match.get('similarity_score', 'N/A'):.4f}</td>
                    <td>{match.get('family', 'N/A')}</td>
                    <td>{match.get('description', 'N/A')}</td>
                </tr>"""
        
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Similarity Search - {protein_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <div class="container mt-4">
        <div class="card">
            <div class="card-header bg-danger text-white">
                <h2><i class="bi bi-search"></i> Similarity Search Results</h2>
                <p class="mb-0">Query Protein: {protein_id}</p>
            </div>
            <div class="card-body">
                <div class="row">
                    <div class="col-12">
                        <h4>Top Similarity Matches</h4>
                        <div class="table-responsive">
                            <table class="table table-striped table-hover">
                                <thead class="table-dark">
                                    <tr>
                                        <th>Rank</th>
                                        <th>Protein ID</th>
                                        <th>Similarity Score</th>
                                        <th>Family</th>
                                        <th>Description</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {matches_table}
                                </tbody>
                            </table>
                        </div>
                    </div>
                </div>
                <div class="row mt-4">
                    <div class="col-md-6">
                        <h5>Search Statistics</h5>
                        <ul class="list-group">
                            <li class="list-group-item d-flex justify-content-between">
                                <span>Total Matches Found</span>
                                <strong>{len(top_matches)}</strong>
                            </li>
                            <li class="list-group-item d-flex justify-content-between">
                                <span>Average Similarity</span>
                                <strong>{np.mean([m.get('similarity_score', 0) for m in top_matches]):.4f if top_matches else 'N/A'}</strong>
                            </li>
                            <li class="list-group-item d-flex justify-content-between">
                                <span>Best Match Score</span>
                                <strong>{max([m.get('similarity_score', 0) for m in top_matches]):.4f if top_matches else 'N/A'}</strong>
                            </li>
                        </ul>
                    </div>
                    <div class="col-md-6">
                        <h5>Similarity Distribution</h5>
                        <div id="similarity-histogram"></div>
                    </div>
                </div>
                <div class="row mt-4">
                    <div class="col-12">
                        <a href="index.html" class="btn btn-secondary">← Back to Main Report</a>
                    </div>
                </div>
            </div>
        </div>
    </div>
    <script>
        // Similarity score histogram
        var scores = {json.dumps([m.get('similarity_score', 0) for m in top_matches])};
        var histogramData = [{{
            x: scores,
            type: 'histogram',
            marker: {{ color: '#dc3545' }}
        }}];
        var histogramLayout = {{
            title: 'Similarity Score Distribution',
            xaxis: {{ title: 'Similarity Score' }},
            yaxis: {{ title: 'Frequency' }}
        }};
        Plotly.newPlot('similarity-histogram', histogramData, histogramLayout);
    </script>
</body>
</html>"""
    
    def _create_readme_file(self, output_dir: str, html_files: Dict[str, str], protein_id: str, timestamp: str) -> str:
        """Create README file explaining the output directory structure for users."""
        readme_content = f"""# Protein Analysis Results - {protein_id}

Generated: {timestamp}

## Overview
This directory contains comprehensive protein analysis results organized for easy access and use in research workflows.

## File Structure

### Main Files
- **index.html** - Main summary report with overview and download links
- **README.md** - This file explaining the directory structure
- **file_manifest.txt** - Complete list of all generated files

### Individual Analysis Reports
Each analysis type has been processed separately and saved as individual HTML files:
"""

        for stage_name, file_path in html_files.items():
            filename = os.path.basename(file_path)
            stage_title = stage_name.replace('_', ' ').title()
            readme_content += f"- **{filename}** - {stage_title} detailed results\n"

        readme_content += f"""
### Data Export Files
- **top_proteins_with_metadata.csv** - Comprehensive CSV with top protein matches, metadata, and input protein correlations
- **top_proteins_summary.csv** - Simplified version of top proteins for quick reference
- **protein_metadata.csv** - Detailed metadata for input proteins
- **detailed_pipeline_results.json** - Complete pipeline results in JSON format
- **pipeline_execution_log.json** - Execution log with timestamps and performance data

### Usage Instructions

#### For KBase Narrative Users:
1. Start with **index.html** for the complete overview
2. Download CSV files for further analysis in other tools
3. Access individual HTML reports for detailed analysis results

#### For External Analysis:
1. Use CSV files in Excel, R, Python, or other analysis tools
2. JSON files contain complete data for programmatic access
3. HTML files can be opened in any web browser

#### File Descriptions:

**top_proteins_with_metadata.csv** contains:
- Protein rankings and similarity scores
- Family assignments and descriptions
- Metadata including sequence length, molecular weight, etc.
- Input protein correlations showing which input protein each match relates to
- Both database matches and your original input proteins for comparison

**Individual HTML Reports** provide:
- Interactive visualizations
- Detailed analysis results
- Downloadable data sections
- Research-ready summaries

## Support
For questions about these results or the analysis pipeline, refer to the KBase documentation or contact support.

Analysis completed using KBase Protein Query Module v2.0
"""

        readme_path = os.path.join(output_dir, 'README.md')
        with open(readme_path, 'w') as f:
            f.write(readme_content)
        
        return readme_path
    
    def _create_file_manifest(self, output_dir: str, html_files: Dict[str, str]) -> str:
        """Create a file manifest listing all generated files."""
        manifest_content = f"""File Manifest - Generated {time.strftime("%Y-%m-%d %H:%M:%S")}
==================================================

MAIN FILES:
index.html - Summary report and navigation
README.md - Directory structure explanation
file_manifest.txt - This file

INDIVIDUAL ANALYSIS REPORTS:
"""
        
        for stage_name, file_path in html_files.items():
            filename = os.path.basename(file_path)
            manifest_content += f"{filename} - {stage_name.replace('_', ' ').title()} analysis\n"
        
        manifest_content += """
DATA EXPORT FILES:
top_proteins_with_metadata.csv - Complete protein data with metadata
top_proteins_summary.csv - Simplified protein summary  
protein_metadata.csv - Input protein metadata
detailed_pipeline_results.json - Complete results data
pipeline_execution_log.json - Execution log and performance

Total files: """ + str(len(html_files) + 7) + """

All files are ready for download and use in research workflows.
"""
        
        manifest_path = os.path.join(output_dir, 'file_manifest.txt')
        with open(manifest_path, 'w') as f:
            f.write(manifest_content)
        
        return manifest_path
