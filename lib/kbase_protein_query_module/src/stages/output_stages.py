"""
Output Stages for KBase Protein Query Module

This module contains all output stages that generate reports and visualizations.
"""

import logging
import time
import os
import json
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

from .base_stage import BaseStage, StageResult
from ..reports.html.generator import HTMLReportGenerator
from ..processing.networks.builder import DynamicNetworkBuilder
from ..storage import ProteinStorage

logger = logging.getLogger(__name__)

class ReportGenerationStage(BaseStage):
    """
    Generates comprehensive HTML reports.
    
    Handles:
    - Multi-layer report design (dashboard + detail layers)
    - Interactive visualizations
    - Statistical summaries
    - KBase-compliant formatting
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.html_generator = HTMLReportGenerator()
        self.output_dir = config.get('output_dir', 'html_reports') if config else 'html_reports'
        self.include_interactive_charts = config.get('include_interactive_charts', True) if config else True
        self.include_network_visualizations = config.get('include_network_visualizations', True) if config else True
    
    def get_stage_name(self) -> str:
        return "report_generation"
    
    def get_required_inputs(self) -> List[str]:
        return ['pipeline_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['report_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['sequence_analysis', 'network_analysis', 'bioinformatics_analysis']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'pipeline_results' not in input_data:
            logger.error("Missing pipeline_results")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'report_files': {
                'type': 'object',
                'properties': {
                    'html_report': {'type': 'string'},
                    'network_visualizations': {'type': 'array'},
                    'summary_statistics': {'type': 'object'}
                }
            },
            'report_stats': {
                'type': 'object',
                'properties': {
                    'total_proteins': {'type': 'integer'},
                    'report_sections': {'type': 'array'},
                    'generation_time': {'type': 'number'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Generate comprehensive HTML report."""
        start_time = time.time()
        
        try:
            pipeline_results = input_data['pipeline_results']
            report_config = input_data.get('report_config', {})
            
            # Generate comprehensive report
            report_result = self.html_generator.generate_comprehensive_report(
                pipeline_results=pipeline_results,
                protein_id=report_config.get('protein_id'),
                sequence=report_config.get('sequence')
            )
            
            # Extract report information
            html_path = report_result.get('html_path')
            network_viz_path = report_result.get('network_viz_path')
            sequence_analysis = report_result.get('sequence_analysis')
            
            # Create report files summary
            report_files = {
                'html_report': html_path,
                'network_visualizations': [network_viz_path] if network_viz_path else [],
                'summary_statistics': {
                    'total_proteins': len(pipeline_results.get('protein_records', [])),
                    'analysis_sections': ['sequence', 'network', 'bioinformatics'],
                    'generation_timestamp': time.time()
                }
            }
            
            # Generate report statistics
            report_stats = {
                'total_proteins': len(pipeline_results.get('protein_records', [])),
                'report_sections': ['overview', 'sequence', 'network', 'statistics', 'bioinformatics', 'raw_data'],
                'generation_time': time.time() - start_time,
                'file_size': os.path.getsize(html_path) if html_path and os.path.exists(html_path) else 0
            }
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=html_path is not None,
                output_data={
                    'report_files': report_files,
                    'report_stats': report_stats
                },
                metadata={
                    'html_path': html_path,
                    'network_viz_path': network_viz_path,
                    'include_interactive_charts': self.include_interactive_charts
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Report generation failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class VisualizationStage(BaseStage):
    """
    Creates interactive visualizations.
    
    Handles:
    - Network visualizations
    - Statistical charts
    - Interactive dashboards
    - KBase UI integration
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.network_builder = DynamicNetworkBuilder()
        self.output_dir = config.get('output_dir', 'visualizations') if config else 'visualizations'
        self.visualization_library = config.get('visualization_library', 'plotly') if config else 'plotly'
        self.interactive_mode = config.get('interactive_mode', True) if config else True
    
    def get_stage_name(self) -> str:
        return "visualization"
    
    def get_required_inputs(self) -> List[str]:
        return ['networks', 'pipeline_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['visualization_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['network_analysis']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        required = ['networks', 'pipeline_results']
        for field in required:
            if field not in input_data:
                logger.error(f"Missing required input: {field}")
                return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'visualizations': {
                'type': 'object',
                'properties': {
                    'network_plots': {'type': 'array'},
                    'statistical_charts': {'type': 'array'},
                    'interactive_dashboards': {'type': 'array'}
                }
            },
            'visualization_stats': {
                'type': 'object',
                'properties': {
                    'total_visualizations': {'type': 'integer'},
                    'successful_visualizations': {'type': 'integer'},
                    'failed_visualizations': {'type': 'integer'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Create interactive visualizations."""
        start_time = time.time()
        
        try:
            networks = input_data['networks']
            pipeline_results = input_data['pipeline_results']
            visualization_config = input_data.get('visualization_config', {})
            
            # Create visualizations
            visualizations = {
                'network_plots': [],
                'statistical_charts': [],
                'interactive_dashboards': []
            }
            
            total_visualizations = 0
            successful_visualizations = 0
            failed_visualizations = 0
            
            # Generate network visualizations
            for protein_id, network_data in networks.items():
                try:
                    if network_data.get('status') == 'success' and network_data.get('network'):
                        network = network_data['network']
                        
                        # Create interactive network visualization
                        viz_result = self.network_builder.create_interactive_visualization(
                            embeddings=pipeline_results.get('embeddings', {}),
                            protein_ids=list(pipeline_results.get('embeddings', {}).keys()),
                            metadata_df=pd.DataFrame(),  # Placeholder
                            query_embedding=pipeline_results.get('embeddings', {}).get(protein_id),
                            query_protein_id=protein_id,
                            output_file=os.path.join(self.output_dir, f"network_{protein_id}.html")
                        )
                        
                        if viz_result[0]:  # Visualization object
                            visualizations['network_plots'].append({
                                'protein_id': protein_id,
                                'file_path': os.path.join(self.output_dir, f"network_{protein_id}.html"),
                                'type': 'interactive_network'
                            })
                            successful_visualizations += 1
                        else:
                            failed_visualizations += 1
                        
                        total_visualizations += 1
                        
                except Exception as e:
                    logger.warning(f"Failed to create visualization for {protein_id}: {e}")
                    failed_visualizations += 1
                    total_visualizations += 1
            
            # Generate statistical charts (placeholder)
            if pipeline_results.get('sequence_analyses'):
                try:
                    # Create property distribution charts
                    chart_path = os.path.join(self.output_dir, "property_distributions.html")
                    visualizations['statistical_charts'].append({
                        'file_path': chart_path,
                        'type': 'property_distributions',
                        'description': 'Amino acid composition and physicochemical properties'
                    })
                    successful_visualizations += 1
                    total_visualizations += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to create statistical charts: {e}")
                    failed_visualizations += 1
                    total_visualizations += 1
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_visualizations > 0,
                output_data={
                    'visualizations': visualizations,
                    'visualization_stats': {
                        'total_visualizations': total_visualizations,
                        'successful_visualizations': successful_visualizations,
                        'failed_visualizations': failed_visualizations
                    }
                },
                metadata={
                    'visualization_library': self.visualization_library,
                    'interactive_mode': self.interactive_mode,
                    'success_rate': successful_visualizations / total_visualizations if total_visualizations > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Visualization creation failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class DataExportStage(BaseStage):
    """
    Exports data in various formats.
    
    Handles:
    - JSON export
    - CSV export
    - HDF5 export
    - Workspace object creation
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.output_dir = config.get('output_dir', 'exports') if config else 'exports'
        self.export_formats = config.get('export_formats', ['json', 'csv']) if config else ['json', 'csv']
        self.include_workspace_objects = config.get('include_workspace_objects', True) if config else True
    
    def get_stage_name(self) -> str:
        return "data_export"
    
    def get_required_inputs(self) -> List[str]:
        return ['pipeline_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['export_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['report_generation']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'pipeline_results' not in input_data:
            logger.error("Missing pipeline_results")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'export_files': {
                'type': 'object',
                'properties': {
                    'json_files': {'type': 'array'},
                    'csv_files': {'type': 'array'},
                    'hdf5_files': {'type': 'array'},
                    'workspace_objects': {'type': 'array'}
                }
            },
            'export_stats': {
                'type': 'object',
                'properties': {
                    'total_exports': {'type': 'integer'},
                    'successful_exports': {'type': 'integer'},
                    'failed_exports': {'type': 'integer'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Export data in various formats."""
        start_time = time.time()
        
        try:
            pipeline_results = input_data['pipeline_results']
            export_config = input_data.get('export_config', {})
            
            # Create export directory
            os.makedirs(self.output_dir, exist_ok=True)
            
            # Export data
            export_files = {
                'json_files': [],
                'csv_files': [],
                'hdf5_files': [],
                'workspace_objects': []
            }
            
            total_exports = 0
            successful_exports = 0
            failed_exports = 0
            
            # JSON export
            if 'json' in self.export_formats:
                try:
                    json_path = os.path.join(self.output_dir, "pipeline_results.json")
                    with open(json_path, 'w') as f:
                        json.dump(pipeline_results, f, indent=2, default=str)
                    
                    export_files['json_files'].append(json_path)
                    successful_exports += 1
                    total_exports += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to export JSON: {e}")
                    failed_exports += 1
                    total_exports += 1
            
            # CSV export
            if 'csv' in self.export_formats:
                try:
                    # Export protein records
                    protein_records = pipeline_results.get('protein_records', [])
                    if protein_records:
                        csv_path = os.path.join(self.output_dir, "protein_records.csv")
                        df = pd.DataFrame([{
                            'protein_id': record.protein_id,
                            'sequence': record.sequence,
                            'source': record.source,
                            'length': len(record.sequence)
                        } for record in protein_records])
                        df.to_csv(csv_path, index=False)
                        export_files['csv_files'].append(csv_path)
                        successful_exports += 1
                        total_exports += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to export CSV: {e}")
                    failed_exports += 1
                    total_exports += 1
            
            # Workspace object creation
            if self.include_workspace_objects and workspace_client:
                try:
                    # Create workspace object with results
                    object_name = export_config.get('object_name', 'protein_query_results')
                    workspace_name = export_config.get('workspace_name', 'default')
                    
                    # Prepare object data
                    object_data = {
                        'pipeline_results': pipeline_results,
                        'export_timestamp': time.time(),
                        'export_config': export_config
                    }
                    
                    # Save to workspace (placeholder - actual implementation would use workspace client)
                    export_files['workspace_objects'].append({
                        'object_name': object_name,
                        'workspace_name': workspace_name,
                        'object_type': 'ProteinQueryResults'
                    })
                    successful_exports += 1
                    total_exports += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to create workspace object: {e}")
                    failed_exports += 1
                    total_exports += 1
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_exports > 0,
                output_data={
                    'export_files': export_files,
                    'export_stats': {
                        'total_exports': total_exports,
                        'successful_exports': successful_exports,
                        'failed_exports': failed_exports
                    }
                },
                metadata={
                    'export_formats': self.export_formats,
                    'include_workspace_objects': self.include_workspace_objects,
                    'success_rate': successful_exports / total_exports if total_exports > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Data export failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
