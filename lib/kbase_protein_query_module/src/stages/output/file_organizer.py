"""
File Organizer for Output Stage

Comprehensive utility to organize, format, and structure all output data files.
Handles JSON metadata, TSV data files, HTML visualizations, and complete output directory structure.
"""

import os
import json
import time
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Union
from pathlib import Path

logger = logging.getLogger(__name__)


class FileOrganizer:
    """
    Comprehensive file organizer for protein analysis outputs.
    Creates structured output directories with proper data formatting.
    
    Output Directory Structure:
    output_dir/
    ├── raw_data/
    │   ├── network/
    │   ├── sequence/
    │   └── embeddings/
    ├── visualizations/
    │   ├── network/
    │   └── other/
    ├── metadata/
    └── pipeline_summary.json
    """
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self._create_output_structure()
        
        self.files_manifest = {
            'created_at': time.strftime('%Y-%m-%d %H:%M:%S'),
            'output_directory': output_dir,
            'data_files': {},
            'visualizations': {},
            'metadata_files': {},
            'directory_structure': self._get_directory_structure()
        }
    
    def _create_output_structure(self):
        """Create the structured output directory."""
        subdirs = [
            'raw_data/network',
            'raw_data/sequence', 
            'raw_data/embeddings',
            'visualizations/network',
            'visualizations/other',
            'metadata'
        ]
        
        for subdir in subdirs:
            os.makedirs(os.path.join(self.output_dir, subdir), exist_ok=True)
        
        logger.info(f"Created output directory structure in: {self.output_dir}")
    
    def _get_directory_structure(self) -> Dict[str, str]:
        """Get the directory structure mapping."""
        return {
            'raw_data_network': os.path.join(self.output_dir, 'raw_data', 'network'),
            'raw_data_sequence': os.path.join(self.output_dir, 'raw_data', 'sequence'),
            'raw_data_embeddings': os.path.join(self.output_dir, 'raw_data', 'embeddings'),
            'visualizations_network': os.path.join(self.output_dir, 'visualizations', 'network'),
            'visualizations_other': os.path.join(self.output_dir, 'visualizations', 'other'),
            'metadata': os.path.join(self.output_dir, 'metadata')
        }
    
    def register_data_file(self, file_type: str, file_path: str, description: str = None, 
                          analysis_type: str = None):
        """Register a data file in the appropriate structured directory."""
        filename = os.path.basename(file_path)
        
        # Determine target directory based on analysis type and file type
        target_dir = self._get_target_directory(file_type, analysis_type)
        
        # Move file to structured location if needed
        if not file_path.startswith(target_dir):
            new_path = os.path.join(target_dir, filename)
            if os.path.exists(file_path) and file_path != new_path:
                os.makedirs(target_dir, exist_ok=True)
                import shutil
                shutil.copy2(file_path, new_path)
                file_path = new_path
        
        self.files_manifest['data_files'][file_type] = {
            'filename': filename,
            'path': file_path,
            'relative_path': os.path.relpath(file_path, self.output_dir),
            'description': description or f"{file_type} data file",
            'analysis_type': analysis_type,
            'size_bytes': os.path.getsize(file_path) if os.path.exists(file_path) else 0
        }
        logger.info(f"Registered {file_type} data file: {filename}")
        return file_path
    
    def _get_target_directory(self, file_type: str, analysis_type: str = None) -> str:
        """Determine the target directory for a file based on type and analysis."""
        if 'network' in file_type.lower() or analysis_type == 'network':
            if file_type.endswith(('.html', '.htm')):
                return self.files_manifest['directory_structure']['visualizations_network']
            else:
                return self.files_manifest['directory_structure']['raw_data_network']
        elif 'sequence' in file_type.lower() or analysis_type == 'sequence':
            return self.files_manifest['directory_structure']['raw_data_sequence']
        elif 'embedding' in file_type.lower() or analysis_type == 'embeddings':
            return self.files_manifest['directory_structure']['raw_data_embeddings']
        elif file_type.endswith(('.html', '.htm')):
            return self.files_manifest['directory_structure']['visualizations_other']
        else:
            return self.files_manifest['directory_structure']['metadata']
    
    def register_visualization(self, viz_type: str, file_path: str, format_type: str = "html"):
        """Register a visualization file in appropriate directory."""
        filename = os.path.basename(file_path)
        
        # Move to appropriate visualization directory
        target_dir = self._get_target_directory(file_path, 'visualization')
        if not file_path.startswith(target_dir):
            new_path = os.path.join(target_dir, filename)
            if os.path.exists(file_path) and file_path != new_path:
                import shutil
                shutil.copy2(file_path, new_path)
                file_path = new_path
        
        self.files_manifest['visualizations'][viz_type] = {
            'filename': filename,
            'path': file_path,
            'relative_path': os.path.relpath(file_path, self.output_dir),
            'format': format_type,
            'size_bytes': os.path.getsize(file_path) if os.path.exists(file_path) else 0
        }
        logger.info(f"Registered {viz_type} visualization: {filename}")
        return file_path
    
    def register_metadata(self, metadata_type: str, file_path: str):
        """Register a metadata file."""
        filename = os.path.basename(file_path)
        self.files_manifest['metadata_files'][metadata_type] = {
            'filename': filename,
            'path': file_path,
            'size_bytes': os.path.getsize(file_path) if os.path.exists(file_path) else 0
        }
        logger.info(f"Registered {metadata_type} metadata: {filename}")
    
    def save_manifest(self) -> str:
        """Save the files manifest as JSON."""
        manifest_path = os.path.join(self.output_dir, 'files_manifest.json')
        with open(manifest_path, 'w') as f:
            json.dump(self.files_manifest, f, indent=2)
        
        logger.info(f"Files manifest saved to: {manifest_path}")
        return manifest_path
    
    # ==========================================
    # DATA FORMATTING METHODS (from data_formatter.py)
    # ==========================================
    
    def format_protein_data(self, protein_data: Dict[str, Any]) -> Dict[str, Any]:
        """Format protein data for consistent output."""
        formatted = {}
        for key, value in protein_data.items():
            if isinstance(value, np.ndarray):
                formatted[key] = value.tolist()
            elif isinstance(value, float):
                formatted[key] = round(value, 3)
            elif isinstance(value, np.float64):
                formatted[key] = round(float(value), 3)
            else:
                formatted[key] = value
        return formatted
    
    def format_similarity_scores(self, similarities: List[float]) -> List[float]:
        """Format similarity scores for display."""
        return [round(float(score), 3) for score in similarities]
    
    def create_comprehensive_pipeline_summary(self, pipeline_results: Dict[str, Any], 
                                            pipeline_metadata: Dict[str, Any] = None) -> str:
        """Create a comprehensive JSON summary of the entire pipeline."""
        summary = {
            'pipeline_info': {
                'generated_at': time.strftime('%Y-%m-%d %H:%M:%S'),
                'version': '2.0.0',
                'total_stages': len(pipeline_results),
                'output_directory': self.output_dir
            },
            'pipeline_metadata': pipeline_metadata or {},
            'analysis_results': {},
            'output_files': {
                'data_files': self.files_manifest['data_files'],
                'visualizations': self.files_manifest['visualizations'],
                'metadata_files': self.files_manifest['metadata_files']
            },
            'directory_structure': self.files_manifest['directory_structure']
        }
        
        # Process each analysis stage
        for stage_name, stage_data in pipeline_results.items():
            summary['analysis_results'][stage_name] = self._summarize_stage_results(stage_name, stage_data)
        
        # Save comprehensive summary
        summary_path = os.path.join(self.output_dir, 'pipeline_summary.json')
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Comprehensive pipeline summary saved to: {summary_path}")
        return summary_path
    
    def _summarize_stage_results(self, stage_name: str, stage_data: Any) -> Dict[str, Any]:
        """Summarize results from a specific analysis stage."""
        if not isinstance(stage_data, dict):
            return {'status': 'completed', 'data_type': str(type(stage_data))}
        
        summary = {
            'status': 'completed',
            'data_keys': list(stage_data.keys()) if isinstance(stage_data, dict) else []
        }
        
        # Stage-specific summaries
        if stage_name == 'network_analysis':
            network_props = stage_data.get('network_properties', {})
            summary.update({
                'network_nodes': network_props.get('num_nodes', 0),
                'network_edges': network_props.get('num_edges', 0),
                'network_density': network_props.get('density', 0),
                'visualization_type': stage_data.get('visualization_type', 'unknown'),
                'csv_files_count': len(stage_data.get('csv_files', {}))
            })
        
        elif stage_name == 'multi_protein_network_analysis':
            individual = stage_data.get('individual_networks', {})
            summary.update({
                'total_query_proteins': len(individual),
                'query_proteins': list(individual.keys()),
                'csv_files_count': len(stage_data.get('csv_data', {}))
            })
        
        elif stage_name == 'sequence_analysis':
            summary.update({
                'total_proteins_analyzed': stage_data.get('total_proteins', 0),
                'successful_analyses': stage_data.get('successful', 0),
                'failed_analyses': stage_data.get('failed', 0),
                'csv_files_count': len(stage_data.get('csv_files', {}))
            })
        
        elif stage_name == 'embedding_generation':
            summary.update({
                'embeddings_shape': stage_data.get('embeddings_shape', []),
                'num_proteins': stage_data.get('num_proteins', 0),
                'h5_file': stage_data.get('h5_path', ''),
                'embedding_model': stage_data.get('model_name', 'unknown')
            })
        
        return summary
    
    def get_file_summary(self) -> Dict[str, Any]:
        """Get a summary of organized files."""
        return {
            'total_data_files': len(self.files_manifest['data_files']),
            'total_visualizations': len(self.files_manifest['visualizations']),
            'total_metadata_files': len(self.files_manifest['metadata_files']),
            'output_directory': self.output_dir,
            'manifest': self.files_manifest
        }


def create_simple_file_list(output_dir: str, files: Dict[str, str]) -> str:
    """Create a simple text file listing all outputs."""
    file_list_path = os.path.join(output_dir, 'output_files.txt')
    
    with open(file_list_path, 'w') as f:
        f.write(f"Protein Query Analysis Output Files\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Output Directory: {output_dir}\n\n")
        
        for file_type, file_path in files.items():
            filename = os.path.basename(file_path)
            f.write(f"{file_type}: {filename}\n")
    
    return file_list_path
