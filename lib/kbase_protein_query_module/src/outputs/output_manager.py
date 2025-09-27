"""
Output Manager for KBase Protein Query Module

This module manages all output generation and organization for the protein query module.
It coordinates with analysis-specific output handlers and creates organized output directories.
"""

import os
import json
import time
import hashlib
import logging
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

@dataclass
class ArtifactRecord:
    """Record for tracking output artifacts."""
    path: str
    kind: str  # 'json', 'csv', 'html', 'image', 'data', 'log'
    description: str = ""
    analysis_type: Optional[str] = None
    size_bytes: int = 0
    checksum_md5: Optional[str] = None

class OutputManager:
    """
    Manages output generation and organization for the protein query module.
    
    This class coordinates with analysis-specific output handlers and creates
    organized output directories with proper metadata and provenance tracking.
    """
    
    def __init__(self, base_output_dir: str, run_id: str, 
                 workspace_name: Optional[str] = None, kb_util=None):
        """
        Initialize the Output Manager.
        
        Args:
            base_output_dir: Base directory for all outputs
            run_id: Unique identifier for this run
            workspace_name: KBase workspace name if applicable
            kb_util: KBase utility library instance for workspace operations
        """
        self.base_output_dir = base_output_dir
        self.run_id = run_id
        self.workspace_name = workspace_name
        self.kb_util = kb_util
        self.timestamp = time.strftime("%Y%m%d_%H%M%S")
        
        # Create main output directory structure
        self.root_dir = os.path.join(
            self.base_output_dir,
            "outputs",
            f"{self.run_id}_{self.timestamp}"
        )
        os.makedirs(self.root_dir, exist_ok=True)
        
        # Initialize tracking
        self.artifacts: List[ArtifactRecord] = []
        self.analysis_outputs: Dict[str, Any] = []
        self.workspace_objects: List[Dict[str, str]] = []  # Track created workspace objects
        
        # Create subdirectories
        self._create_directory_structure()
        
        logger.info(f"OutputManager initialized with root directory: {self.root_dir}")
    
    def _create_directory_structure(self):
        """Create the standard output directory structure."""
        directories = [
            "metadata",
            "process_info", 
            "analysis",
            "logs"
        ]
        
        for directory in directories:
            os.makedirs(os.path.join(self.root_dir, directory), exist_ok=True)
    
    def get_root_dir(self) -> str:
        """Get the root output directory."""
        return self.root_dir
    
    def get_analysis_dir(self, analysis_name: str) -> str:
        """
        Get or create directory for a specific analysis.
        
        Args:
            analysis_name: Name of the analysis
            
        Returns:
            Path to the analysis directory
        """
        analysis_dir = os.path.join(self.root_dir, "analysis", analysis_name)
        os.makedirs(analysis_dir, exist_ok=True)
        return analysis_dir
    
    def _compute_md5(self, filepath: str) -> str:
        """Compute MD5 checksum for a file."""
        md5 = hashlib.md5()
        try:
            with open(filepath, 'rb') as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b''):
                    md5.update(chunk)
        except Exception as e:
            logger.warning(f"Could not compute MD5 for {filepath}: {e}")
        return md5.hexdigest()
    
    def record_artifact(self, path: str, kind: str, description: str = "", 
                       analysis_type: Optional[str] = None) -> str:
        """
        Record an output artifact.
        
        Args:
            path: Path to the artifact file
            kind: Type of artifact (json, csv, html, etc.)
            description: Description of the artifact
            analysis_type: Type of analysis that generated this artifact
            
        Returns:
            Path to the recorded artifact
        """
        try:
            size = os.path.getsize(path) if os.path.exists(path) else 0
            checksum = self._compute_md5(path) if os.path.exists(path) else None
        except Exception as e:
            logger.warning(f"Could not get file info for {path}: {e}")
            size, checksum = 0, None
            
        record = ArtifactRecord(
            path=path,
            kind=kind,
            description=description,
            analysis_type=analysis_type,
            size_bytes=size,
            checksum_md5=checksum
        )
        self.artifacts.append(record)
        return path
    
    def write_json(self, rel_dir: str, filename: str, data: Dict[str, Any], 
                   analysis_type: Optional[str] = None, description: str = "") -> str:
        """
        Write JSON data to a file.
        
        Args:
            rel_dir: Relative directory within output structure
            filename: Name of the file
            data: Data to write as JSON
            analysis_type: Type of analysis generating this output
            description: Description of the file
            
        Returns:
            Path to the written file
        """
        dir_path = os.path.join(self.root_dir, rel_dir)
        os.makedirs(dir_path, exist_ok=True)
        path = os.path.join(dir_path, filename)
        
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, default=str)
        
        return self.record_artifact(
            path, 'json', 
            description=description or filename, 
            analysis_type=analysis_type
        )
    
    def write_text(self, rel_dir: str, filename: str, text: str,
                   analysis_type: Optional[str] = None, description: str = "") -> str:
        """
        Write text data to a file.
        
        Args:
            rel_dir: Relative directory within output structure
            filename: Name of the file
            text: Text content to write
            analysis_type: Type of analysis generating this output
            description: Description of the file
            
        Returns:
            Path to the written file
        """
        dir_path = os.path.join(self.root_dir, rel_dir)
        os.makedirs(dir_path, exist_ok=True)
        path = os.path.join(dir_path, filename)
        
        with open(path, 'w', encoding='utf-8') as f:
            f.write(text)
        
        return self.record_artifact(
            path, 'data',
            description=description or filename,
            analysis_type=analysis_type
        )
    
    def write_csv(self, rel_dir: str, filename: str, csv_text: str,
                  analysis_type: Optional[str] = None, description: str = "") -> str:
        """
        Write CSV data to a file.
        
        Args:
            rel_dir: Relative directory within output structure
            filename: Name of the file
            csv_text: CSV content to write
            analysis_type: Type of analysis generating this output
            description: Description of the file
            
        Returns:
            Path to the written file
        """
        dir_path = os.path.join(self.root_dir, rel_dir)
        os.makedirs(dir_path, exist_ok=True)
        path = os.path.join(dir_path, filename)
        
        with open(path, 'w', encoding='utf-8') as f:
            f.write(csv_text)
        
        return self.record_artifact(
            path, 'csv',
            description=description or filename,
            analysis_type=analysis_type
        )
    
    def save_analysis_output(self, analysis_name: str, result: Dict[str, Any], 
                           output_dir: str) -> Dict[str, Any]:
        """
        Save output from a specific analysis.
        
        Args:
            analysis_name: Name of the analysis
            result: Analysis results to save
            output_dir: Directory to save outputs to
            
        Returns:
            Dictionary of saved file paths
        """
        analysis_dir = self.get_analysis_dir(analysis_name)
        saved_files = {}
        
        try:
            # Save main results as JSON
            results_file = self.write_json(
                f"analysis/{analysis_name}",
                "results.json",
                result,
                analysis_type=analysis_name,
                description=f"Main results from {analysis_name}"
            )
            saved_files["results"] = results_file
            
            # Save summary if available
            if "summary" in result:
                summary_file = self.write_text(
                    f"analysis/{analysis_name}",
                    "summary.txt",
                    str(result["summary"]),
                    analysis_type=analysis_name,
                    description=f"Summary from {analysis_name}"
                )
                saved_files["summary"] = summary_file
            
            # Save any additional data files
            if "data_files" in result:
                for file_info in result["data_files"]:
                    if "content" in file_info and "filename" in file_info:
                        file_path = self.write_text(
                            f"analysis/{analysis_name}",
                            file_info["filename"],
                            file_info["content"],
                            analysis_type=analysis_name,
                            description=file_info.get("description", f"Data file from {analysis_name}")
                        )
                        saved_files[file_info["filename"]] = file_path
            
            # Store analysis output info
            self.analysis_outputs[analysis_name] = {
                "analysis_name": analysis_name,
                "output_dir": analysis_dir,
                "saved_files": saved_files,
                "timestamp": time.time()
            }
            
            logger.info(f"Saved output for analysis {analysis_name} to {analysis_dir}")
            
        except Exception as e:
            logger.error(f"Error saving output for analysis {analysis_name}: {e}")
            raise
        
        return saved_files
    
    def save_metadata(self, config: Dict[str, Any], analyses_run: List[str], 
                     summary: str = "") -> str:
        """
        Save metadata about the run.
        
        Args:
            config: Configuration used for the run
            analyses_run: List of analyses that were run
            summary: Summary of the run
            
        Returns:
            Path to the metadata file
        """
        metadata = {
            "run_id": self.run_id,
            "workspace_name": self.workspace_name,
            "timestamp": self.timestamp,
            "summary": summary,
            "analyses_run": analyses_run,
            "config": config,
            "total_artifacts": len(self.artifacts),
            "analysis_outputs": self.analysis_outputs
        }
        
        return self.write_json(
            "metadata",
            "run_metadata.json",
            metadata,
            description="Run metadata and configuration"
        )
    
    def save_process_info(self, stages_completed: List[str], 
                         execution_time: float, 
                         memory_usage: Optional[Dict[str, Any]] = None) -> str:
        """
        Save process execution information.
        
        Args:
            stages_completed: List of completed stages
            execution_time: Total execution time
            memory_usage: Memory usage information if available
            
        Returns:
            Path to the process info file
        """
        process_info = {
            "run_id": self.run_id,
            "timestamp": self.timestamp,
            "stages_completed": stages_completed,
            "execution_time_seconds": execution_time,
            "memory_usage": memory_usage or {},
            "artifacts_generated": len(self.artifacts)
        }
        
        return self.write_json(
            "process_info",
            "process_info.json",
            process_info,
            description="Process execution information"
        )
    
    def finalize_manifest(self) -> str:
        """
        Create final manifest of all outputs.
        
        Returns:
            Path to the manifest file
        """
        manifest = {
            "run_id": self.run_id,
            "workspace_name": self.workspace_name,
            "timestamp": self.timestamp,
            "root_directory": self.root_dir,
            "total_artifacts": len(self.artifacts),
            "artifacts": [
                {
                    "path": artifact.path,
                    "kind": artifact.kind,
                    "description": artifact.description,
                    "analysis_type": artifact.analysis_type,
                    "size_bytes": artifact.size_bytes,
                    "checksum_md5": artifact.checksum_md5
                }
                for artifact in self.artifacts
            ],
            "analysis_outputs": self.analysis_outputs
        }
        
        return self.write_json(
            "",
            "manifest.json",
            manifest,
            description="Complete manifest of all outputs"
        )
    
    def get_artifact_summary(self) -> Dict[str, Any]:
        """
        Get a summary of all artifacts.
        
        Returns:
            Dictionary containing artifact summary
        """
        summary = {
            "total_artifacts": len(self.artifacts),
            "by_kind": {},
            "by_analysis": {},
            "total_size_bytes": 0
        }
        
        for artifact in self.artifacts:
            # Count by kind
            kind = artifact.kind
            summary["by_kind"][kind] = summary["by_kind"].get(kind, 0) + 1
            
            # Count by analysis
            analysis = artifact.analysis_type or "unknown"
            summary["by_analysis"][analysis] = summary["by_analysis"].get(analysis, 0) + 1
            
            # Sum total size
            summary["total_size_bytes"] += artifact.size_bytes
        
        return summary
    
    def create_workspace_objects(self, analysis_results: Dict[str, Any]) -> List[Dict[str, str]]:
        """
        Create KBase workspace objects from analysis results.
        
        Args:
            analysis_results: Results from all analyses
            
        Returns:
            List of created workspace object references
        """
        created_objects = []
        
        if not self.kb_util or not self.workspace_name:
            logger.warning("No KBUtilLib or workspace name available for workspace object creation")
            return created_objects
        
        try:
            # Create main analysis results object
            main_results_data = {
                'run_id': self.run_id,
                'timestamp': self.timestamp,
                'workspace_name': self.workspace_name,
                'analyses_completed': list(analysis_results.keys()),
                'summary': self._generate_analysis_summary(analysis_results),
                'analysis_results': analysis_results,
                'total_artifacts': len(self.artifacts),
                'output_directory': self.root_dir
            }
            
            main_object_ref = self._save_workspace_object(
                f"{self.run_id}_protein_analysis_results",
                'KBaseProteinQueryModule.ProteinAnalysisResults',
                main_results_data,
                "Complete protein analysis results from all stages"
            )
            
            if main_object_ref:
                created_objects.append({
                    'ref': main_object_ref,
                    'name': f"{self.run_id}_protein_analysis_results",
                    'type': 'KBaseProteinQueryModule.ProteinAnalysisResults',
                    'description': 'Complete protein analysis results from all stages'
                })
            
            # Create individual analysis result objects
            for analysis_name, result_data in analysis_results.items():
                if result_data and 'error' not in result_data:
                    analysis_object_data = {
                        'analysis_name': analysis_name,
                        'run_id': self.run_id,
                        'timestamp': self.timestamp,
                        'results': result_data,
                        'metadata': {
                            'output_directory': self.root_dir,
                            'artifacts_count': len([a for a in self.artifacts if a.analysis_type == analysis_name])
                        }
                    }
                    
                    analysis_object_ref = self._save_workspace_object(
                        f"{self.run_id}_{analysis_name}_results",
                        'KBaseProteinQueryModule.AnalysisResult',
                        analysis_object_data,
                        f"Results from {analysis_name} analysis"
                    )
                    
                    if analysis_object_ref:
                        created_objects.append({
                            'ref': analysis_object_ref,
                            'name': f"{self.run_id}_{analysis_name}_results",
                            'type': 'KBaseProteinQueryModule.AnalysisResult',
                            'description': f"Results from {analysis_name} analysis"
                        })
            
            # Create data export objects for CSV and JSON files
            csv_files = [a for a in self.artifacts if a.kind == 'csv']
            json_files = [a for a in self.artifacts if a.kind == 'json']
            
            if csv_files:
                csv_data = {
                    'run_id': self.run_id,
                    'timestamp': self.timestamp,
                    'files': [{'path': a.path, 'description': a.description, 'size_bytes': a.size_bytes} for a in csv_files],
                    'total_files': len(csv_files)
                }
                
                csv_object_ref = self._save_workspace_object(
                    f"{self.run_id}_data_exports",
                    'KBaseProteinQueryModule.DataExports',
                    csv_data,
                    "CSV data exports from protein analysis"
                )
                
                if csv_object_ref:
                    created_objects.append({
                        'ref': csv_object_ref,
                        'name': f"{self.run_id}_data_exports",
                        'type': 'KBaseProteinQueryModule.DataExports',
                        'description': "CSV data exports from protein analysis"
                    })
            
            if json_files:
                json_data = {
                    'run_id': self.run_id,
                    'timestamp': self.timestamp,
                    'files': [{'path': a.path, 'description': a.description, 'size_bytes': a.size_bytes} for a in json_files],
                    'total_files': len(json_files)
                }
                
                json_object_ref = self._save_workspace_object(
                    f"{self.run_id}_analysis_metadata",
                    'KBaseProteinQueryModule.AnalysisMetadata',
                    json_data,
                    "JSON metadata and configuration files from protein analysis"
                )
                
                if json_object_ref:
                    created_objects.append({
                        'ref': json_object_ref,
                        'name': f"{self.run_id}_analysis_metadata",
                        'type': 'KBaseProteinQueryModule.AnalysisMetadata',
                        'description': "JSON metadata and configuration files from protein analysis"
                    })
            
            # Store created objects
            self.workspace_objects = created_objects
            
            logger.info(f"Created {len(created_objects)} workspace objects")
            return created_objects
            
        except Exception as e:
            logger.error(f"Failed to create workspace objects: {e}")
            return created_objects
    
    def _save_workspace_object(self, object_name: str, object_type: str, 
                              object_data: Dict[str, Any], description: str = "") -> Optional[str]:
        """
        Save a single workspace object using KBUtilLib.
        
        Args:
            object_name: Name for the workspace object
            object_type: KBase object type
            object_data: Data to save
            description: Description of the object
            
        Returns:
            Object reference if successful, None otherwise
        """
        try:
            if hasattr(self.kb_util, 'save_workspace_object'):
                object_ref = self.kb_util.save_workspace_object(
                    self.workspace_name,
                    object_name,
                    object_type,
                    object_data
                )
                logger.info(f"Saved workspace object: {object_name} ({object_type})")
                return object_ref
            else:
                logger.warning("KBUtilLib save_workspace_object method not available")
                return None
                
        except Exception as e:
            logger.error(f"Failed to save workspace object {object_name}: {e}")
            return None
    
    def _generate_analysis_summary(self, analysis_results: Dict[str, Any]) -> str:
        """
        Generate a summary of analysis results for workspace object.
        
        Args:
            analysis_results: Results from all analyses
            
        Returns:
            Summary string
        """
        summary_parts = [
            f"Protein Query Analysis Run {self.run_id}",
            f"Completed {len(analysis_results)} analyses:",
        ]
        
        for analysis_name, result in analysis_results.items():
            if isinstance(result, dict) and 'error' in result:
                summary_parts.append(f"  - {analysis_name}: FAILED ({result['error']})")
            else:
                summary_parts.append(f"  - {analysis_name}: SUCCESS")
        
        summary_parts.extend([
            f"Total artifacts generated: {len(self.artifacts)}",
            f"Output directory: {self.root_dir}",
            f"Generated: {self.timestamp}"
        ])
        
        return "\n".join(summary_parts)
    
    def get_workspace_objects_summary(self) -> Dict[str, Any]:
        """
        Get summary of all created workspace objects.
        
        Returns:
            Dictionary with workspace objects summary
        """
        return {
            'total_objects': len(self.workspace_objects),
            'objects': self.workspace_objects,
            'workspace_name': self.workspace_name,
            'run_id': self.run_id
        }