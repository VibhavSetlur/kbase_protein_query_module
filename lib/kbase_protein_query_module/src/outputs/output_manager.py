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
                 workspace_name: Optional[str] = None):
        """
        Initialize the Output Manager.
        
        Args:
            base_output_dir: Base directory for all outputs
            run_id: Unique identifier for this run
            workspace_name: KBase workspace name if applicable
        """
        self.base_output_dir = base_output_dir
        self.run_id = run_id
        self.workspace_name = workspace_name
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
        self.analysis_outputs: Dict[str, Any] = {}
        
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