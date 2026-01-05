"""
Output Manager for KBase Protein Query Module

Handles output directory structure and provides simple file writing utilities.
"""

import os
import json
import time
import logging
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)

class OutputManager:
    """
    Manages output directory structure for the protein query module.
    
    This class handles:
    - Creating standardized output directory structures
    - Saving analysis results, metadata, and process logs
    - Writing JSON and text files safely
    """
    
    def __init__(
        self, 
        base_output_dir: str, 
        run_id: str, 
        workspace_name: Optional[str] = None, 
        kb_util: Optional[Any] = None
    ):
        """
        Initialize the Output Manager.
        
        Args:
            base_output_dir: Base directory for outputs
            run_id: Unique identifier for the run
            workspace_name: Optional KBase workspace name
            kb_util: Optional KBase utility instance
        """
        logger.debug(f"Initializing OutputManager with base_dir={base_output_dir}, run_id={run_id}")
        self.base_output_dir = base_output_dir
        self.run_id = run_id
        self.workspace_name = workspace_name
        self.kb_util = kb_util
        self.timestamp = time.strftime("%Y%m%d_%H%M%S")
        
        # Create main output directory structure
        self.root_dir = os.path.join(self.base_output_dir, "outputs")
        logger.debug(f"Creating root directory: {self.root_dir}")
        os.makedirs(self.root_dir, exist_ok=True)
        
        # Create subdirectories
        self.metadata_dir = os.path.join(self.root_dir, "metadata")
        self.analysis_dir = os.path.join(self.root_dir, "analysis")
        logger.debug(f"Creating subdirectories: metadata={self.metadata_dir}, analysis={self.analysis_dir}")
        os.makedirs(self.metadata_dir, exist_ok=True)
        os.makedirs(self.analysis_dir, exist_ok=True)
        
        logger.info(f"OutputManager initialized with root directory: {self.root_dir}")

    def get_root_dir(self) -> str:
        """Get the root output directory."""
        return self.root_dir
    
    def get_analysis_dir(self, analysis_name: str) -> str:
        """
        Get or create directory for a specific analysis.
        
        Args:
            analysis_name: Name of the analysis module
            
        Returns:
            Path to the analysis directory
        """
        logger.debug(f"Getting analysis directory for: {analysis_name}")
        analysis_path = os.path.join(self.analysis_dir, analysis_name)
        os.makedirs(analysis_path, exist_ok=True)
        logger.debug(f"Analysis directory: {analysis_path}")
        return analysis_path
    
    def write_json(self, rel_path: str, filename: str, data: Dict[str, Any]) -> str:
        """
        Write JSON data to a file.
        
        Args:
            rel_path: Relative path within root directory (can be empty)
            filename: Name of the file
            data: Dictionary data to write
            
        Returns:
            Absolute path to the written file
        """
        logger.debug(f"Writing JSON file: {filename} in {rel_path}")
        dir_path = os.path.join(self.root_dir, rel_path) if rel_path else self.root_dir
        os.makedirs(dir_path, exist_ok=True)
        path = os.path.join(dir_path, filename)
        
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, default=str)
        
        logger.debug(f"JSON file written: {path}")
        return path
    
    def write_text(self, rel_path: str, filename: str, text: str) -> str:
        """
        Write text data to a file.
        
        Args:
            rel_path: Relative path within root directory
            filename: Name of the file
            text: Text content to write
            
        Returns:
            Absolute path to the written file
        """
        dir_path = os.path.join(self.root_dir, rel_path) if rel_path else self.root_dir
        os.makedirs(dir_path, exist_ok=True)
        path = os.path.join(dir_path, filename)
        
        with open(path, 'w', encoding='utf-8') as f:
            f.write(text)
        
        return path
    
    def save_analysis_output(
        self, 
        analysis_name: str, 
        result: Dict[str, Any], 
        output_dir: str
    ) -> Dict[str, Any]:
        """
        Save output from a specific analysis.
        
        Args:
            analysis_name: Name of the analysis
            result: Analysis results to save
            output_dir: Directory to save outputs to (unused, kept for compatibility)
            
        Returns:
            Dictionary of saved file paths
        """
        logger.debug(f"Saving analysis output for: {analysis_name}")
        analysis_dir = self.get_analysis_dir(analysis_name)
        saved_files = {}
        
        try:
            # Save main results as JSON
            logger.debug(f"Writing results.json for {analysis_name}")
            results_file = self.write_json(
                f"analysis/{analysis_name}",
                "results.json",
                result
            )
            saved_files["results"] = results_file
            
            # Save summary if available
            if "summary" in result:
                logger.debug(f"Writing summary.txt for {analysis_name}")
                summary_file = self.write_text(
                    f"analysis/{analysis_name}",
                    "summary.txt",
                    str(result["summary"])
                )
                saved_files["summary"] = summary_file
            
            logger.info(f"Saved output for analysis {analysis_name} to {analysis_dir}")
            logger.debug(f"Saved files: {list(saved_files.keys())}")
            
        except Exception as e:
            logger.error(f"Error saving output for analysis {analysis_name}: {e}", exc_info=True)
            raise
        
        return saved_files
    
    def save_metadata(
        self, 
        config: Dict[str, Any], 
        analyses_run: List[str], 
        summary: str = "", 
        process_log: Optional[List[Dict[str, Any]]] = None
    ) -> str:
        """
        Save metadata about the run.
        
        Args:
            config: Configuration dictionary used for the run
            analyses_run: List of analyses executed
            summary: Summary string of the run
            process_log: Optional log of process steps
            
        Returns:
            Path to the saved metadata file
        """
        metadata = {
            "run_id": self.run_id,
            "workspace_name": self.workspace_name,
            "timestamp": self.timestamp,
            "summary": summary,
            "analyses_run": analyses_run,
            "config": config
        }
        
        if process_log is not None:
            metadata["process_log"] = process_log
        
        return self.write_json("metadata", "run_metadata.json", metadata)
    
    def save_process_info(
        self, 
        stages_completed: List[str], 
        execution_time: float
    ) -> str:
        """
        Save process execution information.
        
        Args:
            stages_completed: List of completed stages
            execution_time: Total execution time in seconds
            
        Returns:
            Path to the saved process info file
        """
        process_info = {
            "run_id": self.run_id,
            "timestamp": self.timestamp,
            "stages_completed": stages_completed,
            "execution_time_seconds": execution_time
        }
        
        return self.write_json("metadata", "process_info.json", process_info)


def main() -> int:
    """
    Test OutputManager.
    - Creates a temp output directory
    - Writes metadata, process info, and an example analysis result
    - Verifies files exist on disk

    Dependencies:
    - tempfile
    - shutil
    - os
    - json
    - time
    - logging
    """
    import tempfile
    import shutil

    ok = True
    tmpdir = tempfile.mkdtemp(prefix="kpqm_output_test_")
    try:
        mgr = OutputManager(base_output_dir=tmpdir, run_id="test_run")

        # Write metadata and process info
        meta_path = mgr.save_metadata(config={"example": True}, analyses_run=["network_analysis"], summary="Test run")
        proc_path = mgr.save_process_info(stages_completed=["network_analysis"], execution_time=0.123)

        # Write example analysis outputs
        analysis_dir = mgr.get_analysis_dir("network_analysis")
        saved = mgr.save_analysis_output("network_analysis", {"summary": "ok", "success": True}, output_dir=analysis_dir)

        # Verify files
        required = [meta_path, proc_path] + list(saved.values())
        missing = [p for p in required if not (isinstance(p, str) and os.path.exists(p))]
        if missing:
            raise RuntimeError(f"Missing expected files: {missing}")

        print(f"OutputManager test: SUCCESS -> root={mgr.get_root_dir()}")
    except Exception as e:
        ok = False
        print(f"OutputManager test: FAILED - {e}")
    finally:
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
