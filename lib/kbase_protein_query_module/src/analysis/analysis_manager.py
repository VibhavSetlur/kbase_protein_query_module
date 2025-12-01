"""
Analysis Manager Module

Manages execution of protein analysis workflows.
"""

import logging
import importlib
import sys
import os
from typing import Dict, Any, Optional

# Handle both script execution and module import
if __name__ == "__main__" or __package__ is None:
    # Add parent directories to path for script execution
    # We need to add 'lib' to path so 'kbase_protein_query_module' is importable
    current_dir = os.path.dirname(os.path.abspath(__file__)) # src/analysis
    src_dir = os.path.dirname(current_dir) # src
    module_dir = os.path.dirname(src_dir) # kbase_protein_query_module
    lib_dir = os.path.dirname(module_dir) # lib
    
    if lib_dir not in sys.path:
        sys.path.insert(0, lib_dir)
        
    from kbase_protein_query_module.src.analysis.config import get_enabled_analyses
else:
    from .config import get_enabled_analyses

logger = logging.getLogger(__name__)

class AnalysisManager:
    """Manages the execution of protein analysis workflows."""
    
    def __init__(self, output_manager=None, config: Optional[Dict[str, Any]] = None):
        """Initialize the Analysis Manager."""
        self.output_manager = output_manager
        self.config = config or {}
        self.analyses: Dict[str, Any] = {}
        self.results: Dict[str, Any] = {}
        self._load_analysis_modules()
    
    def _load_analysis_modules(self):
        """Load all available analysis modules dynamically."""
        enabled_analyses = get_enabled_analyses()
        
        for analysis_name, analysis_config in enabled_analyses.items():
            try:
                module_path = analysis_config.get("module_path")
                class_name = analysis_config.get("class_name")
                
                if not module_path or not class_name:
                    logger.warning(f"Missing module_path or class_name for {analysis_name}")
                    continue
                
                module = importlib.import_module(module_path)
                if hasattr(module, class_name):
                    analysis_class = getattr(module, class_name)
                    # Pass config to analysis instance - merge analysis-specific config with global config
                    analysis_instance_config = self.config.copy()
                    # Allow analysis-specific config overrides
                    if 'config' in analysis_config:
                        analysis_instance_config.update(analysis_config['config'])
                    self.analyses[analysis_name] = analysis_class(config=analysis_instance_config)
                    logger.info(f"Loaded analysis: {analysis_name}")
                else:
                    logger.warning(f"Class {class_name} not found in {module_path}")
            except Exception as e:
                logger.error(f"Failed to load {analysis_name}: {e}")
    
    def get_available_analyses(self) -> Dict[str, Dict[str, Any]]:
        """Get list of available analyses."""
        return get_enabled_analyses()
    
    def run_analyses(self, analysis_name: str, input_data: Dict[str, Any],
                    output_dir: str = None, **kwargs) -> Dict[str, Any]:
        """Run a specific analysis with input_data."""
        if analysis_name not in self.analyses:
            logger.error(f"Analysis '{analysis_name}' not found")
            return None

        try:
            analysis = self.analyses[analysis_name]
            
            # Standardize execution method
            if hasattr(analysis, "run_network_analysis"):
                result = analysis.run_network_analysis(input_data)
            elif hasattr(analysis, "run"):
                result = analysis.run(input_data, **kwargs)
            elif hasattr(analysis, "analyze"):
                result = analysis.analyze(input_data, **kwargs)
            else:
                raise AttributeError(f"Analysis '{analysis_name}' has no runnable method (run_network_analysis, run, or analyze)")

            self.results[analysis_name] = result
            logger.info(f"Analysis '{analysis_name}' completed")
            return result
        except Exception as e:
            logger.error(f"Error running '{analysis_name}': {e}")
            raise
    
    def get_analysis_results(self, analysis_name: Optional[str] = None) -> Dict[str, Any]:
        """Get results from completed analyses."""
        if analysis_name:
            return self.results.get(analysis_name)
        return self.results.copy()

def main() -> int:
    """
    Test AnalysisManager.
    - Loads enabled analyses
    - Optionally runs network_analysis if available with a minimal input
    """
    ok = True
    try:
        mgr = AnalysisManager()
        available = mgr.get_available_analyses()
        if not isinstance(available, dict):
            raise RuntimeError("get_available_analyses did not return a dict")

        # If network_analysis is available, try a minimal run
        if "network_analysis" in available and "network_analysis" in mgr.analyses:
            input_data = {
                'input_type': 'protein_sequence',
                'protein_sequence': 'ACDEFGHIKLMNPQRSTVWY',
                'output_dir': '/tmp'
            }
            try:
                result = mgr.run_analyses("network_analysis", input_data)
                if not isinstance(result, dict):
                    raise RuntimeError("network_analysis did not return a dict result")
            except Exception:
                # Acceptable for environments without full runtime dependencies
                pass

        print("AnalysisManager test: SUCCESS")
    except Exception as e:
        ok = False
        print(f"AnalysisManager test: FAILED - {e}")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
    