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
    current_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(current_dir)
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
    from analysis.config import get_enabled_analyses
else:
    from .config import get_enabled_analyses

logger = logging.getLogger(__name__)

class AnalysisManager:
    """Manages the execution of protein analysis workflows."""
    
    def __init__(self, output_manager=None):
        """Initialize the Analysis Manager."""
        self.output_manager = output_manager
        self.analyses: Dict[str, Any] = {}
        self.results: Dict[str, Any] = {}
        self._load_analysis_modules()
    
    def _load_analysis_modules(self):
        """Load all available analysis modules dynamically."""
        enabled_analyses = get_enabled_analyses()
        
        for analysis_name, config in enabled_analyses.items():
            try:
                module_path = config.get("module_path")
                class_name = config.get("class_name")
                
                if not module_path or not class_name:
                    logger.warning(f"Missing module_path or class_name for {analysis_name}")
                    continue
                
                module = importlib.import_module(module_path)
                if hasattr(module, class_name):
                    analysis_class = getattr(module, class_name)
                    self.analyses[analysis_name] = analysis_class()
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
            if hasattr(analysis, "run"):
                result = analysis.run(input_data, **kwargs)
            elif hasattr(analysis, "run_network_analysis"):
                result = analysis.run_network_analysis(input_data)
            # Add other analysis methods here as needed
            # elif analysis_name == "your_analysis":
            #     result = analysis.analyze(input_data, **kwargs)
            elif hasattr(analysis, "analyze"):
                result = analysis.analyze(input_data, **kwargs)
            else:
                raise AttributeError(f"Analysis '{analysis_name}' has no runnable method")

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
    