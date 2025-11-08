"""
Analysis Configuration Module

Defines available analyses and their runtime enablement.
"""

from typing import Dict, Any, List
import logging

logger = logging.getLogger(__name__)

def _check_deps(requires_deps: List[str]) -> bool:
    """Check if required dependencies are available."""
    logger.debug(f"Checking dependencies: {requires_deps}")
    for dep in requires_deps:
        try:
            __import__(dep)
            logger.debug(f"Dependency '{dep}' available")
        except ImportError:
            logger.warning(f"Dependency '{dep}' not available")
            return False
    return True

# Analysis metadata
_ANALYSIS_BASE: Dict[str, Dict[str, Any]] = {
    "network_analysis": {
        "name": "Network Analysis",
        "description": "Protein similarity network analysis with interactive visualizations",
        "category": "network",
        "module_path": "kbase_protein_query_module.src.analysis.network_analysis.network_analysis",
        "class_name": "NetworkAnalysis"
#        "requires_deps": ["networkx", "sklearn"],
    }
}

def get_enabled_analyses() -> Dict[str, Dict[str, Any]]:
    """Return enabled analyses based on dependency availability."""
    logger.debug("Getting enabled analyses")
    enabled: Dict[str, Dict[str, Any]] = {}

    for name, meta in _ANALYSIS_BASE.items():
        requires_deps = meta.get("requires_deps", [])
        if requires_deps:
            if not _check_deps(requires_deps):
                logger.info(f"Analysis '{name}' disabled due to missing dependencies")
                continue
        cfg = dict(meta)
        cfg["enabled"] = True
        enabled[name] = cfg
        logger.debug(f"Analysis '{name}' enabled")
    
    logger.info(f"Found {len(enabled)} enabled analyses: {list(enabled.keys())}")
    return enabled



def main() -> int:
    """
    Test config module.
    - Lists enabled analyses
    - Validates returned structure
    """
    ok = True
    try:
        enabled = get_enabled_analyses()
        if not isinstance(enabled, dict):
            raise RuntimeError("get_enabled_analyses did not return a dict")
        # Validate each entry has required keys
        for name, cfg in enabled.items():
            for key in ("module_path", "class_name", "name"):
                if key not in cfg:
                    raise RuntimeError(f"Enabled analysis '{name}' missing key '{key}'")
        print(f"Config test: SUCCESS - enabled={[k for k in enabled.keys()]}")
    except Exception as e:
        ok = False
        print(f"Config test: FAILED - {e}")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
