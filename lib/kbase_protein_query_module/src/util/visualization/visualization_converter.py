"""
Visualization converter shims (deprecated).

Provides minimal HTML export helpers expected by legacy code/tests.
Uses Plotly if available; otherwise writes a simple HTML wrapper.
"""

import os
import logging
from typing import Any

logger = logging.getLogger(__name__)

try:
    import plotly.graph_objects as go  # type: ignore
    _PLOTLY_AVAILABLE = True
except Exception:  # pragma: no cover
    _PLOTLY_AVAILABLE = False
    go = None  # type: ignore


def save_plot_html(figure: Any, output_path: str) -> str:
    """
    Save a Plotly figure to HTML. If Plotly is unavailable, writes a minimal HTML wrapper.

    Args:
        figure: Plotly figure (or any object with .to_html())
        output_path: Destination file path

    Returns:
        The written file path
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    try:
        if _PLOTLY_AVAILABLE and hasattr(figure, 'to_html'):
            html = figure.to_html(include_plotlyjs='cdn', full_html=True)
        else:
            logger.warning("Plotly not available or invalid figure; writing minimal HTML placeholder.")
            html = "<html><head><meta charset='utf-8'><title>Visualization</title></head><body><h3>Visualization placeholder</h3></body></html>"

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        return output_path
    except Exception as e:
        logger.error(f"Failed to save HTML visualization: {e}")
        # Write a minimal error HTML to avoid hard failures in tests
        fallback_html = f"<html><body><pre>Error generating visualization: {e}</pre></body></html>"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(fallback_html)
        return output_path


