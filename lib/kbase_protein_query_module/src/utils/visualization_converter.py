"""
Visualization Converter Utility

Simple utility to convert visualizations to standalone HTML files.
No fancy reports - just raw visualization outputs.
"""

import os
import logging
import base64
from io import BytesIO
from typing import Any

logger = logging.getLogger(__name__)

# Import libraries
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


class VisualizationConverter:
    """Simple utility to convert visualizations to standalone HTML files."""
    
    def detect_visualization_type(self, visualization_obj: Any) -> str:
        """Detect the type of visualization object."""
        if hasattr(visualization_obj, 'to_html') and hasattr(visualization_obj, 'data'):
            return 'plotly'
        elif hasattr(visualization_obj, 'savefig'):
            return 'matplotlib'
        else:
            return 'unknown'
    
    def save_visualization(self, visualization_obj: Any, output_path: str) -> str:
        """Save any visualization to a standalone HTML file."""
        viz_type = self.detect_visualization_type(visualization_obj)
        
        logger.info(f"Saving {viz_type} visualization to: {output_path}")
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        if viz_type == 'plotly' and PLOTLY_AVAILABLE:
            visualization_obj.write_html(output_path, include_plotlyjs='cdn')
        elif viz_type == 'matplotlib' and MATPLOTLIB_AVAILABLE:
            self._save_matplotlib_html(visualization_obj, output_path)
        else:
            logger.warning(f"Unsupported visualization type: {viz_type}")
        
        return output_path
    
    def _save_matplotlib_html(self, fig, output_path: str):
        """Save matplotlib figure as standalone HTML."""
        buffer = BytesIO()
        fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        img_base64 = base64.b64encode(buffer.getvalue()).decode()
        plt.close(fig)
        
        html_content = f'''<!DOCTYPE html>
<html>
<head><title>Visualization</title></head>
<body style="margin:0; padding:20px; text-align:center;">
    <img src="data:image/png;base64,{img_base64}" style="max-width:100%; height:auto;">
</body>
</html>'''
        
        with open(output_path, 'w') as f:
            f.write(html_content)


def save_plot_html(plot_obj: Any, output_path: str) -> str:
    """Quick function to save any plot as HTML."""
    converter = VisualizationConverter()
    return converter.save_visualization(plot_obj, output_path)
