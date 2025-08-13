"""
Chart Generator Component

This module generates chart components for HTML reports.
"""

import logging
from typing import Dict, Any, List

logger = logging.getLogger(__name__)

class ChartGenerator:
    """Generates chart components for HTML reports."""
    
    def __init__(self):
        """Initialize the chart generator."""
        pass
    
    def generate_similarity_chart(self, similarities: List[float]) -> str:
        """Generate similarity distribution chart."""
        return f"""
        <div class="chart-container">
            <h3>Similarity Distribution</h3>
            <div id="similarity-chart" class="chart">
                <p>Similarity scores: {len(similarities)} proteins</p>
                <p>Average similarity: {sum(similarities)/len(similarities):.3f}</p>
            </div>
        </div>
        """
    
    def generate_family_chart(self, family_data: Dict[str, int]) -> str:
        """Generate family distribution chart."""
        return f"""
        <div class="chart-container">
            <h3>Family Distribution</h3>
            <div id="family-chart" class="chart">
                <p>Families found: {len(family_data)}</p>
                <p>Total proteins: {sum(family_data.values())}</p>
            </div>
        </div>
        """
