"""
Protein Details Generator Component

This module generates protein details components for HTML reports.
"""

import logging
from typing import Dict, Any, List

logger = logging.getLogger(__name__)

class ProteinDetailsGenerator:
    """Generates protein details components for HTML reports."""
    
    def __init__(self):
        """Initialize the protein details generator."""
        pass
    
    def generate_protein_details(self, protein_data: Dict[str, Any]) -> str:
        """Generate protein details HTML content."""
        return f"""
        <div class="protein-details">
            <h3>Protein Details</h3>
            <div class="protein-info">
                <p><strong>Protein ID:</strong> {protein_data.get('protein_id', 'Unknown')}</p>
                <p><strong>Sequence Length:</strong> {protein_data.get('sequence_length', 0)}</p>
                <p><strong>Family:</strong> {protein_data.get('family', 'Unknown')}</p>
                <p><strong>Similarity Score:</strong> {protein_data.get('similarity_score', 0):.3f}</p>
            </div>
        </div>
        """
    
    def generate_protein_table(self, proteins: List[Dict[str, Any]]) -> str:
        """Generate protein table HTML content."""
        table_rows = ""
        for protein in proteins:
            table_rows += f"""
            <tr>
                <td>{protein.get('protein_id', 'Unknown')}</td>
                <td>{protein.get('family', 'Unknown')}</td>
                <td>{protein.get('similarity_score', 0):.3f}</td>
                <td>{protein.get('sequence_length', 0)}</td>
            </tr>
            """
        
        return f"""
        <div class="protein-table">
            <h3>Protein Analysis Results</h3>
            <table class="protein-results-table">
                <thead>
                    <tr>
                        <th>Protein ID</th>
                        <th>Family</th>
                        <th>Similarity Score</th>
                        <th>Sequence Length</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>
        """
