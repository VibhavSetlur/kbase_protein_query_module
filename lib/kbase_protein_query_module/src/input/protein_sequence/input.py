"""
Protein Sequence Input Handler
"""

import logging
import time
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)

class ProteinSequenceProcessor:
    """Handles protein sequence input processing."""
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    
    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Process protein sequence input."""
        start_time = time.time()
        
        try:
            logger.info("Processing protein sequence input")
            
            sequence = input_data.get('protein_sequence', '')
            if not sequence:
                raise ValueError("No protein sequence provided")
            
            # Parse and validate
            proteins = self._parse(sequence)
            
            if not proteins:
                raise ValueError("No valid protein sequences found")
            
            processing_time = time.time() - start_time
            logger.info(f"Processed {len(proteins)} protein sequences")
            
            return {
                'success': True,
                'proteins': proteins,
                'processing_time': processing_time
            }
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"Protein sequence processing failed: {e}")
            return {
                'success': False,
                'proteins': [],
                'processing_time': processing_time,
                'error_message': str(e)
            }
    
    def _parse(self, sequence: str) -> List[Dict[str, Any]]:
        """Parse protein sequence."""
        proteins = []
        
        # Check if FASTA format
        if sequence.strip().startswith('>'):
            # FASTA format
            lines = sequence.strip().split('\n')
            current_id = None
            current_seq = []
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        proteins.append(self._create_protein(current_id, ''.join(current_seq)))
                    current_id = line[1:].strip()
                    current_seq = []
                elif line:
                    current_seq.append(line)
            
            if current_id and current_seq:
                proteins.append(self._create_protein(current_id, ''.join(current_seq)))
        else:
            # Direct sequence
            cleaned = self._clean(sequence)
            if cleaned:
                proteins.append(self._create_protein('protein_1', cleaned))
        
        return proteins
    
    def _create_protein(self, protein_id: str, sequence: str) -> Dict[str, Any]:
        """Create protein record."""
        return {
            'protein_id': protein_id,
            'sequence': sequence,
            'source': 'protein_sequence'
        }
    
    def _clean(self, sequence: str) -> str:
        """Clean and validate protein sequence."""
        cleaned = ''.join(char for char in sequence.upper() if char in self.valid_amino_acids)
        return cleaned
    
    def standardize_input(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Standardize input_data format."""
        input_type = input_data.get('input_type', 'protein_sequence')
        
        standardized = {
            'input_type': input_type,
            'protein_sequence': input_data.get('protein_sequence', ''),
            'uniprot_id': [],  # Empty for protein_sequence
            'analysis_name': input_data.get('analysis_name', '')
        }
        
        return standardized
