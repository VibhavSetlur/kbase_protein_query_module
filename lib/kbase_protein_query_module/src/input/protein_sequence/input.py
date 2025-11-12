"""
Protein Sequence Input Handler
"""

import logging
import time
from typing import Dict, Any, List, Optional
from io import StringIO

from Bio.Seq import Seq
from Bio.SeqUtils import IUPACData
from Bio import SeqIO

logger = logging.getLogger(__name__)

class ProteinSequenceProcessor:
    """Handles protein sequence input processing."""
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
    
    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Process protein sequence input."""
        start_time = time.time()
        
        try:
            logger.info("Processing protein sequence input")
            
            sequence = input_data.get('protein_sequence', '')
            if not sequence:
                raise ValueError("No protein sequence provided")
            
            proteins = self._parse_and_validate(sequence)
            
            if not proteins:
                raise ValueError("No valid protein sequences found")
            
            processing_time = time.time() - start_time
            logger.info(f"Processed {len(proteins)} protein sequences")
            
            return {
                'success': True,
                'proteins': proteins,
                'protein_count': len(proteins),
                'processing_time': processing_time
            }
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"Protein sequence processing failed: {e}")
            return {
                'success': False,
                'proteins': {},
                'protein_count': 0,
                'processing_time': processing_time,
                'error_message': str(e)
            }
    
    def _parse_and_validate(self, sequence: str) -> Dict[str, Dict[str, Any]]:
        """Parse and validate protein sequences using BioSeq.
        
        Returns nested dictionary: {'protein_0': {...}, 'protein_1': {...}, ...}
        """
        proteins = {}
        
        sequence = sequence.strip()
        if not sequence:
            return proteins
        
        try:
            if sequence.startswith('>'):
                records = list(SeqIO.parse(StringIO(sequence), "fasta"))
                for idx, record in enumerate(records):
                    validated_seq = self._validate_sequence(str(record.seq))
                    if validated_seq:
                        protein_id = f"protein_{idx}"
                        proteins[protein_id] = {
                            'sequence': validated_seq,
                            'source': 'protein_sequence',
                            'original_id': record.id or protein_id
                        }
            else:
                sequences = self._parse_multiple_formats(sequence)
                for idx, seq in enumerate(sequences):
                    validated_seq = self._validate_sequence(seq)
                    if validated_seq:
                        protein_id = f"protein_{idx}"
                        proteins[protein_id] = {
                            'sequence': validated_seq,
                            'source': 'protein_sequence',
                            'original_id': protein_id
                        }
        except Exception as e:
            logger.warning(f"BioSeq parsing failed, attempting fallback: {e}")
            sequences = self._parse_multiple_formats(sequence.replace('\n', ' '))
            for idx, seq in enumerate(sequences):
                validated_seq = self._validate_sequence(seq)
                if validated_seq:
                    protein_id = f"protein_{idx}"
                    proteins[protein_id] = {
                        'sequence': validated_seq,
                        'source': 'protein_sequence',
                        'original_id': protein_id
                    }
        
        return proteins
    
    def _parse_multiple_formats(self, sequence: str) -> List[str]:
        """Parse sequence from multiple formats: comma-separated, tab-separated, or newline-separated."""
        sequences = []
        
        if not sequence:
            return sequences
        
        sequence = sequence.strip()
        
        if ',' in sequence:
            sequences = [s.strip() for s in sequence.split(',') if s.strip()]
        elif '\t' in sequence:
            sequences = [s.strip() for s in sequence.split('\t') if s.strip()]
        elif '\n' in sequence:
            sequences = [s.strip() for s in sequence.split('\n') if s.strip()]
        else:
            sequences = [sequence]
        
        return sequences
    
    def _validate_sequence(self, sequence: str) -> Optional[str]:
        """Validate and clean protein sequence using BioSeq."""
        if not sequence:
            return None
        
        try:
            seq_obj = Seq(sequence.upper())
            
            valid_chars = set(IUPACData.protein_letters)
            valid_chars.update(['-', '*', 'X', 'B', 'Z', 'J', 'U', 'O'])
            
            cleaned_chars = []
            for char in str(seq_obj):
                if char in valid_chars or char.isspace():
                    if not char.isspace():
                        cleaned_chars.append(char)
            
            cleaned = ''.join(cleaned_chars)
            
            if not cleaned:
                return None
            
            if len(cleaned) < 1:
                return None
            
            return cleaned
        except Exception as e:
            logger.warning(f"Sequence validation error: {e}")
            return None
    
    
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
