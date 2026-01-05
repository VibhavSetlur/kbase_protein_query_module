"""
UniProt ID Input Handler
"""

import logging
import time
import requests
from typing import Dict, Any, List, Optional
import re

from Bio.Seq import Seq
from Bio.SeqUtils import IUPACData

logger = logging.getLogger(__name__)

class UniProtIdProcessor:
    """Handles UniProt ID input processing."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the UniProt ID Processor.
        
        Args:
            config: Configuration dictionary.
        """
        self.config = config or {}
        self.max_retries = self.config.get('max_retries', 3)
        self.timeout = self.config.get('timeout', 30)
        self.api_base_url = "https://rest.uniprot.org/uniprotkb"
    
    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process UniProt ID input.
        
        Args:
            input_data: Raw input dictionary containing 'uniprot_id'.
            
        Returns:
            Dictionary containing processed proteins and metadata.
        """
        start_time = time.time()
        
        try:
            logger.info("Processing UniProt ID input")
            
            uniprot_id = input_data.get('uniprot_id', [])
            if not uniprot_id:
                raise ValueError("No UniProt ID provided")
            
            # Normalize to list
            if isinstance(uniprot_id, str):
                uniprot_id = [uniprot_id.strip()]
            
            # Validate and fetch
            proteins = {}
            idx = 0
            for uid in uniprot_id:
                if self._validate(uid):
                    protein_data = self._fetch(uid)
                    if protein_data:
                        protein_id = f"protein_{idx}"
                        proteins[protein_id] = {
                            'sequence': protein_data['sequence'],
                            'source': 'uniprot',
                            'original_id': protein_data['protein_id']
                        }
                        idx += 1
                else:
                    logger.warning(f"Invalid UniProt ID: {uid}")
            
            if not proteins:
                raise ValueError("No valid UniProt proteins found")
            
            processing_time = time.time() - start_time
            logger.info(f"Processed {len(proteins)} UniProt proteins")
            
            return {
                'success': True,
                'proteins': proteins,
                'protein_count': len(proteins),
                'processing_time': processing_time
            }
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"UniProt ID processing failed: {e}")
            return {
                'success': False,
                'proteins': {},
                'protein_count': 0,
                'processing_time': processing_time,
                'error_message': str(e)
            }
    
    def _validate(self, uniprot_id: str) -> bool:
        """
        Validate UniProt ID format.
        
        Args:
            uniprot_id: UniProt ID string.
            
        Returns:
            True if valid, False otherwise.
        """
        if not uniprot_id or not isinstance(uniprot_id, str):
            return False
        uid = uniprot_id.strip().upper()
        return bool(re.match(r'^[A-Z][0-9A-Z]{5,9}$', uid))
    
    def _fetch(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Fetch protein data from UniProt API.
        
        Args:
            uniprot_id: UniProt ID string.
            
        Returns:
            Dictionary with protein data or None if failed.
        """
        for attempt in range(self.max_retries):
            try:
                url = f"{self.api_base_url}/{uniprot_id}"
                response = requests.get(url, timeout=self.timeout)
                response.raise_for_status()
                data = response.json()
                
                sequence_data = data.get('sequence', {})
                sequence = sequence_data.get('value', '')
                protein_id = data.get('primaryAccession', '')
                
                if not protein_id or not sequence:
                    return None
                
                validated_sequence = self._validate_sequence(sequence)
                if not validated_sequence:
                    logger.warning(f"Invalid sequence for UniProt ID {uniprot_id}")
                    return None
                
                return {
                    'protein_id': protein_id,
                    'sequence': validated_sequence
                }
                
            except requests.exceptions.HTTPError as e:
                if hasattr(e, 'response') and e.response and e.response.status_code == 404:
                    logger.info(f"UniProt ID {uniprot_id} not found (404)")
                    return None
                if attempt == self.max_retries - 1:
                    return None
                time.sleep(1)
            except Exception as e:
                logger.warning(f"Error fetching {uniprot_id}: {e}")
                if attempt == self.max_retries - 1:
                    return None
                time.sleep(1)
        
        return None
    
    def _validate_sequence(self, sequence: str) -> Optional[str]:
        """
        Validate and clean protein sequence using BioSeq.
        
        Args:
            sequence: Input sequence string.
            
        Returns:
            Cleaned sequence string or None if invalid.
        """
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
        """
        Standardize input_data format.
        
        Args:
            input_data: Raw input dictionary.
            
        Returns:
            Standardized input dictionary.
        """
        input_type = input_data.get('input_type', 'uniprot_id')
        
        standardized = {
            'input_type': input_type,
            'uniprot_id': input_data.get('uniprot_id', []),
            'protein_sequence': '',  # Empty for uniprot_id
            'analysis_name': input_data.get('analysis_name', '')
        }
        
        return standardized
