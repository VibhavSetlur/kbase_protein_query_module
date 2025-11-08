"""
UniProt ID Input Handler
"""

import logging
import time
import requests
from typing import Dict, Any, List, Optional
import re

logger = logging.getLogger(__name__)

class UniProtIdProcessor:
    """Handles UniProt ID input processing."""
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.max_retries = self.config.get('max_retries', 3)
        self.timeout = self.config.get('timeout', 30)
        self.api_base_url = "https://rest.uniprot.org/uniprotkb"
    
    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Process UniProt ID input."""
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
            proteins = []
            for uid in uniprot_id:
                if self._validate(uid):
                    protein = self._fetch(uid)
                    if protein:
                        proteins.append(protein)
                else:
                    logger.warning(f"Invalid UniProt ID: {uid}")
            
            if not proteins:
                raise ValueError("No valid UniProt proteins found")
            
            processing_time = time.time() - start_time
            logger.info(f"Processed {len(proteins)} UniProt proteins")
            
            return {
                'success': True,
                'proteins': proteins,
                'processing_time': processing_time
            }
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"UniProt ID processing failed: {e}")
            return {
                'success': False,
                'proteins': [],
                'processing_time': processing_time,
                'error_message': str(e)
            }
    
    def _validate(self, uniprot_id: str) -> bool:
        """Validate UniProt ID format."""
        if not uniprot_id or not isinstance(uniprot_id, str):
            return False
        uid = uniprot_id.strip().upper()
        return bool(re.match(r'^[A-Z][0-9A-Z]{5,9}$', uid))
    
    def _fetch(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """Fetch protein data from UniProt API."""
        for attempt in range(self.max_retries):
            try:
                url = f"{self.api_base_url}/{uniprot_id}"
                response = requests.get(url, timeout=self.timeout)
                response.raise_for_status()
                data = response.json()
                
                # Extract sequence
                sequence_data = data.get('sequence', {})
                sequence = sequence_data.get('value', '')
                protein_id = data.get('primaryAccession', '')
                
                if not protein_id or not sequence:
                    return None
                
                return {
                    'protein_id': protein_id,
                    'sequence': sequence,
                    'source': 'uniprot'
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
    
    def standardize_input(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Standardize input_data format."""
        input_type = input_data.get('input_type', 'uniprot_id')
        
        standardized = {
            'input_type': input_type,
            'uniprot_id': input_data.get('uniprot_id', []),
            'protein_sequence': '',  # Empty for uniprot_id
            'analysis_name': input_data.get('analysis_name', '')
        }
        
        return standardized
