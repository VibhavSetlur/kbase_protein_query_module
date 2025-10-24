"""
UniProt IDs Input Handler

Handles UniProt ID input processing with API integration and validation.
"""

import logging
import time
import requests
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class UniProtData:
    """Container for UniProt data."""
    protein_id: str
    sequence: str
    source: str = "uniprot"
    metadata: Dict[str, Any] = None

class UniProtIdsProcessor:
    """
    Handles UniProt ID input processing.
    
    Supports single IDs, lists, validation, and API data retrieval.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        # API configuration
        self.max_retries = self.config.get('max_retries', 3)
        self.timeout = self.config.get('timeout', 30)
        self.batch_size = self.config.get('batch_size', 100)
        self.api_base_url = "https://rest.uniprot.org/uniprotkb"
    
    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process UniProt ID input.
        
        Args:
            input_data: Input data containing UniProt IDs
            
        Returns:
            Processed data with proteins
        """
        start_time = time.time()
        
        try:
            logger.info("Processing UniProt ID input")
            
            uniprot_ids = input_data.get('uniprot_ids', [])
            if not uniprot_ids:
                raise ValueError("No UniProt IDs provided")
            
            # Normalize input to list
            if isinstance(uniprot_ids, str):
                # Handle comma-separated or newline-separated IDs
                if ',' in uniprot_ids:
                    uniprot_ids = [id.strip() for id in uniprot_ids.split(',') if id.strip()]
                elif '\n' in uniprot_ids:
                    uniprot_ids = [id.strip() for id in uniprot_ids.split('\n') if id.strip()]
                else:
                    uniprot_ids = [uniprot_ids.strip()]
            
            # Validate all IDs
            valid_ids = []
            invalid_ids = []
            
            for uniprot_id in uniprot_ids:
                if self._validate_uniprot_id(uniprot_id):
                    valid_ids.append(uniprot_id)
                else:
                    invalid_ids.append(uniprot_id)
            
            if not valid_ids:
                raise ValueError(f"No valid UniProt IDs found. Invalid IDs: {invalid_ids}")
            
            # Fetch data for valid IDs
            protein_records = []
            fetch_stats = {
                'total_requested': len(valid_ids),
                'successfully_fetched': 0,
                'failed_fetches': 0,
                'errors': []
            }
            
            # Process in batches
            for i in range(0, len(valid_ids), self.batch_size):
                batch = valid_ids[i:i + self.batch_size]
                
                for uniprot_id in batch:
                    try:
                        uniprot_data = self._fetch_uniprot_data(uniprot_id)
                        if uniprot_data:
                            protein_record = self._process_uniprot_data(uniprot_data)
                            if protein_record:
                                protein_records.append(protein_record)
                                fetch_stats['successfully_fetched'] += 1
                            else:
                                fetch_stats['failed_fetches'] += 1
                                fetch_stats['errors'].append(f"Failed to process data for {uniprot_id}")
                        else:
                            fetch_stats['failed_fetches'] += 1
                            fetch_stats['errors'].append(f"No data found for {uniprot_id}")
                    except Exception as e:
                        fetch_stats['failed_fetches'] += 1
                        fetch_stats['errors'].append(f"Failed to fetch {uniprot_id}: {str(e)}")
            
            processing_time = time.time() - start_time
            
            result = {
                'success': True,
                'input_type': 'uniprot_ids',
                'proteins': protein_records,
                'workspace_info': {},
                'processing_time': processing_time,
                'metadata': {
                    'fetch_stats': fetch_stats,
                    'invalid_ids': invalid_ids,
                    'valid_ids_count': len(valid_ids),
                    'processing_time': processing_time,
                    'source': 'uniprot'
                }
            }
            
            logger.info(f"Processed {len(protein_records)} proteins from UniProt")
            return result
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"UniProt ID processing failed: {e}")
            
            return {
                'success': False,
                'input_type': 'uniprot_ids',
                'proteins': [],
                'workspace_info': {},
                'error_message': str(e),
                'processing_time': processing_time,
                'metadata': {
                    'error': str(e),
                    'processing_time': processing_time
                }
            }
    
    def _validate_uniprot_id(self, uniprot_id: str) -> bool:
        """Validate UniProt ID format."""
        if not uniprot_id or not isinstance(uniprot_id, str):
            return False
        
        uniprot_id = uniprot_id.strip().upper()
        
        # Basic UniProt ID validation patterns
        valid_patterns = [
            r'^[A-Z][0-9A-Z]{5}$',  # UniProtKB/Swiss-Prot format (P12345)
            r'^[A-Z][0-9A-Z]{9}$',  # UniProtKB/TrEMBL format (A0A00000000)
        ]
        
        import re
        for pattern in valid_patterns:
            if re.match(pattern, uniprot_id):
                return True
        
        return False
    
    def _fetch_uniprot_data(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """Fetch data from UniProt API."""
        for attempt in range(self.max_retries):
            try:
                url = f"{self.api_base_url}/{uniprot_id}"
                response = requests.get(url, timeout=self.timeout)
                response.raise_for_status()
                
                # Return raw UniProt API data
                return response.json()
                
            except requests.exceptions.HTTPError as e:
                if hasattr(e, 'response') and e.response and e.response.status_code == 404:
                    logger.info(f"UniProt ID {uniprot_id} not found (404)")
                    return None
                logger.warning(f"HTTP error fetching UniProt ID {uniprot_id}: {e}")
                if attempt == self.max_retries - 1:
                    return None
                time.sleep(1)
            except requests.exceptions.RequestException as e:
                logger.warning(f"Request error fetching UniProt ID {uniprot_id}: {e}")
                if attempt == self.max_retries - 1:
                    return None
                time.sleep(1)
            except Exception as e:
                logger.error(f"Error processing UniProt data for {uniprot_id}: {str(e)}")
                return None
        
        return None
    
    def _extract_protein_name(self, data: Dict[str, Any]) -> str:
        """Extract protein name from UniProt data."""
        try:
            protein_description = data.get('proteinDescription', {})
            recommended_name = protein_description.get('recommendedName', {})
            full_name = recommended_name.get('fullName', {})
            if full_name.get('value'):
                return full_name.get('value')
            
            # Fallback to alternative names
            alternative_names = protein_description.get('alternativeNames', [])
            if alternative_names:
                for alt_name in alternative_names:
                    alt_full_name = alt_name.get('fullName', {})
                    if alt_full_name.get('value'):
                        return alt_full_name.get('value')
            
            return 'Unknown Protein'
        except Exception:
            return 'Unknown Protein'
    
    def _extract_organism(self, data: Dict[str, Any]) -> str:
        """Extract organism name from UniProt data."""
        try:
            organism = data.get('organism', {})
            return organism.get('scientificName', 'Unknown')
        except Exception:
            return 'Unknown'
    
    def _process_uniprot_data(self, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Process raw UniProt data into protein record format."""
        try:
            protein_id = data.get('primaryAccession', '')
            sequence_data = data.get('sequence', {})
            sequence = sequence_data.get('value', '')
            
            if not protein_id or not sequence:
                return None
            
            protein_name = self._extract_protein_name(data)
            organism = self._extract_organism(data)
            
            return {
                'protein_id': protein_id,
                'sequence': sequence,
                'source': 'uniprot',
                'metadata': {
                    'protein_name': protein_name,
                    'organism': organism,
                    'sequence_length': sequence_data.get('length', len(sequence)),
                    'molecular_weight': sequence_data.get('molWeight', 0),
                    'crc64': sequence_data.get('crc64', ''),
                    'uniprot_accession': protein_id
                }
            }
        except Exception:
            return None