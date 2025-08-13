"""
Data Extraction Stage for KBase Protein Query Module

This stage extracts protein data from various sources and formats.
"""

import logging
import time
import os
import requests
from typing import Dict, Any, List, Optional, Union

from ..base_stage import BaseStage, StageResult
from ...utils import input_parser
from ...storage import ProteinExistenceChecker

logger = logging.getLogger(__name__)

class DataExtractionStage(BaseStage):
    """
    Extracts protein data from various sources.
    
    Handles:
    - FASTA file/string parsing
    - UniProt data retrieval
    - Workspace object extraction
    - Data format conversion
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.max_retries = config.get('max_retries', 3) if config else 3
        self.timeout = config.get('timeout', 30) if config else 30
        self.batch_size = config.get('batch_size', 100) if config else 100
    
    def get_stage_name(self) -> str:
        return "data_extraction"
    
    def get_required_inputs(self) -> List[str]:
        return ['validated_input']
    
    def get_optional_inputs(self) -> List[str]:
        return ['workspace_client', 'extraction_config']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input parameters."""
        if 'validated_input' not in input_data:
            logger.error("Missing validated_input")
            return False
        
        validated_input = input_data['validated_input']
        if not hasattr(validated_input, 'input_type'):
            logger.error("Invalid validated_input format")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'extracted_data': {
                'type': 'object',
                'properties': {
                    'protein_records': {'type': 'array'},
                    'extraction_stats': {'type': 'object'},
                    'metadata': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Extract protein data from validated input."""
        start_time = time.time()
        
        try:
            validated_input = input_data['validated_input']
            extraction_config = input_data.get('extraction_config', {})
            
            # Extract data based on input type
            if validated_input.input_type == 'uniprot_identifiers':
                extracted_data = self._extract_uniprot_data(validated_input, extraction_config)
            elif validated_input.input_type in ['fasta', 'fasta_file']:
                extracted_data = self._extract_fasta_data(validated_input, extraction_config)
            elif validated_input.input_type == 'workspace_object':
                extracted_data = self._extract_workspace_data(validated_input, extraction_config, workspace_client)
            else:
                extracted_data = self._extract_generic_data(validated_input, extraction_config)
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=True,
                output_data={
                    'extracted_data': extracted_data
                },
                metadata={
                    'input_type': validated_input.input_type,
                    'extraction_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Data extraction failed: {str(e)}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _extract_uniprot_data(self, validated_input, extraction_config: Dict[str, Any]) -> Dict[str, Any]:
        """Extract data from UniProt identifiers."""
        protein_records = []
        extraction_stats = {
            'total_requested': len(validated_input.protein_records),
            'successfully_extracted': 0,
            'failed_extractions': 0,
            'errors': []
        }
        
        # Process in batches
        for i in range(0, len(validated_input.protein_records), self.batch_size):
            batch = validated_input.protein_records[i:i + self.batch_size]
            
            for record in batch:
                try:
                    # Extract UniProt data
                    uniprot_data = self._fetch_uniprot_data(record.id)
                    if uniprot_data:
                        protein_records.append(uniprot_data)
                        extraction_stats['successfully_extracted'] += 1
                    else:
                        extraction_stats['failed_extractions'] += 1
                        extraction_stats['errors'].append(f"No data found for {record.id}")
                        
                except Exception as e:
                    extraction_stats['failed_extractions'] += 1
                    extraction_stats['errors'].append(f"Failed to extract {record.id}: {str(e)}")
        
        return {
            'protein_records': protein_records,
            'extraction_stats': extraction_stats,
            'metadata': {
                'source': 'uniprot',
                'extraction_method': 'api'
            }
        }
    
    def _extract_fasta_data(self, validated_input, extraction_config: Dict[str, Any]) -> Dict[str, Any]:
        """Extract data from FASTA format."""
        protein_records = validated_input.protein_records
        extraction_stats = {
            'total_requested': len(protein_records),
            'successfully_extracted': len(protein_records),
            'failed_extractions': 0,
            'errors': []
        }
        
        return {
            'protein_records': protein_records,
            'extraction_stats': extraction_stats,
            'metadata': {
                'source': validated_input.input_source,
                'extraction_method': 'fasta_parsing'
            }
        }
    
    def _extract_workspace_data(self, validated_input, extraction_config: Dict[str, Any], workspace_client=None) -> Dict[str, Any]:
        """Extract data from workspace object."""
        if not workspace_client:
            raise ValueError("Workspace client required for workspace data extraction")
        
        try:
            # Extract data from workspace object
            object_data = workspace_client.get_objects2({
                'objects': [{'ref': validated_input.input_source}]
            })
            
            # Parse workspace object data
            protein_records = self._parse_workspace_object(object_data['data'][0])
            
            extraction_stats = {
                'total_requested': 1,
                'successfully_extracted': len(protein_records),
                'failed_extractions': 0,
                'errors': []
            }
            
            return {
                'protein_records': protein_records,
                'extraction_stats': extraction_stats,
                'metadata': {
                    'source': validated_input.input_source,
                    'extraction_method': 'workspace_object'
                }
            }
            
        except Exception as e:
            raise ValueError(f"Failed to extract workspace data: {str(e)}")
    
    def _extract_generic_data(self, validated_input, extraction_config: Dict[str, Any]) -> Dict[str, Any]:
        """Extract generic data."""
        extraction_stats = {
            'total_requested': 1,
            'successfully_extracted': 1,
            'failed_extractions': 0,
            'errors': []
        }
        
        return {
            'protein_records': validated_input.protein_records,
            'extraction_stats': extraction_stats,
            'metadata': {
                'source': validated_input.input_source,
                'extraction_method': 'generic'
            }
        }
    
    def _fetch_uniprot_data(self, uniprot_id: str) -> Optional[input_parser.ProteinRecord]:
        """Fetch data from UniProt API."""
        for attempt in range(self.max_retries):
            try:
                url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
                response = requests.get(url, timeout=self.timeout)
                response.raise_for_status()
                
                data = response.json()
                
                # Extract sequence
                sequence = data.get('sequence', {}).get('value', '')
                
                # Extract metadata
                metadata = {
                    'uniprot_id': uniprot_id,
                    'entry_name': data.get('entryName', ''),
                    'protein_name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'organism': data.get('organism', {}).get('scientificName', ''),
                    'length': len(sequence)
                }
                
                return input_parser.ProteinRecord(
                    id=uniprot_id,
                    sequence=sequence,
                    source='uniprot',
                    metadata=metadata
                )
                
            except requests.RequestException as e:
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed to fetch UniProt data for {uniprot_id}: {str(e)}")
                    return None
                time.sleep(1)  # Brief delay before retry
        
        return None
    
    def _parse_workspace_object(self, object_data: Dict[str, Any]) -> List[input_parser.ProteinRecord]:
        """Parse workspace object data into protein records."""
        protein_records = []
        
        # Handle different workspace object types
        object_type = object_data.get('info', [{}])[0].get('type_name', '')
        
        if 'ProteinSequenceSet' in object_type:
            # Parse protein sequence set
            data = object_data.get('data', {})
            sequences = data.get('sequences', [])
            
            for seq_data in sequences:
                protein_records.append(input_parser.ProteinRecord(
                    id=seq_data.get('id', ''),
                    sequence=seq_data.get('sequence', ''),
                    source='workspace',
                    metadata={'object_type': object_type}
                ))
        
        elif 'Genome' in object_type:
            # Parse genome object
            data = object_data.get('data', {})
            features = data.get('features', [])
            
            for feature in features:
                if feature.get('type') == 'CDS':
                    protein_id = feature.get('id', '')
                    sequence = feature.get('protein_translation', '')
                    
                    if sequence:
                        protein_records.append(input_parser.ProteinRecord(
                            id=protein_id,
                            sequence=sequence,
                            source='workspace',
                            metadata={'object_type': object_type, 'feature_type': 'CDS'}
                        ))
        
        else:
            # Generic parsing
            logger.warning(f"Unknown workspace object type: {object_type}")
        
        return protein_records
