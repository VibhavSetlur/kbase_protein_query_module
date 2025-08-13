"""
Input Validation Stage for KBase Protein Query Module

This stage validates and preprocesses input data from various sources.
"""

import logging
import time
import os
import json
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

from ..base_stage import BaseStage, StageResult
from ...utils import input_parser
from ...storage import ProteinExistenceChecker

logger = logging.getLogger(__name__)

@dataclass
class InputData:
    """Container for validated input data."""
    protein_records: List[input_parser.ProteinRecord]
    input_type: str
    input_source: str
    validation_summary: Dict[str, Any]
    metadata: Dict[str, Any] = None

class InputValidationStage(BaseStage):
    """
    Validates and preprocesses input data from various sources.
    
    Handles:
    - Input type detection
    - Basic validation
    - Format standardization
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.max_input_size = config.get('max_input_size', 1000) if config else 1000
        self.allowed_input_types = config.get('allowed_input_types', [
            'fasta', 'fasta_file', 'fasta_url', 'uniprot_identifiers',
            'protein_sequence_set', 'genome_reference', 'single_protein',
            'feature_set', 'genome_set', 'workspace_object'
        ]) if config else [
            'fasta', 'fasta_file', 'fasta_url', 'uniprot_identifiers',
            'protein_sequence_set', 'genome_reference', 'single_protein',
            'feature_set', 'genome_set', 'workspace_object'
        ]
    
    def get_stage_name(self) -> str:
        return "input_validation"
    
    def get_required_inputs(self) -> List[str]:
        return ['input_type', 'input_data']
    
    def get_optional_inputs(self) -> List[str]:
        return ['workspace_client', 'validation_config']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input parameters."""
        required = self.get_required_inputs()
        for field in required:
            if field not in input_data:
                logger.error(f"Missing required input: {field}")
                return False
        
        input_type = input_data.get('input_type')
        if input_type not in self.allowed_input_types:
            logger.error(f"Invalid input type: {input_type}")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'validated_input': {
                'type': 'object',
                'properties': {
                    'input_type': {'type': 'string'},
                    'input_source': {'type': 'string'},
                    'validation_summary': {'type': 'object'},
                    'metadata': {'type': 'object'}
                }
            },
            'validation_stats': {
                'type': 'object',
                'properties': {
                    'total_inputs': {'type': 'integer'},
                    'valid_inputs': {'type': 'integer'},
                    'invalid_inputs': {'type': 'integer'},
                    'warnings': {'type': 'array'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Validate input data."""
        start_time = time.time()
        
        try:
            input_type = input_data['input_type']
            input_data_raw = input_data['input_data']
            validation_config = input_data.get('validation_config', {})
            
            # Initialize validation summary
            validation_summary = {
                'total_inputs': 0,
                'valid_inputs': 0,
                'invalid_inputs': 0,
                'warnings': [],
                'errors': []
            }
            
            # Validate based on input type
            if input_type == 'uniprot_identifiers':
                validated_data = self._validate_uniprot_identifiers(input_data_raw, validation_summary)
            elif input_type == 'fasta':
                validated_data = self._validate_fasta_data(input_data_raw, validation_summary)
            elif input_type == 'fasta_file':
                validated_data = self._validate_fasta_file(input_data_raw, validation_summary)
            elif input_type == 'workspace_object':
                validated_data = self._validate_workspace_object(input_data_raw, validation_summary, workspace_client)
            else:
                validated_data = self._validate_generic_input(input_data_raw, input_type, validation_summary)
            
            # Check input size limits
            if validation_summary['total_inputs'] > self.max_input_size:
                validation_summary['warnings'].append(
                    f"Input size ({validation_summary['total_inputs']}) exceeds limit ({self.max_input_size})"
                )
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=True,
                output_data={
                    'validated_input': validated_data,
                    'validation_stats': validation_summary
                },
                metadata={
                    'input_type': input_type,
                    'validation_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Input validation failed: {str(e)}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _validate_uniprot_identifiers(self, identifiers: List[str], summary: Dict[str, Any]) -> InputData:
        """Validate UniProt identifiers."""
        summary['total_inputs'] = len(identifiers)
        
        valid_identifiers = []
        for identifier in identifiers:
            if self._is_valid_uniprot_id(identifier):
                valid_identifiers.append(identifier)
                summary['valid_inputs'] += 1
            else:
                summary['invalid_inputs'] += 1
                summary['warnings'].append(f"Invalid UniProt identifier: {identifier}")
        
        return InputData(
            protein_records=[input_parser.ProteinRecord(id=id, sequence="", source="uniprot") for id in valid_identifiers],
            input_type="uniprot_identifiers",
            input_source="uniprot",
            validation_summary=summary
        )
    
    def _validate_fasta_data(self, fasta_data: str, summary: Dict[str, Any]) -> InputData:
        """Validate FASTA data."""
        try:
            protein_records = input_parser.parse_fasta_string(fasta_data)
            summary['total_inputs'] = len(protein_records)
            summary['valid_inputs'] = len(protein_records)
            
            return InputData(
                protein_records=protein_records,
                input_type="fasta",
                input_source="fasta_string",
                validation_summary=summary
            )
        except Exception as e:
            summary['errors'].append(f"FASTA parsing failed: {str(e)}")
            raise
    
    def _validate_fasta_file(self, file_path: str, summary: Dict[str, Any]) -> InputData:
        """Validate FASTA file."""
        try:
            protein_records = input_parser.parse_fasta_file(file_path)
            summary['total_inputs'] = len(protein_records)
            summary['valid_inputs'] = len(protein_records)
            
            return InputData(
                protein_records=protein_records,
                input_type="fasta_file",
                input_source=file_path,
                validation_summary=summary
            )
        except Exception as e:
            summary['errors'].append(f"FASTA file parsing failed: {str(e)}")
            raise
    
    def _validate_workspace_object(self, object_ref: str, summary: Dict[str, Any], workspace_client=None) -> InputData:
        """Validate workspace object."""
        if not workspace_client:
            raise ValueError("Workspace client required for workspace object validation")
        
        try:
            # Validate object reference format
            if not self._is_valid_workspace_ref(object_ref):
                raise ValueError(f"Invalid workspace reference format: {object_ref}")
            
            summary['total_inputs'] = 1
            summary['valid_inputs'] = 1
            
            return InputData(
                protein_records=[],
                input_type="workspace_object",
                input_source=object_ref,
                validation_summary=summary
            )
        except Exception as e:
            summary['errors'].append(f"Workspace object validation failed: {str(e)}")
            raise
    
    def _validate_generic_input(self, data: Any, input_type: str, summary: Dict[str, Any]) -> InputData:
        """Validate generic input data."""
        summary['total_inputs'] = 1
        summary['valid_inputs'] = 1
        
        return InputData(
            protein_records=[],
            input_type=input_type,
            input_source="generic",
            validation_summary=summary
        )
    
    def _is_valid_uniprot_id(self, identifier: str) -> bool:
        """Check if identifier is a valid UniProt ID."""
        # Basic UniProt ID validation
        if not identifier or len(identifier) < 3:
            return False
        
        # Check for common UniProt ID patterns
        valid_patterns = [
            r'^[A-Z][0-9A-Z]{5}$',  # UniProtKB/Swiss-Prot format
            r'^[A-Z][0-9A-Z]{9}$',  # UniProtKB/TrEMBL format
            r'^[A-Z][0-9A-Z]{4}$',  # Short format
        ]
        
        import re
        for pattern in valid_patterns:
            if re.match(pattern, identifier):
                return True
        
        return False
    
    def _is_valid_workspace_ref(self, ref: str) -> bool:
        """Check if reference is a valid workspace reference."""
        # Basic workspace reference validation
        if not ref or '/' not in ref:
            return False
        
        parts = ref.split('/')
        if len(parts) != 2:
            return False
        
        try:
            int(parts[0])  # Workspace ID should be numeric
            int(parts[1])  # Object ID should be numeric
            return True
        except ValueError:
            return False
