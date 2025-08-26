"""
Input Stages for KBase Protein Query Module

This module contains all input stages that handle data ingestion and validation.
Integrates with DataFileUtil for workspace object operations and file management.
"""

import logging
import time
import os
import json
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

from .base_stage import BaseStage, StageResult
from ..utils import input_parser
from ..storage import ProteinExistenceChecker

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
            raw_input_data = input_data['input_data']
            validation_config = input_data.get('validation_config', {})
            
            # Basic validation
            validation_summary = {
                'input_type': input_type,
                'input_size': len(str(raw_input_data)) if raw_input_data else 0,
                'validation_passed': True,
                'warnings': []
            }
            
            # Check input size
            if validation_summary['input_size'] > self.max_input_size * 1000:  # Rough size estimate
                validation_summary['warnings'].append(f"Input size may be too large: {validation_summary['input_size']} characters")
            
            # Type-specific validation
            if input_type == 'fasta':
                if not isinstance(raw_input_data, str) or not raw_input_data.strip():
                    validation_summary['validation_passed'] = False
                    validation_summary['warnings'].append("FASTA input must be a non-empty string")
            
            elif input_type == 'uniprot_identifiers':
                if not isinstance(raw_input_data, (str, list)):
                    validation_summary['validation_passed'] = False
                    validation_summary['warnings'].append("UniProt identifiers must be string or list")
            
            elif input_type == 'workspace_object':
                if not workspace_client:
                    validation_summary['warnings'].append("Workspace client not provided for workspace object validation")
            
            # Create validated input container
            validated_input = InputData(
                protein_records=[],  # Will be populated by DataExtractionStage
                input_type=input_type,
                input_source=str(raw_input_data)[:100] + "..." if len(str(raw_input_data)) > 100 else str(raw_input_data),
                validation_summary=validation_summary,
                metadata={
                    'validation_timestamp': time.time(),
                    'validation_config': validation_config
                }
            )
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=validation_summary['validation_passed'],
                output_data={
                    'validated_input': validated_input,
                    'validation_stats': {
                        'total_inputs': 1,
                        'valid_inputs': 1 if validation_summary['validation_passed'] else 0,
                        'invalid_inputs': 0 if validation_summary['validation_passed'] else 1,
                        'warnings': validation_summary['warnings']
                    }
                },
                metadata={
                    'input_type': input_type,
                    'validation_summary': validation_summary
                },
                execution_time=execution_time,
                warnings=validation_summary['warnings']
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Input validation failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class DataExtractionStage(BaseStage):
    """
    Extracts protein data from validated input sources.
    
    Handles:
    - FASTA parsing
    - UniProt ID resolution
    - Workspace object extraction
    - Data normalization
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.input_parser = input_parser.InputParser()
        self.max_proteins = config.get('max_proteins', 1000) if config else 1000
    
    def get_stage_name(self) -> str:
        return "data_extraction"
    
    def get_required_inputs(self) -> List[str]:
        return ['validated_input']
    
    def get_optional_inputs(self) -> List[str]:
        return ['workspace_client', 'extraction_config']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['input_validation']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'validated_input' not in input_data:
            logger.error("Missing validated_input")
            return False
        
        validated_input = input_data['validated_input']
        if not isinstance(validated_input, InputData):
            logger.error("validated_input must be InputData instance")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'protein_records': {
                'type': 'array',
                'items': {
                    'type': 'object',
                    'properties': {
                        'protein_id': {'type': 'string'},
                        'source': {'type': 'string'},
                        'sequence': {'type': 'string'},
                        'metadata': {'type': 'object'}
                    }
                }
            },
            'extraction_stats': {
                'type': 'object',
                'properties': {
                    'total_extracted': {'type': 'integer'},
                    'successful_extractions': {'type': 'integer'},
                    'failed_extractions': {'type': 'integer'},
                    'sequence_lengths': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Extract protein data from validated input."""
        start_time = time.time()
        
        try:
            validated_input = input_data['validated_input']
            extraction_config = input_data.get('extraction_config', {})
            
            # Set workspace client for parser
            if workspace_client:
                self.input_parser.workspace_client = workspace_client
            
            # Extract protein records
            protein_records = self.input_parser.parse_input(
                validated_input.input_type,
                validated_input.input_source
            )
            
            # Apply limits
            if len(protein_records) > self.max_proteins:
                logger.warning(f"Limiting extraction to {self.max_proteins} proteins from {len(protein_records)} found")
                protein_records = protein_records[:self.max_proteins]
            
            # Validate extracted records
            validation_results = self.input_parser.validate_records(protein_records)
            
            # Calculate statistics
            sequence_lengths = {}
            successful_extractions = 0
            failed_extractions = 0
            
            for record in protein_records:
                if record.sequence and len(record.sequence.strip()) > 0:
                    successful_extractions += 1
                    length = len(record.sequence)
                    sequence_lengths[record.protein_id] = length
                else:
                    failed_extractions += 1
            
            # Update validated input with extracted records
            validated_input.protein_records = protein_records
            validated_input.metadata.update({
                'extraction_timestamp': time.time(),
                'extraction_config': extraction_config,
                'validation_results': validation_results
            })
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_extractions > 0,
                output_data={
                    'protein_records': protein_records,
                    'validated_input': validated_input,
                    'extraction_stats': {
                        'total_extracted': len(protein_records),
                        'successful_extractions': successful_extractions,
                        'failed_extractions': failed_extractions,
                        'sequence_lengths': sequence_lengths,
                        'validation_results': validation_results
                    }
                },
                metadata={
                    'input_type': validated_input.input_type,
                    'extraction_summary': validation_results
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Data extraction failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class WorkspaceObjectStage(BaseStage):
    """
    Handles workspace object operations and data retrieval.
    
    Handles:
    - Workspace object validation
    - Data retrieval from workspace using DataFileUtil
    - Object type detection
    - Metadata extraction
    - File operations for protein sequence sets and genomes
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.supported_object_types = config.get('supported_object_types', [
            'KBaseSequences.ProteinSequenceSet',
            'KBaseGenomes.Genome',
            'KBaseSets.GenomeSet',
            'KBaseSets.FeatureSet',
            'KBaseGenomes.GenomeAnnotation',
            'KBaseSets.AssemblySet'
        ]) if config else [
            'KBaseSequences.ProteinSequenceSet',
            'KBaseGenomes.Genome',
            'KBaseSets.GenomeSet',
            'KBaseSets.FeatureSet',
            'KBaseGenomes.GenomeAnnotation',
            'KBaseSets.AssemblySet'
        ]
    
    def get_stage_name(self) -> str:
        return "workspace_object"
    
    def get_required_inputs(self) -> List[str]:
        return ['workspace_ref', 'workspace_client']
    
    def get_optional_inputs(self) -> List[str]:
        return ['object_config']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['input_validation']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate workspace object input."""
        if 'workspace_ref' not in input_data:
            logger.error("Missing workspace_ref")
            return False
        
        if 'workspace_client' not in input_data:
            logger.error("Missing workspace_client")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'workspace_object': {
                'type': 'object',
                'properties': {
                    'object_info': {'type': 'object'},
                    'object_data': {'type': 'object'},
                    'object_type': {'type': 'string'},
                    'metadata': {'type': 'object'}
                }
            },
            'object_stats': {
                'type': 'object',
                'properties': {
                    'object_size': {'type': 'integer'},
                    'object_type': {'type': 'string'},
                    'extraction_success': {'type': 'boolean'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None, datafileutil_client=None) -> StageResult:
        """Process workspace object using DataFileUtil for file operations."""
        start_time = time.time()
        
        try:
            workspace_ref = input_data['workspace_ref']
            ws_client = input_data.get('workspace_client', workspace_client)
            dfu_client = input_data.get('datafileutil_client', datafileutil_client)
            object_config = input_data.get('object_config', {})
            
            if not ws_client:
                raise ValueError("Workspace client is required")
            
            # Get object info
            object_info = ws_client.get_object_info3({
                'objects': [{'ref': workspace_ref}]
            })['infos'][0]
            
            object_type = object_info[2]
            
            # Validate object type
            if object_type not in self.supported_object_types:
                logger.warning(f"Object type {object_type} not in supported types: {self.supported_object_types}")
            
            # Get object data
            object_data = ws_client.get_objects2({
                'objects': [{'ref': workspace_ref}]
            })['data'][0]
            
            # Process object data based on type
            processed_data = self._process_object_by_type(
                object_data['data'], 
                object_type, 
                dfu_client,
                object_config
            )
            
            # Extract metadata
            metadata = {
                'object_name': object_info[1],
                'object_type': object_type,
                'object_version': object_info[4],
                'object_size': len(str(object_data['data'])),
                'extraction_timestamp': time.time(),
                'object_config': object_config,
                'processed_data': processed_data
            }
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=True,
                output_data={
                    'workspace_object': {
                        'object_info': object_info,
                        'object_data': object_data['data'],
                        'object_type': object_type,
                        'metadata': metadata,
                        'processed_data': processed_data
                    },
                    'object_stats': {
                        'object_size': metadata['object_size'],
                        'object_type': object_type,
                        'extraction_success': True,
                        'processing_summary': processed_data.get('summary', {})
                    }
                },
                metadata={
                    'workspace_ref': workspace_ref,
                    'object_type': object_type
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Workspace object processing failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _process_object_by_type(self, object_data: Dict[str, Any], object_type: str, 
                               dfu_client=None, config: Dict[str, Any] = None) -> Dict[str, Any]:
        """Process workspace object data based on its type."""
        config = config or {}
        
        if object_type == 'KBaseSequences.ProteinSequenceSet':
            return self._process_protein_sequence_set(object_data, dfu_client, config)
        elif object_type == 'KBaseGenomes.Genome':
            return self._process_genome(object_data, dfu_client, config)
        elif object_type == 'KBaseSets.GenomeSet':
            return self._process_genome_set(object_data, dfu_client, config)
        elif object_type == 'KBaseSets.FeatureSet':
            return self._process_feature_set(object_data, dfu_client, config)
        else:
            # Generic processing for unsupported types
            return {
                'type': 'generic',
                'data': object_data,
                'summary': {
                    'total_items': len(object_data) if isinstance(object_data, (list, dict)) else 1,
                    'processing_method': 'generic'
                }
            }
    
    def _process_protein_sequence_set(self, data: Dict[str, Any], dfu_client=None, 
                                    config: Dict[str, Any] = None) -> Dict[str, Any]:
        """Process ProteinSequenceSet object."""
        sequences = data.get('sequences', [])
        protein_records = []
        
        for seq in sequences:
            protein_records.append({
                'protein_id': seq.get('id', 'unknown'),
                'sequence': seq.get('sequence', ''),
                'source': 'ProteinSequenceSet',
                'metadata': {
                    'description': seq.get('description', ''),
                    'md5': seq.get('md5', ''),
                    'length': len(seq.get('sequence', ''))
                }
            })
        
        return {
            'type': 'protein_sequence_set',
            'protein_records': protein_records,
            'summary': {
                'total_proteins': len(protein_records),
                'processing_method': 'protein_sequence_set'
            }
        }
    
    def _process_genome(self, data: Dict[str, Any], dfu_client=None, 
                       config: Dict[str, Any] = None) -> Dict[str, Any]:
        """Process Genome object."""
        features = data.get('features', [])
        protein_records = []
        
        for feature in features:
            if feature.get('type') == 'CDS':  # Only process coding sequences
                protein_seq = feature.get('protein_translation', '')
                if protein_seq:
                    protein_records.append({
                        'protein_id': feature.get('id', 'unknown'),
                        'sequence': protein_seq,
                        'source': 'Genome',
                        'metadata': {
                            'feature_type': feature.get('type', ''),
                            'location': feature.get('location', []),
                            'function': feature.get('function', ''),
                            'length': len(protein_seq)
                        }
                    })
        
        return {
            'type': 'genome',
            'protein_records': protein_records,
            'summary': {
                'total_proteins': len(protein_records),
                'total_features': len(features),
                'processing_method': 'genome'
            }
        }
    
    def _process_genome_set(self, data: Dict[str, Any], dfu_client=None, 
                           config: Dict[str, Any] = None) -> Dict[str, Any]:
        """Process GenomeSet object."""
        genome_refs = data.get('genome_refs', [])
        all_protein_records = []
        
        # Note: This would require additional workspace calls to process each genome
        # For now, return summary information
        return {
            'type': 'genome_set',
            'genome_refs': genome_refs,
            'summary': {
                'total_genomes': len(genome_refs),
                'processing_method': 'genome_set',
                'note': 'Individual genomes need to be processed separately'
            }
        }
    
    def _process_feature_set(self, data: Dict[str, Any], dfu_client=None, 
                            config: Dict[str, Any] = None) -> Dict[str, Any]:
        """Process FeatureSet object."""
        features = data.get('elements', [])
        protein_records = []
        
        for feature in features:
            if feature.get('type') == 'CDS':
                protein_seq = feature.get('protein_translation', '')
                if protein_seq:
                    protein_records.append({
                        'protein_id': feature.get('id', 'unknown'),
                        'sequence': protein_seq,
                        'source': 'FeatureSet',
                        'metadata': {
                            'feature_type': feature.get('type', ''),
                            'function': feature.get('function', ''),
                            'length': len(protein_seq)
                        }
                    })
        
        return {
            'type': 'feature_set',
            'protein_records': protein_records,
            'summary': {
                'total_proteins': len(protein_records),
                'total_features': len(features),
                'processing_method': 'feature_set'
            }
        }
