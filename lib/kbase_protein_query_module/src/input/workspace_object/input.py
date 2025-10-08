"""
Workspace Object Input Handler

This module handles workspace object input processing for KBase narrative integration.
"""

import logging
import time
import os
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class WorkspaceData:
    """Container for workspace object data."""
    object_ref: str
    object_type: str
    protein_records: List[Dict[str, Any]]
    metadata: Dict[str, Any]

class WorkspaceObjectProcessor:
    """
    Handles workspace object retrieval and processing.
    
    Handles:
    - Workspace object validation
    - Object data retrieval
    - Data format conversion
    """
    
    def __init__(self, config: Dict[str, Any] = None, kb_util=None):
        self.config = config or {}
        self.kb_util = kb_util
        self.max_retries = config.get('max_retries', 3) if config else 3
        self.timeout = config.get('timeout', 30) if config else 30
    
    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process workspace object input.
        
        Args:
            input_data: Input data containing workspace object reference
            
        Returns:
            Processed data with proteins and workspace info
        """
        start_time = time.time()
        
        try:
            logger.info("Processing workspace object input")
            
            workspace_ref = input_data.get('workspace_object', '')
            if not workspace_ref:
                raise ValueError("workspace_object is required for workspace object input type")
            
            # Validate workspace reference format
            if not self._is_valid_workspace_ref(workspace_ref):
                raise ValueError(f"Invalid workspace reference format: {workspace_ref}")
            
            # Get workspace client
            workspace_client = self._get_workspace_client()
            
            # Retrieve object data
            object_data = self._retrieve_object_data(workspace_client, workspace_ref)
            
            # Parse object data into proteins
            protein_records = self._parse_workspace_object(object_data)
            
            processing_time = time.time() - start_time
            
            result = {
                'success': True,
                'input_type': 'workspace_object',
                'proteins': protein_records,
                'workspace_info': {
                    'object_ref': workspace_ref,
                    'object_type': self._extract_object_type(object_data),
                    'object_name': self._extract_object_name(object_data),
                    'workspace_id': self._extract_workspace_id(object_data),
                    'object_id': self._extract_object_id(object_data)
                },
                'processing_time': processing_time,
                'metadata': {
                    'total_proteins': len(protein_records),
                    'processing_time': processing_time,
                    'source': 'workspace_object',
                    'workspace_info': {
                        'object_ref': workspace_ref,
                        'object_type': self._extract_object_type(object_data),
                        'object_name': self._extract_object_name(object_data),
                        'workspace_id': self._extract_workspace_id(object_data),
                        'object_id': self._extract_object_id(object_data)
                    }
                }
            }
            
            logger.info(f"Processed workspace object with {len(protein_records)} proteins")
            return result
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"Workspace object processing failed: {e}")
            
            return {
                'success': False,
                'input_type': 'workspace_object',
                'proteins': [],
                'workspace_info': {},
                'error_message': str(e),
                'processing_time': processing_time,
                'metadata': {
                    'error': str(e),
                    'processing_time': processing_time
                }
            }
    
    def _is_valid_workspace_ref(self, ref: str) -> bool:
        """Check if reference is a valid workspace reference."""
        if not ref or '/' not in ref:
            return False
        
        parts = ref.split('/')
        if len(parts) != 3:
            return False
        
        try:
            # For numeric references, validate they're integers
            if parts[0].isdigit() and parts[1].isdigit():
                int(parts[0])  # Workspace ID should be numeric
                int(parts[1])  # Object ID should be numeric
            # For string references, just validate they're not empty
            elif parts[0] and parts[1]:
                pass  # String references are valid
            else:
                return False
            return True
        except ValueError:
            return False
    
    def _get_workspace_client(self):
        """Get workspace client from KBUtilLib or create new one."""
        if self.kb_util and hasattr(self.kb_util, 'get_client'):
            try:
                return self.kb_util.get_client()
            except Exception as e:
                logger.warning(f"Could not get workspace client from KBUtilLib: {e}")
        
        # Fallback to creating new workspace client
        callback_url = os.environ.get('SDK_CALLBACK_URL')
        if not callback_url:
            raise RuntimeError("SDK_CALLBACK_URL not found in environment. Cannot create Workspace client.")
        
        try:
            from installed_clients.WorkspaceClient import Workspace
            return Workspace(callback_url)
        except Exception as e:
            raise RuntimeError(f"Could not create workspace client: {e}")
    
    def _retrieve_object_data(self, workspace_client, workspace_ref: str) -> Dict[str, Any]:
        """Retrieve object data from workspace."""
        try:
            object_data = workspace_client.get_objects2({
                'objects': [{'ref': workspace_ref}]
            })
            
            if not object_data or 'data' not in object_data or not object_data['data']:
                raise ValueError(f"Object not found: {workspace_ref}")
            
            return object_data
            
        except Exception as e:
            raise ValueError(f"Failed to retrieve workspace object {workspace_ref}: {str(e)}")
    
    def _parse_workspace_object(self, object_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse workspace object data into protein records."""
        protein_records = []
        
        # Check if object_data contains multiple objects (array of objects)
        if 'data' in object_data and isinstance(object_data['data'], list):
            # Process multiple objects
            for obj in object_data['data']:
                protein_records.extend(self._parse_single_object(obj))
        else:
            # Process single object
            protein_records.extend(self._parse_single_object(object_data))
        
        return protein_records
    
    def _parse_single_object(self, object_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse a single workspace object into protein records."""
        protein_records = []
        
        # Handle different workspace object types
        info = object_data.get('info', [])
        if info and len(info) > 2:
            object_type = info[2]  # Type name is at index 2 in KBase info
        else:
            object_type = ''
        data = object_data.get('data', {}) or {}
        
        if 'ProteinSequenceSet' in object_type:
            # Parse protein sequence set
            sequences = data.get('sequences', [])
            
            for seq_data in sequences:
                protein_records.append({
                    'protein_id': seq_data.get('id', ''),
                    'sequence': seq_data.get('sequence', ''),
                    'source': 'workspace_object',
                    'metadata': {
                        'object_type': object_type,
                        'original_id': seq_data.get('id', ''),
                        'description': seq_data.get('description', '')
                    }
                })
        
        elif 'Genome' in object_type:
            # Parse genome object - handle both features and proteins
            if 'features' in data:
                features = data.get('features', [])
                
                for feature in features:
                    protein_id = feature.get('id', '')
                    sequence = feature.get('protein_translation', '')
                    
                    # Process if it's a CDS feature or has protein translation
                    if sequence and (feature.get('type') == 'CDS' or 'protein_translation' in feature):
                        protein_records.append({
                            'protein_id': protein_id,
                            'sequence': sequence,
                            'source': 'workspace_object',
                            'metadata': {
                                'object_type': object_type,
                                'feature_type': feature.get('type', 'unknown'),
                                'feature_id': protein_id,
                                'location': feature.get('location', [])
                            }
                        })
            elif 'proteins' in data:
                # Handle direct proteins structure
                proteins = data.get('proteins', [])
                
                for protein in proteins:
                    protein_id = protein.get('id', '')
                    sequence = protein.get('sequence', '')
                    
                    if sequence:
                        protein_records.append({
                            'protein_id': protein_id,
                            'sequence': sequence,
                            'source': 'workspace_object',
                            'metadata': {
                                'object_type': object_type,
                                'original_id': protein_id
                            }
                        })
        
        elif 'GenomeSet' in object_type:
            # Parse genome set
            genome_refs = data.get('elements', {}).get('genomes', [])
            
            for genome_ref in genome_refs:
                try:
                    workspace_client = self._get_workspace_client()
                    genome_data = workspace_client.get_objects2({
                        'objects': [{'ref': genome_ref}]
                    })
                    
                    if genome_data and 'data' in genome_data and genome_data['data']:
                        genome_obj = genome_data['data'][0]
                        genome_proteins = self._parse_workspace_object(genome_obj)
                        
                        # Add genome context to metadata
                        for protein in genome_proteins:
                            protein['metadata']['genome_ref'] = genome_ref
                            protein['metadata']['genome_set'] = True
                        
                        protein_records.extend(genome_proteins)
                        
                except Exception as e:
                    logger.warning(f"Failed to parse genome {genome_ref} from genome set: {e}")
                    continue
        
        else:
            logger.warning(f"Unknown workspace object type: {object_type}")
        
        return protein_records
    
    def _extract_object_type(self, object_data: Dict[str, Any]) -> str:
        """Extract object type from object data."""
        if 'data' in object_data and isinstance(object_data['data'], list):
            # Multiple objects - return type of first object
            first_obj = object_data['data'][0] if object_data['data'] else {}
            info = first_obj.get('info', [])
        else:
            # Single object
            info = object_data.get('info', [])
        
        return info[2] if len(info) > 2 else ''
    
    def _extract_object_name(self, object_data: Dict[str, Any]) -> str:
        """Extract object name from object data."""
        if 'data' in object_data and isinstance(object_data['data'], list):
            # Multiple objects - return name of first object
            first_obj = object_data['data'][0] if object_data['data'] else {}
            info = first_obj.get('info', [])
        else:
            # Single object
            info = object_data.get('info', [])
        
        return info[1] if len(info) > 1 else ''
    
    def _extract_workspace_id(self, object_data: Dict[str, Any]) -> str:
        """Extract workspace ID from object data."""
        if 'data' in object_data and isinstance(object_data['data'], list):
            # Multiple objects - return workspace ID of first object
            first_obj = object_data['data'][0] if object_data['data'] else {}
            info = first_obj.get('info', [])
        else:
            # Single object
            info = object_data.get('info', [])
        
        return info[6] if len(info) > 6 else ''
    
    def _extract_object_id(self, object_data: Dict[str, Any]) -> str:
        """Extract object ID from object data."""
        if 'data' in object_data and isinstance(object_data['data'], list):
            # Multiple objects - return ID of first object
            first_obj = object_data['data'][0] if object_data['data'] else {}
            info = first_obj.get('info', [])
        else:
            # Single object
            info = object_data.get('info', [])
        
        return info[0] if len(info) > 0 else ''