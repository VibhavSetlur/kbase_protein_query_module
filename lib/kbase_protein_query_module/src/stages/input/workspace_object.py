"""
Workspace Object Stage for KBase Protein Query Module

This stage handles workspace object retrieval and processing.
"""

import logging
import time
from typing import Dict, Any, List, Optional, Union

from ..base_stage import BaseStage, StageResult
from ...utils import input_parser

logger = logging.getLogger(__name__)

class WorkspaceObjectStage(BaseStage):
    """
    Handles workspace object retrieval and processing.
    
    Handles:
    - Workspace object validation
    - Object data retrieval
    - Data format conversion
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.max_retries = config.get('max_retries', 3) if config else 3
        self.timeout = config.get('timeout', 30) if config else 30
    
    def get_stage_name(self) -> str:
        return "workspace_object"
    
    def get_required_inputs(self) -> List[str]:
        return ['validated_input']
    
    def get_optional_inputs(self) -> List[str]:
        return ['workspace_client', 'object_config']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input parameters."""
        if 'validated_input' not in input_data:
            logger.error("Missing validated_input")
            return False
        
        validated_input = input_data['validated_input']
        if validated_input.input_type != 'workspace_object':
            logger.error("Input type must be workspace_object")
            return False
        
        if not validated_input.input_source:
            logger.error("Missing workspace object reference")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'workspace_data': {
                'type': 'object',
                'properties': {
                    'object_info': {'type': 'object'},
                    'object_data': {'type': 'object'},
                    'extraction_stats': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Process workspace object."""
        start_time = time.time()
        
        try:
            validated_input = input_data['validated_input']
            object_config = input_data.get('object_config', {})
            
            if not workspace_client:
                raise ValueError("Workspace client required for workspace object processing")
            
            # Retrieve workspace object
            workspace_data = self._retrieve_workspace_object(
                validated_input.input_source, 
                workspace_client, 
                object_config
            )
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=True,
                output_data={
                    'workspace_data': workspace_data
                },
                metadata={
                    'object_ref': validated_input.input_source,
                    'processing_time': execution_time
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Workspace object processing failed: {str(e)}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _retrieve_workspace_object(self, object_ref: str, workspace_client, object_config: Dict[str, Any]) -> Dict[str, Any]:
        """Retrieve workspace object data."""
        for attempt in range(self.max_retries):
            try:
                # Get object info
                object_info = workspace_client.get_object_info3({
                    'objects': [{'ref': object_ref}]
                })
                
                # Get object data
                object_data = workspace_client.get_objects2({
                    'objects': [{'ref': object_ref}]
                })
                
                extraction_stats = {
                    'object_retrieved': True,
                    'object_type': object_info['infos'][0][2],
                    'object_name': object_info['infos'][0][1],
                    'object_size': len(str(object_data['data'][0])),
                    'retrieval_attempts': attempt + 1
                }
                
                return {
                    'object_info': object_info['infos'][0],
                    'object_data': object_data['data'][0],
                    'extraction_stats': extraction_stats
                }
                
            except Exception as e:
                if attempt == self.max_retries - 1:
                    raise ValueError(f"Failed to retrieve workspace object {object_ref}: {str(e)}")
                time.sleep(1)  # Brief delay before retry
        
        raise ValueError(f"Failed to retrieve workspace object {object_ref} after {self.max_retries} attempts")
