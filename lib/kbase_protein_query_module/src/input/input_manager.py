"""
Input Manager for KBase Protein Query Module

Centralized input handling for protein sequences and UniProt IDs.
"""

import logging
import time
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

from .workspace_object.input import WorkspaceObjectProcessor
from .protein_sequence.input import ProteinSequenceProcessor
from .uniprot_ids.input import UniProtIdsProcessor

logger = logging.getLogger(__name__)

@dataclass
class InputResult:
    """Result container for input processing."""
    success: bool
    proteins: List[Dict[str, Any]]
    workspace_info: Dict[str, Any]
    input_type: str
    processing_time: float
    error_message: Optional[str] = None
    metadata: Dict[str, Any] = None

class InputManager:
    """Manages all input processing for the protein query module."""
    
    def __init__(self, config: Dict[str, Any], kb_util=None):
        """Initialize the Input Manager."""
        self.config = config
        self.kb_util = kb_util
        
        # Initialize input processors
        self.workspace_processor = WorkspaceObjectProcessor(config, kb_util)
        self.protein_sequence_processor = ProteinSequenceProcessor(config)
        self.uniprot_processor = UniProtIdsProcessor(config)
        
        # Configure enabled input types
        self.enabled_input_types = config.get('enabled_input_types', [
            'protein_input', 'uniprot_ids'
        ])
        
        # Track last processed data for access by other components
        self.last_processed_data = {}
        
        logger.info("InputManager initialized")
    
    def process_input(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Process input data through the complete input pipeline."""
        start_time = time.time()
        
        try:
            logger.info("Starting input processing")
            
            input_type = input_data.get('input_type', 'protein_input')
            
            # Validate input type
            if input_type not in self.enabled_input_types:
                raise ValueError(f"Input type '{input_type}' is not enabled. Available types: {self.enabled_input_types}")
            
            # Route to appropriate processor
            if input_type == 'protein_input':
                result = ProteinSequenceProcessor(self.config).process(input_data)
            elif input_type == 'uniprot_ids':
                result = UniProtIdsProcessor(self.config).process(input_data)
            else:
                raise ValueError(f"Unsupported input type: {input_type}")
            
            # Add processing metadata
            processing_time = result.get('processing_time', time.time() - start_time) if isinstance(result, dict) else time.time() - start_time
            result['processing_time'] = processing_time
            result['input_type'] = input_type
            
            # Store processed data for access by other components
            self.last_processed_data = result.copy()
            
            logger.info(f"Input processing completed in {processing_time:.2f} seconds")
            logger.info(f"Processed {len(result.get('proteins', []))} proteins")
            
            return result
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"Input processing failed: {e}")
            
            return {
                'success': False,
                'proteins': [],
                'workspace_info': {},
                'input_type': input_data.get('input_type', 'unknown'),
                'processing_time': processing_time,
                'error_message': str(e),
                'metadata': {}
            }
    
    
    def get_supported_input_types(self) -> List[str]:
        """Get list of supported input types."""
        return self.enabled_input_types.copy()
    
    def is_input_type_enabled(self, input_type: str) -> bool:
        """Check if an input type is enabled."""
        return input_type in self.enabled_input_types
    
    def enable_input_type(self, input_type: str):
        """Enable an input type."""
        if input_type not in self.enabled_input_types:
            self.enabled_input_types.append(input_type)
            logger.info(f"Enabled input type: {input_type}")
    
    def disable_input_type(self, input_type: str):
        """Disable an input type."""
        if input_type in self.enabled_input_types:
            self.enabled_input_types.remove(input_type)
            logger.info(f"Disabled input type: {input_type}")
    
    def validate_input_data(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data without processing it."""
        try:
            input_type = input_data.get('input_type')
            if not input_type:
                logger.error("input_type is required")
                return False
            
            if input_type not in self.enabled_input_types:
                logger.error(f"Input type '{input_type}' is not enabled")
                return False
            
            # Check required fields based on input type
            if input_type == 'protein_input' and not input_data.get('protein_input'):
                logger.error("protein_input is required for protein_input type")
                return False
            elif input_type == 'uniprot_ids' and not input_data.get('uniprot_ids'):
                logger.error("uniprot_ids is required for uniprot_ids type")
                return False
            
            return True
            
        except Exception as e:
            logger.error(f"Input validation failed: {e}")
            return False
