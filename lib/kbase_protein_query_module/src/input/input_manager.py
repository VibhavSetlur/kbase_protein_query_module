"""
Input Manager for KBase Protein Query Module

Centralized input handling for protein sequences and UniProt IDs.
"""

import logging
import time
import sys
import os
from typing import Dict, Any, List, Optional

# Handle both script execution and module import
if __name__ == "__main__" or __package__ is None:
    # Add parent directories to path for script execution
    current_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(current_dir)
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
    from input.protein_sequence.input import ProteinSequenceProcessor
    from input.uniprot_id.input import UniProtIdProcessor
else:
    from .protein_sequence.input import ProteinSequenceProcessor
    from .uniprot_id.input import UniProtIdProcessor

logger = logging.getLogger(__name__)

class InputManager:
    """Manages all input processing for the protein query module."""
    
    def __init__(self, config: Dict[str, Any], kb_util=None):
        """Initialize the Input Manager."""
        self.config = config or {}
        self.kb_util = kb_util
        
        logger.info("InputManager initialized")
    
    def normalize_input_params(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Normalize input parameters from test/narrative format to internal format.
        
        Handles:
        - uniprot_ids -> uniprot_id
        - protein_input -> protein_sequence
        """
        normalized = input_data.copy()
        
        # Map uniprot_ids to uniprot_id
        if 'uniprot_ids' in normalized and 'uniprot_id' not in normalized:
            normalized['uniprot_id'] = normalized.pop('uniprot_ids')
            if normalized.get('input_type') == 'uniprot_ids':
                normalized['input_type'] = 'uniprot_id'
        
        # Map protein_input to protein_sequence
        if 'protein_input' in normalized and 'protein_sequence' not in normalized:
            protein_input = normalized.pop('protein_input')
            # Handle list or string
            if isinstance(protein_input, list):
                # Join multiple sequences or take first
                normalized['protein_sequence'] = '\n'.join(protein_input) if len(protein_input) > 1 else protein_input[0]
            else:
                normalized['protein_sequence'] = protein_input
            if normalized.get('input_type') == 'protein_input':
                normalized['input_type'] = 'protein_sequence'
        
        return normalized
    
    def process_input(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """Process input data through the complete input pipeline."""
        start_time = time.time()
        
        try:
            logger.info("Starting input processing")
            
            # Normalize input parameters first
            normalized_input = self.normalize_input_params(input_data)
            input_type = normalized_input.get('input_type', 'protein_sequence')
            
            # Route to appropriate processor
            if input_type == 'protein_sequence':
                processor = ProteinSequenceProcessor(self.config)
                result = processor.process(normalized_input)
                standardized = processor.standardize_input(normalized_input)
            elif input_type == 'uniprot_id':
                processor = UniProtIdProcessor(self.config)
                result = processor.process(normalized_input)
                standardized = processor.standardize_input(normalized_input)
            else:
                raise ValueError(f"Unsupported input type: {input_type}")
            
            # Merge processing result with standardized input
            output = {
                **standardized,
                **result
            }
            
            processing_time = time.time() - start_time
            output['processing_time'] = processing_time
            
            logger.info(f"Input processing completed in {processing_time:.2f} seconds")
            logger.info(f"Processed {len(output.get('proteins', []))} proteins")
            
            return output
            
        except Exception as e:
            processing_time = time.time() - start_time
            logger.error(f"Input processing failed: {e}")
            
            return {
                'success': False,
                'input_type': input_data.get('input_type', 'unknown'),
                'protein_sequence': '',
                'uniprot_id': [],
                'analysis_name': '',
                'proteins': [],
                'processing_time': processing_time,
                'error_message': str(e)
            }


def main() -> int:
    """
    Test InputManager.
    
    Dependencies:
    - time
    - logging
    - typing
    - protein_sequence.input
    - uniprot_id.input
    """
    ok = True
    try:
        mgr = InputManager(config={})
        input_data = {
            'input_type': 'protein_sequence',
            'protein_sequence': 'ACDEFGHIKLMNPQRSTVWY'
        }
        result = mgr.process_input(input_data)
        if not isinstance(result, dict) or result.get('processing_time') is None:
            raise RuntimeError('Unexpected result structure from process_input')
        print("InputManager test: SUCCESS")
    except Exception as e:
        ok = False
        print(f"InputManager test: FAILED - {e}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
