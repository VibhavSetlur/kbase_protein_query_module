"""
Input handling module for KBase Protein Query Module.
"""

from .input_manager import InputManager
from .protein_sequence import ProteinSequenceProcessor
from .uniprot_id import UniProtIdProcessor

__all__ = [
    'InputManager',
    'ProteinSequenceProcessor',
    'UniProtIdProcessor'
]
