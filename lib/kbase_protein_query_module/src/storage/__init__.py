"""
Storage Module for Protein Data Management

This module provides comprehensive storage solutions for protein data including:
- Hierarchical storage for large-scale protein datasets
- Family assignment and management
- Protein existence checking
- Metadata storage and retrieval
"""

from .protein_storage import (
    ProteinStorage, 
    CompressedMetadataStorage, 
    MemoryEfficientLoader, 
    ProteinIDsIndex,
    create_storage_structure
)
from .protein_family_assigner import ProteinFamilyAssigner
from .protein_existence_checker import ProteinExistenceChecker

__all__ = [
    'ProteinStorage',
    'CompressedMetadataStorage', 
    'MemoryEfficientLoader',
    'ProteinIDsIndex',
    'ProteinFamilyAssigner',
    'ProteinExistenceChecker',
    'create_storage_structure'
]
