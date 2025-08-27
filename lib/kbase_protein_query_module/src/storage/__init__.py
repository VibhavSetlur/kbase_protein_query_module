"""
Storage Module for Protein Data Management

This module provides comprehensive storage solutions for protein data including:
- Hierarchical storage for large-scale protein datasets
- Family assignment and management  
- Protein existence checking
- Metadata storage and retrieval
- Extensible indexing strategies for different use cases
- Performance-optimized storage backends
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
from .indexing_strategy import (
    IndexingStrategy, 
    IndexingConfig, 
    IndexingRegistry,
    FAISSIndexingStrategy,
    get_indexing_registry, 
    register_indexing_strategy
)

__all__ = [
    # Core storage components
    'ProteinStorage',
    'CompressedMetadataStorage', 
    'MemoryEfficientLoader',
    'ProteinIDsIndex',
    'ProteinFamilyAssigner',
    'ProteinExistenceChecker',
    'create_storage_structure',
    
    # Extensible indexing framework
    'IndexingStrategy',
    'IndexingConfig',
    'IndexingRegistry', 
    'FAISSIndexingStrategy',
    'get_indexing_registry',
    'register_indexing_strategy'
]
