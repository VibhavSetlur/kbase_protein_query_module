"""
Processing Module for KBase Protein Query Module

This module contains all processing logic for protein data:
- Embedding generation using ESM-2 models
- Similarity search and indexing
- Network construction and analysis
- Data transformation and processing pipelines
"""

from .embeddings.generator import ProteinEmbeddingGenerator, generate_embeddings_from_fasta
from .similarity.hierarchical_index import HierarchicalIndex, StreamingIndex, create_optimized_indexes
from .networks.builder import (
    DynamicNetworkBuilder, 
    create_localized_network,
    save_network,
    load_network,
    visualize_interactive_protein_network
)

__all__ = [
    # Embedding generation
    'ProteinEmbeddingGenerator',
    'generate_embeddings_from_fasta',
    
    # Similarity search
    'HierarchicalIndex',
    'StreamingIndex', 
    'create_optimized_indexes',
    
    # Network analysis
    'DynamicNetworkBuilder',
    'create_localized_network',
    'save_network',
    'load_network',
    'visualize_interactive_protein_network'
]
