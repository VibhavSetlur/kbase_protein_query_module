"""
Test configuration and fixtures for KBase Protein Query Module

This module provides test fixtures and configuration for the protein query module tests.
"""

import pytest
import tempfile
import os
import sys
import numpy as np
import torch
from unittest.mock import MagicMock
from pathlib import Path

# Ensure both module root and lib are on sys.path so tests can import 'lib.*' and 'kbase_protein_query_module.*'
module_root = str(Path(__file__).parent.parent)
lib_dir = str(Path(module_root) / 'lib')
if module_root not in sys.path:
    sys.path.insert(0, module_root)
if lib_dir not in sys.path:
    sys.path.insert(0, lib_dir)

# Import core modules for testing
from kbase_protein_query_module.src.core import BaseStage, StageResult, PipelineConfig
from kbase_protein_query_module.src.stages import (
    STAGE_REGISTRY, STAGE_DEPENDENCIES, get_stage_class, get_stage_dependencies
)
from kbase_protein_query_module.src.workflows import ProteinQueryWorkflow, WorkflowResult

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    import shutil
    shutil.rmtree(temp_dir)

@pytest.fixture
def sample_protein_data():
    """Sample protein data for testing."""
    return {
        'proteins': [
            {
                'protein_id': 'P12345',
                'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'
            },
            {
                'protein_id': 'P67890', 
                'sequence': 'MKTGFLVKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'
            }
        ]
    }

@pytest.fixture
def sample_config():
    """Sample pipeline configuration for testing."""
    return PipelineConfig(
        input_proteins=['P12345', 'P67890'],
        embedding_model='esm2_t6_8M_UR50D',
        similarity_threshold=0.8,
        max_similar_proteins=50
    )

@pytest.fixture
def sample_workflow(sample_config):
    """Sample workflow instance for testing."""
    return ProteinQueryWorkflow(sample_config)

@pytest.fixture(scope="session")
def mock_embedding_generator():
    """Create a mock embedding generator for testing purposes."""
    
    class MockESMModel:
        def __init__(self):
            self.config = type('Config', (), {'hidden_size': 320})()
        
        def __call__(self, **kwargs):
            # Return mock embeddings
            batch_size = kwargs.get('input_ids', torch.tensor([[0]])).shape[0]
            return type('Outputs', (), {
                'last_hidden_state': torch.randn(batch_size, 10, 320)
            })()
        
        def eval(self):
            pass
        
        def to(self, device):
            return self
    
    class MockTokenizer:
        def __init__(self):
            self.vocab_size = 33
        
        def __call__(self, text, **kwargs):
            # Return mock tokens
            return {
                'input_ids': torch.tensor([[1, 2, 3, 4, 5]]),
                'attention_mask': torch.tensor([[1, 1, 1, 1, 1]])
            }
    
    # Create mock generator
    mock_gen = MagicMock()
    mock_gen.model = MockESMModel()
    mock_gen.tokenizer = MockTokenizer()
    mock_gen.device = 'cpu'
    mock_gen.embedding_dim = 320
    mock_gen.generate_embedding = lambda seq, protein_id=None: np.random.randn(320).astype(np.float32)
    mock_gen.generate_embeddings_batch = lambda seqs, ids, batch_size=8: {id_: np.random.randn(320).astype(np.float32) for id_ in ids}
    
    return mock_gen

@pytest.fixture(scope="session")
def mock_pipeline_config():
    """Create a mock pipeline config for testing."""
    return PipelineConfig(
        input_type="FASTA",
        input_data="test.fasta",
        enabled_stages=["input_validation", "embedding_generation"],
        similarity_threshold=0.8,
        max_similar_proteins=10,
        embedding_model="esm2_t6_8M_UR50D",
        output_format="html"
    )


