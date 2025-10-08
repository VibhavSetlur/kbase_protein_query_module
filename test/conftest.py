import sys
import types

# Stub visualization module paths referenced by src.util.__init__
_vis_stub = types.ModuleType('visualization')
sys.modules['kbase_protein_query_module.src.util.visualization'] = _vis_stub
sys.modules['lib.kbase_protein_query_module.src.util.visualization'] = _vis_stub

"""
Test configuration and fixtures for KBase Protein Query Module

This module provides test fixtures and configuration for the protein query module tests.
"""

import pytest
import tempfile
import os
import sys
import numpy as np
from unittest.mock import MagicMock
from pathlib import Path

# Ensure both module root and lib are on sys.path so tests can import 'lib.*' and 'kbase_protein_query_module.*'
module_root = str(Path(__file__).parent.parent)
lib_dir = str(Path(module_root) / 'lib')
if module_root not in sys.path:
    sys.path.insert(0, module_root)
if lib_dir not in sys.path:
    sys.path.insert(0, lib_dir)

# Import core modules for testing (with error handling)
try:
    from kbase_protein_query_module.src.core import BaseStage, StageResult, PipelineConfig
    
    from kbase_protein_query_module.src.workflows import ProteinQueryWorkflow, WorkflowResult
    CORE_MODULES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import core modules: {e}")
    CORE_MODULES_AVAILABLE = False
    # Create mock classes for basic testing
    class BaseStage:
        pass
    class StageResult:
        pass
    class PipelineConfig:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)
    class ProteinQueryWorkflow:
        def __init__(self, config):
            self.config = config
    class WorkflowResult:
        pass

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
            # Return mock embeddings without torch dependency
            batch_size = 1  # Default batch size
            if 'input_ids' in kwargs:
                # Try to get batch size from input_ids if it's a numpy array or list
                input_ids = kwargs['input_ids']
                if hasattr(input_ids, 'shape'):
                    batch_size = input_ids.shape[0]
                elif isinstance(input_ids, (list, tuple)):
                    batch_size = len(input_ids)
            
            return type('Outputs', (), {
                'last_hidden_state': np.random.randn(batch_size, 10, 320)
            })()
        
        def eval(self):
            pass
        
        def to(self, device):
            return self
    
    class MockTokenizer:
        def __init__(self):
            self.vocab_size = 33
        
        def __call__(self, text, **kwargs):
            # Return mock tokens without torch dependency
            return {
                'input_ids': np.array([[1, 2, 3, 4, 5]]),
                'attention_mask': np.array([[1, 1, 1, 1, 1]])
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


