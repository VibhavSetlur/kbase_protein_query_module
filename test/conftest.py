"""
Test configuration and fixtures for KBase Protein Query Module

This module provides test fixtures and configuration for the protein query module tests.
"""

import pytest
import tempfile
import os
import sys
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


