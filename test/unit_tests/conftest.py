"""
Pytest configuration and shared fixtures for unit tests.

This module provides common test fixtures, mock objects, and test utilities
that are shared across all unit tests in the KBase Protein Query Module.
"""

import os
import sys
import tempfile
import shutil
import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, MagicMock, patch
from typing import Dict, Any, List
import logging

# Add the lib directory to the Python path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../lib'))

# Test configuration
TEST_CONFIG = {
    'enabled_input_types': ['protein_input', 'uniprot_ids', 'workspace_object'],
    'enabled_analyses': ['network_analysis'],
    'output_dir': '/tmp/test_outputs',
    'workspace_name': 'test_workspace',
    'storage_config': {
        'use_compression': True,
        'chunk_size': 100,
        'max_family_size': 1000
    },
    'similarity_config': {
        'similarity_threshold': 0.7,
        'top_k_matches': 10,
        'use_family_assignment': True
    }
}

@pytest.fixture
def test_config():
    """Provide test configuration."""
    import copy
    return copy.deepcopy(TEST_CONFIG)

@pytest.fixture
def temp_dir():
    """Provide a temporary directory for test files."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)

@pytest.fixture
def mock_workspace_client():
    """Mock KBase workspace client."""
    mock_client = Mock()
    mock_client.get_workspace_info.return_value = [1, 'test_workspace', 'test_user', '2024-01-01T00:00:00+0000', 1, 'n', 'n', 'test_user']
    mock_client.get_objects2.return_value = {
        'data': [{
            'data': {
                'features': [
                    {
                        'id': 'P12345',
                        'type': 'CDS',
                        'protein_translation': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
                        'location': [[0, 1, '+']]
                    }
                ]
            },
            'info': [1, 'test_object', 'KBaseGenomes.Genome-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]
        }]
    }
    mock_client.save_objects.return_value = [[1, 'test_workspace', 'test_object', 'KBaseReport.Report-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]]
    return mock_client

@pytest.fixture
def mock_kb_util():
    """Mock KBase utility library."""
    mock_util = Mock()
    mock_util.log.return_value = None
    return mock_util

@pytest.fixture
def sample_protein_sequences():
    """Sample protein sequences for testing."""
    return [
        'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
        'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
        'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
    ]

@pytest.fixture
def sample_uniprot_ids():
    """Sample UniProt IDs for testing."""
    return ['P12345', 'P67890', 'Q13547']

@pytest.fixture
def sample_embeddings():
    """Sample protein embeddings for testing."""
    return np.random.rand(3, 320).astype(np.float32)

@pytest.fixture
def sample_protein_ids():
    """Sample protein IDs for testing."""
    return ['P12345', 'P67890', 'Q13547']

@pytest.fixture
def sample_metadata_df():
    """Sample protein metadata DataFrame for testing."""
    return pd.DataFrame({
        'Entry': ['P12345', 'P67890', 'Q13547'],
        'Protein names': ['Protein 1', 'Protein 2', 'Protein 3'],
        'Organism': ['Homo sapiens', 'Mus musculus', 'Escherichia coli'],
        'EC number': ['1.1.1.1', '2.3.4.5', '3.6.7.8'],
        'Function [CC]': ['Function 1', 'Function 2', 'Function 3'],
        'Protein families': ['Family A', 'Family B', 'Family C'],
        'Reviewed': ['Swiss-Prot', 'Swiss-Prot', 'TrEMBL']
    })

@pytest.fixture
def mock_embedding_generator():
    """Mock protein embedding generator."""
    mock_gen = Mock()
    mock_gen.generate_embedding.return_value = np.random.rand(320).astype(np.float32)
    mock_gen.generate_embeddings_batch.return_value = {
        'P12345': np.random.rand(320).astype(np.float32),
        'P67890': np.random.rand(320).astype(np.float32)
    }
    mock_gen.get_embedding_dimension.return_value = 320
    mock_gen.get_model_info.return_value = {'model': 'esm2_t6_8M_UR50D', 'dimension': 320}
    return mock_gen

@pytest.fixture
def mock_similarity_search_results():
    """Mock similarity search results."""
    return [
        {'protein_id': 'P12345', 'similarity': 0.95, 'family': 'Family A'},
        {'protein_id': 'P67890', 'similarity': 0.87, 'family': 'Family B'},
        {'protein_id': 'Q13547', 'similarity': 0.78, 'family': 'Family C'}
    ]

@pytest.fixture
def mock_network_analysis_results():
    """Mock network analysis results."""
    return {
        'network_graph': {
            'nodes': ['P12345', 'P67890', 'Q13547'],
            'edges': [{'source': 'P12345', 'target': 'P67890', 'weight': 0.85}]
        },
        'node_centrality': {
            'P12345': {'degree': 2, 'betweenness': 0.1, 'closeness': 0.5},
            'P67890': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.33},
            'Q13547': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.33}
        },
        'community_structure': {
            'communities': [['P12345', 'P67890'], ['Q13547']],
            'modularity': 0.25
        }
    }

@pytest.fixture
def mock_family_assignments():
    """Mock protein family assignments."""
    return {
        'P12345': {'family_id': 'Family_A', 'confidence': 0.95, 'similarity': 0.92},
        'P67890': {'family_id': 'Family_B', 'confidence': 0.87, 'similarity': 0.85},
        'Q13547': {'family_id': 'Family_C', 'confidence': 0.78, 'similarity': 0.75}
    }

@pytest.fixture
def sample_workspace_object_ref():
    """Sample workspace object reference."""
    return '1/2/3'

@pytest.fixture
def mock_context():
    """Mock execution context for KBase methods."""
    return {
        'token': 'test_token',
        'provenance': [{'ws_name': 'test_workspace'}]
    }

@pytest.fixture(autouse=True)
def setup_logging():
    """Setup logging for tests."""
    logging.basicConfig(level=logging.WARNING)
    yield
    # Cleanup after test

@pytest.fixture
def disable_external_calls():
    """Disable external API calls during testing."""
    with patch('requests.get'), \
         patch('requests.post'), \
         patch('urllib.request.urlopen'):
        yield

class TestDataGenerator:
    """Helper class for generating test data."""
    
    @staticmethod
    def create_protein_data(count: int = 5) -> List[Dict[str, Any]]:
        """Create sample protein data."""
        return [
            {
                'protein_id': f'P{10000 + i}',
                'sequence': f'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG{i}',
                'source': 'test',
                'metadata': {'organism': 'Test organism', 'function': 'Test function'}
            }
            for i in range(count)
        ]
    
    @staticmethod
    def create_embeddings(count: int = 5, dimension: int = 320) -> np.ndarray:
        """Create sample embeddings."""
        return np.random.rand(count, dimension).astype(np.float32)
    
    @staticmethod
    def create_similarity_matrix(size: int = 5) -> np.ndarray:
        """Create sample similarity matrix."""
        matrix = np.random.rand(size, size)
        # Make it symmetric
        matrix = (matrix + matrix.T) / 2
        # Set diagonal to 1
        np.fill_diagonal(matrix, 1.0)
        return matrix

@pytest.fixture
def test_data_generator():
    """Provide test data generator."""
    return TestDataGenerator()

