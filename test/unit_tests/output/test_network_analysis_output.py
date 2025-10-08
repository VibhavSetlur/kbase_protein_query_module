"""
Unit tests for NetworkAnalysisOutput.

Tests network analysis output generation and visualization.
"""

import pytest
import sys
import os
import tempfile
import shutil
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.output.analysis.network_analysis.output import NetworkAnalysisOutput, AnalysisOutputResult


class TestNetworkAnalysisOutput:
    """Test cases for NetworkAnalysisOutput."""
    
    def test_initialization(self):
        """Test NetworkAnalysisOutput initialization."""
        output_handler = NetworkAnalysisOutput()
        
        assert output_handler.output_manager is None
        
        # Test with output manager
        mock_output_manager = Mock()
        output_handler = NetworkAnalysisOutput(mock_output_manager)
        
        assert output_handler.output_manager == mock_output_manager
    
    def test_generate_outputs_success(self, temp_dir):
        """Test successful output generation."""
        output_handler = NetworkAnalysisOutput()
        
        # Mock analysis data
        analysis_data = {
            'network_graph': {
                'nodes': ['P12345', 'P67890', 'Q13547'],
                'edges': [
                    {'source': 'P12345', 'target': 'P67890', 'weight': 0.85},
                    {'source': 'P67890', 'target': 'Q13547', 'weight': 0.72}
                ]
            },
            'node_centrality': {
                'P12345': {'degree': 2, 'betweenness': 0.1, 'closeness': 0.5},
                'P67890': {'degree': 2, 'betweenness': 0.2, 'closeness': 0.6},
                'Q13547': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.4}
            },
            'community_structure': {
                'communities': [['P12345', 'P67890'], ['Q13547']],
                'modularity': 0.25
            },
            'pathway_enrichment': [
                {'pathway': 'Pathway A', 'p_value': 0.01, 'enrichment_score': 2.5},
                {'pathway': 'Pathway B', 'p_value': 0.05, 'enrichment_score': 1.8}
            ]
        }
        
        # Mock the output manager
        mock_output_manager = Mock()
        mock_output_manager.write_json.return_value = os.path.join(temp_dir, 'network_analysis.json')
        mock_output_manager.write_html.return_value = os.path.join(temp_dir, 'network_analysis.html')
        
        output_handler.output_manager = mock_output_manager
        
        # Mock plotly availability
        with patch('kbase_protein_query_module.src.output.analysis.network_analysis.output.PLOTLY_AVAILABLE', True):
            result = output_handler.generate_outputs(analysis_data, 'network_analysis')
        
        assert isinstance(result, AnalysisOutputResult)
        assert result.success is True
        assert len(result.output_files) > 0
        assert 'nodes' in result.metadata
        assert 'edges' in result.metadata
        assert result.summary is not None
    
    def test_generate_outputs_without_plotly(self, temp_dir):
        """Test output generation without plotly."""
        output_handler = NetworkAnalysisOutput()
        
        analysis_data = {
            'network_graph': {
                'nodes': ['P12345', 'P67890'],
                'edges': [{'source': 'P12345', 'target': 'P67890', 'weight': 0.85}]
            },
            'node_centrality': {
                'P12345': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.5},
                'P67890': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.5}
            },
            'community_structure': {
                'communities': [['P12345', 'P67890']],
                'modularity': 0.0
            },
            'pathway_enrichment': []
        }
        
        # Mock the output manager
        mock_output_manager = Mock()
        mock_output_manager.write_json.return_value = os.path.join(temp_dir, 'network_analysis.json')
        
        output_handler.output_manager = mock_output_manager
        
        # Mock plotly not available
        with patch('kbase_protein_query_module.src.output.analysis.network_analysis.output.PLOTLY_AVAILABLE', False):
            result = output_handler.generate_outputs(analysis_data, 'network_analysis')
        
        assert isinstance(result, AnalysisOutputResult)
        assert result.success is True
        assert len(result.output_files) > 0
        # Should not include HTML visualization
        html_files = [f for f in result.output_files if f.endswith('.html')]
        assert len(html_files) == 0
    
    def test_generate_outputs_minimal_data(self, temp_dir):
        """Test output generation with minimal data."""
        output_handler = NetworkAnalysisOutput()
        
        # Minimal analysis data
        analysis_data = {
            'network_graph': {
                'nodes': ['P12345'],
                'edges': []
            },
            'node_centrality': {
                'P12345': {'degree': 0, 'betweenness': 0.0, 'closeness': 0.0}
            },
            'community_structure': {
                'communities': [['P12345']],
                'modularity': 0.0
            },
            'pathway_enrichment': []
        }
        
        # Mock the output manager
        mock_output_manager = Mock()
        mock_output_manager.write_json.return_value = os.path.join(temp_dir, 'network_analysis.json')
        
        output_handler.output_manager = mock_output_manager
        
        result = output_handler.generate_outputs(analysis_data, 'network_analysis')
        
        assert isinstance(result, AnalysisOutputResult)
        assert result.success is True
        assert result.metadata['nodes'] == 1
        assert result.metadata['edges'] == 0
        assert result.metadata['communities'] == 1
    
    def test_generate_outputs_error_handling(self, temp_dir):
        """Test error handling in output generation."""
        output_handler = NetworkAnalysisOutput()
        
        # Invalid analysis data
        analysis_data = {
            'network_graph': None,  # Invalid data
            'node_centrality': {},
            'community_structure': {},
            'pathway_enrichment': []
        }
        
        # Mock the output manager to raise exception
        mock_output_manager = Mock()
        mock_output_manager.write_json.side_effect = Exception("Write error")
        
        output_handler.output_manager = mock_output_manager
        
        result = output_handler.generate_outputs(analysis_data, 'network_analysis')
        
        assert isinstance(result, AnalysisOutputResult)
        assert result.success is False
        assert result.error_message is not None
        assert 'Write error' in result.error_message
    
    def test_generate_network_statistics(self, temp_dir):
        """Test network statistics generation."""
        output_handler = NetworkAnalysisOutput()
        
        network_graph = {
            'nodes': ['P12345', 'P67890', 'Q13547'],
            'edges': [
                {'source': 'P12345', 'target': 'P67890', 'weight': 0.85},
                {'source': 'P67890', 'target': 'Q13547', 'weight': 0.72}
            ]
        }
        
        node_centrality = {
            'P12345': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.33},
            'P67890': {'degree': 2, 'betweenness': 1.0, 'closeness': 1.0},
            'Q13547': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.33}
        }
        
        community_structure = {
            'communities': [['P12345', 'P67890'], ['Q13547']],
            'modularity': 0.25
        }
        
        stats = output_handler._generate_network_statistics(
            network_graph, node_centrality, community_structure
        )
        
        assert isinstance(stats, dict)
        assert 'network_size' in stats
        assert 'centrality_summary' in stats
        assert 'community_summary' in stats
        
        network_size = stats['network_size']
        assert network_size['nodes'] == 3
        assert network_size['edges'] == 2
        
        community_summary = stats['community_summary']
        assert community_summary['communities'] == 2
        assert community_summary['modularity'] == 0.25
    
    def test_get_output_description(self):
        """Test getting output description."""
        output_handler = NetworkAnalysisOutput()
        
        description = output_handler.get_output_description()
        
        assert isinstance(description, str)
        assert len(description) > 0
        assert 'network' in description.lower()
    
    def test_get_supported_formats(self):
        """Test getting supported output formats."""
        output_handler = NetworkAnalysisOutput()
        
        formats = output_handler.get_supported_formats()
        
        assert isinstance(formats, list)
        assert 'json' in formats
        assert 'html' in formats
    
    def test_validate_analysis_data_valid(self):
        """Test validation with valid analysis data."""
        output_handler = NetworkAnalysisOutput()
        
        valid_data = {
            'network_graph': {
                'nodes': ['P12345', 'P67890'],
                'edges': [{'source': 'P12345', 'target': 'P67890', 'weight': 0.85}]
            },
            'node_centrality': {
                'P12345': {'degree': 1},
                'P67890': {'degree': 1}
            },
            'community_structure': {
                'communities': [['P12345', 'P67890']],
                'modularity': 0.0
            }
        }
        
        assert output_handler.validate_analysis_data(valid_data) is True
    
    def test_validate_analysis_data_invalid(self):
        """Test validation with invalid analysis data."""
        output_handler = NetworkAnalysisOutput()
        
        invalid_data = {
            'node_centrality': {},
            'community_structure': {}
            # Missing network_graph key entirely
        }
        
        assert output_handler.validate_analysis_data(invalid_data) is False
        
        # Test with missing keys - validation only checks for network_graph presence
        incomplete_data = {
            'network_graph': {'nodes': [], 'edges': []}
            # Missing node_centrality and community_structure but has network_graph
        }
        
        # Since validation only checks for network_graph presence, this should be True
        assert output_handler.validate_analysis_data(incomplete_data) is True
    
    def test_generate_visualization_html(self, temp_dir):
        """Test HTML visualization generation."""
        output_handler = NetworkAnalysisOutput()
        
        network_graph = {
            'nodes': ['P12345', 'P67890'],
            'edges': [{'source': 'P12345', 'target': 'P67890', 'weight': 0.85}]
        }
        
        node_centrality = {
            'P12345': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.5},
            'P67890': {'degree': 1, 'betweenness': 0.0, 'closeness': 0.5}
        }
        
        community_structure = {
            'communities': [['P12345', 'P67890']],
            'modularity': 0.0
        }
        
        # Mock the output manager
        mock_output_manager = Mock()
        mock_output_manager.write_html.return_value = os.path.join(temp_dir, 'network_analysis.html')
        
        output_handler.output_manager = mock_output_manager
        
        # Mock plotly availability
        with patch('kbase_protein_query_module.src.output.analysis.network_analysis.output.PLOTLY_AVAILABLE', True):
            html_path = output_handler._generate_visualization_html(
                network_graph, node_centrality, community_structure, 'network_analysis'
            )
        
        assert html_path is not None
        # Note: html_path is mocked, so we don't check file existence
        assert html_path == os.path.join(temp_dir, 'network_analysis.html')
    
    def test_analysis_output_result_dataclass(self):
        """Test AnalysisOutputResult dataclass."""
        result = AnalysisOutputResult(
            success=True,
            output_files=['file1.json', 'file2.html'],
            metadata={'nodes': 3, 'edges': 2},
            summary='Analysis completed successfully',
            error_message=''
        )
        
        assert result.success is True
        assert len(result.output_files) == 2
        assert result.metadata['nodes'] == 3
        assert result.summary == 'Analysis completed successfully'
        assert result.error_message == ''
    
    def test_empty_network_handling(self, temp_dir):
        """Test handling of empty network."""
        output_handler = NetworkAnalysisOutput()
        
        empty_analysis_data = {
            'network_graph': {
                'nodes': [],
                'edges': []
            },
            'node_centrality': {},
            'community_structure': {
                'communities': [],
                'modularity': 0.0
            },
            'pathway_enrichment': []
        }
        
        # Mock the output manager
        mock_output_manager = Mock()
        mock_output_manager.write_json.return_value = os.path.join(temp_dir, 'empty_network.json')
        
        output_handler.output_manager = mock_output_manager
        
        result = output_handler.generate_outputs(empty_analysis_data, 'network_analysis')
        
        assert isinstance(result, AnalysisOutputResult)
        assert result.success is True
        assert result.metadata['nodes'] == 0
        assert result.metadata['edges'] == 0
        assert result.metadata['communities'] == 0
    
    def test_large_network_handling(self, temp_dir):
        """Test handling of large network."""
        output_handler = NetworkAnalysisOutput()
        
        # Create large network data
        nodes = [f'P{i:05d}' for i in range(1000)]
        edges = [{'source': f'P{i:05d}', 'target': f'P{(i+1)%1000:05d}', 'weight': 0.5} for i in range(1000)]
        
        large_analysis_data = {
            'network_graph': {
                'nodes': nodes,
                'edges': edges
            },
            'node_centrality': {
                node: {'degree': 2, 'betweenness': 0.001, 'closeness': 0.5}
                for node in nodes
            },
            'community_structure': {
                'communities': [nodes[i:i+100] for i in range(0, 1000, 100)],
                'modularity': 0.8
            },
            'pathway_enrichment': []
        }
        
        # Mock the output manager
        mock_output_manager = Mock()
        mock_output_manager.write_json.return_value = os.path.join(temp_dir, 'large_network.json')
        
        output_handler.output_manager = mock_output_manager
        
        result = output_handler.generate_outputs(large_analysis_data, 'network_analysis')
        
        assert isinstance(result, AnalysisOutputResult)
        assert result.success is True
        assert result.metadata['nodes'] == 1000
        assert result.metadata['edges'] == 1000
        assert result.metadata['communities'] == 10

