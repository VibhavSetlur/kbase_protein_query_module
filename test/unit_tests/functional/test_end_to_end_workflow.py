"""
Functional tests for end-to-end workflow execution.

Tests complete workflows from input through analysis to output.
"""

import pytest
import sys
import os
import numpy as np
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.core.workflow_orchestrator import WorkflowOrchestrator, PipelineConfig
from kbase_protein_query_module.src.input.input_manager import InputManager
from kbase_protein_query_module.src.analysis.analysis_manager import AnalysisManager
from kbase_protein_query_module.src.output.output_manager import OutputManager


class TestEndToEndWorkflow:
    """Test cases for end-to-end workflow execution."""
    
    def test_complete_protein_input_workflow(self, test_config, mock_kb_util, temp_dir, sample_protein_sequences):
        """Test complete workflow with protein input."""
        # Setup workflow orchestrator
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        
        # Initialize components
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        # Prepare input data
        input_data = {
            'input_type': 'protein_input',
            'protein_input': sample_protein_sequences,
            'workspace_name': 'test_workspace',
            'analysis_name': 'test_analysis'
        }
        
        # Mock all the components
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': 3, 'edges': 2},
                            'node_centrality': {},
                            'community_structure': {'communities': [['protein_0', 'protein_1'], ['protein_2']]}
                        }
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['network_analysis.json', 'network_analysis.html'],
                        'summary': 'Analysis completed successfully'
                    }
                    
                    # Run workflow
                    result = orchestrator.run_workflow(
                        input_data,
                        temp_dir,
                        'test_workspace',
                        selected_analyses=['network_analysis']
                    )
                    
                    # Verify result
                    assert result.success is True
                    assert 'network_analysis' in result.analyses_completed
                    assert len(result.analysis_results) > 0
                    assert result.execution_time > 0
                    assert result.output_directory == temp_dir
    
    def test_complete_uniprot_workflow(self, test_config, mock_kb_util, temp_dir, sample_uniprot_ids):
        """Test complete workflow with UniProt input."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        input_data = {
            'input_type': 'uniprot_ids',
            'uniprot_ids': sample_uniprot_ids,
            'workspace_name': 'test_workspace',
            'analysis_name': 'uniprot_analysis'
        }
        
        # Mock components
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': uniprot_id, 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG', 'source': 'uniprot'}
                    for uniprot_id in sample_uniprot_ids
                ],
                'input_type': 'uniprot_ids',
                'processing_time': 0.2
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': len(sample_uniprot_ids), 'edges': 1},
                            'node_centrality': {},
                            'community_structure': {'communities': [[uniprot_id] for uniprot_id in sample_uniprot_ids]}
                        }
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['network_analysis.json'],
                        'summary': 'UniProt analysis completed'
                    }
                    
                    result = orchestrator.run_workflow(
                        input_data,
                        temp_dir,
                        'test_workspace',
                        selected_analyses=['network_analysis']
                    )
                    
                    assert result.success is True
                    assert 'network_analysis' in result.analyses_completed
    
    def test_complete_workspace_workflow(self, test_config, mock_kb_util, temp_dir, mock_workspace_client, sample_workspace_object_ref):
        """Test complete workflow with workspace object input."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        input_data = {
            'input_type': 'workspace_object',
            'workspace_object': sample_workspace_object_ref,
            'workspace_name': 'test_workspace',
            'analysis_name': 'workspace_analysis'
        }
        
        # Mock components
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': 'P12345', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG', 'source': 'workspace_object'}
                ],
                'input_type': 'workspace_object',
                'processing_time': 0.3
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': 1, 'edges': 0},
                            'node_centrality': {},
                            'community_structure': {'communities': [['P12345']]}
                        }
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['network_analysis.json'],
                        'summary': 'Workspace analysis completed'
                    }
                    
                    result = orchestrator.run_workflow(
                        input_data,
                        temp_dir,
                        'test_workspace',
                        selected_analyses=['network_analysis']
                    )
                    
                    assert result.success is True
                    assert 'network_analysis' in result.analyses_completed
    
    def test_workflow_with_multiple_analyses(self, test_config, mock_kb_util, temp_dir, sample_protein_sequences):
        """Test workflow with multiple analyses."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': sample_protein_sequences,
            'workspace_name': 'test_workspace',
            'analysis_name': 'multi_analysis'
        }
        
        # Mock components
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': 3, 'edges': 2},
                            'node_centrality': {},
                            'community_structure': {'communities': [['protein_0', 'protein_1'], ['protein_2']]}
                        }
                    },
                    'sequence_analysis': {
                        'success': True,
                        'analysis_type': 'sequence_analysis',
                        'results': {
                            'similarity_matrix': np.random.rand(3, 3).tolist(),
                            'pairwise_similarities': []
                        }
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['analysis_output.json'],
                        'summary': 'Multiple analyses completed'
                    }
                    
                    result = orchestrator.run_workflow(
                        input_data,
                        temp_dir,
                        'test_workspace',
                        selected_analyses=['network_analysis', 'sequence_analysis']
                    )
                    
                    assert result.success is True
                    assert len(result.analyses_completed) == 2
                    assert 'network_analysis' in result.analyses_completed
                    assert 'sequence_analysis' in result.analyses_completed
    
    def test_workflow_error_handling(self, test_config, mock_kb_util, temp_dir):
        """Test workflow error handling and recovery."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'],
            'workspace_name': 'test_workspace',
            'analysis_name': 'error_test'
        }
        
        # Mock input processing to fail
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': False,
                'error_message': 'Input processing failed'
            }
            
            result = orchestrator.run_workflow(
                input_data,
                temp_dir,
                'test_workspace'
            )
            
            assert result.success is False
            assert result.error_message is not None
            assert 'Input processing failed' in result.error_message
    
    def test_workflow_partial_analysis_failure(self, test_config, mock_kb_util, temp_dir, sample_protein_sequences):
        """Test workflow with partial analysis failure."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': sample_protein_sequences,
            'workspace_name': 'test_workspace',
            'analysis_name': 'partial_failure_test'
        }
        
        # Mock components
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                # One analysis succeeds, one fails
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': 3, 'edges': 2},
                            'node_centrality': {},
                            'community_structure': {'communities': [['protein_0', 'protein_1'], ['protein_2']]}
                        }
                    },
                    'sequence_analysis': {
                        'success': False,
                        'error': 'Sequence analysis failed'
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['network_analysis.json'],
                        'summary': 'Partial analysis completed'
                    }
                    
                    result = orchestrator.run_workflow(
                        input_data,
                        temp_dir,
                        'test_workspace',
                        selected_analyses=['network_analysis', 'sequence_analysis']
                    )
                    
                    # Should still succeed overall, but with partial results
                    assert result.success is True
                    assert len(result.analyses_completed) == 1
                    assert 'network_analysis' in result.analyses_completed
                    assert 'sequence_analysis' not in result.analyses_completed
    
    def test_workflow_output_generation(self, test_config, mock_kb_util, temp_dir, sample_protein_sequences):
        """Test workflow output generation."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': sample_protein_sequences,
            'workspace_name': 'test_workspace',
            'analysis_name': 'output_test'
        }
        
        # Mock components
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': 3, 'edges': 2},
                            'node_centrality': {},
                            'community_structure': {'communities': [['protein_0', 'protein_1'], ['protein_2']]}
                        }
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['network_analysis.json', 'network_analysis.html', 'summary.txt'],
                        'summary': 'Output generation completed'
                    }
                    
                    with patch.object(orchestrator.output_manager, 'save_metadata') as mock_metadata:
                        mock_metadata.return_value = 'metadata.json'
                        
                        with patch.object(orchestrator.output_manager, 'save_process_info') as mock_process:
                            mock_process.return_value = 'process_info.json'
                            
                            result = orchestrator.run_workflow(
                                input_data,
                                temp_dir,
                                'test_workspace',
                                selected_analyses=['network_analysis']
                            )
                            
                            assert result.success is True
                            assert len(result.final_output.get('output_files', [])) > 0
                            mock_output.assert_called()
                            mock_metadata.assert_called()
                            mock_process.assert_called()
    
    def test_workflow_performance_metrics(self, test_config, mock_kb_util, temp_dir, sample_protein_sequences):
        """Test workflow performance metrics collection."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': sample_protein_sequences,
            'workspace_name': 'test_workspace',
            'analysis_name': 'performance_test'
        }
        
        # Mock components with timing
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.15
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': 3, 'edges': 2},
                            'node_centrality': {},
                            'community_structure': {'communities': [['protein_0', 'protein_1'], ['protein_2']]}
                        },
                        'processing_time': 0.25
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['network_analysis.json'],
                        'summary': 'Performance test completed'
                    }
                    
                    result = orchestrator.run_workflow(
                        input_data,
                        temp_dir,
                        'test_workspace',
                        selected_analyses=['network_analysis']
                    )
                    
                    assert result.success is True
                    assert result.execution_time > 0
                    assert result.execution_time >= 0.4  # Input + Analysis time
                    assert 'stages_completed' in result.final_output
    
    def test_workflow_configuration_validation(self, test_config, mock_kb_util, temp_dir):
        """Test workflow configuration validation."""
        # Test with invalid configuration
        invalid_config = {
            'input_type': 'invalid_type',
            'data': 'some_data'
        }
        
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        result = orchestrator.run_workflow(
            invalid_config,
            temp_dir,
            'test_workspace'
        )
        
        assert result.success is False
        assert result.error_message is not None
    
    def test_workflow_cleanup(self, test_config, mock_kb_util, temp_dir, sample_protein_sequences):
        """Test workflow cleanup and resource management."""
        orchestrator = WorkflowOrchestrator(test_config, mock_kb_util)
        orchestrator.initialize_components(temp_dir, 'test_workspace')
        
        # Mock successful workflow
        with patch.object(orchestrator.input_manager, 'process_input') as mock_input:
            mock_input.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            
            with patch.object(orchestrator.analysis_manager, 'run_multiple_analyses') as mock_analysis:
                mock_analysis.return_value = {
                    'network_analysis': {
                        'success': True,
                        'analysis_type': 'network_analysis',
                        'results': {
                            'network_graph': {'nodes': 3, 'edges': 2},
                            'node_centrality': {},
                            'community_structure': {'communities': [['protein_0', 'protein_1'], ['protein_2']]}
                        }
                    }
                }
                
                with patch.object(orchestrator.output_manager, 'save_analysis_output') as mock_output:
                    mock_output.return_value = {
                        'success': True,
                        'output_files': ['network_analysis.json'],
                        'summary': 'Cleanup test completed'
                    }
                    
                    result = orchestrator.run_workflow(
                        {
                            'input_type': 'protein_input',
                            'protein_input': sample_protein_sequences,
                            'workspace_name': 'test_workspace',
                            'analysis_name': 'cleanup_test'
                        },
                        temp_dir,
                        'test_workspace',
                        selected_analyses=['network_analysis']
                    )
                    
                    assert result.success is True
        
        # Test cleanup
        orchestrator.cleanup()
        
        # Verify cleanup was called (if implemented)
        # This test ensures the cleanup method exists and can be called
        assert hasattr(orchestrator, 'cleanup')

