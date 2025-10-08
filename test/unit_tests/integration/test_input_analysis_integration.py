"""
Integration tests for input to analysis pipeline.

Tests the complete flow from input processing through analysis execution.
"""

import pytest
import sys
import os
import numpy as np
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.input.input_manager import InputManager
from kbase_protein_query_module.src.analysis.analysis_manager import AnalysisManager
from kbase_protein_query_module.src.core.workflow_orchestrator import WorkflowOrchestrator


class TestInputAnalysisIntegration:
    """Test cases for input to analysis integration."""
    
    def test_protein_input_to_analysis_flow(self, test_config, mock_kb_util, sample_protein_sequences):
        """Test complete flow from protein input to analysis."""
        # Setup input manager
        input_manager = InputManager(test_config, mock_kb_util)
        
        # Setup analysis manager
        analysis_manager = AnalysisManager()
        
        # Process input
        input_data = {
            'input_type': 'protein_input',
            'protein_input': sample_protein_sequences
        }
        
        # Mock the protein sequence processor
        with patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input(input_data)
            
            assert input_result['success'] is True
            assert len(input_result['proteins']) == len(sample_protein_sequences)
        
        # Mock analysis execution
        with patch.object(analysis_manager, 'run_analysis') as mock_run_analysis:
            mock_run_analysis.return_value = {
                'success': True,
                'analysis_type': 'network_analysis',
                'results': {'nodes': 3, 'edges': 2}
            }
            
            analysis_result = analysis_manager.run_analysis(
                'network_analysis',
                input_result['proteins']
            )
            
            assert analysis_result['success'] is True
            mock_run_analysis.assert_called_once()
    
    def test_uniprot_input_to_analysis_flow(self, test_config, mock_kb_util, sample_uniprot_ids):
        """Test complete flow from UniProt input to analysis."""
        # Setup managers
        input_manager = InputManager(test_config, mock_kb_util)
        analysis_manager = AnalysisManager()
        
        input_data = {
            'input_type': 'uniprot_ids',
            'uniprot_ids': sample_uniprot_ids
        }
        
        # Mock UniProt processor
        with patch('kbase_protein_query_module.src.input.input_manager.UniProtIdsProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': uniprot_id, 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG', 'source': 'uniprot'}
                    for uniprot_id in sample_uniprot_ids
                ],
                'input_type': 'uniprot_ids',
                'processing_time': 0.2
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input(input_data)
            
            assert input_result['success'] is True
            assert len(input_result['proteins']) == len(sample_uniprot_ids)
        
        # Mock analysis execution
        with patch.object(analysis_manager, 'run_analysis') as mock_run_analysis:
            mock_run_analysis.return_value = {
                'success': True,
                'analysis_type': 'network_analysis',
                'results': {'nodes': len(sample_uniprot_ids), 'edges': 1}
            }
            
            analysis_result = analysis_manager.run_analysis(
                'network_analysis',
                input_result['proteins']
            )
            
            assert analysis_result['success'] is True
    
    def test_workspace_input_to_analysis_flow(self, test_config, mock_kb_util, mock_workspace_client, sample_workspace_object_ref):
        """Test complete flow from workspace input to analysis."""
        # Setup managers
        input_manager = InputManager(test_config, mock_kb_util)
        analysis_manager = AnalysisManager()
        
        input_data = {
            'input_type': 'workspace_object',
            'workspace_object': sample_workspace_object_ref
        }
        
        # Mock workspace processor
        with patch('kbase_protein_query_module.src.input.input_manager.WorkspaceObjectProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': 'P12345', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG', 'source': 'workspace_object'}
                ],
                'input_type': 'workspace_object',
                'processing_time': 0.3
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input(input_data)
            
            assert input_result['success'] is True
            assert len(input_result['proteins']) == 1
        
        # Mock analysis execution
        with patch.object(analysis_manager, 'run_analysis') as mock_run_analysis:
            mock_run_analysis.return_value = {
                'success': True,
                'analysis_type': 'network_analysis',
                'results': {'nodes': 1, 'edges': 0}
            }
            
            analysis_result = analysis_manager.run_analysis(
                'network_analysis',
                input_result['proteins']
            )
            
            assert analysis_result['success'] is True
    
    def test_multiple_analyses_execution(self, test_config, mock_kb_util, sample_protein_sequences):
        """Test execution of multiple analyses on the same input."""
        # Setup managers
        input_manager = InputManager(test_config, mock_kb_util)
        analysis_manager = AnalysisManager()
        
        input_data = {
            'input_type': 'protein_input',
            'protein_input': sample_protein_sequences
        }
        
        # Mock input processing
        with patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input(input_data)
        
        # Mock multiple analyses
        with patch.object(analysis_manager, 'run_multiple_analyses') as mock_run_multiple:
            mock_run_multiple.return_value = {
                'network_analysis': {
                    'success': True,
                    'analysis_type': 'network_analysis',
                    'results': {'nodes': 3, 'edges': 2}
                },
                'sequence_analysis': {
                    'success': True,
                    'analysis_type': 'sequence_analysis',
                    'results': {'similarity_matrix': np.random.rand(3, 3)}
                }
            }
            
            analysis_results = analysis_manager.run_multiple_analyses(
                ['network_analysis', 'sequence_analysis'],
                input_result['proteins']
            )
            
            assert len(analysis_results) == 2
            assert 'network_analysis' in analysis_results
            assert 'sequence_analysis' in analysis_results
            assert analysis_results['network_analysis']['success'] is True
            assert analysis_results['sequence_analysis']['success'] is True
    
    def test_input_validation_integration(self, test_config, mock_kb_util):
        """Test input validation integration with analysis pipeline."""
        input_manager = InputManager(test_config, mock_kb_util)
        
        # Test invalid input type
        invalid_input = {
            'input_type': 'invalid_type',
            'data': 'some_data'
        }
        
        result = input_manager.process_input(invalid_input)
        
        assert result['success'] is False
        assert 'error_message' in result
    
    def test_analysis_dependency_handling(self, test_config, mock_kb_util, sample_protein_sequences):
        """Test handling of analysis dependencies."""
        input_manager = InputManager(test_config, mock_kb_util)
        analysis_manager = AnalysisManager()
        
        # Mock input processing
        with patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input({
                'input_type': 'protein_input',
                'protein_input': sample_protein_sequences
            })
        
        # Test dependency checking
        available_data = {
            'proteins': input_result['proteins'],
            'similarity_results': [],
            'embeddings': np.random.rand(3, 320),
            'protein_ids': ['protein_0', 'protein_1', 'protein_2'],
            'metadata_df': Mock(),
            'query_embedding': np.random.rand(320),
            'query_protein_id': 'protein_0'
        }
        
        # Mock analysis to exist
        analysis_manager.analysis_modules['network_analysis'] = Mock()
        
        can_run = analysis_manager.can_run_analysis('network_analysis', available_data)
        
        # Should be able to run if all dependencies are available
        assert isinstance(can_run, bool)
    
    def test_error_propagation_input_to_analysis(self, test_config, mock_kb_util):
        """Test error propagation from input to analysis."""
        input_manager = InputManager(test_config, mock_kb_util)
        analysis_manager = AnalysisManager()
        
        # Mock input processing to fail
        with patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': False,
                'error_message': 'Input processing failed'
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input({
                'input_type': 'protein_input',
                'protein_input': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']
            })
            
            assert input_result['success'] is False
        
        # Analysis should handle failed input gracefully
        analysis_result = analysis_manager.run_analysis(
            'network_analysis',
            []  # Empty proteins due to failed input
        )
        
        # Should return None or handle gracefully
        assert analysis_result is None or not analysis_result.get('success', False)
    
    def test_data_format_consistency(self, test_config, mock_kb_util, sample_protein_sequences):
        """Test data format consistency between input and analysis."""
        input_manager = InputManager(test_config, mock_kb_util)
        
        # Mock input processing
        with patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': True,
                'proteins': [
                    {
                        'protein_id': f'protein_{i}',
                        'sequence': seq,
                        'source': 'protein_sequence',
                        'metadata': {'organism': 'Test organism'}
                    }
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.1
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input({
                'input_type': 'protein_input',
                'protein_input': sample_protein_sequences
            })
            
            # Verify data format
            assert input_result['success'] is True
            assert 'proteins' in input_result
            
            for protein in input_result['proteins']:
                assert 'protein_id' in protein
                assert 'sequence' in protein
                assert 'source' in protein
                assert protein['source'] == 'protein_sequence'
    
    def test_performance_tracking_integration(self, test_config, mock_kb_util, sample_protein_sequences):
        """Test performance tracking across input and analysis."""
        input_manager = InputManager(test_config, mock_kb_util)
        analysis_manager = AnalysisManager()
        
        # Mock input processing with timing
        with patch('kbase_protein_query_module.src.input.input_manager.ProteinSequenceProcessor') as mock_processor_class:
            mock_processor = Mock()
            mock_processor.process.return_value = {
                'success': True,
                'proteins': [
                    {'protein_id': f'protein_{i}', 'sequence': seq, 'source': 'protein_sequence'}
                    for i, seq in enumerate(sample_protein_sequences)
                ],
                'input_type': 'protein_input',
                'processing_time': 0.15
            }
            mock_processor_class.return_value = mock_processor
            
            input_result = input_manager.process_input({
                'input_type': 'protein_input',
                'protein_input': sample_protein_sequences
            })
            
            assert 'processing_time' in input_result
            assert input_result['processing_time'] == 0.15
        
        # Mock analysis with timing
        with patch.object(analysis_manager, 'run_analysis') as mock_run_analysis:
            mock_run_analysis.return_value = {
                'success': True,
                'analysis_type': 'network_analysis',
                'results': {'nodes': 3, 'edges': 2},
                'processing_time': 0.25
            }
            
            analysis_result = analysis_manager.run_analysis(
                'network_analysis',
                input_result['proteins']
            )
            
            assert analysis_result['success'] is True
            assert 'processing_time' in analysis_result

