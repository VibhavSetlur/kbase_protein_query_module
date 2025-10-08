"""
Unit tests for AnalysisManager.

Tests analysis coordination, dependency management, and execution.
"""

import pytest
import sys
import os
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.analysis.analysis_manager import AnalysisManager


@pytest.fixture
def mock_analysis_manager():
    """Fixture that provides AnalysisManager with mocked module loading."""
    with patch.object(AnalysisManager, '_load_analysis_modules'):
        return AnalysisManager()


class TestAnalysisManager:
    """Test cases for AnalysisManager."""
    
    def test_initialization(self, test_config, mock_analysis_manager):
        """Test AnalysisManager initialization."""
        manager = mock_analysis_manager
        
        assert manager.output_manager is None
        assert isinstance(manager.analysis_modules, dict)
        assert isinstance(manager.results, dict)
        assert isinstance(manager.analyses, dict)
    
    def test_initialization_with_output_manager(self, test_config):
        """Test initialization with output manager."""
        mock_output_manager = Mock()
        manager = AnalysisManager(output_manager=mock_output_manager)
        
        assert manager.output_manager == mock_output_manager
    
    @patch('kbase_protein_query_module.src.analysis.analysis_manager.get_enabled_analyses')
    def test_load_analysis_modules(self, mock_get_analyses, test_config):
        """Test loading analysis modules."""
        mock_get_analyses.return_value = {
            'network_analysis': {
                'enabled': True,
                'module_path': 'analysis.network_analysis.network_analysis',
                'class_name': 'NetworkAnalysis'
            }
        }
        
        # Mock the analysis class
        mock_analysis_class = Mock()
        
        with patch('importlib.import_module') as mock_import:
            mock_module = Mock()
            mock_module.NetworkAnalysis = mock_analysis_class
            mock_import.return_value = mock_module
            
            manager = AnalysisManager()
            
            assert 'network_analysis' in manager.analysis_modules
            assert 'network_analysis' in manager.analyses
    
    def test_get_available_analyses(self, test_config):
        """Test getting available analyses."""
        manager = AnalysisManager()
        
        analyses = manager.get_available_analyses()
        
        assert isinstance(analyses, dict)
    
    def test_get_analysis_dependencies(self, test_config):
        """Test getting analysis dependencies."""
        manager = AnalysisManager()
        
        dependencies = manager.get_analysis_dependencies('network_analysis')
        
        assert isinstance(dependencies, list)
    
    def test_can_run_analysis_with_available_data(self, test_config):
        """Test checking if analysis can run with available data."""
        manager = AnalysisManager()
        
        available_data = {
            'similarity_results': [],
            'embeddings': Mock(),
            'protein_ids': [],
            'metadata_df': Mock(),
            'query_embedding': Mock(),
            'query_protein_id': 'test'
        }
        
        # Mock the analysis to exist
        manager.analysis_modules['network_analysis'] = Mock()
        
        can_run = manager.can_run_analysis('network_analysis', available_data)
        
        # Should return True if analysis exists and dependencies are met
        assert isinstance(can_run, bool)
    
    def test_can_run_analysis_disabled(self, test_config):
        """Test checking disabled analysis."""
        manager = AnalysisManager()
        
        available_data = {'similarity_results': []}
        
        can_run = manager.can_run_analysis('disabled_analysis', available_data)
        
        assert can_run is False
    
    def test_can_run_analysis_missing_dependencies(self, test_config, mock_analysis_manager):
        """Test checking analysis with missing dependencies."""
        manager = mock_analysis_manager
        
        available_data = {}  # Missing required dependencies
        
        # Mock the analysis to exist
        manager.analysis_modules['network_analysis'] = Mock()
        
        # Mock get_analysis_dependencies to return some dependencies
        with patch('kbase_protein_query_module.src.analysis.analysis_manager.get_analysis_dependencies') as mock_get_deps:
            mock_get_deps.return_value = ['similarity_results', 'protein_data']
            
            can_run = manager.can_run_analysis('network_analysis', available_data)
            
            assert can_run is False
    
    def test_run_analysis_success(self, test_config):
        """Test running analysis successfully."""
        manager = AnalysisManager()
        
        # Mock analysis class
        mock_analysis = Mock()
        mock_analysis.analyze.return_value = {'success': True, 'results': 'test'}
        
        manager.analyses['test_analysis'] = mock_analysis
        
        result = manager.run_analysis(
            'test_analysis',
            ['protein1', 'protein2'],
            output_dir='/tmp'
        )
        
        assert result is not None
        assert result['success'] is True
        assert 'test_analysis' in manager.results
        mock_analysis.analyze.assert_called_once()
    
    def test_run_analysis_not_found(self, test_config):
        """Test running non-existent analysis."""
        manager = AnalysisManager()
        
        result = manager.run_analysis(
            'nonexistent_analysis',
            ['protein1', 'protein2']
        )
        
        assert result is None
    
    def test_run_analysis_error(self, test_config, mock_analysis_manager):
        """Test running analysis with error."""
        manager = mock_analysis_manager
        
        # Mock analysis class that raises exception
        mock_analysis = Mock()
        mock_analysis.analyze.side_effect = Exception("Analysis error")
        
        manager.analyses['test_analysis'] = mock_analysis
        
        with pytest.raises(Exception, match="Analysis error"):
            manager.run_analysis(
                'test_analysis',
                ['protein1', 'protein2']
            )
    
    def test_run_multiple_analyses(self, test_config):
        """Test running multiple analyses."""
        manager = AnalysisManager()
        
        # Mock analysis classes
        mock_analysis1 = Mock()
        mock_analysis1.analyze.return_value = {'success': True, 'results': 'test1'}
        
        mock_analysis2 = Mock()
        mock_analysis2.analyze.return_value = {'success': True, 'results': 'test2'}
        
        manager.analyses['analysis1'] = mock_analysis1
        manager.analyses['analysis2'] = mock_analysis2
        
        results = manager.run_multiple_analyses(
            ['analysis1', 'analysis2'],
            ['protein1', 'protein2']
        )
        
        assert 'analysis1' in results
        assert 'analysis2' in results
        assert results['analysis1']['success'] is True
        assert results['analysis2']['success'] is True
    
    def test_run_multiple_analyses_with_error(self, test_config, mock_analysis_manager):
        """Test running multiple analyses with one failing."""
        manager = mock_analysis_manager
        
        # Mock analysis classes - one successful, one failing
        mock_analysis1 = Mock()
        mock_analysis1.analyze.return_value = {'success': True, 'results': 'test1'}
        
        mock_analysis2 = Mock()
        mock_analysis2.analyze.side_effect = Exception("Analysis error")
        
        manager.analyses['analysis1'] = mock_analysis1
        manager.analyses['analysis2'] = mock_analysis2
        
        results = manager.run_multiple_analyses(
            ['analysis1', 'analysis2'],
            ['protein1', 'protein2']
        )
        
        assert 'analysis1' in results
        assert 'analysis2' in results
        assert results['analysis1']['success'] is True
        assert 'error' in results['analysis2']
    
    def test_order_analyses_by_dependencies(self, test_config):
        """Test ordering analyses by dependencies."""
        manager = AnalysisManager()
        
        # Mock dependency function
        def mock_get_dependencies(analysis_name):
            if analysis_name == 'analysis2':
                return ['analysis1']
            return []
        
        with patch.object(manager, 'get_analysis_dependencies', side_effect=mock_get_dependencies):
            ordered = manager._order_analyses_by_dependencies(['analysis2', 'analysis1'])
            
            # analysis1 should come before analysis2
            assert ordered.index('analysis1') < ordered.index('analysis2')
    
    def test_get_analysis_results_all(self, test_config):
        """Test getting all analysis results."""
        manager = AnalysisManager()
        
        # Add some results
        manager.results['analysis1'] = {'success': True, 'data': 'test1'}
        manager.results['analysis2'] = {'success': True, 'data': 'test2'}
        
        results = manager.get_analysis_results()
        
        assert len(results) == 2
        assert 'analysis1' in results
        assert 'analysis2' in results
    
    def test_get_analysis_results_specific(self, test_config):
        """Test getting specific analysis results."""
        manager = AnalysisManager()
        
        # Add some results
        manager.results['analysis1'] = {'success': True, 'data': 'test1'}
        manager.results['analysis2'] = {'success': True, 'data': 'test2'}
        
        result = manager.get_analysis_results('analysis1')
        
        assert result == {'success': True, 'data': 'test1'}
    
    def test_get_analysis_results_not_found(self, test_config):
        """Test getting results for non-existent analysis."""
        manager = AnalysisManager()
        
        result = manager.get_analysis_results('nonexistent')
        
        assert result is None
    
    def test_clear_results(self, test_config):
        """Test clearing analysis results."""
        manager = AnalysisManager()
        
        # Add some results
        manager.results['analysis1'] = {'success': True, 'data': 'test1'}
        manager.results['analysis2'] = {'success': True, 'data': 'test2'}
        
        assert len(manager.results) == 2
        
        manager.clear_results()
        
        assert len(manager.results) == 0
    
    def test_register_analysis(self, test_config):
        """Test registering analysis class."""
        manager = AnalysisManager()
        
        mock_analysis_class = Mock()
        
        manager.register_analysis('custom_analysis', mock_analysis_class)
        
        assert 'custom_analysis' in manager.analysis_modules
        assert 'custom_analysis' in manager.analyses
        assert manager.analyses['custom_analysis'] == mock_analysis_class
    
    def test_get_analysis(self, test_config):
        """Test getting analysis class by name."""
        manager = AnalysisManager()
        
        mock_analysis_class = Mock()
        manager.analyses['test_analysis'] = mock_analysis_class
        
        result = manager.get_analysis('test_analysis')
        
        assert result == mock_analysis_class
    
    def test_get_analysis_not_found(self, test_config):
        """Test getting non-existent analysis class."""
        manager = AnalysisManager()
        
        result = manager.get_analysis('nonexistent')
        
        assert result is None
    
    def test_list_analyses(self, test_config):
        """Test listing all analyses."""
        manager = AnalysisManager()
        
        # Add some analyses
        manager.analyses['analysis1'] = Mock()
        manager.analyses['analysis2'] = Mock()
        
        analyses = manager.list_analyses()
        
        assert len(analyses) == 2
        assert 'analysis1' in analyses
        assert 'analysis2' in analyses
    
    def test_analysis_module_loading_error(self, test_config):
        """Test handling of analysis module loading errors."""
        manager = AnalysisManager()
        
        # Mock get_enabled_analyses to return problematic config
        with patch('kbase_protein_query_module.src.analysis.analysis_manager.get_enabled_analyses') as mock_get:
            mock_get.return_value = {
                'bad_analysis': {
                    'enabled': True,
                    'module_path': 'nonexistent.module',
                    'class_name': 'NonExistentClass'
                }
            }
            
            # Should not raise exception, just log error
            manager._load_analysis_modules()
            
            # Bad analysis should not be loaded
            assert 'bad_analysis' not in manager.analysis_modules
    
    def test_analysis_class_not_found(self, test_config):
        """Test handling when analysis class is not found in module."""
        manager = AnalysisManager()
        
        with patch('kbase_protein_query_module.src.analysis.analysis_manager.get_enabled_analyses') as mock_get:
            mock_get.return_value = {
                'bad_analysis': {
                    'enabled': True,
                    'module_path': 'builtins',  # Use builtin module
                    'class_name': 'NonExistentClass'
                }
            }
            
            # Should not raise exception, just log error
            manager._load_analysis_modules()
            
            # Bad analysis should not be loaded
            assert 'bad_analysis' not in manager.analysis_modules
    
    def test_configuration_validation_failure(self, test_config):
        """Test handling configuration validation failure."""
        with patch('kbase_protein_query_module.src.analysis.analysis_manager.validate_analysis_config', return_value=False):
            # Should not raise exception, just log warning
            manager = AnalysisManager()
            
            assert manager.analysis_modules is not None
    
    def test_analysis_with_kwargs(self, test_config):
        """Test running analysis with additional keyword arguments."""
        manager = AnalysisManager()
        
        # Mock analysis class
        mock_analysis = Mock()
        mock_analysis.analyze.return_value = {'success': True, 'results': 'test'}
        
        manager.analyses['test_analysis'] = mock_analysis
        
        result = manager.run_analysis(
            'test_analysis',
            ['protein1', 'protein2'],
            output_dir='/tmp',
            custom_param='custom_value',
            another_param=123
        )
        
        assert result is not None
        # Verify that analyze was called with kwargs
        mock_analysis.analyze.assert_called_once()
        call_args = mock_analysis.analyze.call_args
        assert 'custom_param' in call_args.kwargs
        assert call_args.kwargs['custom_param'] == 'custom_value'
        assert call_args.kwargs['another_param'] == 123

