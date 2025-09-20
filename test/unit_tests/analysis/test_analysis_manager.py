"""
Unit tests for AnalysisManager
"""
import pytest
from unittest.mock import Mock, patch
from lib.kbase_protein_query_module.src.analysis.analysis_manager import AnalysisManager


class TestAnalysisManager:
    """Test cases for AnalysisManager"""
    
    def test_analysis_manager_initialization(self):
        """Test AnalysisManager initializes correctly"""
        manager = AnalysisManager()
        assert manager.analyses == {}
        assert manager.results == {}
    
    def test_register_analysis(self):
        """Test registering an analysis"""
        manager = AnalysisManager()
        
        mock_analysis = Mock()
        mock_analysis.get_metadata.return_value = Mock(name="test_analysis")
        
        manager.register_analysis("test_analysis", mock_analysis)
        assert "test_analysis" in manager.analyses
        assert manager.analyses["test_analysis"] == mock_analysis
    
    def test_get_analysis(self):
        """Test getting a registered analysis"""
        manager = AnalysisManager()
        
        mock_analysis = Mock()
        manager.analyses["test_analysis"] = mock_analysis
        
        result = manager.get_analysis("test_analysis")
        assert result == mock_analysis
        
        # Test non-existent analysis
        result = manager.get_analysis("non_existent")
        assert result is None
    
    def test_list_analyses(self):
        """Test listing available analyses"""
        manager = AnalysisManager()
        
        mock_analysis1 = Mock()
        mock_analysis1.get_metadata.return_value = Mock(name="analysis1")
        mock_analysis2 = Mock()
        mock_analysis2.get_metadata.return_value = Mock(name="analysis2")
        
        manager.analyses["analysis1"] = mock_analysis1
        manager.analyses["analysis2"] = mock_analysis2
        
        analyses = manager.list_analyses()
        assert len(analyses) == 2
        assert "analysis1" in analyses
        assert "analysis2" in analyses
    
    def test_run_analysis(self):
        """Test running a single analysis"""
        manager = AnalysisManager()
        
        mock_analysis = Mock()
        mock_analysis.analyze.return_value = {"result": "success"}
        manager.analyses["test_analysis"] = mock_analysis
        
        proteins = [{"id": "P12345", "sequence": "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"}]
        
        result = manager.run_analysis("test_analysis", proteins)
        assert result == {"result": "success"}
        mock_analysis.analyze.assert_called_once_with(proteins)
    
    def test_run_analysis_not_found(self):
        """Test running non-existent analysis"""
        manager = AnalysisManager()
        
        proteins = [{"id": "P12345", "sequence": "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"}]
        
        result = manager.run_analysis("non_existent", proteins)
        assert result is None
    
    def test_run_analysis_with_error(self):
        """Test running analysis that raises exception"""
        manager = AnalysisManager()
        
        mock_analysis = Mock()
        mock_analysis.analyze.side_effect = Exception("Analysis failed")
        manager.analyses["test_analysis"] = mock_analysis
        
        proteins = [{"id": "P12345", "sequence": "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"}]
        
        result = manager.run_analysis("test_analysis", proteins)
        assert result is None
    
    def test_run_multiple_analyses(self):
        """Test running multiple analyses"""
        manager = AnalysisManager()
        
        mock_analysis1 = Mock()
        mock_analysis1.analyze.return_value = {"result1": "success"}
        mock_analysis2 = Mock()
        mock_analysis2.analyze.return_value = {"result2": "success"}
        
        manager.analyses["analysis1"] = mock_analysis1
        manager.analyses["analysis2"] = mock_analysis2
        
        proteins = [{"id": "P12345", "sequence": "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"}]
        
        results = manager.run_multiple_analyses(["analysis1", "analysis2"], proteins)
        assert len(results) == 2
        assert "analysis1" in results
        assert "analysis2" in results
        assert results["analysis1"] == {"result1": "success"}
        assert results["analysis2"] == {"result2": "success"}
    
    def test_get_analysis_results(self):
        """Test getting analysis results"""
        manager = AnalysisManager()
        
        manager.results["test_analysis"] = {"result": "success"}
        
        result = manager.get_analysis_results("test_analysis")
        assert result == {"result": "success"}
        
        # Test non-existent result
        result = manager.get_analysis_results("non_existent")
        assert result is None
    
    def test_clear_results(self):
        """Test clearing analysis results"""
        manager = AnalysisManager()
        
        manager.results["test_analysis"] = {"result": "success"}
        assert len(manager.results) == 1
        
        manager.clear_results()
        assert len(manager.results) == 0
