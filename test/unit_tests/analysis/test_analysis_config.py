"""
Unit tests for analysis configuration.

Tests configuration management, validation, and analysis discovery.
"""

import pytest
import sys
import os
from unittest.mock import patch

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.analysis.config import (
    get_enabled_analyses,
    get_analysis_by_category,
    is_analysis_enabled,
    get_analysis_dependencies,
    validate_analysis_config,
    ANALYSIS_CONFIG
)


class TestAnalysisConfig:
    """Test cases for analysis configuration."""
    
    def test_analysis_config_structure(self):
        """Test that ANALYSIS_CONFIG has proper structure."""
        assert hasattr(ANALYSIS_CONFIG, '__getitem__')  # It should be dict-like
        
        # Check that network_analysis is configured
        assert 'network_analysis' in ANALYSIS_CONFIG
        
        network_config = ANALYSIS_CONFIG['network_analysis']
        assert 'enabled' in network_config
        assert 'name' in network_config
        assert 'description' in network_config
        assert 'category' in network_config
        assert 'dependencies' in network_config
        assert 'output_type' in network_config
        assert 'module_path' in network_config
        assert 'class_name' in network_config
    
    def test_get_enabled_analyses(self):
        """Test getting enabled analyses."""
        enabled_analyses = get_enabled_analyses()
        
        assert isinstance(enabled_analyses, dict)
        
        # Should only return enabled analyses
        for analysis_name, config in enabled_analyses.items():
            assert config.get('enabled', False) is True
    
    def test_get_enabled_analyses_empty_when_none_enabled(self):
        """Test getting enabled analyses when none are enabled."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {}):
            enabled_analyses = get_enabled_analyses()
            assert enabled_analyses == {}
    
    def test_get_analysis_by_category(self):
        """Test getting analyses by category."""
        network_analyses = get_analysis_by_category('network')
        
        assert isinstance(network_analyses, dict)
        
        # Should only return analyses in the specified category
        for analysis_name, config in network_analyses.items():
            assert config.get('category') == 'network'
    
    def test_get_analysis_by_category_nonexistent(self):
        """Test getting analyses by non-existent category."""
        nonexistent_analyses = get_analysis_by_category('nonexistent')
        
        assert isinstance(nonexistent_analyses, dict)
        assert len(nonexistent_analyses) == 0
    
    def test_is_analysis_enabled_true(self):
        """Test checking if enabled analysis is enabled."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'test_analysis': {'enabled': True}
        }):
            assert is_analysis_enabled('test_analysis') is True
    
    def test_is_analysis_enabled_false(self):
        """Test checking if disabled analysis is enabled."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'test_analysis': {'enabled': False}
        }):
            assert is_analysis_enabled('test_analysis') is False
    
    def test_is_analysis_enabled_nonexistent(self):
        """Test checking if non-existent analysis is enabled."""
        assert is_analysis_enabled('nonexistent_analysis') is False
    
    def test_get_analysis_dependencies(self):
        """Test getting analysis dependencies."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'test_analysis': {'dependencies': ['dep1', 'dep2']}
        }):
            dependencies = get_analysis_dependencies('test_analysis')
            assert dependencies == ['dep1', 'dep2']
    
    def test_get_analysis_dependencies_empty(self):
        """Test getting dependencies for analysis with no dependencies."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'test_analysis': {}
        }):
            dependencies = get_analysis_dependencies('test_analysis')
            assert dependencies == []
    
    def test_get_analysis_dependencies_nonexistent(self):
        """Test getting dependencies for non-existent analysis."""
        dependencies = get_analysis_dependencies('nonexistent_analysis')
        assert dependencies == []
    
    def test_validate_analysis_config_valid(self):
        """Test validating valid analysis configuration."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'analysis1': {
                'enabled': True,
                'dependencies': []
            },
            'analysis2': {
                'enabled': True,
                'dependencies': ['analysis1']
            }
        }):
            assert validate_analysis_config() is True
    
    def test_validate_analysis_config_missing_dependency(self):
        """Test validating configuration with missing dependency."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'analysis1': {
                'enabled': True,
                'dependencies': ['nonexistent_analysis']
            }
        }):
            # Should return True but log warning
            assert validate_analysis_config() is True
    
    def test_validate_analysis_config_disabled_dependency(self):
        """Test validating configuration with disabled dependency."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'analysis1': {
                'enabled': True,
                'dependencies': ['disabled_analysis']
            },
            'disabled_analysis': {
                'enabled': False,
                'dependencies': []
            }
        }):
            # Should return True but log warning
            assert validate_analysis_config() is True
    
    def test_validate_analysis_config_exception(self):
        """Test validating configuration with exception."""
        with patch('kbase_protein_query_module.src.analysis.config.get_enabled_analyses', side_effect=Exception("Config error")):
            assert validate_analysis_config() is False
    
    def test_network_analysis_configuration(self):
        """Test network analysis specific configuration."""
        network_config = ANALYSIS_CONFIG.get('network_analysis', {})
        
        # Check if network analysis is available (may be disabled due to missing dependencies)
        if network_config.get('enabled', False):
            assert network_config.get('enabled') is True
            assert network_config.get('name') == 'Network Analysis'
            assert network_config.get('category') == 'network'
            assert network_config.get('output_type') == 'interactive_html'
            assert 'NetworkAnalysis' in network_config.get('class_name', '')
        else:
            # If disabled due to missing dependencies, verify the structure is still correct
            assert network_config.get('name') == 'Network Analysis'
            assert network_config.get('category') == 'network'
            assert network_config.get('output_type') == 'interactive_html'
            assert 'NetworkAnalysis' in network_config.get('class_name', '')
            # The enabled flag should be False when dependencies are missing
            assert network_config.get('enabled') is False
    
    def test_analysis_config_immutability(self):
        """Test that analysis configuration cannot be accidentally modified."""
        original_config = ANALYSIS_CONFIG.copy()
        
        # Try to modify the config
        try:
            ANALYSIS_CONFIG['test_analysis'] = {'enabled': True}
        except Exception:
            pass  # Expected if config is read-only
        
        # Config should remain unchanged
        assert ANALYSIS_CONFIG == original_config
    
    def test_configuration_types(self):
        """Test that configuration values have correct types."""
        for analysis_name, config in ANALYSIS_CONFIG.items():
            assert isinstance(config, dict)
            assert isinstance(config.get('enabled'), bool)
            assert isinstance(config.get('name'), str)
            assert isinstance(config.get('description'), str)
            assert isinstance(config.get('category'), str)
            assert isinstance(config.get('dependencies'), list)
            assert isinstance(config.get('output_type'), str)
            assert isinstance(config.get('module_path'), str)
            assert isinstance(config.get('class_name'), str)
    
    def test_dependency_cycles_detection(self):
        """Test that dependency cycles are handled gracefully."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {
            'analysis1': {
                'enabled': True,
                'dependencies': ['analysis2']
            },
            'analysis2': {
                'enabled': True,
                'dependencies': ['analysis1']  # Circular dependency
            }
        }):
            # Should not crash, just return True
            assert validate_analysis_config() is True
    
    def test_empty_configuration(self):
        """Test behavior with empty configuration."""
        with patch('kbase_protein_query_module.src.analysis.config.ANALYSIS_CONFIG', {}):
            assert get_enabled_analyses() == {}
            assert get_analysis_by_category('any') == {}
            assert is_analysis_enabled('any') is False
            assert get_analysis_dependencies('any') == []
            assert validate_analysis_config() is True
    
    def test_configuration_consistency(self):
        """Test that configuration is internally consistent."""
        enabled_analyses = get_enabled_analyses()
        
        for analysis_name, config in enabled_analyses.items():
            # Check that all required fields are present
            required_fields = ['name', 'description', 'category', 'dependencies', 'output_type', 'module_path', 'class_name']
            for field in required_fields:
                assert field in config, f"Missing field '{field}' in analysis '{analysis_name}'"
            
            # Check that dependencies are valid
            dependencies = config.get('dependencies', [])
            assert isinstance(dependencies, list)
            
            # Check that category is not empty
            assert config.get('category', '').strip() != ''
            
            # Check that class name is not empty
            assert config.get('class_name', '').strip() != ''

