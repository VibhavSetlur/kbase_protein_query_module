"""
Unit tests for Analysis Configuration
"""
import pytest
from lib.kbase_protein_query_module.src.analysis.config import (
    get_enabled_analyses, get_analysis_by_category, is_analysis_enabled,
    get_analysis_dependencies, validate_analysis_config, ANALYSIS_CONFIG, OUTPUT_CONFIG
)


class TestAnalysisConfig:
    """Test cases for Analysis Configuration"""
    
    def test_get_enabled_analyses(self):
        """Test getting enabled analyses"""
        enabled = get_enabled_analyses()
        assert isinstance(enabled, dict)
        assert len(enabled) > 0
        
        # All returned analyses should be enabled
        for name, config in enabled.items():
            assert config.get("enabled", False) is True
    
    def test_get_analysis_by_category(self):
        """Test getting analyses by category"""
        # Test with existing category
        network_analyses = get_analysis_by_category("network")
        assert isinstance(network_analyses, dict)
        
        # Test with non-existent category
        empty_analyses = get_analysis_by_category("non_existent")
        assert isinstance(empty_analyses, dict)
        assert len(empty_analyses) == 0
    
    def test_is_analysis_enabled(self):
        """Test checking if analysis is enabled"""
        # Test with enabled analysis
        assert is_analysis_enabled("network_analysis") is True
        assert is_analysis_enabled("sequence_analysis") is True
        
        # Test with non-existent analysis
        assert is_analysis_enabled("non_existent") is False
    
    def test_get_analysis_dependencies(self):
        """Test getting analysis dependencies"""
        # Test with analysis that has dependencies
        deps = get_analysis_dependencies("network_analysis")
        assert isinstance(deps, list)
        assert "similarity_search" in deps
        assert "embeddings" in deps
        
        # Test with analysis that has no dependencies
        deps = get_analysis_dependencies("sequence_analysis")
        assert isinstance(deps, list)
        assert "embeddings" in deps
        
        # Test with non-existent analysis
        deps = get_analysis_dependencies("non_existent")
        assert isinstance(deps, list)
        assert len(deps) == 0
    
    def test_validate_analysis_config(self):
        """Test validating analysis configuration"""
        # Should return True for valid config
        assert validate_analysis_config() is True
    
    def test_analysis_config_structure(self):
        """Test analysis configuration structure"""
        assert isinstance(ANALYSIS_CONFIG, dict)
        assert len(ANALYSIS_CONFIG) > 0
        
        # Check that each analysis has required fields
        for name, config in ANALYSIS_CONFIG.items():
            assert "enabled" in config
            assert "name" in config
            assert "description" in config
            assert "category" in config
            assert "dependencies" in config
            assert "output_type" in config
    
    def test_output_config_structure(self):
        """Test output configuration structure"""
        assert isinstance(OUTPUT_CONFIG, dict)
        assert "base_output_dir" in OUTPUT_CONFIG
        assert "create_analysis_folders" in OUTPUT_CONFIG
        assert "include_metadata" in OUTPUT_CONFIG
        assert "include_process_info" in OUTPUT_CONFIG
        assert "formats" in OUTPUT_CONFIG
        assert "compression" in OUTPUT_CONFIG
    
    def test_analysis_categories(self):
        """Test that all analyses have valid categories"""
        categories = set()
        for config in ANALYSIS_CONFIG.values():
            category = config.get("category")
            assert category is not None
            categories.add(category)
        
        # Should have multiple categories
        assert len(categories) > 1
    
    def test_analysis_dependencies_consistency(self):
        """Test that dependencies reference valid analyses"""
        analysis_names = set(ANALYSIS_CONFIG.keys())
        
        for name, config in ANALYSIS_CONFIG.items():
            dependencies = config.get("dependencies", [])
            for dep in dependencies:
                assert dep in analysis_names, f"Analysis '{name}' has invalid dependency '{dep}'"
