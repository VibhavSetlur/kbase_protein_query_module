"""
Integration tests for directory creation and file system operations.

These tests ensure that the module properly creates output directories
and handles file system operations correctly, preventing runtime errors
that would only be discovered in the KBase narrative environment.
"""

import os
import tempfile
import shutil
import pytest
from unittest.mock import Mock

from kbase_protein_query_module import kbase_protein_query_moduleImpl as impl_module


class TestDirectoryCreation:
    """Test directory creation and file system operations."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp(prefix='protein_query_dir_test_')
        self.mock_ctx = Mock()
        self.mock_ctx.token = 'mock_token'
        self.mock_ctx.provenance = [{'service': 'kbase_protein_query_module', 'method': 'run_protein_query_analysis'}]
    
    def teardown_method(self):
        """Clean up test environment."""
        try:
            shutil.rmtree(self.temp_dir)
        except Exception:
            pass  # Ignore cleanup errors in tests
    
    def test_output_directory_creation_protein_sequence(self):
        """Test that output directories are created properly for protein sequence input."""
        service = impl_module.kbase_protein_query_module({})
        
        params = {
            'workspace_name': 'test_workspace',
            'analysis_name': 'dir_test_protein',
            'input_type': 'protein_input',
            'protein_input': ['MKFLVNVALVFMVVYISYIYA'],
            'analysis_stages': ['network_analysis'],
            'output_dir': self.temp_dir
        }
        
        # Run the analysis
        result = service.run_protein_query_analysis(self.mock_ctx, params)
        
        # Check that output directory was created
        expected_output_dir = os.path.join(self.temp_dir, 'outputs')
        assert os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}"
        
        # Check that run-specific subdirectory was created
        contents = os.listdir(expected_output_dir)
        assert len(contents) > 0, f"No run directories found in {expected_output_dir}"
        
        run_dir = os.path.join(expected_output_dir, contents[0])
        assert os.path.isdir(run_dir), f"Run directory not created properly: {run_dir}"
        
        # Check that expected subdirectories exist
        expected_subdirs = ['metadata', 'process_info', 'analysis', 'logs']
        for subdir in expected_subdirs:
            subdir_path = os.path.join(run_dir, subdir)
            assert os.path.exists(subdir_path), f"Expected subdirectory not created: {subdir_path}"
        
        # Verify the result doesn't contain directory creation errors
        assert result and len(result) > 0
        output = result[0]
        
        # Should not have directory creation errors
        summary = output.get('summary', '')
        assert 'No such file or directory' not in summary, f"Directory creation error found: {summary}"
        assert 'outputs' not in summary or 'A url is required' in summary, f"Unexpected directory error: {summary}"
    
    def test_output_directory_creation_uniprot_ids(self):
        """Test that output directories are created properly for UniProt IDs input."""
        service = impl_module.kbase_protein_query_module({})
        
        params = {
            'workspace_name': 'test_workspace',
            'analysis_name': 'dir_test_uniprot',
            'input_type': 'uniprot_ids',
            'uniprot_ids': ['P12345', 'Q67890'],
            'analysis_stages': ['network_analysis'],
            'output_dir': self.temp_dir
        }
        
        # Run the analysis
        result = service.run_protein_query_analysis(self.mock_ctx, params)
        
        # Check that output directory was created
        expected_output_dir = os.path.join(self.temp_dir, 'outputs')
        assert os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}"
        
        # Check that run-specific subdirectory was created
        contents = os.listdir(expected_output_dir)
        assert len(contents) > 0, f"No run directories found in {expected_output_dir}"
        
        run_dir = os.path.join(expected_output_dir, contents[0])
        assert os.path.isdir(run_dir), f"Run directory not created properly: {run_dir}"
        
        # Check that expected subdirectories exist
        expected_subdirs = ['metadata', 'process_info', 'analysis', 'logs']
        for subdir in expected_subdirs:
            subdir_path = os.path.join(run_dir, subdir)
            assert os.path.exists(subdir_path), f"Expected subdirectory not created: {subdir_path}"
        
        # Verify the result doesn't contain directory creation errors
        assert result and len(result) > 0
        output = result[0]
        
        # Should not have directory creation errors
        summary = output.get('summary', '')
        assert 'No such file or directory' not in summary, f"Directory creation error found: {summary}"
        assert 'outputs' not in summary or 'A url is required' in summary, f"Unexpected directory error: {summary}"
    
    def test_output_directory_creation_nonexistent_parent(self):
        """Test that output directories are created even when parent doesn't exist."""
        # Create a nested directory that doesn't exist yet
        nested_dir = os.path.join(self.temp_dir, 'nested', 'deep', 'directory')
        
        service = impl_module.kbase_protein_query_module({})
        
        params = {
            'workspace_name': 'test_workspace',
            'analysis_name': 'dir_test_nested',
            'input_type': 'protein_input',
            'protein_input': ['MKFLVNVALVFMVVYISYIYA'],
            'analysis_stages': ['network_analysis'],
            'output_dir': nested_dir
        }
        
        # Run the analysis
        result = service.run_protein_query_analysis(self.mock_ctx, params)
        
        # Check that the nested directory was created
        assert os.path.exists(nested_dir), f"Nested directory not created: {nested_dir}"
        
        # Check that output directory was created
        expected_output_dir = os.path.join(nested_dir, 'outputs')
        assert os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}"
        
        # Check that run-specific subdirectory was created
        contents = os.listdir(expected_output_dir)
        assert len(contents) > 0, f"No run directories found in {expected_output_dir}"
        
        run_dir = os.path.join(expected_output_dir, contents[0])
        assert os.path.isdir(run_dir), f"Run directory not created properly: {run_dir}"
        
        # Verify the result doesn't contain directory creation errors
        assert result and len(result) > 0
        output = result[0]
        
        # Should not have directory creation errors
        summary = output.get('summary', '')
        assert 'No such file or directory' not in summary, f"Directory creation error found: {summary}"
        assert 'outputs' not in summary or 'A url is required' in summary, f"Unexpected directory error: {summary}"
    
    def test_workflow_orchestrator_directory_creation(self):
        """Test workflow orchestrator directory creation directly."""
        from kbase_protein_query_module.src.core.workflow_orchestrator import WorkflowOrchestrator
        from kbase_protein_query_module.src.core.pipeline_config import PipelineConfig
        
        # Create a nested directory that doesn't exist yet
        nested_dir = os.path.join(self.temp_dir, 'workflow', 'test')
        
        config = {
            'workspace_name': 'test_workspace',
            'analysis_name': 'workflow_dir_test',
            'input_type': 'protein_input',
            'protein_input': ['MKFLVNVALVFMVVYISYIYA'],
            'analysis_stages': ['network_analysis']
        }
        
        pipeline_config = PipelineConfig(config)
        pipeline_config.output_dir = nested_dir
        
        # Create orchestrator
        orchestrator = WorkflowOrchestrator(config=pipeline_config)
        
        # Execute workflow
        result = orchestrator.execute(pipeline_config)
        
        # Check that the nested directory was created
        assert os.path.exists(nested_dir), f"Nested directory not created: {nested_dir}"
        
        # Check that output directory was created
        expected_output_dir = os.path.join(nested_dir, 'outputs')
        assert os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}"
        
        # Check that run-specific subdirectory was created
        contents = os.listdir(expected_output_dir)
        assert len(contents) > 0, f"No run directories found in {expected_output_dir}"
        
        run_dir = os.path.join(expected_output_dir, contents[0])
        assert os.path.isdir(run_dir), f"Run directory not created properly: {run_dir}"
        
        # Check that expected subdirectories exist
        expected_subdirs = ['metadata', 'process_info', 'analysis', 'logs']
        for subdir in expected_subdirs:
            subdir_path = os.path.join(run_dir, subdir)
            assert os.path.exists(subdir_path), f"Expected subdirectory not created: {subdir_path}"
    
    def test_output_manager_directory_creation(self):
        """Test output manager directory creation directly."""
        from kbase_protein_query_module.src.output.output_manager import OutputManager
        
        # Create a nested directory that doesn't exist yet
        nested_dir = os.path.join(self.temp_dir, 'output_manager', 'test')
        
        # Create output manager
        output_manager = OutputManager(
            base_output_dir=nested_dir,
            run_id='test_run_123',
            workspace_name='test_workspace'
        )
        
        # Check that the nested directory was created
        assert os.path.exists(nested_dir), f"Nested directory not created: {nested_dir}"
        
        # Check that output directory was created
        expected_output_dir = os.path.join(nested_dir, 'outputs')
        assert os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}"
        
        # Check that run-specific subdirectory was created
        contents = os.listdir(expected_output_dir)
        assert len(contents) > 0, f"No run directories found in {expected_output_dir}"
        
        run_dir = os.path.join(expected_output_dir, contents[0])
        assert os.path.isdir(run_dir), f"Run directory not created properly: {run_dir}"
        
        # Check that expected subdirectories exist
        expected_subdirs = ['metadata', 'process_info', 'analysis', 'logs']
        for subdir in expected_subdirs:
            subdir_path = os.path.join(run_dir, subdir)
            assert os.path.exists(subdir_path), f"Expected subdirectory not created: {subdir_path}"
    
    def test_error_handling_missing_directory_permissions(self):
        """Test error handling when directory creation fails due to permissions."""
        # Create a directory with restricted permissions (if possible)
        restricted_dir = os.path.join(self.temp_dir, 'restricted')
        os.makedirs(restricted_dir)
        
        # Try to make it read-only (this might not work on all systems)
        try:
            os.chmod(restricted_dir, 0o444)  # Read-only
            
            service = impl_module.kbase_protein_query_module({})
            
            params = {
                'workspace_name': 'test_workspace',
                'analysis_name': 'permission_test',
                'input_type': 'protein_input',
                'protein_input': ['MKFLVNVALVFMVVYISYIYA'],
                'analysis_stages': ['network_analysis'],
                'output_dir': os.path.join(restricted_dir, 'outputs')
            }
            
            # This should either work or fail gracefully
            result = service.run_protein_query_analysis(self.mock_ctx, params)
            
            # If it works, that's fine
            if result and len(result) > 0:
                output = result[0]
                # Should not crash the application
                assert 'job_id' in output
                assert 'analysis_result_ref' in output
            
        except PermissionError:
            # Permission denied is expected behavior
            pass
        except Exception as e:
            # Other exceptions should be handled gracefully
            assert 'No such file or directory' not in str(e), f"Unexpected directory error: {e}"
        finally:
            # Restore permissions for cleanup
            try:
                os.chmod(restricted_dir, 0o755)
            except Exception:
                pass
