"""
Unit tests for OutputManager.

Tests output coordination, file management, and artifact tracking.
"""

import pytest
import sys
import os
import tempfile
import shutil
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.src.output.output_manager import OutputManager, ArtifactRecord


class TestOutputManager:
    """Test cases for OutputManager."""
    
    def test_initialization(self, temp_dir, mock_kb_util):
        """Test OutputManager initialization."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        assert output_manager.base_output_dir == temp_dir
        assert output_manager.run_id == 'test_run_id'
        assert output_manager.workspace_name == 'test_workspace'
        assert output_manager.kb_util == mock_kb_util
        assert isinstance(output_manager.artifacts, list)
        assert isinstance(output_manager.workspace_objects, list)
    
    def test_get_root_dir(self, temp_dir, mock_kb_util):
        """Test getting root output directory."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        root_dir = output_manager.get_root_dir()
        
        # The actual implementation creates outputs/{run_id}_{timestamp} structure
        assert root_dir.startswith(os.path.join(temp_dir, 'outputs', 'test_run_id_'))
        assert os.path.exists(root_dir)
    
    def test_get_analysis_dir(self, temp_dir, mock_kb_util):
        """Test getting analysis-specific directory."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        analysis_dir = output_manager.get_analysis_dir('network_analysis')
        
        # The actual implementation creates outputs/{run_id}_{timestamp}/analysis/{analysis_name} structure
        assert analysis_dir.startswith(os.path.join(temp_dir, 'outputs', 'test_run_id_'))
        assert analysis_dir.endswith(os.path.join('analysis', 'network_analysis'))
        assert os.path.exists(analysis_dir)
    
    def test_record_artifact(self, temp_dir, mock_kb_util):
        """Test recording artifacts."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        # Create a test file
        test_file = os.path.join(temp_dir, 'test_file.txt')
        with open(test_file, 'w') as f:
            f.write('test content')
        
        artifact_path = output_manager.record_artifact(
            test_file,
            'text',
            'Test artifact',
            'network_analysis'
        )
        
        assert artifact_path is not None
        assert len(output_manager.artifacts) == 1
        
        artifact = output_manager.artifacts[0]
        assert artifact.path == test_file
        assert artifact.kind == 'text'
        assert artifact.description == 'Test artifact'
        assert artifact.analysis_type == 'network_analysis'
        assert artifact.size_bytes > 0
    
    def test_write_json(self, temp_dir, mock_kb_util):
        """Test writing JSON data."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        test_data = {'test': 'data', 'number': 123}
        
        file_path = output_manager.write_json(
            'test_subdir',
            'test.json',
            test_data,
            'network_analysis',
            'Test JSON file'
        )
        
        assert file_path is not None
        assert os.path.exists(file_path)
        
        # Verify content
        import json
        with open(file_path, 'r') as f:
            loaded_data = json.load(f)
        
        assert loaded_data == test_data
        
        # Verify artifact was recorded
        assert len(output_manager.artifacts) == 1
        assert output_manager.artifacts[0].kind == 'json'
    
    def test_write_text(self, temp_dir, mock_kb_util):
        """Test writing text data."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        test_text = 'This is test text content'
        
        file_path = output_manager.write_text(
            'test_subdir',
            'test.txt',
            test_text,
            'network_analysis',
            'Test text file'
        )
        
        assert file_path is not None
        assert os.path.exists(file_path)
        
        # Verify content
        with open(file_path, 'r') as f:
            loaded_text = f.read()
        
        assert loaded_text == test_text
        
        # Verify artifact was recorded
        assert len(output_manager.artifacts) == 1
        assert output_manager.artifacts[0].kind == 'data'
    
    def test_write_csv(self, temp_dir, mock_kb_util):
        """Test writing CSV data."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        csv_data = 'name,value\nitem1,100\nitem2,200'
        
        file_path = output_manager.write_csv(
            'test_subdir',
            'test.csv',
            csv_data,
            'network_analysis',
            'Test CSV file'
        )
        
        assert file_path is not None
        assert os.path.exists(file_path)
        
        # Verify content
        with open(file_path, 'r') as f:
            loaded_csv = f.read()
        
        assert loaded_csv == csv_data
        
        # Verify artifact was recorded
        assert len(output_manager.artifacts) == 1
        assert output_manager.artifacts[0].kind == 'csv'
    
    def test_write_html(self, temp_dir, mock_kb_util):
        """Test writing HTML data."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        html_content = '<html><body><h1>Test HTML</h1></body></html>'
        
        file_path = output_manager.write_html(
            'test_subdir',
            'test.html',
            html_content,
            'network_analysis',
            'Test HTML file'
        )
        
        assert file_path is not None
        assert os.path.exists(file_path)
        
        # Verify content
        with open(file_path, 'r') as f:
            loaded_html = f.read()
        
        assert loaded_html == html_content
        
        # Verify artifact was recorded
        assert len(output_manager.artifacts) == 1
        assert output_manager.artifacts[0].kind == 'html'
    
    def test_save_analysis_output(self, temp_dir, mock_kb_util):
        """Test saving analysis output."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        analysis_result = {
            'success': True,
            'analysis_type': 'network_analysis',
            'results': {
                'network_graph': {'nodes': 3, 'edges': 2},
                'statistics': {'density': 0.5}
            },
            'metadata': {'processing_time': 1.5}
        }
        
        result = output_manager.save_analysis_output(
            'network_analysis',
            analysis_result,
            os.path.join(temp_dir, 'analysis_output')
        )
        
        assert isinstance(result, dict)
        assert 'results' in result
        assert len(result['results']) > 0
    
    def test_save_metadata(self, temp_dir, mock_kb_util):
        """Test saving metadata."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        config = {
            'input_type': 'protein_input',
            'analysis_name': 'test_analysis'
        }
        
        analyses_run = ['network_analysis', 'sequence_analysis']
        summary = 'Test analysis completed successfully'
        
        file_path = output_manager.save_metadata(config, analyses_run, summary)
        
        assert file_path is not None
        assert os.path.exists(file_path)
        
        # Verify content
        import json
        with open(file_path, 'r') as f:
            metadata = json.load(f)
        
        assert 'config' in metadata
        assert 'analyses_run' in metadata
        assert 'summary' in metadata
        assert 'timestamp' in metadata
        assert metadata['config'] == config
        assert metadata['analyses_run'] == analyses_run
        assert metadata['summary'] == summary
    
    def test_save_process_info(self, temp_dir, mock_kb_util):
        """Test saving process information."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        stages_completed = ['input_processing', 'analysis', 'output_generation']
        execution_time = 15.5
        memory_usage = {'peak_memory': '512MB', 'average_memory': '256MB'}
        
        file_path = output_manager.save_process_info(stages_completed, execution_time, memory_usage)
        
        assert file_path is not None
        assert os.path.exists(file_path)
        
        # Verify content
        import json
        with open(file_path, 'r') as f:
            process_info = json.load(f)
        
        assert 'stages_completed' in process_info
        assert 'execution_time_seconds' in process_info
        assert 'memory_usage' in process_info
        assert 'timestamp' in process_info
        assert process_info['stages_completed'] == stages_completed
        assert process_info['execution_time_seconds'] == execution_time
        assert process_info['memory_usage'] == memory_usage
    
    def test_finalize_manifest(self, temp_dir, mock_kb_util):
        """Test finalizing the manifest."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        # Add some artifacts
        output_manager.record_artifact(
            os.path.join(temp_dir, 'test1.txt'),
            'text',
            'Test artifact 1'
        )
        output_manager.record_artifact(
            os.path.join(temp_dir, 'test2.json'),
            'json',
            'Test artifact 2'
        )
        
        manifest_path = output_manager.finalize_manifest()
        
        assert manifest_path is not None
        assert os.path.exists(manifest_path)
        
        # Verify manifest content
        import json
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
        
        assert 'artifacts' in manifest
        assert 'analysis_outputs' in manifest
        assert 'run_id' in manifest
        assert len(manifest['artifacts']) == 2
    
    def test_get_artifact_summary(self, temp_dir, mock_kb_util):
        """Test getting artifact summary."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        # Add artifacts of different types
        output_manager.record_artifact(
            os.path.join(temp_dir, 'test1.txt'),
            'text',
            'Test artifact 1'
        )
        output_manager.record_artifact(
            os.path.join(temp_dir, 'test2.json'),
            'json',
            'Test artifact 2'
        )
        output_manager.record_artifact(
            os.path.join(temp_dir, 'test3.html'),
            'html',
            'Test artifact 3'
        )
        
        summary = output_manager.get_artifact_summary()
        
        assert isinstance(summary, dict)
        assert 'total_artifacts' in summary
        assert 'by_kind' in summary
        assert 'by_analysis' in summary
        assert 'total_size_bytes' in summary
        
        assert summary['total_artifacts'] == 3
        assert summary['by_kind']['text'] == 1
        assert summary['by_kind']['json'] == 1
        assert summary['by_kind']['html'] == 1
    
    def test_create_workspace_objects(self, temp_dir, mock_kb_util, mock_workspace_client):
        """Test creating workspace objects."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        analysis_results = {
            'network_analysis': {
                'success': True,
                'results': {
                    'network_graph': {'nodes': 3, 'edges': 2},
                    'statistics': {'density': 0.5}
                }
            }
        }
        
        # Mock workspace client
        output_manager.wsClient = mock_workspace_client
        
        objects = output_manager.create_workspace_objects(analysis_results)
        
        assert isinstance(objects, list)
        # Should create objects for successful analyses
        assert len(objects) >= 1
    
    def test_zip_and_upload_outputs(self, temp_dir, mock_kb_util):
        """Test zipping and uploading outputs."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        # Add some artifacts
        test_file = os.path.join(temp_dir, 'test.txt')
        with open(test_file, 'w') as f:
            f.write('test content')
        
        output_manager.record_artifact(test_file, 'text', 'Test artifact')
        
        # Mock the DataFileUtil client
        mock_dfu = Mock()
        mock_dfu.file_to_shock.return_value = {
            'shock_id': 'test_shock_id',
            'node_file_name': 'test_archive.zip',
            'shock_url': 'https://shock.test/node/test_shock_id'
        }
        
        with patch('installed_clients.DataFileUtilClient.DataFileUtil', return_value=mock_dfu):
            result = output_manager.zip_and_upload_outputs()
        
        assert isinstance(result, dict)
        assert 'shock_id' in result
        assert 'node_file_name' in result
        assert 'shock_url' in result
        assert result['shock_id'] == 'test_shock_id'
    
    def test_artifact_record_dataclass(self):
        """Test ArtifactRecord dataclass."""
        artifact = ArtifactRecord(
            path='/test/path/file.txt',
            kind='text',
            description='Test artifact',
            analysis_type='network_analysis',
            size_bytes=1024,
            checksum_md5='test_checksum'
        )
        
        assert artifact.path == '/test/path/file.txt'
        assert artifact.kind == 'text'
        assert artifact.description == 'Test artifact'
        assert artifact.analysis_type == 'network_analysis'
        assert artifact.size_bytes == 1024
        assert artifact.checksum_md5 == 'test_checksum'
    
    def test_error_handling_invalid_path(self, temp_dir, mock_kb_util):
        """Test error handling with invalid paths."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        # Test with invalid file path
        invalid_path = '/invalid/path/file.txt'
        
        # Should handle gracefully
        artifact_path = output_manager.record_artifact(
            invalid_path,
            'text',
            'Invalid artifact'
        )
        
        # Should return None or handle error gracefully
        assert artifact_path is None or isinstance(artifact_path, str)
    
    def test_directory_creation_automatic(self, temp_dir, mock_kb_util):
        """Test that directories are created automatically."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        # Try to write to a nested directory that doesn't exist
        nested_dir = 'deeply/nested/directory'
        
        file_path = output_manager.write_text(
            nested_dir,
            'test.txt',
            'test content'
        )
        
        assert file_path is not None
        assert os.path.exists(file_path)
        
        # Verify the nested directory was created
        expected_dir = os.path.join(output_manager.root_dir, nested_dir)
        assert os.path.exists(expected_dir)
    
    def test_multiple_analyses_output(self, temp_dir, mock_kb_util):
        """Test handling outputs from multiple analyses."""
        output_manager = OutputManager(temp_dir, 'test_run_id', 'test_workspace', mock_kb_util)
        
        # Save outputs from different analyses
        network_result = {
            'success': True,
            'analysis_type': 'network_analysis',
            'results': {'network_graph': {'nodes': 3, 'edges': 2}}
        }
        
        sequence_result = {
            'success': True,
            'analysis_type': 'sequence_analysis',
            'results': {'similarity_matrix': [[1.0, 0.5], [0.5, 1.0]]}
        }
        
        output1 = output_manager.save_analysis_output(
            'network_analysis',
            network_result,
            os.path.join(temp_dir, 'network_output')
        )
        
        output2 = output_manager.save_analysis_output(
            'sequence_analysis',
            sequence_result,
            os.path.join(temp_dir, 'sequence_output')
        )
        
        assert isinstance(output1, dict)
        assert isinstance(output2, dict)
        assert 'results' in output1
        assert 'results' in output2
        
        # Should have artifacts from both analyses
        assert len(output_manager.artifacts) >= 2

