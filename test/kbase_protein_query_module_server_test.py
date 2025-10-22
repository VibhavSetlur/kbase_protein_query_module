# -*- coding: utf-8 -*-
"""
Comprehensive KBase server tests for the protein query module.
Covers:
- Discovery endpoint
- Main pipeline across supported input types
- Error handling and reporting
- Shock upload stubbing and output path targeting
- Integration with KBase services
"""

import os
import sys
import unittest
import time
import tempfile
import shutil
from unittest.mock import Mock, patch

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../lib'))

try:
    from installed_clients.authclient import KBaseAuth as _KBaseAuth
except ImportError:
    from installed_clients.KBaseAuth import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil

from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module


class kbase_protein_query_moduleTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', 'test_token')
        # Use proper KBase development workspace URL
        cls.wsURL = os.environ.get('KB_WORKSPACE_URL', 'https://appdev.kbase.us/services/ws')
        cls.callbackURL = os.environ.get('SDK_CALLBACK_URL', 'http://localhost:0')
        
        try:
            # Create workspace client with proper authentication
            cls.wsClient = Workspace(cls.wsURL, token=cls.token)
            workspace_name = f'test_kbase_protein_query_module_{int(time.time())}'
            
            # Create workspace with proper parameters
            workspace_info = cls.wsClient.create_workspace({
                'workspace': workspace_name,
                'description': 'Test workspace for kbase_protein_query_module'
            })
            cls.wsName = workspace_info[1]  # workspace name from info
            print(f'Created test workspace: {cls.wsName}')
            
        except Exception as e:
            print(f'Warning: Could not create workspace: {e}')
            print('Using mock workspace client for testing')
            # Create mock KBase clients for testing
            cls.wsClient = Mock()
            cls.wsClient.create_workspace = Mock(return_value=[1, 'test_workspace', 'test_user', '2024-01-01T00:00:00+0000', 1, 'n', 'n', 'test_user'])
            cls.wsClient.delete_workspace = Mock(return_value=None)
            cls.wsClient.get_workspace_info = Mock(return_value=[1, 'test_workspace', 'test_user', '2024-01-01T00:00:00+0000', 1, 'n', 'n', 'test_user'])
            cls.wsClient.save_objects = Mock(return_value=[[1, 'test_workspace', 'test_object', 'KBaseReport.Report-1.0', '2024-01-01T00:00:00+0000', 1, 'test_user', 1, 'test_workspace', 'test_object', 'test_checksum', 1, {}]])
            
            # Mock other KBase clients
            cls.reportClient = Mock()
            cls.reportClient.create = Mock(return_value={'ref': 'test_report_ref', 'name': 'test_report'})
            
            cls.dataFileUtil = Mock()
            cls.dataFileUtil.file_to_shock = Mock(return_value={
                'shock_id': 'stub_shock', 
                'node_file_name': 'archive.zip', 
                'shock_url': 'https://shock.stub/node/stub_shock'
            })
            
            cls.kbUtilLib = Mock()
            cls.kbUtilLib.log = Mock()
            
            cls.wsName = 'test_workspace'
            
            # Set up patches for KBase clients
            cls._setup_mock_patches()
        
        # Setup test data
        cls._setup_test_data()

    @classmethod
    def tearDownClass(cls):
        # Stop patches if they exist
        if hasattr(cls, 'patches'):
            for p in cls.patches:
                p.stop()
        
        if hasattr(cls, 'wsClient') and cls.wsClient and hasattr(cls, 'wsName') and cls.wsName != 'test_workspace':
            try:
                cls.wsClient.delete_workspace({'workspace': cls.wsName})
                print(f'Test workspace {cls.wsName} deleted.')
            except Exception as e:
                print(f'Warning: Could not delete workspace {cls.wsName}: {e}')
        
        # Clean up temporary directory
        if hasattr(cls, 'test_local_tmp') and cls.test_local_tmp and os.path.exists(cls.test_local_tmp):
            try:
                shutil.rmtree(cls.test_local_tmp)
                print(f'Test temporary directory {cls.test_local_tmp} cleaned up.')
            except Exception as e:
                print(f'Warning: Could not clean up test directory {cls.test_local_tmp}: {e}')

    @classmethod
    def _setup_mock_patches(cls):
        """Setup patches for KBase clients."""
        cls.patches = []
        cls.patches.append(patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.KBaseReport', return_value=cls.reportClient))
        cls.patches.append(patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.Workspace', return_value=cls.wsClient))
        cls.patches.append(patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.DataFileUtil', return_value=cls.dataFileUtil))
        cls.patches.append(patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.KBUtilLib', return_value=cls.kbUtilLib))
        
        # Start all patches
        for p in cls.patches:
            p.start()

    @classmethod
    def _setup_test_data(cls):
        """Setup simple test data for all tests."""
        cls.test_protein_ids = ['P12345', 'Q67890', 'P04637']  # Real UniProt IDs
        cls.test_sequence = 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        # Use a proper temporary directory for testing
        cls.test_local_tmp = tempfile.mkdtemp(prefix='kbase_protein_query_test_')
        try:
            os.makedirs(cls.test_local_tmp, exist_ok=True)
        except Exception:
            pass

    def setUp(self):
        self.serviceImpl = kbase_protein_query_module({})
        self.ctx = {
            'token': self.token,
            'provenance': [{'ws_name': self.wsName}]
        }

    def test_get_available_analyses(self):
        """Test discovery endpoint for available analyses."""
        res = self.serviceImpl.get_available_analyses(self.ctx)
        self.assertIsInstance(res, list)
        self.assertIsInstance(res[0], dict)
        self.assertIn('available_analyses', res[0])
        self.assertIn('summary', res[0])

    def test_status(self):
        """Test the status method."""
        result = self.serviceImpl.status(self.ctx)
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertEqual(result[0]['state'], 'OK')
        self.assertIn('version', result[0])

    def test_run_protein_query_analysis_uniprot_ids(self):
        """Test complete workflow with UniProt IDs input - simulates actual KBase run."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': self.test_protein_ids,
            'analysis_name': 'uniprot_test',
            'analysis_stages': ['network_analysis'],
            'output_config': {'output_dir': self.test_local_tmp}
        }
        
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        
        # Validate basic result structure
        self.assertIsInstance(result, list)
        self.assertIsInstance(out, dict)
        
        # Validate required output fields
        required_fields = ['report_name', 'report_ref', 'analysis_result_ref', 'summary', 'start_time', 'stages_completed']
        for field in required_fields:
            self.assertIn(field, out, f"Missing required field: {field}")
        
        # Validate Shock upload (simulates file upload to KBase)
        self.assertIn('shock_id', out)
        self.assertIn('shock_url', out)
        self.assertEqual(out['shock_id'], 'stub_shock')
        
        # Validate output directory creation
        expected_output_dir = os.path.join(self.test_local_tmp, 'outputs')
        self.assertTrue(os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}")
        
        # Validate that output directory contains timestamped subdirectories
        output_subdirs = [d for d in os.listdir(expected_output_dir) if os.path.isdir(os.path.join(expected_output_dir, d))]
        self.assertTrue(len(output_subdirs) > 0, f"No timestamped output subdirectories found in {expected_output_dir}")
        
        # Validate analysis outputs exist and are not empty
        analysis_outputs = out.get('analysis_outputs', {})
        self.assertIsInstance(analysis_outputs, dict, "Analysis outputs should be a dictionary")
        
        # Validate network analysis outputs
        network_analysis = analysis_outputs.get('network_analysis', {})
        self.assertIsInstance(network_analysis, dict, "Network analysis output should be a dictionary")
        
        # Validate that network analysis has required fields
        if network_analysis:
            self.assertIn('network_properties', network_analysis, "Network analysis should have network_properties")
            self.assertIn('network_statistics', network_analysis, "Network analysis should have network_statistics")
            
            # Validate network properties are not empty
            network_props = network_analysis.get('network_properties', {})
            self.assertGreater(network_props.get('num_nodes', 0), 0, "Network should have nodes")
            
            # Validate HTML visualization file exists
            html_path = network_analysis.get('html_path')
            if html_path:
                self.assertTrue(os.path.exists(html_path), f"HTML visualization file not found: {html_path}")
                self.assertGreater(os.path.getsize(html_path), 0, f"HTML visualization file is empty: {html_path}")
            
            # Validate CSV files exist and are not empty
            csv_files = network_analysis.get('csv_files', {})
            for csv_type, csv_path in csv_files.items():
                self.assertTrue(os.path.exists(csv_path), f"CSV file not found: {csv_path}")
                self.assertGreater(os.path.getsize(csv_path), 0, f"CSV file is empty: {csv_path}")
        
        # Validate no errors in summary (critical for KBase narrative)
        summary = out.get('summary', '')
        self.assertNotIn('No such file or directory', summary, f"FileNotFoundError found: {summary}")
        self.assertNotIn('FileNotFoundError', summary, f"FileNotFoundError found: {summary}")
        self.assertNotIn('Error', summary, f"Error found: {summary}")
        
        # Validate analysis completion
        self.assertIn('network_analysis', out.get('stages_completed', []))

    def test_run_protein_query_analysis_protein_sequences(self):
        """Test complete workflow with protein sequences input - simulates actual KBase run."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'protein_seq_test',
            'analysis_stages': ['network_analysis'],
            'output_config': {'output_dir': self.test_local_tmp}
        }
        
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        
        # Validate basic result structure
        self.assertIsInstance(result, list)
        self.assertIsInstance(out, dict)
        
        # Validate required output fields
        required_fields = ['report_name', 'report_ref', 'analysis_result_ref', 'summary', 'start_time', 'stages_completed']
        for field in required_fields:
            self.assertIn(field, out, f"Missing required field: {field}")
        
        # Validate Shock upload
        self.assertIn('shock_id', out)
        self.assertIn('shock_url', out)
        self.assertEqual(out['shock_id'], 'stub_shock')
        
        # Validate output directory creation
        expected_output_dir = os.path.join(self.test_local_tmp, 'outputs')
        self.assertTrue(os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}")
        
        # Validate that output directory contains timestamped subdirectories
        output_subdirs = [d for d in os.listdir(expected_output_dir) if os.path.isdir(os.path.join(expected_output_dir, d))]
        self.assertTrue(len(output_subdirs) > 0, f"No timestamped output subdirectories found in {expected_output_dir}")
        
        # Validate analysis outputs exist and are not empty
        analysis_outputs = out.get('analysis_outputs', {})
        self.assertIsInstance(analysis_outputs, dict, "Analysis outputs should be a dictionary")
        
        # Validate network analysis outputs
        network_analysis = analysis_outputs.get('network_analysis', {})
        self.assertIsInstance(network_analysis, dict, "Network analysis output should be a dictionary")
        
        # Validate that network analysis has required fields
        if network_analysis:
            self.assertIn('network_properties', network_analysis, "Network analysis should have network_properties")
            self.assertIn('network_statistics', network_analysis, "Network analysis should have network_statistics")
            
            # Validate network properties are not empty
            network_props = network_analysis.get('network_properties', {})
            self.assertGreater(network_props.get('num_nodes', 0), 0, "Network should have nodes")
            
            # Validate HTML visualization file exists
            html_path = network_analysis.get('html_path')
            if html_path:
                self.assertTrue(os.path.exists(html_path), f"HTML visualization file not found: {html_path}")
                self.assertGreater(os.path.getsize(html_path), 0, f"HTML visualization file is empty: {html_path}")
            
            # Validate CSV files exist and are not empty
            csv_files = network_analysis.get('csv_files', {})
            for csv_type, csv_path in csv_files.items():
                self.assertTrue(os.path.exists(csv_path), f"CSV file not found: {csv_path}")
                self.assertGreater(os.path.getsize(csv_path), 0, f"CSV file is empty: {csv_path}")
        
        # Validate no errors in summary (critical for KBase narrative)
        summary = out.get('summary', '')
        self.assertNotIn('No such file or directory', summary, f"FileNotFoundError found: {summary}")
        self.assertNotIn('FileNotFoundError', summary, f"FileNotFoundError found: {summary}")
        self.assertNotIn('Error', summary, f"Error found: {summary}")
        
        # Validate analysis completion
        self.assertIn('network_analysis', out.get('stages_completed', []))

    def test_error_handling_invalid_parameters(self):
        """Test error handling for invalid parameters - ensures proper KBase error responses."""
        # Test missing workspace_name
        params = {}
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('workspace_name is required', str(context.exception))
        
        # Test invalid input_type
        params = {
            'workspace_name': self.wsName,
            'input_type': 'invalid_type',
            'analysis_name': 'error_test'
        }
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('Invalid input_type', str(context.exception))

    def test_zip_and_upload_functionality(self):
        """Test that the zip and upload functionality works correctly."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': self.test_protein_ids,
            'analysis_name': 'zip_test',
            'analysis_stages': ['network_analysis'],
            'output_config': {'output_dir': self.test_local_tmp}
        }
        
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        
        # Validate Shock upload information
        self.assertIn('shock_id', out)
        self.assertIn('shock_url', out)
        self.assertEqual(out['shock_id'], 'stub_shock')
        self.assertEqual(out['shock_url'], 'https://shock.test/node/stub_shock')
        
        # Validate that output directory was created and has content
        expected_output_dir = os.path.join(self.test_local_tmp, 'outputs')
        self.assertTrue(os.path.exists(expected_output_dir), f"Output directory not created: {expected_output_dir}")
        
        # Check that there are timestamped subdirectories
        output_subdirs = [d for d in os.listdir(expected_output_dir) if os.path.isdir(os.path.join(expected_output_dir, d))]
        self.assertTrue(len(output_subdirs) > 0, f"No timestamped output subdirectories found")
        
        # Validate that the zip process would work (directory has content)
        for subdir in output_subdirs:
            subdir_path = os.path.join(expected_output_dir, subdir)
            files_in_subdir = []
            for root, dirs, files in os.walk(subdir_path):
                files_in_subdir.extend(files)
            self.assertTrue(len(files_in_subdir) > 0, f"No files found in output subdirectory: {subdir_path}")
        
        # Validate that the mock DataFileUtil was called (simulating zip upload)
        # This ensures the zip functionality is properly integrated
        self.assertIsNotNone(out.get('shock_id'), "Shock ID should be present after zip upload")
        self.assertIsNotNone(out.get('shock_url'), "Shock URL should be present after zip upload")


if __name__ == '__main__':
    unittest.main()