# -*- coding: utf-8 -*-
"""
KBase server tests for the protein query module.
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

from installed_clients.WorkspaceClient import Workspace
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil

from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module


class kbase_protein_query_moduleTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', 'test_token')
        cls.wsURL = os.environ.get('KB_WORKSPACE_URL', 'https://appdev.kbase.us/services/ws')
        cls.callbackURL = os.environ.get('SDK_CALLBACK_URL', 'http://localhost:0')
        
        # Create mock KBase clients for testing
        cls.wsClient = Mock()
        cls.wsClient.create_workspace = Mock(return_value=[1, 'test_workspace', 'test_user', '2024-01-01T00:00:00+0000', 1, 'n', 'n', 'test_user'])
        cls.wsClient.delete_workspace = Mock(return_value=None)
        
        cls.reportClient = Mock()
        cls.reportClient.create_extended_report = Mock(return_value={'name': 'test_report', 'ref': 'test_report_ref'})
        
        cls.dataFileUtil = Mock()
        cls.dataFileUtil.file_to_shock = Mock(return_value={
            'shock_id': 'test_shock_id_12345', 
            'shock_url': 'https://shock.test/node/test_shock_id_12345',
            'handle': {
                'file_name': 'output.zip',
                'id': 'test_handle_id'
            }
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
        # Note: DataFileUtil is now mocked at the instance level in setUp() to avoid cross-test contamination
        cls.patches.append(patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.KBUtilLib', return_value=cls.kbUtilLib))
        
        # Mock WorkflowOrchestrator
        cls.workflow_mock = Mock()
        cls.workflow_mock.run_workflow.return_value = Mock(
            final_output={
                'summary': 'Test analysis completed successfully',
                'protein_count': 3
            },
            analyses_completed=['network_analysis']
        )
        cls.patches.append(patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.WorkflowOrchestrator', return_value=cls.workflow_mock))
        
        # Mock PipelineConfig
        cls.pipeline_config_mock = Mock()
        cls.patches.append(patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.PipelineConfig', return_value=cls.pipeline_config_mock))
        
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
        # Create a fresh DataFileUtil mock for each test to avoid cross-test contamination
        self.dataFileUtil = Mock()
        self.dataFileUtil.file_to_shock = Mock(return_value={
            'shock_id': 'test_shock_id_12345', 
            'shock_url': 'https://shock.test/node/test_shock_id_12345',
            'handle': {
                'file_name': 'output.zip',
                'id': 'test_handle_id'
            }
        })
        # Apply the fresh mock to the service instance - this should override any class-level patches
        self.serviceImpl.dfu = self.dataFileUtil
        
        # Reset workflow mock
        self.workflow_mock.reset_mock()
        self.workflow_mock.run_workflow.return_value = Mock(
            final_output={
                'summary': 'Test analysis completed successfully',
                'protein_count': 3
            },
            analyses_completed=['network_analysis']
        )
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

    def test_mock_setup_verification(self):
        """Test that mocks are properly set up for each test."""
        # Verify DataFileUtil mock is working
        result = self.dataFileUtil.file_to_shock({'file_path': '/test.zip', 'make_handle': 1})
        self.assertEqual(result['shock_id'], 'test_shock_id_12345')
        self.assertEqual(result['shock_url'], 'https://shock.test/node/test_shock_id_12345')
        
        # Verify no side effects are set
        self.assertIsNone(self.dataFileUtil.file_to_shock.side_effect)

    def test_run_protein_query_analysis_uniprot_ids(self):
        """Test workflow with UniProt IDs input."""
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
        
        # Validate required output fields (only report references)
        required_fields = ['report_name', 'report_ref']
        for field in required_fields:
            self.assertIn(field, out, f"Missing required field: {field}")
        
        # Validate report references
        self.assertEqual(out['report_name'], 'test_report')
        self.assertEqual(out['report_ref'], 'test_report_ref')
        
        # Verify DataFileUtil was called correctly
        self.dataFileUtil.file_to_shock.assert_called_once()
        call_args = self.dataFileUtil.file_to_shock.call_args[0][0]
        self.assertEqual(call_args.get('make_handle'), 1)
        self.assertIn('file_path', call_args)
        self.assertTrue(str(call_args['file_path']).endswith('.zip'))

    def test_run_protein_query_analysis_protein_sequences(self):
        """Test workflow with protein sequences input."""
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
        
        # Validate required output fields (only report references)
        required_fields = ['report_name', 'report_ref']
        for field in required_fields:
            self.assertIn(field, out, f"Missing required field: {field}")
        
        # Validate report references
        self.assertEqual(out['report_name'], 'test_report')
        self.assertEqual(out['report_ref'], 'test_report_ref')

    def test_error_handling_invalid_parameters(self):
        """Test error handling for invalid parameters."""
        # Test missing workspace_name
        params = {}
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('workspace_name is required', str(context.exception))

    def test_zip_and_upload_functionality(self):
        """Test zip and upload functionality."""
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
        
        # Validate report references
        self.assertIn('report_name', out)
        self.assertIn('report_ref', out)
        self.assertEqual(out['report_name'], 'test_report')
        self.assertEqual(out['report_ref'], 'test_report_ref')
        
        # Verify DataFileUtil was called with correct parameters
        self.dataFileUtil.file_to_shock.assert_called_once()
        call_args = self.dataFileUtil.file_to_shock.call_args[0][0]
        self.assertEqual(call_args.get('make_handle'), 1)
        self.assertIn('file_path', call_args)
        self.assertTrue(str(call_args['file_path']).endswith('.zip'))

    def test_shock_upload_failure_handling(self):
        """Test that Shock upload failures are handled properly."""
        # Mock DataFileUtil to raise an exception
        self.dataFileUtil.file_to_shock.side_effect = Exception("Shock upload failed")
        
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': self.test_protein_ids,
            'analysis_name': 'shock_failure_test',
            'analysis_stages': ['network_analysis'],
            'output_config': {'output_dir': self.test_local_tmp}
        }
        
        # The method should now handle Shock upload failures gracefully
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        
        # Validate that the analysis still completes successfully
        self.assertIsInstance(result, list)
        self.assertIsInstance(out, dict)
        
        # Validate required output fields are present (only report references)
        required_fields = ['report_name', 'report_ref']
        for field in required_fields:
            self.assertIn(field, out, f"Missing required field: {field}")
        
        # Validate report references
        self.assertEqual(out['report_name'], 'test_report')
        self.assertEqual(out['report_ref'], 'test_report_ref')
        
        # Verify DataFileUtil was called and failed as expected
        self.dataFileUtil.file_to_shock.assert_called_once()



if __name__ == '__main__':
    unittest.main()