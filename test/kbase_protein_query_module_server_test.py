# -*- coding: utf-8 -*-
"""
Comprehensive KBase server tests for the upgraded backend.
Covers:
- Discovery endpoint
- Main pipeline across supported input types
- Error handling and reporting
- Shock upload stubbing and output path targeting tmp@test_local when specified
- Integration with new test suite infrastructure
"""

import os
import sys
import unittest
import time
import pytest
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


@pytest.mark.kbase
@pytest.mark.server
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
            cls.reportClient.create_extended_report = Mock(return_value={'ref': 'test_report_ref', 'name': 'test_report'})
            
            cls.dataFileUtil = Mock()
            cls.dataFileUtil.save_objects = Mock(return_value=[{'ref': 'test_data_ref'}])
            
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
        cls.test_protein_ids = ['P00001', 'P00002', 'P00003']
        cls.test_sequence = 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'
        # Use test-local tmp dir to mimic container behavior
        cls.test_local_tmp = '/kb/module/work/tmp@test_local/'
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
        # Stub DFU for Shock upload for all tests
        class _DFUStub:
            def __init__(self, *args, **kwargs):
                pass
            def file_to_shock(self, params):
                # Ensure zip of an existing directory
                path = params.get('file_path')
                assert path and os.path.isdir(path)
                return {'shock_id': 'stub_shock', 'node_file_name': 'archive.zip', 'shock_url': 'https://shock.stub/node/stub_shock'}
        import sys as _sys
        _sys.modules['installed_clients.DataFileUtilClient'] = type('m', (), {'DataFileUtil': _DFUStub})

    @pytest.mark.kbase
    def test_get_available_analyses(self):
        """Test discovery endpoint for available analyses."""
        res = self.serviceImpl.get_available_analyses(self.ctx)
        self.assertIsInstance(res, list)
        self.assertIsInstance(res[0], dict)
        self.assertIn('available_analyses', res[0])
        self.assertIn('summary', res[0])

    @pytest.mark.integration
    @pytest.mark.kbase
    def test_run_pipeline_uniprot_ids_with_tmp_output_and_shock(self):
        """End-to-end: uniprot_ids input, Shock stubbed, outputs under tmp@test_local when requested."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': self.test_protein_ids,
            'analysis_name': 'test_analysis',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        self.assertIsInstance(result, list)
        self.assertIsInstance(out, dict)
        # Report fields
        self.assertTrue(len(out.get('report_name', '')))
        self.assertTrue(len(out.get('report_ref', '')))
        # Shock fields present from stub
        self.assertEqual(out.get('shock_id'), 'stub_shock')
        # Stages list (may be empty if no analyses enabled but key should exist)
        self.assertIn('stages_completed', out)

    def test_run_pipeline_protein_sequences(self):
        """End-to-end: protein_input with direct sequences."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'seq_analysis',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        self.assertIsInstance(out, dict)
        self.assertIn('report_name', out)
        self.assertIn('report_ref', out)

    def test_run_pipeline_protein_sequence_string(self):
        """End-to-end: protein_input with single string sequence (like MMMM)."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': 'MMMM',
            'analysis_name': 'string_seq_analysis',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        self.assertIsInstance(out, dict)
        self.assertIn('report_name', out)
        self.assertIn('report_ref', out)

    def test_uniprot_ids_input(self):
        """End-to-end: uniprot_ids input."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': self.test_protein_ids,
            'analysis_name': 'uniprot_test',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        out = result[0]
        self.assertIsInstance(out, dict)
        self.assertIn('report_name', out)

    @pytest.mark.kbase
    def test_status(self):
        """Test the status method."""
        result = self.serviceImpl.status(self.ctx)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertEqual(result[0]['state'], 'OK')
        self.assertIn('version', result[0])

    def test_workspace_connection(self):
        """Test workspace connectivity."""
        try:
            ws_info = self.wsClient.get_workspace_info({'workspace': self.wsName})
            self.assertIsNotNone(ws_info)
        except Exception as e:
            # Fallback to list_workspace_info
            try:
                ws_list = self.wsClient.list_workspace_info({'perm': 'a'})
                self.assertIsNotNone(ws_list)
            except Exception as e2:
                self.fail(f"Workspace connection failed: {e2}")

    def test_error_handling(self):
        """Test error handling for invalid parameters."""
        # Missing workspace_name and required fields - should raise ValueError
        params = {}
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        
        # Check that the error message is appropriate
        self.assertIn('workspace_name is required', str(context.exception))

    def test_parameter_validation(self):
        """Test parameter validation with unsupported input_type and missing fields."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'unsupported_type',
            'analysis_name': 'x'
        }
        # Now validation errors properly raise exceptions instead of returning error results
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        
        # Check that the error message is appropriate
        self.assertIn('Invalid input_type', str(context.exception))
        self.assertIn('unsupported_type', str(context.exception))

    def test_all_supported_input_types(self):
        """Test supported input types: protein_input, uniprot_ids."""
        # protein_input
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'pi',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        self.assertIsInstance(self.serviceImpl.run_protein_query_analysis(self.ctx, params), list)

        # uniprot_ids
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': self.test_protein_ids,
            'analysis_name': 'ui',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        self.assertIsInstance(self.serviceImpl.run_protein_query_analysis(self.ctx, params), list)

    def test_individual_methods_absent(self):
        """Deprecated endpoints removed; assert absence for unified API."""
        self.assertFalse(hasattr(self.serviceImpl, 'check_protein_existence'))
        self.assertFalse(hasattr(self.serviceImpl, 'generate_protein_embedding'))
        self.assertFalse(hasattr(self.serviceImpl, 'assign_family_fast'))
        self.assertFalse(hasattr(self.serviceImpl, 'find_top_matches_from_embedding'))
        self.assertFalse(hasattr(self.serviceImpl, 'summarize_and_visualize_results'))

    def test_output_minimal_contract(self):
        """Validate minimal output contract for upgraded API."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'output_test',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        output = result[0]
        # Required minimal fields
        for field in ['report_name', 'report_ref', 'analysis_result_ref', 'summary', 'start_time', 'stages_completed']:
            self.assertIn(field, output)
        # Optional Shock fields (present due to stub)
        self.assertEqual(output.get('shock_id'), 'stub_shock')

    def test_error_scenarios(self):
        """Test various error scenarios."""
        # Invalid input type - should raise ValueError
        params = {
            'workspace_name': self.wsName,
            'input_type': 'invalid_type',
            'analysis_name': 'error_test'
        }
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('Invalid input_type', str(context.exception))
        
        # Missing workspace - should raise ValueError
        params = {
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'no_workspace'
        }
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('workspace_name is required', str(context.exception))
        
        # Empty analysis name - should raise ValueError
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': ''
        }
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('analysis_name cannot be empty', str(context.exception))

    def test_legacy_compatibility(self):
        """Legacy alias removed in revamped API."""
        self.assertFalse(hasattr(self.serviceImpl, 'run_kbase_protein_query_module'))

    def test_performance_smoke(self):
        """Basic smoke: pipeline completes within reasonable time."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence] * 5,
            'analysis_name': 'perf_test',
            'output_config': {'output_dir': self.test_local_tmp}
        }
        start_time = time.time()
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        execution_time = time.time() - start_time
        self.assertIsInstance(result, list)
        self.assertLess(execution_time, 60)



if __name__ == '__main__':
    unittest.main()
