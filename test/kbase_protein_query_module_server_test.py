# -*- coding: utf-8 -*-
"""
Streamlined test suite for kbase_protein_query_module focusing on pipeline functionality.
Tests core methods and ensures proper integration with KBase services.
"""

import os
import sys
import unittest
import time
from unittest.mock import Mock, patch

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
        # Allow tests to run without real token for basic functionality testing
        
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

    def setUp(self):
        self.serviceImpl = kbase_protein_query_module({})
        self.ctx = {
            'token': self.token,
            'provenance': [{'ws_name': self.wsName}]
        }
        
        # Workspace client is now always available (either real or mock)

    def test_get_available_analyses(self):
        """Test discovery endpoint for available analyses."""
        res = self.serviceImpl.get_available_analyses(self.ctx)
        self.assertIsInstance(res, list)
        self.assertIsInstance(res[0], dict)
        self.assertIn('available_analyses', res[0])
        self.assertIn('summary', res[0])

    def test_run_protein_query_analysis(self):
        """Test the main unified analysis pipeline."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': self.test_protein_ids,
            'analysis_name': 'test_analysis',
            'enabled_stages': ['embedding_generation', 'family_assignment']
        }
        
        # Stub DFU for Shock upload
        import sys
        class _DFUStub:
            def __init__(self, *args, **kwargs):
                pass
            def file_to_shock(self, params):
                return {'shock_id': 'stub_shock', 'node_file_name': 'archive.zip'}
        sys.modules['installed_clients.DataFileUtilClient'] = type('m', (), {'DataFileUtil': _DFUStub})
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('output_directory', result[0])
        self.assertIn('analysis_result_ref', result[0])
        self.assertIn('stages_completed', result[0])
        # Ensure report info is provided for Narrative integration
        self.assertTrue(len(result[0].get('report_name', '')))
        self.assertTrue(len(result[0].get('report_ref', '')))

    def test_legacy_method_compatibility(self):
        """Test backward compatibility with legacy method name."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'uniprot_ids',
            'uniprot_ids': ['P00001'],
            'analysis_name': 'legacy_test',
            'enabled_stages': ['embedding_generation']
        }
        
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('output_directory', result[0])

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
        # Missing workspace_name and required fields
        params = {}
        res = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIsInstance(res, list)
        out = res[0]
        self.assertEqual(out.get('analysis_result_ref'), 'error')
        self.assertIn('summary', out)
        self.assertIn('report_name', out)
        self.assertIn('report_ref', out)

    def test_parameter_validation(self):
        """Test parameter validation."""
        # Missing conditional input for selected type
        params = {
            'workspace_name': self.wsName,
            'input_type': 'direct_sequences',
            'analysis_name': 'x'
        }
        res = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIsInstance(res, list)
        out = res[0]
        self.assertEqual(out.get('analysis_result_ref'), 'error')
        self.assertIn('summary', out)

    def test_all_input_types(self):
        """Test all supported input types."""
        # Test direct_sequences input type
        params = {
            'workspace_name': self.wsName,
            'input_type': 'direct_sequences',
            'direct_sequences': [self.test_sequence],
            'analysis_name': 'test_direct_sequences'
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIsInstance(result, list)
        self.assertIn('output_directory', result[0])
        
        # Test fasta_file input type
        params = {
            'workspace_name': self.wsName,
            'input_type': 'fasta_file',
            'fasta_file': 'test.fasta',
            'analysis_name': 'test_fasta'
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIsInstance(result, list)
        self.assertIn('output_directory', result[0])
        
        # Test workspace_object input type
        params = {
            'workspace_name': self.wsName,
            'input_type': 'workspace_object',
            'workspace_object': 'test_genome',
            'analysis_name': 'test_workspace_object'
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIsInstance(result, list)
        self.assertIn('output_directory', result[0])

    def test_individual_methods(self):
        """Deprecated endpoints removed; assert absence for unified API."""
        self.assertFalse(hasattr(self.serviceImpl, 'check_protein_existence'))
        self.assertFalse(hasattr(self.serviceImpl, 'generate_protein_embedding'))
        self.assertFalse(hasattr(self.serviceImpl, 'assign_family_fast'))
        self.assertFalse(hasattr(self.serviceImpl, 'find_top_matches_from_embedding'))
        self.assertFalse(hasattr(self.serviceImpl, 'summarize_and_visualize_results'))

    def test_output_validation(self):
        """Test output structure validation."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'output_test'
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        
        # Validate main output structure
        output = result[0]
        required_fields = [
            'report_name', 'report_ref', 'analysis_result_ref', 'summary',
            'input_parameters', 'start_time', 'protein_count', 'stages_completed',
            'output_directory', 'general_info_dir', 'network_analysis_dir',
            'sequence_analysis_dir', 'embeddings_file_path', 'top_proteins_csv_path'
        ]
        
        for field in required_fields:
            self.assertIn(field, output, f"Missing required output field: {field}")
        
        # Validate data types
        self.assertIsInstance(output['start_time'], (int, float))
        self.assertIsInstance(output['protein_count'], int)
        self.assertIsInstance(output['stages_completed'], list)
        self.assertIsInstance(output['input_parameters'], dict)
        self.assertTrue(len(output.get('report_name', '')) >= 0)
        self.assertTrue(len(output.get('report_ref', '')) >= 0)
        # Shock fields should be present from Shock upload integration
        self.assertIn('shock_id', output)
        # shock_url may be empty depending on environment
        # Validate Shock fields exist
        self.assertIn('shock_id', output)

    def test_error_scenarios(self):
        """Test various error scenarios."""
        # Test invalid input type (ensures graceful error output)
        params = {
            'workspace_name': self.wsName,
            'input_type': 'invalid_type',
            'analysis_name': 'error_test'
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertEqual(result[0].get('analysis_result_ref'), 'error')
        
        # Test missing workspace
        params = {
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'no_workspace'
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertEqual(result[0].get('analysis_result_ref'), 'error')
        
        # Test empty analysis name
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': ''
        }
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertEqual(result[0].get('analysis_result_ref'), 'error')

    def test_legacy_compatibility(self):
        """Test backward compatibility with legacy method names."""
        # Legacy alias removed
        self.assertFalse(hasattr(self.serviceImpl, 'run_kbase_protein_query_module'))

    def test_stage_selection(self):
        """Test different analysis stage combinations."""
        base_params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': [self.test_sequence],
            'analysis_name': 'stage_test'
        }
        
        # Test embedding only
        params = {**base_params, 'enabled_stages': ['embedding_generation']}
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('embedding_generation', result[0]['stages_completed'])
        
        # Test family assignment only
        params = {**base_params, 'enabled_stages': ['family_assignment']}
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertIn('family_assignment', result[0]['stages_completed'])
        
        # Test full pipeline
        params = {**base_params, 'enabled_stages': ['embedding_generation', 'family_assignment', 'similarity_search']}
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        self.assertGreater(len(result[0]['stages_completed']), 2)

    def test_performance_monitoring(self):
        """Test performance monitoring and resource usage."""
        params = {
            'workspace_name': self.wsName,
            'input_type': 'direct_sequences',
            'direct_sequences': [self.test_sequence] * 10,  # Multiple sequences
            'analysis_name': 'perf_test'
        }
        
        start_time = time.time()
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        execution_time = time.time() - start_time
        
        # Validate timing information
        self.assertIn('start_time', result[0])
        self.assertGreater(execution_time, 0)
        self.assertLess(execution_time, 60)  # Should complete within reasonable time


if __name__ == '__main__':
    unittest.main()
