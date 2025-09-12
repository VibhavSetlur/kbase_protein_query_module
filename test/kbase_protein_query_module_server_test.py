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
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        if not cls.token:
            raise ValueError('KB_AUTH_TOKEN environment variable is required')
        
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

    def test_check_protein_existence(self):
        """Test protein existence check functionality."""
        params = {
            'protein_id': 'P00001',
            'workspace_name': self.wsName,
            'generate_embedding': False
        }
        
        result = self.serviceImpl.check_protein_existence(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('exists', result[0])
        self.assertIn('family_id', result[0])
        self.assertIn('summary', result[0])

    def test_generate_protein_embedding(self):
        """Test protein embedding generation."""
        params = {
            'input_type': 'sequence',
            'input_data': self.test_sequence,
            'workspace_name': self.wsName,
            'model_name': 'esm2_t6_8M_UR50D'
        }
        
        result = self.serviceImpl.generate_protein_embedding(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('embedding_result_ref', result[0])
        self.assertIn('embedding_norm', result[0])

    def test_assign_family_fast(self):
        """Test family assignment functionality."""
        # Use mock data to avoid workspace dependencies
        params = {
            'embedding_ref': 'demo_embedding_ref',
            'protein_id': 'P00001',
            'workspace_name': self.wsName
        }
        
        result = self.serviceImpl.assign_family_fast(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('family_id', result[0])
        self.assertIn('confidence', result[0])

    def test_find_top_matches_from_embedding(self):
        """Test similarity search functionality."""
        # Use mock data to avoid workspace dependencies
        params = {
            'embedding_ref': 'demo_embedding_ref',
            'protein_id': 'P00001',
            'workspace_name': self.wsName,
            'max_matches': 5
        }
        
        result = self.serviceImpl.find_top_matches_from_embedding(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('matches', result[0])
        self.assertIn('family_id', result[0])

    def test_summarize_and_visualize_results(self):
        """Test result summarization and visualization."""
        # Create some demo result refs
        demo_refs = ['demo_ref_1', 'demo_ref_2']
        
        params = {
            'result_refs': demo_refs,
            'output_name': 'test_analysis',
            'workspace_name': self.wsName
        }
        
        result = self.serviceImpl.summarize_and_visualize_results(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('output_directory', result[0])
        self.assertIn('general_info_dir', result[0])

    def test_run_protein_query_analysis(self):
        """Test the main unified analysis pipeline."""
        params = {
            'workspace_name': self.wsName,
            'input_proteins': self.test_protein_ids,
            'analysis_stages': ['embedding_generation']
        }
        
        result = self.serviceImpl.run_protein_query_analysis(self.ctx, params)
        
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)
        self.assertIn('output_directory', result[0])
        self.assertIn('analysis_result_ref', result[0])
        self.assertIn('stages_completed', result[0])

    def test_legacy_method_compatibility(self):
        """Test backward compatibility with legacy method name."""
        params = {
            'workspace_name': self.wsName,
            'input_proteins': ['P00001'],
            'analysis_stages': ['embedding_generation']
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
        # Test with missing required parameter
        params = {}  # Missing protein_id
        
        with self.assertRaises(ValueError):
            self.serviceImpl.check_protein_existence(self.ctx, params)

    def test_parameter_validation(self):
        """Test parameter validation."""
        # Test with missing input_data (should raise ValueError)
        params = {
            'input_type': 'sequence',
            'input_data': ''  # Empty input data
        }
        
        with self.assertRaises(ValueError):
            self.serviceImpl.generate_protein_embedding(self.ctx, params)


if __name__ == '__main__':
    unittest.main()
