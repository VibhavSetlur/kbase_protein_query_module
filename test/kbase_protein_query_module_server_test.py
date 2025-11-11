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
        
        # Do NOT mock WorkflowOrchestrator to allow real analysis and outputs
        # No PipelineConfig class - using simple dict configs now
        
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
        
        # No orchestrator mock; real workflow will run
        self.ctx = {
            'token': self.token,
            'provenance': [{'ws_name': self.wsName}]
        }


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
        zip_path = call_args['file_path']
        self.assertTrue(str(zip_path).endswith('.zip'))
        
        # Verify that the zip file was actually created
        self.assertTrue(os.path.exists(zip_path), f"Output zip file does not exist: {zip_path}")
        self.assertGreater(os.path.getsize(zip_path), 0, "Output zip file is empty")
        
        # Verify that output directory structure exists
        output_dir = os.path.join(self.test_local_tmp, 'outputs')
        self.assertTrue(os.path.exists(output_dir), f"Output directory does not exist: {output_dir}")
        
        # CRITICAL: Verify that analysis actually ran and created outputs
        # Check if analysis output directory exists - if not, analysis didn't run
        analysis_dir = os.path.join(output_dir, 'analysis', 'network_analysis')
        
        if not os.path.exists(analysis_dir):
            # Analysis directory doesn't exist - check metadata to see what happened
            metadata_dir = os.path.join(output_dir, 'metadata')
            if os.path.exists(metadata_dir):
                import json
                metadata_file = os.path.join(metadata_dir, 'final_results.json')
                if os.path.exists(metadata_file):
                    with open(metadata_file, 'r') as f:
                        workflow_result = json.load(f)
                    analysis_results = workflow_result.get('analysis_results', {})
                    if 'network_analysis' in analysis_results:
                        na_result = analysis_results['network_analysis']
                        if na_result is None or (isinstance(na_result, dict) and na_result.get('success') is False):
                            error_msg = na_result.get('error_message', 'Analysis not available') if isinstance(na_result, dict) else 'Analysis returned None'
                            self.fail(f"❌ Network analysis was requested but did NOT run: {error_msg}\n"
                                     f"Required dependencies: transformers (for embeddings), networkx and sklearn (for network analysis)\n"
                                     f"This test REQUIRES these dependencies to verify real output generation.\n"
                                     f"Install dependencies: pip install torch transformers networkx scikit-learn plotly")
            
            # If we can't determine why, fail with clear message
            self.fail(f"❌ Network analysis output directory does not exist: {analysis_dir}\n"
                     f"This means network analysis did NOT run.\n"
                     f"Required dependencies: transformers, networkx, sklearn\n"
                     f"Install: pip install torch transformers networkx scikit-learn plotly\n"
                     f"Coverage report shows network_analysis.py has only 20% coverage - analysis code is NOT executing!")
        
        # Analysis directory exists - verify it contains real outputs
        results_json = os.path.join(analysis_dir, 'results.json')
        self.assertTrue(os.path.exists(results_json), f"Results JSON does not exist: {results_json}")
        
        import json
        with open(results_json, 'r') as f:
            results = json.load(f)
        self.assertIn('success', results, "Results JSON missing 'success' field")
        self.assertTrue(results.get('success'), 
                       f"❌ Network analysis failed: {results.get('error_message', 'Unknown error')}")
        
        # Verify output files were created
        output_files = results.get('output_files', [])
        self.assertGreater(len(output_files), 0, 
                          "❌ No output files were created by network analysis. Expected: HTML visualization, TSV files.")
        
        # Verify each output file exists and is not empty
        for file_path in output_files:
            if isinstance(file_path, str):
                # Try absolute path first
                full_path = file_path if os.path.isabs(file_path) else os.path.join(output_dir, file_path)
                if not os.path.exists(full_path):
                    # Try relative to analysis_dir
                    full_path = os.path.join(analysis_dir, os.path.basename(file_path))
                self.assertTrue(os.path.exists(full_path), 
                              f"❌ Output file does not exist: {file_path} (tried: {full_path})")
                file_size = os.path.getsize(full_path)
                self.assertGreater(file_size, 0, 
                                 f"❌ Output file is empty (0 bytes): {full_path}")
        
        # Verify specific file types were created
        file_names = [os.path.basename(f) for f in output_files if isinstance(f, str)]
        has_html = any('html' in f.lower() for f in file_names)
        has_tsv = any('tsv' in f.lower() for f in file_names)
        
        if not has_html:
            self.fail("❌ Network analysis did not create HTML visualization file. Expected: network_visualization_*.html")
        if not has_tsv:
            self.fail("❌ Network analysis did not create TSV files. Expected: network_statistics_*.tsv and/or network_edges_*.tsv")
        
        # If we get here, analysis ran and created real outputs!
        print(f"✅ Network analysis ran successfully and created {len(output_files)} output files")
        print(f"   Files: {', '.join(file_names)}")

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
        
        # Verify DataFileUtil was called correctly
        self.dataFileUtil.file_to_shock.assert_called_once()
        call_args = self.dataFileUtil.file_to_shock.call_args[0][0]
        self.assertEqual(call_args.get('make_handle'), 1)
        self.assertIn('file_path', call_args)
        zip_path = call_args['file_path']
        self.assertTrue(str(zip_path).endswith('.zip'))
        
        # Verify that the zip file was actually created
        self.assertTrue(os.path.exists(zip_path), f"Output zip file does not exist: {zip_path}")
        self.assertGreater(os.path.getsize(zip_path), 0, "Output zip file is empty")
        
        # Verify that output directory structure exists
        output_dir = os.path.join(self.test_local_tmp, 'outputs')
        self.assertTrue(os.path.exists(output_dir), f"Output directory does not exist: {output_dir}")
        
        # CRITICAL: Verify that analysis actually ran and created outputs
        # Check if analysis output directory exists - if not, analysis didn't run
        analysis_dir = os.path.join(output_dir, 'analysis', 'network_analysis')
        
        if not os.path.exists(analysis_dir):
            # Analysis directory doesn't exist - check metadata to see what happened
            metadata_dir = os.path.join(output_dir, 'metadata')
            if os.path.exists(metadata_dir):
                import json
                metadata_file = os.path.join(metadata_dir, 'final_results.json')
                if os.path.exists(metadata_file):
                    with open(metadata_file, 'r') as f:
                        workflow_result = json.load(f)
                    analysis_results = workflow_result.get('analysis_results', {})
                    if 'network_analysis' in analysis_results:
                        na_result = analysis_results['network_analysis']
                        if na_result is None or (isinstance(na_result, dict) and na_result.get('success') is False):
                            error_msg = na_result.get('error_message', 'Analysis not available') if isinstance(na_result, dict) else 'Analysis returned None'
                            self.fail(f"❌ Network analysis was requested but did NOT run: {error_msg}\n"
                                     f"Required dependencies: transformers (for embeddings), networkx and sklearn (for network analysis)\n"
                                     f"This test REQUIRES these dependencies to verify real output generation.\n"
                                     f"Install dependencies: pip install torch transformers networkx scikit-learn plotly")
            
            # If we can't determine why, fail with clear message
            self.fail(f"❌ Network analysis output directory does not exist: {analysis_dir}\n"
                     f"This means network analysis did NOT run.\n"
                     f"Required dependencies: transformers, networkx, sklearn\n"
                     f"Install: pip install torch transformers networkx scikit-learn plotly\n"
                     f"Coverage report shows network_analysis.py has only 20% coverage - analysis code is NOT executing!")
        
        # Analysis directory exists - verify it contains real outputs
        results_json = os.path.join(analysis_dir, 'results.json')
        self.assertTrue(os.path.exists(results_json), f"Results JSON does not exist: {results_json}")
        
        import json
        with open(results_json, 'r') as f:
            results = json.load(f)
        self.assertIn('success', results, "Results JSON missing 'success' field")
        self.assertTrue(results.get('success'), 
                       f"❌ Network analysis failed: {results.get('error_message', 'Unknown error')}")
        
        # Verify output files were created
        output_files = results.get('output_files', [])
        self.assertGreater(len(output_files), 0, 
                          "❌ No output files were created by network analysis. Expected: HTML visualization, TSV files.")
        
        # Verify each output file exists and is not empty
        for file_path in output_files:
            if isinstance(file_path, str):
                # Try absolute path first
                full_path = file_path if os.path.isabs(file_path) else os.path.join(output_dir, file_path)
                if not os.path.exists(full_path):
                    # Try relative to analysis_dir
                    full_path = os.path.join(analysis_dir, os.path.basename(file_path))
                self.assertTrue(os.path.exists(full_path), 
                              f"❌ Output file does not exist: {file_path} (tried: {full_path})")
                file_size = os.path.getsize(full_path)
                self.assertGreater(file_size, 0, 
                                 f"❌ Output file is empty (0 bytes): {full_path}")
        
        # Verify specific file types were created
        file_names = [os.path.basename(f) for f in output_files if isinstance(f, str)]
        has_html = any('html' in f.lower() for f in file_names)
        has_tsv = any('tsv' in f.lower() for f in file_names)
        
        if not has_html:
            self.fail("❌ Network analysis did not create HTML visualization file. Expected: network_visualization_*.html")
        if not has_tsv:
            self.fail("❌ Network analysis did not create TSV files. Expected: network_statistics_*.tsv and/or network_edges_*.tsv")
        
        # If we get here, analysis ran and created real outputs!
        print(f"✅ Network analysis ran successfully and created {len(output_files)} output files")
        print(f"   Files: {', '.join(file_names)}")




if __name__ == '__main__':
    unittest.main()