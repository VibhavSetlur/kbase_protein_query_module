import os
import time
import unittest
from unittest.mock import Mock, patch

from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module


class PublicApiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', 'test_token')
        cls.wsName = 'test_workspace'

        # Mocks for KBase clients
        cls.wsClient = Mock()
        # Minimal responses to satisfy save/get operations used by Impl
        cls.wsClient.create_workspace = Mock(return_value=[1, cls.wsName, 'test_user', '2024-01-01T00:00:00+0000', 1, 'n', 'n', 'test_user'])
        cls.wsClient.delete_workspace = Mock(return_value=None)
        cls.wsClient.get_workspace_info = Mock(return_value=[1, cls.wsName, 'test_user', '2024-01-01T00:00:00+0000', 1, 'n', 'n', 'test_user'])
        cls.wsClient.save_objects = Mock(return_value=[[1, cls.wsName, 'obj', 'KBaseReport.Report-1.0', '2024-01-01T00:00:00+0000', 1, f"{cls.wsName}/1/1", 1, cls.wsName, 'obj', 'checksum', 1, {}]])
        cls.wsClient.get_objects2 = Mock(return_value={'data': [{'data': {'type': 'test'}}]})

        cls.reportClient = Mock()
        cls.reportClient.create_extended_report = Mock(return_value={'ref': '1/1/1', 'name': 'test_report'})

        cls.dataFileUtil = Mock()
        cls.dataFileUtil.save_objects = Mock(return_value=[{'ref': '1/2/1'}])

        cls.kbUtilLib = Mock()
        cls.kbUtilLib.log = Mock()
        cls.kbUtilLib.get_workspace_client = Mock(return_value=cls.wsClient)
        cls.kbUtilLib.save_workspace_object = Mock(return_value='1/3/1')
        cls.kbUtilLib.create_report = Mock(return_value={'ref': '1/4/1', 'name': 'kb_util_report'})

        # Patch installed clients used inside Impl to return our mocks
        cls.patches = [
            patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.KBaseReport', return_value=cls.reportClient),
            patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.Workspace', return_value=cls.wsClient),
            patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.DataFileUtil', return_value=cls.dataFileUtil),
            patch('kbase_protein_query_module.kbase_protein_query_moduleImpl.KBUtilLib', return_value=cls.kbUtilLib),
        ]
        for p in cls.patches:
            p.start()

    @classmethod
    def tearDownClass(cls):
        for p in getattr(cls, 'patches', []):
            p.stop()

    def setUp(self):
        self.impl = kbase_protein_query_module({})
        self.ctx = {
            'token': self.token,
            'provenance': [{'ws_name': self.wsName}]
        }

    def test_status(self):
        res = self.impl.status(self.ctx)
        self.assertIsInstance(res, list)
        self.assertEqual(res[0]['state'], 'OK')

    def test_check_protein_existence(self):
        # Deprecated direct endpoint: ensure attribute not present in unified API
        self.assertFalse(hasattr(self.impl, 'check_protein_existence'))

    def test_generate_protein_embedding(self):
        # Deprecated direct endpoint: ensure attribute not present in unified API
        self.assertFalse(hasattr(self.impl, 'generate_protein_embedding'))

    def test_assign_family_fast(self):
        # Deprecated direct endpoint: ensure attribute not present in unified API
        self.assertFalse(hasattr(self.impl, 'assign_family_fast'))

    def test_find_top_matches_from_embedding(self):
        # Deprecated direct endpoint: ensure attribute not present in unified API
        self.assertFalse(hasattr(self.impl, 'find_top_matches_from_embedding'))

    def test_summarize_and_visualize_results(self):
        # Deprecated direct endpoint: ensure attribute not present in unified API
        self.assertFalse(hasattr(self.impl, 'summarize_and_visualize_results'))

    def test_run_protein_query_analysis_sequence_input(self):
        params = {
            'workspace_name': self.wsName,
            'input_type': 'protein_input',
            'protein_input': ['M' * 80],
            'analysis_name': 'test_analysis',
            # No legacy stage names; allow orchestrator defaults
        }
        res = self.impl.run_protein_query_analysis(self.ctx, params)
        out = res[0]
        # Shock integration fields
        self.assertIn('shock_id', out)
        self.assertIn('stages_completed', out)
        self.assertIsInstance(out['stages_completed'], list)

    def test_get_available_analyses(self):
        res = self.impl.get_available_analyses(self.ctx)
        out = res[0]
        self.assertIn('available_analyses', out)
        self.assertIn('summary', out)


if __name__ == '__main__':
    unittest.main()



