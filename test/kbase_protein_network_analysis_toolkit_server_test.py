# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from kbase_protein_network_analysis_toolkit.kbase_protein_network_analysis_toolkitImpl import kbase_protein_network_analysis_toolkit
from kbase_protein_network_analysis_toolkit.kbase_protein_network_analysis_toolkitServer import MethodContext
from kbase_protein_network_analysis_toolkit.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class kbase_protein_network_analysis_toolkitTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kbase_protein_network_analysis_toolkit'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kbase_protein_network_analysis_toolkit',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kbase_protein_network_analysis_toolkit(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_run_kbase_protein_network_analysis_toolkit(self):
        params = {'workspace_name': self.wsName, 'parameter_1': 'Hello World!'}
        ret = self.serviceImpl.run_kbase_protein_network_analysis_toolkit(self.ctx, params)
        self.assertIsInstance(ret, list)
        self.assertGreaterEqual(len(ret), 1)
        output = ret[0]
        self.assertIn('report_name', output)
        self.assertIn('report_ref', output)
        self.assertIn('input_parameters', output)
        self.assertEqual(output['input_parameters'], params)
        self.assertIn('summary', output)
        self.assertIn('start_time', output)

    def test_run_complete_workflow(self):
        # This test will only check that the workflow runs and returns a report, not the full biological correctness
        params = {
            'workspace_name': self.wsName,
            'query_sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ',
            'query_protein_id': 'TEST_PROT',
            'k_similar': 2,
            'network_method': 'mutual_knn',
            'save_results': False
        }
        ret = self.serviceImpl.run_complete_workflow(self.ctx, params)
        self.assertIsInstance(ret, list)
        self.assertGreaterEqual(len(ret), 1)
        output = ret[0]
        self.assertIn('report_name', output)
        self.assertIn('report_ref', output)
        self.assertIn('workflow_status', output)
        self.assertIn(output['workflow_status'], ['success', 'error'])
        self.assertIn('summary', output)
        self.assertIn('input_parameters', output)
        self.assertEqual(output['input_parameters'], params)
