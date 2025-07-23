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
        if 'KB_AUTH_TOKEN' not in os.environ or not os.environ['KB_AUTH_TOKEN']:
            os.environ['KB_AUTH_TOKEN'] = 'DUMMY_TOKEN'
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
        if 'SDK_CALLBACK_URL' not in os.environ:
            os.environ['SDK_CALLBACK_URL'] = 'http://dummy-callback-url'
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

    # --- Unit tests ---
    def test_assign_protein_family(self):
        import numpy as np
        import tempfile
        import os
        from kbase_protein_network_analysis_toolkit.assign_protein_family import AssignProteinFamily
        temp_dir = tempfile.mkdtemp()
        try:
            centroid_file = os.path.join(temp_dir, 'family_centroids.npz')
            family_ids = np.array(['FAM1', 'FAM2'])
            centroids = np.array([[1.0, 0.0], [0.0, 1.0]])
            eigenprotein_ids = np.array(['P1', 'P2'])
            np.savez(centroid_file, family_ids=family_ids, centroids=centroids, eigenprotein_ids=eigenprotein_ids)
            assigner = AssignProteinFamily()
            assigner.load_family_centroids(centroid_file)
            embedding = [0.99, 0.01]
            ret = assigner.assign_family(embedding)
            self.assertEqual(ret['family_id'], 'FAM1')
            self.assertAlmostEqual(ret['confidence'], 0.99, places=1)
            self.assertEqual(ret['eigenprotein_id'], 'P1')
            embedding = [0.01, 0.99]
            ret = assigner.assign_family(embedding)
            self.assertEqual(ret['family_id'], 'FAM2')
            self.assertAlmostEqual(ret['confidence'], 0.99, places=1)
            self.assertEqual(ret['eigenprotein_id'], 'P2')
        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_check_existence(self):
        import tempfile
        import shutil
        import os
        import pandas as pd
        import numpy as np
        from kbase_protein_network_analysis_toolkit.check_existence import ProteinExistenceChecker
        from kbase_protein_network_analysis_toolkit.storage import ProteinStorage, CompressedMetadataStorage
        temp_dir = tempfile.mkdtemp()
        try:
            storage = ProteinStorage(base_dir=temp_dir)
            family_id = 'FAM1'
            protein_ids = ['P00001', 'P00002', 'P00003']
            embeddings = np.random.randn(3, 8).astype(np.float32)
            metadata = pd.DataFrame({
                'Protein names': ['Alpha', 'Beta', 'Gamma'],
                'Organism': ['E. coli', 'E. coli', 'E. coli']
            }, index=protein_ids)
            storage.store_family_embeddings(family_id, embeddings, protein_ids, metadata)
            metadata_storage = CompressedMetadataStorage(metadata_dir=str(storage.metadata_dir))
            metadata_storage.store_metadata(metadata, family_id=family_id)
            checker = ProteinExistenceChecker(storage=storage)
            # protein exists
            result = checker.check_protein_existence('P00001')
            self.assertTrue(result['exists'])
            self.assertEqual(result['family_id'], family_id)
            self.assertIsInstance(result['metadata'], dict)
            self.assertEqual(result['metadata']['Protein names'], 'Alpha')
            # protein does not exist
            result = checker.check_protein_existence('P99999')
            self.assertFalse(result['exists'])
            self.assertIsNone(result['family_id'])
            self.assertIsNone(result['metadata'])
            # empty protein id
            with self.assertRaises(ValueError):
                checker.check_protein_existence('')
            # corrupted metadata
            meta_file = os.path.join(storage.metadata_dir, f'family_{family_id}_metadata.parquet')
            os.remove(meta_file)
            result = checker.check_protein_existence('P00002')
            self.assertTrue(result['exists'])
            # performance batch
            import time
            ids = protein_ids * 100
            start = time.time()
            for pid in ids:
                checker.check_protein_existence(pid)
            elapsed = time.time() - start
            self.assertLess(elapsed, 2.0)
        finally:
            shutil.rmtree(temp_dir)

    def test_embedding_generator(self):
        import numpy as np
        from kbase_protein_network_analysis_toolkit.embedding_generator import ProteinEmbeddingGenerator
        generator = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D", device="cpu")
        seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        emb = generator.generate_embedding(seq, pooling_method="mean")
        self.assertEqual(emb.shape[0], generator.embedding_dim)
        emb = generator.generate_embedding(seq, pooling_method="cls")
        self.assertEqual(emb.shape[0], generator.embedding_dim)
        emb = generator.generate_embedding(seq, pooling_method="max")
        self.assertEqual(emb.shape[0], generator.embedding_dim)
        seqs = [seq, seq[:20], seq[:10]]
        ids = ["A", "B", "C"]
        embs = generator.generate_embeddings_batch(seqs, ids, pooling_method="mean", batch_size=2)
        self.assertEqual(len(embs), 3)
        for e in embs.values():
            self.assertEqual(e.shape[0], generator.embedding_dim)
        with self.assertRaises(ValueError):
            generator.generate_embedding(seq, pooling_method="invalid")

    def test_network_builder(self):
        import numpy as np
        import pandas as pd
        import os
        from kbase_protein_network_analysis_toolkit.network_builder import DynamicNetworkBuilder, visualize_interactive_protein_network
        embeddings = np.random.randn(10, 8).astype(np.float32)
        protein_ids = [f"P{i:03d}" for i in range(10)]
        metadata = pd.DataFrame({
            'Protein names': [f'Prot{i}' for i in range(10)],
            'Organism': ['E. coli']*10
        }, index=protein_ids)
        query_emb = np.random.randn(1, 8).astype(np.float32)
        query_id = "QUERY"
        builder = DynamicNetworkBuilder(k_neighbors=3, similarity_threshold=0.1)
        G = builder.build_mutual_knn_network(embeddings, protein_ids, query_emb, query_id)
        self.assertGreaterEqual(len(G.nodes), 11)
        self.assertGreater(len(G.edges), 0)
        G = builder.build_threshold_network(embeddings, protein_ids, query_emb, query_id)
        self.assertGreaterEqual(len(G.nodes), 11)
        G = builder.build_hybrid_network(embeddings, protein_ids, query_emb, query_id)
        self.assertGreaterEqual(len(G.nodes), 11)
        G = builder.build_mutual_knn_network(embeddings, protein_ids)
        props = builder.analyze_network_properties(G)
        self.assertIn('density', props)
        self.assertIn('average_degree', props)
        out = 'test_network.html'
        fig, G = visualize_interactive_protein_network(embeddings, protein_ids, metadata, query_embedding=query_emb, query_protein_id=query_id, output_file=out)
        self.assertIsNotNone(fig)
        self.assertTrue(os.path.exists(out))
        os.remove(out)

    def test_similarity_index(self):
        import tempfile
        import shutil
        import numpy as np
        import os
        from kbase_protein_network_analysis_toolkit.similarity_index import HierarchicalIndex, StreamingIndex
        temp_dir = tempfile.mkdtemp()
        try:
            embeddings = np.random.randn(10, 8).astype(np.float32)
            protein_ids = [f"P{i:03d}" for i in range(10)]
            index = HierarchicalIndex(base_dir=temp_dir, index_type='faiss', quantization='none', cache_size=2)
            family_id = 'FAMX'
            index.create_family_index(family_id, embeddings, protein_ids)
            query = np.random.randn(8).astype(np.float32)
            sims, ids = index.search_family(family_id, query, top_k=5)
            self.assertEqual(len(ids), 5)
            self.assertEqual(len(sims), 5)
            for i in range(5):
                fam = f'FAM{i}'
                index.create_family_index(fam, embeddings, protein_ids)
                index.search_family(fam, np.random.randn(8).astype(np.float32), top_k=1)
            self.assertLessEqual(len(index._family_cache), 2)
            streaming = StreamingIndex(storage_dir=temp_dir, batch_size=5)
            fam_assign = {pid: family_id for pid in protein_ids}
            emb_file = os.path.join(temp_dir, 'embeddings.h5')
            import h5py
            with h5py.File(emb_file, 'w') as f:
                f.create_dataset('embeddings', data=embeddings)
                f.create_dataset('protein_ids', data=protein_ids, dtype=h5py.special_dtype(vlen=str))
            streaming_file = streaming.create_streaming_index(emb_file, protein_ids, fam_assign)
            query = np.random.randn(8).astype(np.float32)
            results = streaming.stream_search(query, emb_file, streaming_file, top_k=3)
            self.assertLessEqual(len(results), 3)
        finally:
            shutil.rmtree(temp_dir)

    def test_storage(self):
        import tempfile
        import shutil
        import os
        import pandas as pd
        import numpy as np
        from kbase_protein_network_analysis_toolkit.storage import ProteinStorage, CompressedMetadataStorage
        temp_dir = tempfile.mkdtemp()
        try:
            storage = ProteinStorage(base_dir=temp_dir, chunk_size=2)
            family_id = 'FAM2'
            protein_ids = ['A', 'B', 'C', 'D']
            embeddings = np.random.randn(4, 4).astype(np.float32)
            metadata = pd.DataFrame({'desc': ['a', 'b', 'c', 'd']}, index=protein_ids)
            storage.store_family_embeddings(family_id, embeddings, protein_ids, metadata)
            emb, ids = storage.load_family_embeddings(family_id)
            self.assertEqual(len(ids), 4)
            np.testing.assert_array_almost_equal(emb, embeddings)
            storage.store_family_embeddings(family_id, embeddings, protein_ids)
            batches = list(storage.stream_family_embeddings(family_id, batch_size=2))
            self.assertEqual(len(batches), 2)
            self.assertEqual(batches[0][0].shape[0], 2)
            meta_storage = CompressedMetadataStorage(metadata_dir=str(storage.metadata_dir))
            meta_storage.store_metadata(metadata, family_id=family_id)
            loaded = meta_storage.load_metadata(family_id=family_id, protein_ids=['A', 'C'])
            self.assertEqual(list(loaded.index), ['A', 'C'])
            from kbase_protein_network_analysis_toolkit.storage import _create_artificial_families
            ids = [f'P{i:05d}' for i in range(7)]
            fams = list(_create_artificial_families(ids, max_family_size=3))
            self.assertEqual(len(fams), 3)
            self.assertEqual(len(fams[0][1]), 3)
            with self.assertRaises(FileNotFoundError):
                storage.load_family_embeddings('NOFAM')
        finally:
            shutil.rmtree(temp_dir)

    def test_workflow_orchestrator(self):
        import tempfile
        import shutil
        import os
        import yaml
        import numpy as np
        from kbase_protein_network_analysis_toolkit.workflow_orchestrator import ProteinNetworkWorkflow
        temp_dir = tempfile.mkdtemp()
        try:
            config_file = os.path.join(temp_dir, 'config.yaml')
            config = {
                'storage': {'optimized_storage_dir': temp_dir},
                'embedding': {'model_name': 'esm2_t6_8M_UR50D', 'device': 'cpu'},
                'logging': {'log_file': os.path.join(temp_dir, 'test.log'), 'level': 'INFO'}
            }
            with open(config_file, 'w') as f:
                yaml.dump(config, f)
            from kbase_protein_network_analysis_toolkit.embedding_generator import ProteinEmbeddingGenerator
            embedding_generator = ProteinEmbeddingGenerator(model_name='esm2_t6_8M_UR50D', device='cpu')
            embedding_dim = embedding_generator.get_embedding_dim() if hasattr(embedding_generator, 'get_embedding_dim') else 320
            from kbase_protein_network_analysis_toolkit.storage import ProteinStorage
            storage = ProteinStorage(base_dir=temp_dir)
            family_id = 'FAMW'
            protein_ids = ['X1', 'X2', 'X3']
            embeddings = np.random.randn(3, embedding_dim).astype(np.float32)
            storage.store_family_embeddings(family_id, embeddings, protein_ids)
            workflow = ProteinNetworkWorkflow(config_file=config_file)
            seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
            result = workflow.run_optimized_workflow(seq, query_protein_id="X1", k_similar=2, network_method="mutual_knn", save_results=False)
            self.assertEqual(result['status'], 'success')
            self.assertIn('query_embedding', result)
            self.assertIn('family_id', result)
            self.assertIn('similar_proteins', result)
            self.assertIn('network', result)
            self.assertIn('network_properties', result)
            self.assertIn('performance_metrics', result)
            with self.assertRaises(FileNotFoundError):
                from kbase_protein_network_analysis_toolkit.workflow_orchestrator import ProteinNetworkWorkflow
                ProteinNetworkWorkflow(config_file='nonexistent.yaml')
        finally:
            shutil.rmtree(temp_dir)
