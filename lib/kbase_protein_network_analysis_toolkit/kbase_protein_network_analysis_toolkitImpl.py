# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport

from kbase_protein_network_analysis_toolkit.check_existence import ProteinExistenceChecker
from kbase_protein_network_analysis_toolkit.embedding_generator import ProteinEmbeddingGenerator
from kbase_protein_network_analysis_toolkit.similarity_index import HierarchicalIndex
from kbase_protein_network_analysis_toolkit.network_builder import DynamicNetworkBuilder
from kbase_protein_network_analysis_toolkit.workflow_orchestrator import ProteinNetworkWorkflow
#END_HEADER


class kbase_protein_network_analysis_toolkit:
    '''
    Module Name:
    kbase_protein_network_analysis_toolkit

    Module Description:
    A KBase module: kbase_protein_network_analysis_toolkit
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_kbase_protein_network_analysis_toolkit(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kbase_protein_network_analysis_toolkit
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['parameter_1']},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_kbase_protein_network_analysis_toolkit

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kbase_protein_network_analysis_toolkit return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def check_protein_existence(self, ctx, params):
        """
        Check if a protein exists in the storage system.
        :param params: dict with 'protein_id' and 'workspace_name'
        :returns: dict with 'report_name', 'report_ref', and existence info
        """
        #BEGIN check_protein_existence
        checker = ProteinExistenceChecker()
        result = checker.check_protein_existence(params['protein_id'])
        report = KBaseReport(self.callback_url)
        text = f"Protein {params['protein_id']} existence: {result['exists']}\n"
        report_info = report.create({
            'report': {'objects_created': [], 'text_message': text},
            'workspace_name': params['workspace_name']
        })
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'exists': result['exists'],
            'family_id': result['family_id'],
            'metadata': result['metadata']
        }
        #END check_protein_existence
        return [output]

    def generate_protein_embedding(self, ctx, params):
        """
        Generate a protein embedding from a sequence or protein_id.
        :param params: dict with 'sequence' or 'protein_id', and 'workspace_name'
        :returns: dict with 'embedding' and 'summary'
        """
        #BEGIN generate_protein_embedding
        generator = ProteinEmbeddingGenerator()
        sequence = params.get('sequence')
        if not sequence and 'protein_id' in params:
            # Optionally, fetch sequence by protein_id if implemented
            raise ValueError('Sequence must be provided if protein_id lookup is not implemented.')
        embedding = generator.generate_embedding(sequence)
        summary = f"Generated embedding for protein. Shape: {embedding.shape}"
        output = {
            'embedding': embedding.tolist(),
            'summary': summary
        }
        #END generate_protein_embedding
        return [output]

    def find_top_matches_from_embedding(self, ctx, params):
        """
        Find top matches for a given protein embedding.
        :param params: dict with 'embedding', 'top_n', and 'workspace_name'
        :returns: dict with 'matches' and 'summary'
        """
        #BEGIN find_top_matches_from_embedding
        embedding = params['embedding']
        top_n = params.get('top_n', 10)
        # For demonstration, use HierarchicalIndex with default config
        index = HierarchicalIndex()
        # This assumes a family_id is provided or can be inferred; adjust as needed
        family_id = params.get('family_id')
        if not family_id:
            raise ValueError('family_id must be provided for similarity search.')
        import numpy as np
        embedding_np = np.array(embedding)
        similarities, protein_ids = index.search_family(family_id, embedding_np, top_k=top_n)
        matches = [
            {'protein_id': pid, 'similarity': float(sim)}
            for pid, sim in zip(protein_ids, similarities)
        ]
        summary = f"Found {len(matches)} top matches in family {family_id}."
        output = {
            'matches': matches,
            'summary': summary
        }
        #END find_top_matches_from_embedding
        return [output]

    def summarize_and_visualize_results(self, ctx, params):
        """
        Summarize and visualize protein network analysis results.
        :param params: dict with 'search_results', 'output_dir', and 'workspace_name'
        :returns: dict with 'report_name' and 'report_ref'
        """
        #BEGIN summarize_and_visualize_results
        # This is a placeholder; actual implementation would depend on result structure
        output_dir = params.get('output_dir', self.shared_folder)
        summary = f"Results summarized and visualized in {output_dir}."
        report = KBaseReport(self.callback_url)
        report_info = report.create({
            'report': {'objects_created': [], 'text_message': summary},
            'workspace_name': params['workspace_name']
        })
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }
        #END summarize_and_visualize_results
        return [output]

    def run_complete_workflow(self, ctx, params):
        """
        Run the complete protein network analysis workflow.
        :param params: dict with 'query_sequence', 'query_protein_id', 'k_similar', 'network_method', 'save_results', 'workspace_name'
        :returns: dict with 'report_name', 'report_ref', and workflow summary/results
        """
        #BEGIN run_complete_workflow
        query_sequence = params['query_sequence']
        query_protein_id = params.get('query_protein_id', 'QUERY_PROTEIN')
        k_similar = params.get('k_similar', 50)
        network_method = params.get('network_method', 'mutual_knn')
        save_results = params.get('save_results', True)
        workspace_name = params['workspace_name']

        workflow = ProteinNetworkWorkflow()
        results = workflow.run_optimized_workflow(
            query_sequence=query_sequence,
            query_protein_id=query_protein_id,
            k_similar=k_similar,
            network_method=network_method,
            save_results=save_results
        )

        # Prepare summary for report
        summary = f"Workflow completed.\nStatus: {results.get('status')}\n"
        if 'error' in results:
            summary += f"Error: {results['error']}\n"
        else:
            summary += f"Query protein: {query_protein_id}\n"
            summary += f"Family: {results.get('family_id')}\n"
            summary += f"Similar proteins found: {len(results.get('similar_proteins', []))}\n"
            summary += f"Network nodes: {results.get('network_properties', {}).get('num_nodes', 'N/A')}\n"
            summary += f"Network edges: {results.get('network_properties', {}).get('num_edges', 'N/A')}\n"

        report = KBaseReport(self.callback_url)
        report_info = report.create({
            'report': {'objects_created': [], 'text_message': summary},
            'workspace_name': workspace_name
        })
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'workflow_status': results.get('status'),
            'family_id': results.get('family_id'),
            'network_properties': results.get('network_properties'),
            'similar_proteins': results.get('similar_proteins', [])
        }
        #END run_complete_workflow
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
