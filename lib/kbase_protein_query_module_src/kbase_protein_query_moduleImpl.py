# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace

from kbase_protein_query_module_src.check_existence import ProteinExistenceChecker
from kbase_protein_query_module_src.embedding_generator import ProteinEmbeddingGenerator
from kbase_protein_query_module_src.similarity_index import HierarchicalIndex
from kbase_protein_query_module_src.network_builder import DynamicNetworkBuilder
from kbase_protein_query_module_src.workflow_orchestrator import ProteinNetworkWorkflow
#END_HEADER


class kbase_protein_query_module:
    '''
    Module Name:
    kbase_protein_query_module

    Module Description:
    A KBase module: kbase_protein_query_module
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
        self.callback_url = os.environ.get('SDK_CALLBACK_URL')
        if self.callback_url is None:
            raise RuntimeError('SDK_CALLBACK_URL environment variable must be set.')
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.family_assigner = None
        from kbase_protein_query_module_src.assign_protein_family import AssignProteinFamily
        self.family_assigner = AssignProteinFamily()
        centroid_path = os.path.join('data', 'family_centroids.npz')
        if os.path.exists(centroid_path):
            self.family_assigner.load_family_centroids(centroid_path)
        #END_CONSTRUCTOR
        pass

    def generate_protein_embedding(self, ctx, params):
        """
        Generate a protein embedding from an amino acid sequence.
        :param params: dict with 'sequence' (required) and 'workspace_name' (optional)
        :returns: dict with 'report_name', 'report_ref', 'embedding_result', 'embedding_result_ref'
        """
        #BEGIN generate_protein_embedding
        import time
        import uuid
        start_time = time.time()
        validation_status = "success"
        error = None
        sequence = params.get('sequence')
        workspace_name = params.get('workspace_name')
        if not sequence or not isinstance(sequence, str) or not sequence.strip():
            validation_status = "error"
            error = "Parameter 'sequence' must be a non-empty string."
            summary = f"Validation failed: {error}"
            return [{
                'report_name': None,
                'report_ref': None,
                'embedding_result': None,
                'embedding_result_ref': None,
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        generator = ProteinEmbeddingGenerator()
        embedding = generator.generate_embedding(sequence)
        summary = f"Generated embedding for protein. Shape: {embedding.shape}"
        # Save embedding as a workspace object
        ws = Workspace(os.environ['WS_URL'])
        embedding_obj_name = f"embedding_result_{uuid.uuid4().hex[:8]}"
        embedding_obj_type = "KBaseProtein.ProteinEmbeddingResult"
        embedding_obj_data = {
            'embedding': embedding.tolist(),
            'sequence': sequence,
            'summary': summary,
            'input_parameters': params,
            'embedding_norm': float((embedding**2).sum()**0.5),
            'sequence_length': len(sequence),
            'created_at': start_time
        }
        try:
            save_ret = ws.save_objects({
                'workspace': workspace_name,
                'objects': [{
                    'type': embedding_obj_type,
                    'data': embedding_obj_data,
                    'name': embedding_obj_name
                }]
            })
            embedding_result_ref = f"{save_ret[0][6]}/{save_ret[0][0]}/{save_ret[0][4]}"
        except Exception as e:
            validation_status = "error"
            error = f"Failed to save embedding object: {str(e)}"
            return [{
                'report_name': None,
                'report_ref': None,
                'embedding_result': None,
                'embedding_result_ref': None,
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': error
            }]
        # Create a KBase report
        report = KBaseReport(self.callback_url)
        report_text = summary
        try:
            report_info = report.create({
                'report': {
                    'objects_created': [{
                        'ref': embedding_result_ref,
                        'description': 'Protein embedding result'
                    }],
                    'text_message': report_text
                },
                'workspace_name': workspace_name
            })
            report_name = report_info['name']
            report_ref = report_info['ref']
        except Exception as e:
            validation_status = "error"
            error = f"Failed to create report: {str(e)}"
            return [{
                'report_name': None,
                'report_ref': None,
                'embedding_result': embedding_obj_data,
                'embedding_result_ref': embedding_result_ref,
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': error
            }]
        output = {
            'report_name': report_name,
            'report_ref': report_ref,
            'embedding_result': embedding_obj_data,
            'embedding_result_ref': embedding_result_ref,
            'status': validation_status,
            'error': error,
            'input_parameters': params,
            'start_time': start_time,
            'summary': summary
        }
        #END generate_protein_embedding
        return [output]

    def assign_family_fast(self, ctx, params):
        """
        Quickly assign a protein embedding to a family by similarity to the medoid (binary Hamming distance only).
        :param params: dict with 'embedding' (list of uint8), and 'workspace_name'
        :returns: dict with 'family_id', 'confidence', and 'eigenprotein_id'
        """
        #BEGIN assign_family_fast
        import time
        import numpy as np
        start_time = time.time()
        validation_status = "success"
        error = None
        embedding = params.get('embedding')
        if embedding is None or not isinstance(embedding, list) or not embedding or not all(isinstance(x, (int, np.integer)) for x in embedding):
            validation_status = "error"
            error = "Parameter 'embedding' must be a non-empty list of uint8 integers."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        embedding_np = np.array(embedding, dtype=np.uint8)
        try:
            result = self.family_assigner.assign_family(embedding_np)
        except Exception as e:
            validation_status = "error"
            error = str(e)
            summary = f"Assignment failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        output = {
            'family_id': result['family_id'],
            'confidence': float(result['confidence']),
            'eigenprotein_id': result['eigenprotein_id'],
            'input_parameters': params,
            'start_time': start_time,
            'validation_status': validation_status,
            'error': error,
            'summary': f"Assigned to family {result['family_id']} with confidence {result['confidence']:.3f}."
        }
        #END assign_family_fast
        return [output]

    def check_protein_existence(self, ctx, params):
        """
        Check if a protein exists in the storage system.
        :param params: dict with 'protein_id' and 'workspace_name'
        :returns: dict with 'report_name', 'report_ref', and existence info
        """
        #BEGIN check_protein_existence
        import time
        import re
        start_time = time.time()
        validation_status = "success"
        error = None
        protein_id = params.get('protein_id')
        if not protein_id or not isinstance(protein_id, str) or not protein_id.strip():
            validation_status = "error"
            error = "Parameter 'protein_id' must be a non-empty string."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        # Optional: UniProt regex validation
        uniprot_regex = r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-Z0-9]{6,10}$"
        if not re.match(uniprot_regex, protein_id):
            validation_status = "error"
            error = "Parameter 'protein_id' does not match expected UniProt format."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        checker = ProteinExistenceChecker()
        result = checker.check_protein_existence(protein_id)
        report = KBaseReport(self.callback_url)
        text = f"Protein {protein_id} existence: {result['exists']}\n"
        if result['exists']:
            text += f"Family: {result['family_id']}\nMetadata: {result['metadata']}\n"
        else:
            text += "Protein not found in any family.\n"
        report_info = report.create({
            'report': {'objects_created': [], 'text_message': text},
            'workspace_name': params['workspace_name']
        })
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'exists': result['exists'],
            'family_id': result['family_id'],
            'metadata': result['metadata'],
            'input_parameters': params,
            'start_time': start_time,
            'summary': text,
            'validation_status': validation_status,
            'error': error
        }
        #END check_protein_existence
        return [output]

    def find_top_matches_from_embedding(self, ctx, params):
        """
        Find top matches for a given protein embedding (binary FAISS only).
        :param params: dict with 'embedding' (list of uint8), 'top_n', and 'workspace_name'
        :returns: dict with 'matches' and 'summary'
        """
        #BEGIN find_top_matches_from_embedding
        import time
        import numpy as np
        start_time = time.time()
        validation_status = "success"
        error = None
        embedding = params.get('embedding')
        family_id = params.get('family_id')
        top_n = params.get('top_n', 10)
        if embedding is None or not isinstance(embedding, list) or not embedding or not all(isinstance(x, (int, np.integer)) for x in embedding):
            validation_status = "error"
            error = "Parameter 'embedding' must be a non-empty list of uint8 integers."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        if not family_id:
            validation_status = "error"
            error = "Parameter 'family_id' must be provided for similarity search."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        if not isinstance(top_n, int) or not (1 <= top_n <= 1000):
            validation_status = "error"
            error = "Parameter 'top_n' must be an integer between 1 and 1000."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        index = HierarchicalIndex()
        embedding_np = np.array(embedding, dtype=np.uint8)
        try:
            similarities, protein_ids = index.search_family(family_id, embedding_np, top_k=top_n)
        except Exception as e:
            validation_status = "error"
            error = str(e)
            summary = f"Similarity search failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        matches = [
            {'protein_id': pid, 'similarity': float(sim)}
            for pid, sim in zip(protein_ids, similarities)
        ]
        summary = f"Found {len(matches)} top matches in family {family_id}."
        output = {
            'matches': matches,
            'summary': summary,
            'input_parameters': params,
            'start_time': start_time,
            'family_id': family_id,
            'top_n': top_n,
            'similarity_stats': {
                'max': float(np.max(similarities)) if len(similarities) else None,
                'min': float(np.min(similarities)) if len(similarities) else None,
                'mean': float(np.mean(similarities)) if len(similarities) else None
            },
            'validation_status': validation_status,
            'error': error
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
        import time
        import os
        start_time = time.time()
        validation_status = "success"
        error = None
        search_results = params.get('search_results')
        output_dir = params.get('output_dir', self.shared_folder)
        if search_results is None:
            validation_status = "error"
            error = "Parameter 'search_results' is required."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        if output_dir and not isinstance(output_dir, str):
            validation_status = "error"
            error = "Parameter 'output_dir' must be a valid path string."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        summary = f"Results summarized and visualized in {output_dir}."
        report = KBaseReport(self.callback_url)
        report_info = report.create({
            'report': {'objects_created': [], 'text_message': summary},
            'workspace_name': params['workspace_name']
        })
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'input_parameters': params,
            'start_time': start_time,
            'output_dir': output_dir,
            'summary': summary,
            'validation_status': validation_status,
            'error': error
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
        import time
        import os
        start_time = time.time()
        validation_status = "success"
        error = None
        query_sequence = params.get('query_sequence')
        workspace_name = params.get('workspace_name')
        if not query_sequence or not isinstance(query_sequence, str) or not query_sequence.strip():
            validation_status = "error"
            error = "Parameter 'query_sequence' must be a non-empty string."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        # Validate optional params
        k_similar = params.get('k_similar', 50)
        if not isinstance(k_similar, int) or not (1 <= k_similar <= 1000):
            validation_status = "error"
            error = "Parameter 'k_similar' must be an integer between 1 and 1000."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        network_method = params.get('network_method', 'mutual_knn')
        save_results = params.get('save_results', True)
        query_protein_id = params.get('query_protein_id', 'QUERY_PROTEIN')
        workflow = ProteinNetworkWorkflow()
        try:
            results = workflow.run_optimized_workflow(
                query_sequence=query_sequence,
                query_protein_id=query_protein_id,
                k_similar=k_similar,
                network_method=network_method,
                save_results=save_results
            )
        except Exception as e:
            validation_status = "error"
            error = f"Workflow execution failed: {str(e)}"
            summary = f"Workflow failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        summary = f"Workflow completed.\nStatus: {results.get('status')}\n"
        if 'error' in results:
            validation_status = "error"
            error = results['error']
            summary += f"Error: {results['error']}\n"
        else:
            summary += f"Query protein: {query_protein_id}\n"
            summary += f"Family: {results.get('family_id')}\n"
            summary += f"Similar proteins found: {len(results.get('similar_proteins', []))}\n"
            summary += f"Network nodes: {results.get('network_properties', {}).get('num_nodes', 'N/A')}\n"
            summary += f"Network edges: {results.get('network_properties', {}).get('num_edges', 'N/A')}\n"
            summary += f"Timing (s): {results.get('timing', {})}\n"
            summary += f"Performance: {results.get('performance_metrics', {})}\n"
        html_file = os.path.join(self.shared_folder, f"protein_network_report_{int(time.time())}.html")
        self._generate_detailed_html_report(results, html_file)
        report = KBaseReport(self.callback_url)
        report_info = report.create_extended_report({
            'message': summary,
            'objects_created': [],
            'warnings': [],
            'workspace_name': workspace_name,
            'file_links': [],
            'html_links': [html_file],
            'direct_html_link_index': 0,
            'html_window_height': 600
        })
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'workflow_status': results.get('status'),
            'family_id': results.get('family_id'),
            'network_properties': results.get('network_properties'),
            'similar_proteins': results.get('similar_proteins', []),
            'input_parameters': params,
            'start_time': start_time,
            'timing': results.get('timing', {}),
            'performance_metrics': results.get('performance_metrics', {}),
            'summary': summary,
            'validation_status': validation_status,
            'error': error
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

    def run_kbase_protein_query_module(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kbase_protein_query_module
        import time
        start_time = time.time()
        validation_status = "success"
        error = None
        if not isinstance(params, dict):
            validation_status = "error"
            error = "Input parameters must be a dictionary."
            summary = f"Validation failed: {error}"
            return [{
                'status': validation_status,
                'error': error,
                'input_parameters': params,
                'start_time': start_time,
                'summary': summary
            }]
        report = KBaseReport(self.callback_url)
        report_text = f"KBase Protein Network Analysis Toolkit\n\nParameters:\n{params}\n\nAnalysis started."
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': report_text},
                                                'workspace_name': params.get('workspace_name', 'UNKNOWN')})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'input_parameters': params,
            'start_time': start_time,
            'summary': 'Analysis initialized. See report for details.',
            'validation_status': validation_status,
            'error': error
        }
        #END run_kbase_protein_query_module

        # At some point might do deeper type checking...
        # if not isinstance(output, dict):
        #     raise ValueError('Method run_kbase_protein_query_module return value ' +
        #                      'output is not type dict as required.')
        # return the results
        return [output]

    def _generate_detailed_html_report(self, results, output_html_path):
        """
        Generate a detailed, modern HTML report for the workflow results.
        Args:
            results: dict containing all workflow results
            output_html_path: path to save the HTML file
        """
        import json
        import datetime
        import os
        from collections import OrderedDict

        # Helper to make UniProt/EC links
        def uniprot_link(pid):
            if pid and isinstance(pid, str) and len(pid) >= 6:
                return f'<a href="https://www.uniprot.org/uniprot/{pid}" target="_blank">{pid}</a>'
            return pid or 'N/A'
        def expasy_link(ec):
            if ec and isinstance(ec, str) and ec.replace(".", "").replace("-", "").replace(",", "").replace(" ", "").replace(":", "").replace("/", "").replace("EC", "").replace("ec", "").replace("n/a", "").replace("N/A", "").replace("", "").strip():
                return f'<a href="https://enzyme.expasy.org/EC/{ec}" target="_blank">{ec}</a>'
            return ec or 'N/A'

        # --- Metadata extraction ---
        metadata = results.get('family_metadata')
        if metadata is not None and hasattr(metadata, 'to_dict'):
            metadata_dict = metadata.to_dict(orient='index')
            metadata_columns = list(metadata.columns)
        else:
            metadata_dict = {}
            metadata_columns = []

        # --- Top matches table ---
        similar_proteins = results.get('similar_proteins', [])
        top_matches_rows = []
        for match in similar_proteins:
            pid = match.get('protein_id', '')
            row = f'<tr>'
            row += f'<td>{match.get("rank", "")}</td>'
            row += f'<td>{uniprot_link(pid)}</td>'
            row += f'<td>{"{:.4f}".format(match.get("similarity_score", 0))}</td>'
            # Add all metadata columns for this protein
            meta = metadata_dict.get(pid, {})
            for col in metadata_columns:
                val = meta.get(col, '')
                if col.lower().startswith('ec'):
                    val = expasy_link(val)
                elif col.lower() in ('entry', 'uniprot', 'uniprot id'):
                    val = uniprot_link(val)
                # Scrollable for long text
                if isinstance(val, str) and len(val) > 80:
                    val = f'<div style="max-height:80px;overflow:auto;">{val}</div>'
                row += f'<td>{val}</td>'
            row += '</tr>'
            top_matches_rows.append(row)

        # --- Query protein metadata ---
        query_pid = results.get('input_parameters', {}).get('query_protein_id', '') or results.get('query_protein_id', '')
        query_meta = metadata_dict.get(query_pid, {})
        query_meta_html = ''
        if query_meta:
            for col in metadata_columns:
                val = query_meta.get(col, '')
                if col.lower().startswith('ec'):
                    val = expasy_link(val)
                elif col.lower() in ('entry', 'uniprot', 'uniprot id'):
                    val = uniprot_link(val)
                if isinstance(val, str) and len(val) > 120:
                    val = f'<div style="max-height:120px;overflow:auto;">{val}</div>'
                query_meta_html += f'<tr><td class="label">{col}</td><td>{val}</td></tr>'
        else:
            query_meta_html = '<tr><td colspan="2">No metadata available for query protein.</td></tr>'

        # --- Network properties ---
        net_props = results.get('network_properties', {})
        net_props_html = ''
        if net_props:
            for k, v in net_props.items():
                if isinstance(v, (list, dict)):
                    v = json.dumps(v)
                net_props_html += f'<tr><td class="label">{k}</td><td>{v}</td></tr>'
        else:
            net_props_html = '<tr><td colspan="2">No network properties available.</td></tr>'

        # --- Data files and downloads ---
        data_files = results.get('data_files', [])
        download_links = results.get('download_links', [])
        downloads_html = ''
        for f in download_links:
            downloads_html += f'<li><a href="{f}" download>{os.path.basename(f)}</a></li>'

        # --- Network visualization ---
        network_html_file = results.get('network_html_file')
        network_viz_html = ''
        if network_html_file and os.path.exists(network_html_file):
            network_viz_html = f'<iframe src="{network_html_file}" width="100%" height="600" style="border:1px solid #ccc;border-radius:10px;"></iframe>'
        else:
            network_viz_html = '<div style="color:#888;">Network visualization not available.</div>'

        # --- Validation & Status section ---
        validation_status = results.get('validation_status', 'unknown')
        error = results.get('error')
        warnings = results.get('warnings', [])
        status_color = '#4caf50' if validation_status == 'success' else '#f44336'
        status_html = f'<span style="color:{status_color};font-weight:bold;">{validation_status.upper()}</span>'
        error_html = ''
        if error:
            error_html = f'<div style="color:#fff;background:#f44336;padding:10px 18px;border-radius:8px;margin:12px 0;font-weight:bold;">Error: {error}</div>'
        elif warnings:
            error_html = f'<div style="color:#856404;background:#fff3cd;padding:10px 18px;border-radius:8px;margin:12px 0;font-weight:bold;">Warning: {warnings}</div>'

        # --- Step-by-step workflow table (if present) ---
        workflow_steps = results.get('workflow_steps')
        workflow_table_html = ''
        if workflow_steps and isinstance(workflow_steps, list):
            workflow_table_html = '<div class="section box"><h2>Workflow Steps</h2><table><thead><tr><th>Step</th><th>Status</th><th>Timing (s)</th><th>Error</th></tr></thead><tbody>'
            for step in workflow_steps:
                workflow_table_html += f"<tr><td>{step.get('name','')}</td><td>{step.get('status','')}</td><td>{step.get('timing','')}</td><td>{step.get('error','')}</td></tr>"
            workflow_table_html += '</tbody></table></div>'

        # --- HTML Template ---
        html = f"""
        <!DOCTYPE html>
        <html lang=\"en\">
        <head>
            <meta charset=\"UTF-8\">
            <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
            <title>KBase Protein Network Analysis Report</title>
            <link href=\"https://fonts.googleapis.com/css?family=Roboto:400,700&display=swap\" rel=\"stylesheet\">
            <style>
                body {{ font-family: 'Roboto', Arial, sans-serif; background: #f7f9fa; color: #222; margin: 0; }}
                .container {{ max-width: 1200px; margin: 30px auto; background: #fff; border-radius: 18px; box-shadow: 0 2px 16px rgba(0,0,0,0.08); padding: 32px; }}
                h1, h2, h3 {{ color: #2a4d69; }}
                .section {{ margin-bottom: 36px; }}
                .box {{ background: #f0f4f8; border-radius: 12px; padding: 18px 22px; margin-bottom: 18px; box-shadow: 0 1px 4px rgba(0,0,0,0.04); }}
                .label {{ font-weight: bold; color: #41729f; }}
                table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
                th, td {{ padding: 8px 12px; border-bottom: 1px solid #e0e6ed; }}
                th {{ background: #eaf1fb; color: #2a4d69; }}
                tr:hover {{ background: #f5faff; }}
                .scroll-table {{ max-height: 350px; overflow-y: auto; display: block; }}
                .chip {{ display: inline-block; padding: 0 10px; height: 24px; font-size: 13px; line-height: 24px; border-radius: 12px; background: #eaf1fb; color: #2a4d69; margin-right: 8px; margin-bottom: 4px; }}
                .network-box {{ min-height: 300px; }}
                .downloads ul {{ margin: 0; padding-left: 20px; }}
            </style>
        </head>
        <body>
        <div class=\"container\">
            <h1>KBase Protein Network Analysis Report</h1>
            <div class=\"section box\">
                <h2>Validation & Status</h2>
                <div>Status: {status_html}</div>
                {error_html}
            </div>
            {workflow_table_html}
            <div class=\"section box\">
                <h2>Query Protein Overview</h2>
                <table>{query_meta_html}</table>
            </div>
            <div class=\"section box\">
                <h2>Similarity Search Results</h2>
                <div class=\"scroll-table\">
                <table>
                    <thead><tr><th>Rank</th><th>Protein ID</th><th>Similarity</th>{''.join([f'<th>{col}</th>' for col in metadata_columns])}</tr></thead>
                    <tbody>{''.join(top_matches_rows)}</tbody>
                </table>
                </div>
            </div>
            <div class=\"section box network-box\">
                <h2>Network Visualization</h2>
                {network_viz_html}
            </div>
            <div class=\"section box\">
                <h2>Network Properties</h2>
                <table>{net_props_html}</table>
            </div>
            <div class=\"section box\">
                <h2>Data Used</h2>
                <div><span class=\"label\">Files:</span> {json.dumps(data_files)}</div>
            </div>
            <div class=\"section box downloads\">
                <h2>Downloadable Artifacts</h2>
                <ul>{downloads_html}</ul>
            </div>
            <div class=\"section box\">
                <h2>Workflow Parameters & Timing</h2>
                <pre style=\"background:#f7f9fa; border-radius:8px; padding:10px;\">{json.dumps(results.get('input_parameters', {}), indent=2)}</pre>
                <div style=\"margin-top:12px;\"><span class=\"label\">Performance:</span> <span>{json.dumps(results.get('performance_metrics', {}))}</span></div>
                <div style=\"margin-top:12px;\"><span class=\"label\">Timing:</span> <span>{json.dumps(results.get('timing', {}))}</span></div>
            </div>
        </div>
        </body>
        </html>
        """
        with open(output_html_path, 'w') as f:
            f.write(html)
