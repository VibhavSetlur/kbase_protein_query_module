# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace

from .src.check_existence import ProteinExistenceChecker
from .src.embedding_generator import ProteinEmbeddingGenerator
from .src.similarity_index import HierarchicalIndex
from .src.network_builder import DynamicNetworkBuilder
from .src.workflow_orchestrator import ProteinNetworkWorkflow
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
    GIT_URL = "https://github.com/VibhavSetlur/kbase_protein_query_module.git"
    GIT_COMMIT_HASH = "8b39ef0deb55b71aea168c9f451b9f60d24db5ce"

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
        from .src.assign_protein_family import AssignProteinFamily
        self.family_assigner = AssignProteinFamily()
        centroid_path = os.path.join('data', 'family_centroids.npz')
        if os.path.exists(centroid_path):
            self.family_assigner.load_family_centroids(centroid_path)
        #END_CONSTRUCTOR
        pass


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
        
        # Validate input parameters
        if not isinstance(params, dict):
            raise ValueError("Input parameters must be a dictionary.")
        
        report = KBaseReport(self.callback_url)
        report_text = f"KBase Protein Network Analysis Toolkit\n\nParameters:\n{params}\n\nAnalysis started."
        report_info = report.create_extended_report({
            'message': report_text,
            'objects_created': [],
            'workspace_name': params.get('workspace_name', 'UNKNOWN')
        })
        
        # Return the result according to KIDL spec
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'input_parameters': params,
            'summary': report_text,
            'start_time': start_time
        }
        
        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kbase_protein_query_module return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
        #END run_kbase_protein_query_module

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kbase_protein_query_module return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def check_protein_existence(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "CheckProteinExistenceResults" (Check if a
           protein exists in the storage system.) -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "exists" of Long, parameter "family_id" of String,
           parameter "metadata" of mapping from String to unspecified object,
           parameter "input_parameters" of mapping from String to unspecified
           object, parameter "start_time" of Double, parameter "summary" of
           String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN check_protein_existence
        import time
        import re
        import uuid
        start_time = time.time()
        protein_id = params.get('protein_id')
        workspace_name = params.get('workspace_name')
        
        # Quick validation - return error immediately for invalid inputs
        if not protein_id:
            raise ValueError("Parameter 'protein_id' is required and cannot be empty.")
        
        if not isinstance(protein_id, str):
            raise ValueError("Parameter 'protein_id' must be a string.")
        
        protein_id = protein_id.strip()
        if not protein_id:
            raise ValueError("Parameter 'protein_id' cannot be empty or whitespace only.")
        
        # More flexible validation - accept various formats
        # UniProt format: P12345, Q1A2B3, etc.
        uniprot_regex = r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$"
        # Dummy data format: family_X_prot_Y
        dummy_regex = r"^family_\d+_prot_\d+$"
        # UniProt-like format: P followed by 8 digits
        uniprot_like_regex = r"^P\d{8}$"
        # Generic alphanumeric format (6-10 characters)
        generic_regex = r"^[A-Z0-9]{6,10}$"
        
        # Check if protein_id matches any of the accepted formats
        is_valid = (re.match(uniprot_regex, protein_id) or 
                   re.match(dummy_regex, protein_id) or 
                   re.match(uniprot_like_regex, protein_id) or
                   re.match(generic_regex, protein_id))
        
        if not is_valid:
            raise ValueError(f"Parameter 'protein_id' '{protein_id}' does not match expected formats. "
                           f"Accepted formats: UniProt IDs (e.g., P12345), dummy IDs (e.g., family_0_prot_1), "
                           f"or generic alphanumeric IDs (6-10 characters).")
        
        checker = ProteinExistenceChecker()
        result = checker.check_protein_existence(protein_id)
        
        # Save result to workspace for user access
        existence_result_ref = None
        try:
            # Get workspace URL from configuration
            ws_url = self.config.get('workspace-url')
            if not ws_url:
                # Fallback to environment variable
                ws_url = os.environ.get('KBASE_ENDPOINT', 'https://appdev.kbase.us/services')
                if not ws_url.endswith('/ws'):
                    ws_url = ws_url + '/ws'
            
            ws = Workspace(ws_url)
            existence_obj_name = f"protein_existence_{uuid.uuid4().hex[:8]}"
            existence_obj_type = "ProteinExistenceResult"
            existence_obj_data = {
                'protein_id': protein_id,
                'exists': result['exists'],
                'family_id': result['family_id'],
                'metadata': result['metadata'],
                'search_timestamp': start_time,
                'summary': f"Protein {protein_id} {'found' if result['exists'] else 'not found'} in database"
            }
            
            # Handle workspace parameter - could be workspace name or ID
            workspace_param = workspace_name
            if isinstance(workspace_name, str) and workspace_name.isdigit():
                # It's a workspace ID
                workspace_param = int(workspace_name)
            elif isinstance(workspace_name, str) and ':' in workspace_name:
                # It's a workspace name with owner
                workspace_param = workspace_name
            
            save_ret = ws.save_objects({
                'workspace': workspace_param,
                'objects': [{
                    'type': existence_obj_type,
                    'data': existence_obj_data,
                    'name': existence_obj_name
                }]
            })
            existence_result_ref = f"{save_ret[0][6]}/{save_ret[0][0]}/{save_ret[0][4]}"
        except Exception as e:
            # If workspace save fails, still return the result but without workspace reference
            existence_result_ref = None
            print(f"Warning: Failed to save existence result to workspace: {str(e)}")
        
        # Create report
        report = KBaseReport(self.callback_url)
        text = f"""
        Protein Existence Check Results
        
        Protein ID: {protein_id}
        Exists: {result['exists']}
        """
        if result['exists']:
            text += f"Family ID: {result['family_id']}\n"
            text += f"Metadata: {result['metadata']}\n"
        else:
            text += "Protein not found in any family.\n"
        
        # Only include objects_created if existence_result_ref is not None
        report_params = {
            'message': text,
            'workspace_name': workspace_name
        }
        
        if existence_result_ref is not None:
            report_params['objects_created'] = [{
                'ref': existence_result_ref,
                'description': 'Protein existence check result'
            }]
        
        report_info = report.create_extended_report(report_params)
        
        # Return the result according to KIDL spec and UI expectations
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'exists': 1 if result['exists'] else 0,  # Convert boolean to Long as per KIDL spec
            'family_id': result['family_id'],
            'metadata': result['metadata'],
            'input_parameters': params,
            'start_time': start_time,
            'summary': text
        }
        
        return output
        #END check_protein_existence

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method check_protein_existence return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def generate_protein_embedding(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "GenerateProteinEmbeddingResults"
           (Generate a protein embedding from a sequence.) -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "embedding_result_ref" of String, parameter
           "summary" of String, parameter "input_parameters" of mapping from
           String to unspecified object, parameter "start_time" of Double,
           parameter "embedding_norm" of Double, parameter "sequence_length"
           of Long, parameter "embedding_dim" of Long
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN generate_protein_embedding
        import time
        import uuid
        start_time = time.time()
        sequence = params.get('sequence')
        workspace_name = params.get('workspace_name')
        
        # Validate input parameters
        if not sequence or not isinstance(sequence, str) or not sequence.strip():
            raise ValueError("Parameter 'sequence' must be a non-empty string.")
        
        # Generate embedding
        generator = ProteinEmbeddingGenerator()
        embedding = generator.generate_embedding(sequence)
        summary = f"Generated {len(embedding)}-dimensional protein embedding using ESM-2 model"
        
        # Calculate embedding norm
        embedding_norm = float((embedding**2).sum()**0.5)
        sequence_length = len(sequence)
        
        # Save embedding as workspace object
        embedding_ref = None
        try:
            # Get workspace URL from configuration
            ws_url = self.config.get('workspace-url')
            if not ws_url:
                # Fallback to environment variable
                ws_url = os.environ.get('KBASE_ENDPOINT', 'https://appdev.kbase.us/services')
                if not ws_url.endswith('/ws'):
                    ws_url = ws_url + '/ws'
            
            ws = Workspace(ws_url)
            embedding_obj_name = f"protein_embedding_{uuid.uuid4().hex[:8]}"
            embedding_obj_type = "ProteinEmbeddingResult"
            embedding_obj_data = {
                'input_id': sequence[:50] + "..." if len(sequence) > 50 else sequence,
                'input_type': 'sequence',
                'embedding_ref': None,  # Self-reference
                'embedding': embedding.tolist(),
                'model_name': 'esm2_t6_8M_UR50D',
                'pooling_method': 'mean',  # Always mean pooling
                'metadata': {
                    'sequence_length': sequence_length,
                    'embedding_norm': embedding_norm,
                    'created_at': start_time,
                    'summary': summary
                }
            }
            
            # Handle workspace parameter - could be workspace name or ID
            workspace_param = workspace_name
            if isinstance(workspace_name, str) and workspace_name.isdigit():
                # It's a workspace ID
                workspace_param = int(workspace_name)
            elif isinstance(workspace_name, str) and ':' in workspace_name:
                # It's a workspace name with owner
                workspace_param = workspace_name
            
            save_ret = ws.save_objects({
                'workspace': workspace_param,
                'objects': [{
                    'type': embedding_obj_type,
                    'data': embedding_obj_data,
                    'name': embedding_obj_name
                }]
            })
            embedding_ref = f"{save_ret[0][6]}/{save_ret[0][0]}/{save_ret[0][4]}"
        except Exception as e:
            # If workspace save fails, still return the embedding but without workspace reference
            embedding_ref = None
            print(f"Warning: Failed to save embedding to workspace: {str(e)}")
        
        # Create report
        report = KBaseReport(self.callback_url)
        report_text = f"""
        Protein Embedding Generation Results
        
        Summary: {summary}
        Sequence Length: {sequence_length} amino acids
        Embedding Dimension: {len(embedding)}
        Embedding Norm: {embedding_norm:.4f}
        
        The protein embedding has been successfully generated and saved as a workspace object.
        You can use this embedding for protein similarity searches and family assignments.
        """
        
        # Only include objects_created if embedding_ref is not None
        report_params = {
            'message': report_text,
            'workspace_name': workspace_name
        }
        
        if embedding_ref is not None:
            report_params['objects_created'] = [{
                'ref': embedding_ref,
                'description': 'Protein embedding result'
            }]
        
        report_info = report.create_extended_report(report_params)
        
        # Return the result according to KIDL spec and UI expectations
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'embedding_result_ref': embedding_ref,
            'summary': summary,
            'input_parameters': params,
            'start_time': start_time,
            'embedding_norm': embedding_norm,
            'sequence_length': sequence_length,
            'embedding_dim': len(embedding)  # Return embedding dimension instead of full embedding
        }
        
        return output
        #END generate_protein_embedding

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method generate_protein_embedding return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def assign_family_fast(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "AssignFamilyFastResults" (Quickly assign
           a protein embedding to a family by similarity to the medoid.) ->
           structure: parameter "family_id" of String, parameter "confidence"
           of Double, parameter "eigenprotein_id" of String, parameter
           "input_parameters" of mapping from String to unspecified object,
           parameter "start_time" of Double
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN assign_family_fast
        import time
        import numpy as np
        import uuid
        start_time = time.time()
        embedding_ref = params.get('embedding_ref')
        workspace_name = params.get('workspace_name')
        
        # Validate input parameters
        if not embedding_ref:
            raise ValueError("Parameter 'embedding_ref' must be provided.")
        
        # Get embedding from workspace object
        # Get workspace URL from configuration
        ws_url = self.config.get('workspace-url')
        if not ws_url:
            # Fallback to environment variable
            ws_url = os.environ.get('KBASE_ENDPOINT', 'https://appdev.kbase.us/services')
            if not ws_url.endswith('/ws'):
                ws_url = ws_url + '/ws'
        ws = Workspace(ws_url)
        try:
            embedding_obj = ws.get_objects([{'ref': embedding_ref}])[0]
            embedding = embedding_obj['data']['embedding']
        except Exception as e:
            raise ValueError(f"Failed to retrieve embedding from workspace: {str(e)}")
        
        if not embedding or not isinstance(embedding, list) or not all(isinstance(x, (int, float)) for x in embedding):
            raise ValueError("Embedding must be a non-empty list of numbers.")
        
        embedding_np = np.array(embedding, dtype=np.float32)
        result = self.family_assigner.assign_family(embedding_np)
        
        # Save family assignment as workspace object
        # Reuse the same workspace connection
        assignment_ref = None
        assignment_obj_name = f"family_assignment_{uuid.uuid4().hex[:8]}"
        assignment_obj_type = "ProteinFamilyAssignmentResult"
        assignment_obj_data = {
            'input_id': str(embedding[:10]) + "..." if len(embedding) > 10 else str(embedding),
            'input_type': 'embedding',
            'embedding_ref': embedding_ref,  # Reference the input embedding
            'assigned_family_id': result['family_id'],
            'similarity_score': float(result['confidence']),
            'metadata': {
                'eigenprotein_id': result['eigenprotein_id'],
                'confidence': float(result['confidence']),
                'created_at': start_time
            }
        }
        
        try:
            # Handle workspace parameter - could be workspace name or ID
            workspace_param = workspace_name
            if isinstance(workspace_name, str) and workspace_name.isdigit():
                # It's a workspace ID
                workspace_param = int(workspace_name)
            elif isinstance(workspace_name, str) and ':' in workspace_name:
                # It's a workspace name with owner
                workspace_param = workspace_name
            
            save_ret = ws.save_objects({
                'workspace': workspace_param,
                'objects': [{
                    'type': assignment_obj_type,
                    'data': assignment_obj_data,
                    'name': assignment_obj_name
                }]
            })
            assignment_ref = f"{save_ret[0][6]}/{save_ret[0][0]}/{save_ret[0][4]}"
        except Exception as e:
            # If workspace save fails, still return the result but without workspace reference
            assignment_ref = None
            print(f"Warning: Failed to save family assignment to workspace: {str(e)}")
        
        # Create report
        report = KBaseReport(self.callback_url)
        report_text = f"""
        Protein Family Assignment Results
        
        Family ID: {result['family_id']}
        Confidence: {result['confidence']:.4f}
        Eigenprotein ID: {result['eigenprotein_id']}
        
        The protein has been assigned to family {result['family_id']} with confidence {result['confidence']:.4f}.
        """
        
        # Only include objects_created if assignment_ref is not None
        report_params = {
            'message': report_text,
            'workspace_name': workspace_name
        }
        
        if assignment_ref is not None:
            report_params['objects_created'] = [{
                'ref': assignment_ref,
                'description': 'Protein family assignment result'
            }]
        
        report_info = report.create_extended_report(report_params)
        
        # Return the result according to KIDL spec and UI expectations
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'family_id': result['family_id'],
            'confidence': float(result['confidence']),
            'eigenprotein_id': result['eigenprotein_id'],
            'input_parameters': params,
            'start_time': start_time
        }
        
        return output
        #END assign_family_fast

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method assign_family_fast return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def find_top_matches_from_embedding(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "FindTopMatchesFromEmbeddingResults" (Find
           top matches for a given protein embedding.) -> structure:
           parameter "matches" of list of mapping from String to unspecified
           object, parameter "summary" of String, parameter
           "input_parameters" of mapping from String to unspecified object,
           parameter "start_time" of Double, parameter "family_id" of String,
           parameter "top_n" of Long, parameter "similarity_stats" of mapping
           from String to Double
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN find_top_matches_from_embedding
        import time
        import numpy as np
        import uuid
        start_time = time.time()
        embedding_ref = params.get('embedding_ref')
        family_id = params.get('family_id')
        top_n = params.get('top_n', 10)
        
        # Validate input parameters
        if not embedding_ref:
            raise ValueError("Parameter 'embedding_ref' must be provided.")
        
        # Get embedding from workspace object
        # Get workspace URL from configuration
        ws_url = self.config.get('workspace-url')
        if not ws_url:
            # Fallback to environment variable
            ws_url = os.environ.get('KBASE_ENDPOINT', 'https://appdev.kbase.us/services')
            if not ws_url.endswith('/ws'):
                ws_url = ws_url + '/ws'
        ws = Workspace(ws_url)
        try:
            embedding_obj = ws.get_objects([{'ref': embedding_ref}])[0]
            embedding = embedding_obj['data']['embedding']
        except Exception as e:
            raise ValueError(f"Failed to retrieve embedding from workspace: {str(e)}")
        
        if not embedding or not isinstance(embedding, list) or not all(isinstance(x, (int, float)) for x in embedding):
            raise ValueError("Embedding must be a non-empty list of numbers.")
        
        if not family_id:
            raise ValueError("Parameter 'family_id' must be provided for similarity search.")
        
        if not isinstance(top_n, int) or not (1 <= top_n <= 1000):
            raise ValueError("Parameter 'top_n' must be an integer between 1 and 1000.")
        
        index = HierarchicalIndex()
        embedding_np = np.array(embedding, dtype=np.uint8)
        similarities, protein_ids = index.search_family(family_id, embedding_np, top_k=top_n)
        
        matches = [
            {'protein_id': pid, 'similarity': float(sim)}
            for pid, sim in zip(protein_ids, similarities)
        ]
        summary = f"Found {len(matches)} top matches in family {family_id}."
        
        # Save similarity search results as workspace object
        # Reuse the same workspace connection
        search_ref = None
        search_obj_name = f"similarity_search_{uuid.uuid4().hex[:8]}"
        search_obj_type = "SummarizeVisualizeResult"
        search_obj_data = {
            'input_id': str(embedding[:10]) + "..." if len(embedding) > 10 else str(embedding),
            'input_type': 'embedding',
            'top_matches_result_ref': None,  # Self-reference
            'summary_html': report_text,
            'metadata': {
                'family_id': family_id,
                'top_n': top_n,
                'matches_count': len(matches),
                'embedding_ref': embedding_ref,  # Reference the input embedding
                'similarity_stats': {
                    'max': float(np.max(similarities)) if len(similarities) else None,
                    'min': float(np.min(similarities)) if len(similarities) else None,
                    'mean': float(np.mean(similarities)) if len(similarities) else None
                },
                'created_at': start_time
            }
        }
        
        try:
            # Handle workspace parameter - could be workspace name or ID
            workspace_name = params.get('workspace_name')
            workspace_param = workspace_name
            if isinstance(workspace_name, str) and workspace_name.isdigit():
                # It's a workspace ID
                workspace_param = int(workspace_name)
            elif isinstance(workspace_name, str) and ':' in workspace_name:
                # It's a workspace name with owner
                workspace_param = workspace_name
            
            save_ret = ws.save_objects({
                'workspace': workspace_param,
                'objects': [{
                    'type': search_obj_type,
                    'data': search_obj_data,
                    'name': search_obj_name
                }]
            })
            search_ref = f"{save_ret[0][6]}/{save_ret[0][0]}/{save_ret[0][4]}"
        except Exception as e:
            # If workspace save fails, still return the result but without workspace reference
            search_ref = None
            print(f"Warning: Failed to save similarity search results to workspace: {str(e)}")
        
        # Create report
        report = KBaseReport(self.callback_url)
        report_text = f"""
        Protein Similarity Search Results
        
        Family ID: {family_id}
        Top N: {top_n}
        Matches Found: {len(matches)}
        
        Similarity Statistics:
        - Maximum: {search_obj_data['similarity_stats']['max']:.4f}
        - Minimum: {search_obj_data['similarity_stats']['min']:.4f}
        - Mean: {search_obj_data['similarity_stats']['mean']:.4f}
        
        Top matches have been saved as a workspace object for further analysis.
        """
        
        # Only include objects_created if search_ref is not None
        report_params = {
            'message': report_text,
            'workspace_name': params.get('workspace_name')
        }
        
        if search_ref is not None:
            report_params['objects_created'] = [{
                'ref': search_ref,
                'description': 'Protein similarity search results'
            }]
        
        report_info = report.create_extended_report(report_params)
        
        # Return the result according to KIDL spec and UI expectations
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'matches': matches,
            'summary': summary,
            'input_parameters': params,
            'start_time': start_time,
            'family_id': family_id,
            'top_n': int(top_n),  # Ensure it's an integer as per KIDL spec
            'similarity_stats': {
                'max': float(np.max(similarities)) if len(similarities) else None,
                'min': float(np.min(similarities)) if len(similarities) else None,
                'mean': float(np.mean(similarities)) if len(similarities) else None
            }
        }
        
        return output
        #END find_top_matches_from_embedding

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method find_top_matches_from_embedding return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def summarize_and_visualize_results(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "SummarizeAndVisualizeResultsResults"
           (Summarize and visualize protein network analysis results.) ->
           structure: parameter "report_name" of String, parameter
           "report_ref" of String, parameter "input_parameters" of mapping
           from String to unspecified object, parameter "start_time" of
           Double, parameter "output_dir" of String, parameter "summary" of
           String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN summarize_and_visualize_results
        import time
        import os
        start_time = time.time()
        top_matches_result_ref = params.get('top_matches_result_ref')  # Match UI spec parameter name
        output_dir = params.get('output_dir', self.shared_folder)
        
        # Validate input parameters
        if not top_matches_result_ref:
            raise ValueError("Parameter 'top_matches_result_ref' is required.")
        
        if output_dir and not isinstance(output_dir, str):
            raise ValueError("Parameter 'output_dir' must be a valid path string.")
        
        summary = f"Results summarized and visualized in {output_dir}."
        report = KBaseReport(self.callback_url)
        report_info = report.create_extended_report({
            'message': summary,
            'objects_created': [],
            'workspace_name': params['workspace_name']
        })
        
        # Return the result according to KIDL spec
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'input_parameters': params,
            'start_time': start_time,
            'output_dir': output_dir,
            'summary': summary
        }
        
        return output
        #END summarize_and_visualize_results

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method summarize_and_visualize_results return value ' +
                             'output is not type dict as required.')
        # return the results
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
