# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import logging
import time
import shutil
from typing import Dict, Any, List, Optional, Union

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBUtilLibClient import KBUtilLib

# Import the new architecture
from .src import WorkflowOrchestrator
from .src.analysis import get_enabled_analyses

logger = logging.getLogger(__name__)
#END_HEADER


class kbase_protein_query_module:
    '''
    Module Name:
    kbase_protein_query_module

    Module Description:
    A KBase module: kbase_protein_query_module

Protein query and analysis module with comprehensive network analysis capabilities.
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "2.0.0"
    GIT_URL = "https://github.com/VibhavSetlur/kbase_protein_query_module.git"
    GIT_COMMIT_HASH = "5898726743862cb2bf4084ca7990e1c4e337d463"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config or {}
        self.callback_url = os.environ.get('SDK_CALLBACK_URL', 'http://localhost:0')
        self.shared_folder = self.config.get('scratch', '/tmp')
        
        # Initialize KBase clients
        self.dfu = DataFileUtil(self.callback_url)
        self.kb_util = KBUtilLib(self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'), 
                               scratch=self.shared_folder)
        #END_CONSTRUCTOR
        pass


    def run_protein_query_analysis(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ProteinQueryAnalysisResults" ->
           structure: parameter "job_id" of String, parameter
           "analysis_result_ref" of String, parameter "summary" of String,
           parameter "input_parameters" of mapping from String to unspecified
           object, parameter "start_time" of Float, parameter
           "protein_count" of Int, parameter "stages_completed" of list of
           String, parameter "report_name" of String, parameter "report_ref"
           of String, parameter "shock_id" of String, parameter "shock_url"
           of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_protein_query_analysis
        start_time = time.time()
        
        # Extract workspace information
        workspace_name = params.get('workspace_name')
        if not workspace_name:
            raise ValueError("workspace_name is required")

        # Extract analysis name
        analysis_name = params.get('analysis_name', 'protein_analysis')
        
        input_type = params.get('input_type')
        if not input_type:
            logger.error("input_type is required")
            raise ValueError("input_type is required")
        
        # Prepare configuration - simple dict, no PipelineConfig class
        output_config = params.get('output_config', {})
        output_dir = output_config.get('output_dir', '/tmp/protein_query_output')
        os.makedirs(output_dir, exist_ok=True)
        
        # Get analyses to run - handle both 'analyses' and 'analysis_stages' parameter names
        selected_analyses = params.get('analyses') or params.get('analysis_stages') or ['network_analysis']
        
        # Resolve module root directory for data paths
        # In Docker: /kb/module, in local: try to find project root
        # Note: Data path resolution in NetworkAnalysis will check for reference data
        # in environment variables (KB_DATA_DIR, KB_REFDATA_DIR) and common mount points
        module_root = os.environ.get('KB_MODULE_DIR', '/kb/module')
        if not os.path.exists(module_root):
            # Try to find module root by looking for lib directory
            # Start from this file and navigate up
            current_file = os.path.abspath(__file__)
            # From: lib/kbase_protein_query_module/kbase_protein_query_moduleImpl.py
            # To: module root (parent of lib/)
            lib_dir = os.path.dirname(os.path.dirname(current_file))
            module_root = os.path.dirname(lib_dir)
        
        logger.info(f"Using module_root: {module_root}")
        
        # Log reference data environment variables if set
        kb_data_dir = os.environ.get('KB_DATA_DIR')
        kb_refdata_dir = os.environ.get('KB_REFDATA_DIR')
        if kb_data_dir:
            logger.info(f"KB_DATA_DIR environment variable set: {kb_data_dir}")
        if kb_refdata_dir:
            logger.info(f"KB_REFDATA_DIR environment variable set: {kb_refdata_dir}")
        
        workflow_config = {
            'output_dir': output_dir,
            'workspace_name': workspace_name,
            'selected_analyses': selected_analyses,
            'module_root': module_root,
            # Data paths - will be resolved by NetworkAnalysis with reference data support
            # NetworkAnalysis will check:
            # 1. Environment variables (KB_DATA_DIR, KB_REFDATA_DIR)
            # 2. Common reference data mount points (/data, /kb/data, /refdata)
            # 3. Module root data directory (for test data in Docker image)
            'embeddings_file': params.get('embeddings_file') or 'data/embeddings/embeddings.tsv',
            'index_path': params.get('index_path') or 'data/indexes',
        }
        
        # Execute workflow through orchestrator
        workflow = WorkflowOrchestrator(config=workflow_config, kb_util=self.kb_util)
        result = workflow.run_workflow(
            input_data=params,
            output_dir=output_dir,
            workspace_name=workspace_name,
            selected_analyses=selected_analyses,
        )
        
        # Check if workflow succeeded - if not, raise an error to fail the job
        if not result.get('success', False):
            error_message = result.get('error_message', 'Workflow execution failed')
            failed_analyses = result.get('failed_analyses', [])
            if failed_analyses:
                error_message = f"Workflow failed: {', '.join(failed_analyses)} analysis(es) failed. {error_message}"
            logger.error(error_message)
            raise RuntimeError(error_message)
        
        # Create a zip archive of the output directory and upload to Shock
        try:
            archive_base = os.path.join(self.shared_folder, "output")
            dir_to_zip = result.get('output_directory') or output_dir
            zip_path = shutil.make_archive(archive_base, 'zip', dir_to_zip)
            shock_info = self.dfu.file_to_shock({'file_path': zip_path, 'make_handle': 1})
            if not shock_info or not isinstance(shock_info, dict) or not shock_info.get('shock_id'):
                logger.warning("DataFileUtil.file_to_shock returned invalid response, creating report without file links")
                shock_info = None
        except Exception as e:
            logger.error(f"Failed to upload results to Shock: {e}")
            logger.warning("Creating report without file links due to Shock upload failure")
            shock_info = None
        
        # Create KBase report with file links 
        report_info = self._create_kbase_report(result, analysis_name, workspace_name, shock_info)
        
        # Return KBase report references for narrative display
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }
        
        #END run_protein_query_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_protein_query_analysis return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def get_available_analyses(self, ctx):
        """Return enabled analyses in a simple structure for the Narrative/tests."""
        # ctx is the context object
        # return variables are: output
        #BEGIN get_available_analyses
        try:
            analyses = list((get_enabled_analyses() or {}).keys())
        except Exception as e:
            logger.error(f"Failed to retrieve available analyses: {e}")
            analyses = []
        output = {"analyses": analyses}
        #END get_available_analyses

        if not isinstance(output, dict):
            raise ValueError('Method get_available_analyses return value output is not type dict as required.')
        return [output]

    
    def status(self, ctx):
        """
        :returns: instance of type "StatusResults" -> structure: parameter
           "state" of String, parameter "message" of String, parameter
           "version" of String, parameter "git_url" of String, parameter
           "git_commit_hash" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN status
        output = {
            'state': 'OK',
            'message': 'Protein Query Module is operational',
            'version': self.VERSION,
            'git_url': self.GIT_URL,
            'git_commit_hash': self.GIT_COMMIT_HASH
        }
        #END status

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method status return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    
    def _create_kbase_report(self, result, analysis_name, workspace_name, shock_info=None):
        """Create KBase report."""
        # Create report using KBaseReport client
        report_client = KBaseReport(self.callback_url)
        
        # Prepare file links with proper None checking
        file_links = []
        if shock_info and isinstance(shock_info, dict) and shock_info.get('shock_id'):
            # Safely extract file name with proper None handling
            handle_info = shock_info.get('handle', {})
            if isinstance(handle_info, dict):
                file_name = handle_info.get('file_name', 'output.zip')
            else:
                file_name = 'output.zip'
            
            file_links.append({
                'shock_id': shock_info['shock_id'],
                'name': file_name,
                'label': 'Analysis Results'
            })
        
        # Create comprehensive report message
        workflow_success = result.get('success', False)
        summary = (result.get('final_output') or {}).get('summary', 'Analysis completed')
        protein_count = (result.get('final_output') or {}).get('protein_count', 0)
        stages_completed = result.get('analyses_completed', [])
        error_message = result.get('error_message', '')
        
        if workflow_success:
            message = f"Protein Query Analysis completed successfully.\n\n"
        else:
            message = f"Protein Query Analysis FAILED.\n\n"
        
        message += f"Summary: {summary}\n"
        message += f"Proteins processed: {protein_count}\n"
        if stages_completed:
            message += f"Stages completed: {', '.join(stages_completed)}\n"
        if error_message:
            message += f"Error: {error_message}\n"
        if file_links:
            message += f"Output files: {len(file_links)} file(s) available for download"
        
        # Create report following KBase documentation
        report_params = {
            'message': message,
            'workspace_name': workspace_name,
            'report_object_name': f"{analysis_name}_report",
            'objects_created': [],
            'warnings': [],
            'file_links': file_links,
        }
        # Only set HTML fields if we actually add HTML links
        # Avoids Narrative error "Report not found for index 0"
        html_links: List[Dict[str, Any]] = []
        if html_links:
            report_params['html_links'] = html_links
            report_params['direct_html_link_index'] = 0
            report_params['html_window_height'] = 333
        
        report_info = report_client.create_extended_report(report_params)
        
        # Return references 
        return {
            'name': report_info['name'],
            'ref': report_info['ref']
        }
    
