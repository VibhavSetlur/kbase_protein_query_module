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
from .src import WorkflowOrchestrator, PipelineConfig

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
        
        # Validate input parameters
        self._validate_input_parameters(params)
        
        # Create pipeline configuration
        pipeline_config = self._create_pipeline_config(params, workspace_name, ctx)
        
        # Execute workflow through orchestrator
        workflow = WorkflowOrchestrator(config=pipeline_config, kb_util=self.kb_util)
        output_dir = getattr(pipeline_config, 'output_dir', '/tmp/protein_query_output')
        selected_analyses = params.get('analysis_stages')
        result = workflow.run_workflow(params, output_dir, selected_analyses=selected_analyses)
        
        # Create a zip archive of the output directory and upload to Shock (align with KBase examples)
        try:
            archive_base = os.path.join(self.shared_folder, f"{analysis_name}_results")
            zip_path = shutil.make_archive(archive_base, 'zip', output_dir)
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
        """
        :returns: instance of type "GetAvailableAnalysesResults" ->
           structure: parameter "available_analyses" of mapping from String
           to unspecified object, parameter "summary" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN get_available_analyses
        
        try:
            # Import the analysis configuration
            from .src.analysis import get_enabled_analyses
            
            # Get available analyses
            available_analyses = get_enabled_analyses()
            
            # Create summary
            analysis_count = len(available_analyses)
            analysis_names = list(available_analyses.keys())
            summary = f"Found {analysis_count} available analyses: {', '.join(analysis_names)}"
            
            output = {
                'available_analyses': available_analyses,
                'summary': summary
            }
            
            logger.info(f"Retrieved {analysis_count} available analyses for UI")
            
        except Exception as e:
            logger.error(f"Error getting available analyses: {e}")
            output = {
                'available_analyses': {},
                'summary': f"Error retrieving available analyses: {str(e)}"
            }
        
        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method get_available_analyses return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
        #END get_available_analyses

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
    
    def _validate_input_parameters(self, params):
        """Validate input parameters and raise exceptions for invalid inputs."""
        input_type = params.get('input_type', 'protein_input')
        
        if input_type == 'protein_input':
            protein_input = params.get('protein_input')
            if not protein_input:
                raise ValueError("protein_input is required for protein_input type")
        
        elif input_type == 'uniprot_ids':
            uniprot_ids = params.get('uniprot_ids', [])
            if not uniprot_ids:
                raise ValueError("uniprot_ids is required for uniprot_ids input type")
        
        else:
            raise ValueError(f"Invalid input_type: {input_type}. Must be one of: protein_input, uniprot_ids")
    
    def _create_pipeline_config(self, params, workspace_name, ctx):
        """Create pipeline configuration from parameters."""
        from .src import PipelineConfig
        
        # Handle output_config parameter from narrative
        output_config = params.get('output_config', {})
        output_dir = output_config.get('output_dir', '/tmp/protein_query_output')
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Get input proteins based on input type
        input_proteins = []
        if params.get('input_type') == 'protein_input':
            input_proteins = params.get('protein_sequences', [])
        elif params.get('input_type') == 'uniprot_ids':
            input_proteins = params.get('uniprot_ids', [])
        
        # Create PipelineConfig with only supported parameters
        return PipelineConfig(
            input_proteins=input_proteins,
            input_type=params.get('input_type', 'protein_input'),
            output_dir=output_dir,
            workspace_name=workspace_name,
            selected_analyses=params.get('analysis_stages', ['network_analysis'])
        )
    
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
        summary = result.final_output.get('summary', 'Analysis completed successfully')
        protein_count = result.final_output.get('protein_count', 0)
        stages_completed = result.analyses_completed
        
        message = f"Protein Query Analysis completed successfully.\n\n"
        message += f"Summary: {summary}\n"
        message += f"Proteins processed: {protein_count}\n"
        message += f"Stages completed: {', '.join(stages_completed)}\n"
        if file_links:
            message += f"Output files: {len(file_links)} file(s) available for download"
        
        # Create report following KBase documentation exactly
        report_params = {
            'message': message,
            'workspace_name': workspace_name,
            'report_object_name': f"{analysis_name}_report",
            'objects_created': [],
            'warnings': [],
            'file_links': file_links,
            'html_links': [],
            'direct_html_link_index': 0,
            'html_window_height': 333
        }
        
        report_info = report_client.create_extended_report(report_params)
        
        # Return references 
        return {
            'name': report_info['name'],
            'ref': report_info['ref']
        }
    
