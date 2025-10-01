# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import logging
import time
import uuid
import numpy as np
from typing import Dict, Any, List, Optional, Union
from pathlib import Path

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBUtilLibClient import KBUtilLib

# Import internal orchestrator and config (lazy via package __init__)
from .src import WorkflowOrchestrator, PipelineConfig

logger = logging.getLogger(__name__)
#END_HEADER


class kbase_protein_query_module:
    '''
    Module Name:
    kbase_protein_query_module

    Module Description:
    A KBase module: kbase_protein_query_module

This module provides comprehensive protein query analysis capabilities using UniProt IDs as the canonical identifier:

COMPREHENSIVE ANALYSIS WORKFLOW:
1. CheckProteinExistence: Verify protein exists using UniProt ID, optionally generate embedding
2. GenerateProteinEmbeddings: Create embeddings from sequence input or protein check results
3. AssignProteinFamily: Assign proteins to families using similarity to centroids
4. FindTopMatches: Perform similarity search within families
5. SummarizeAndVisualize: Generate comprehensive HTML reports with network analysis
6. RunProteinQueryAnalysis: Unified pipeline for comprehensive protein analysis

ADVANCED CAPABILITIES:
- UniProt ID canonical identifier system (exact match only)
- ESM-2 protein language model for embedding generation
- Efficient FAISS-based similarity search and clustering
- Family assignment using binary centroid similarity
- Comprehensive metadata storage and retrieval
- HTML report generation with network visualization
- Workspace object management for downstream analysis
- Bioinformatics integration with protein databases
- Network analysis and protein relationship mapping
- Advanced similarity metrics and statistical analysis
- Modular pipeline architecture with configurable stages
- Real-time performance monitoring and error handling

Authors: Vibhav Setlur
Contact: https://kbase.us/contact-us/
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "2.0.0"
    GIT_URL = "https://github.com/VibhavSetlur/kbase_protein_query_module.git"
    GIT_COMMIT_HASH = "36203034384319fef4abcbb4d82c8a8e3a07f512"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config or {}
        self.callback_url = os.environ.get('SDK_CALLBACK_URL', 'http://localhost:0')
        self.shared_folder = self.config.get('scratch', '/tmp')
        
        # Initialize basic logging
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        
        # Initialize scalability components (inline to avoid compilation overwrites)
        try:
            from .src.core.resource_manager import ResourceManager, ResourceLimits
            from .src.core.parallel_processor import ParallelProcessor
            
            # Initialize resource manager with server-aware percentage limits
            resource_limits = ResourceLimits(
                max_memory_percent=60.0,    # Use max 60% of server memory
                max_cpu_percent=70.0,       # Use max 70% of server CPU  
                max_disk_percent=80.0,      # Use max 80% of server disk
                batch_size_proteins=500,    # Conservative batch sizes for shared servers
                max_concurrent_tasks=2,     # Limited concurrency for server environments
                server_safety_margin=0.2    # 20% safety margin for other processes
            )
            
            self.resource_manager = ResourceManager(resource_limits)
            self.parallel_processor = ParallelProcessor(
                max_workers=4,
                resource_manager=self.resource_manager
            )
            
            logger.info("Scalability components initialized successfully")
            
        except Exception as e:
            logger.warning(f"Could not initialize scalability components: {e}")
            self.resource_manager = None
            self.parallel_processor = None
        
        # Initialize KBase clients with KBUtilLib
        try:
            self.dfu = DataFileUtil(self.callback_url)
            self.kb_util = KBUtilLib(self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'), 
                                   scratch=self.shared_folder)
            logger.info("KBUtilLib client initialized successfully")
        except Exception as e:
            logger.warning(f"Could not initialize clients: {e}")
            self.dfu = None
            self.kb_util = None
        
        # Initialize components (simplified)
        self.existence_checker = None
        self.family_assigner = None
        self.embedding_generator = None
        self.hierarchical_index = None
        self.network_builder = None
        self.html_report_generator = None
        self.workflow_orchestrator = None
        
        #END_CONSTRUCTOR
        pass

    def _create_mock_pipeline_results(self, input_proteins, analysis_stages):
        """Create mock pipeline results for testing purposes."""
        pipeline_results = {}
        
        for stage in analysis_stages:
            if stage == 'embedding_generation':
                pipeline_results[stage] = {
                    'status': 'completed',
                    'results': [{'protein_id': protein_id, 'embedding': [0.1] * 1280} for protein_id in input_proteins]
                }
            elif stage == 'family_assignment':
                pipeline_results[stage] = {
                    'status': 'completed',
                    'results': [{'protein_id': protein_id, 'family_id': f'family_{i}', 'confidence': 0.8} for i, protein_id in enumerate(input_proteins)]
                }
            elif stage == 'similarity_search':
                pipeline_results[stage] = {
                    'status': 'completed',
                    'results': [{'protein_id': protein_id, 'matches': [{'id': f'match_{i}_{j}', 'similarity': 0.9} for j in range(5)]} for i, protein_id in enumerate(input_proteins)]
                }
            else:
                pipeline_results[stage] = {
                    'status': 'completed',
                    'results': [{'protein_id': protein_id, 'result': f'mock_{stage}_result'} for protein_id in input_proteins]
                }
        
        return pipeline_results

    def _log_with_kbutillib(self, level, message, context=None):
        """
        Log message using KBUtilLib if available, otherwise use standard logging.
        """
        try:
            if hasattr(self, 'kb_util') and self.kb_util and hasattr(self.kb_util, 'log'):
                self.kb_util.log(level, message, context)
            else:
                # Fallback to standard logging
                if level == 'ERROR':
                    logger.error(message)
                elif level == 'WARNING':
                    logger.warning(message)
                elif level == 'INFO':
                    logger.info(message)
                else:
                    logger.debug(message)
        except Exception as e:
            # Fallback to standard logging if KBUtilLib logging fails
            logger.error(f"KBUtilLib logging failed: {e}, message: {message}")

    def _deprecated_delegate(self, ctx, params, origin_method: str):
        raise NotImplementedError("Legacy methods are removed. Use run_protein_query_analysis only.")

    def check_protein_existence(self, ctx, params):
        raise NotImplementedError("Legacy endpoint removed. Use run_protein_query_analysis.")

    def generate_protein_embedding(self, ctx, params):
        raise NotImplementedError("Legacy endpoint removed. Use run_protein_query_analysis.")

    def assign_family_fast(self, ctx, params):
        raise NotImplementedError("Legacy endpoint removed. Use run_protein_query_analysis.")

    def find_top_matches_from_embedding(self, ctx, params):
        raise NotImplementedError("Legacy endpoint removed. Use run_protein_query_analysis.")

    def summarize_and_visualize_results(self, ctx, params):
        raise NotImplementedError("Legacy endpoint removed. Use run_protein_query_analysis.")

    def run_protein_query_analysis(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ProteinQueryAnalysisResults" (Unified
           Protein Query Analysis Pipeline This method provides a single
           entry point for comprehensive protein analysis, supporting
           multiple input types and configurable analysis stages.) ->
           structure: parameter "report_name" of String, parameter
           "report_ref" of String, parameter "analysis_result_ref" of String,
           parameter "summary" of String, parameter "input_parameters" of
           mapping from String to unspecified object, parameter "start_time"
           of Double, parameter "html_report_path" of String, parameter
           "protein_count" of Long, parameter "stages_completed" of list of
           String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_protein_query_analysis
        start_time = time.time()
        # Initialize safe defaults for exception paths
        analysis_name = (params or {}).get('analysis_name') or 'protein_analysis'
        workspace_name = (params or {}).get('workspace_name')
        try:
            try:
                workspace_name = params.get('workspace_name') or (
                    ctx.get('provenance', [{}])[0].get('ws_name') if isinstance(ctx.get('provenance'), list) and ctx.get('provenance') else ctx.get('ws_name')
                )
            except Exception:
                workspace_name = params.get('workspace_name')
            if not workspace_name:
                raise ValueError("workspace_name is required")

            input_type = params.get('input_type')
            input_kwargs = {}
            
            # Handle conditional input parameters based on input_type
            if input_type == 'protein_input':
                protein_input = params.get('protein_input')
                if not protein_input:
                    raise ValueError("Protein input is required when input_type is 'protein_input'")
                # Parse protein sequences from input
                input_proteins = []
                if isinstance(protein_input, str):
                    # Handle both direct sequences and FASTA format
                    lines = protein_input.strip().split('\n')
                    for line in lines:
                        line = line.strip()
                        if line and not line.startswith('>'):
                            input_proteins.append(line)
                    if not input_proteins:
                        input_proteins = [protein_input]  # Single sequence
                elif isinstance(protein_input, list):
                    input_proteins = protein_input
                else:
                    raise ValueError("Protein input must be a string or list")
                input_kwargs['input_proteins'] = input_proteins
                
            elif input_type == 'uniprot_ids':
                uniprot_ids = params.get('uniprot_ids', [])
                if not uniprot_ids:
                    raise ValueError("UniProt IDs are required when input_type is 'uniprot_ids'")
                # Handle both single ID and list of IDs
                if isinstance(uniprot_ids, str):
                    uniprot_ids = [uniprot_ids.strip() for uniprot_ids in uniprot_ids.split(',') if uniprot_ids.strip()]
                input_kwargs['input_proteins'] = uniprot_ids
                
            elif input_type == 'workspace_object':
                workspace_object = params.get('workspace_object')
                if not workspace_object:
                    raise ValueError("Workspace object is required when input_type is 'workspace_object'")
                input_kwargs['workspace_object_ref'] = workspace_object
            elif input_type == 'direct_sequences':
                # Backward compatibility: map to new schema
                sequences = params.get('direct_sequences')
                if not sequences:
                    raise ValueError("direct_sequences requires sequence input")
                if isinstance(sequences, str):
                    sequences = [sequences]
                input_kwargs['input_proteins'] = sequences
                input_type = 'protein_input'
                
            else:
                raise ValueError(f"Invalid input_type: {input_type}. Must be one of: protein_input, uniprot_ids, workspace_object")

            # Get analysis name for output
            analysis_name = params.get('analysis_name')
            if not analysis_name:
                raise ValueError("analysis_name is required")
            
            # Stage selection is now managed by the orchestrator/config. Do not hardcode legacy stage names.
            enabled_stages = params.get('enabled_stages') or []
            analysis_config = params.get('analysis_config') or {}
            storage_config = params.get('storage_config') or {}
            output_config = params.get('output_config') or {}

            config = PipelineConfig(
                enabled_stages=enabled_stages,
                stage_configs=analysis_config,
                storage_config=storage_config,
                similarity_config=storage_config,
                **input_kwargs
            )
            setattr(config, 'output_dir', output_config.get('output_dir') or self.shared_folder)
            if 'generate_html_report' in output_config:
                config.generate_html_report = bool(output_config.get('generate_html_report'))
            if 'generate_network_visualization' in output_config:
                config.generate_network_visualization = bool(output_config.get('generate_network_visualization'))

            auth_token = ctx.get('token') if isinstance(ctx, dict) else os.environ.get('KB_AUTH_TOKEN')
            config.workspace_url = self.callback_url
            config.auth_token = auth_token
            try:
                workspace_client = self.kb_util.get_workspace_client() if (self.kb_util and hasattr(self.kb_util, 'get_workspace_client')) else Workspace(self.callback_url)
            except Exception:
                workspace_client = None
            setattr(config, 'workspace_client', workspace_client)
            setattr(config, 'workspace_name', workspace_name)

            workflow = WorkflowOrchestrator(config=config, kb_util=self.kb_util)
            wf_result = workflow.execute()
            final = wf_result.final_output or {}
            stages_completed = wf_result.stages_completed or final.get('stages_completed', [])
            
            # Shock upload info from orchestrator
            shock_info = final.get('shock', {}) or {}
            shock_id = shock_info.get('shock_id') or shock_info.get('id') or ''
            shock_url = shock_info.get('download_url') or ''

            # Ensure report fields are always present for Narrative integration
            report_name = final.get('report_name') or final.get('report', {}).get('name') or f"{analysis_name}_report"
            report_ref = final.get('report_ref') or final.get('report', {}).get('ref') or ''

            # Always create a concise report referencing the Shock archive
            if not report_ref:
                try:
                    report_client = KBaseReport(self.callback_url)
                    report_message = final.get('summary') or 'Protein query analysis completed'
                    if shock_id:
                        report_message += f"\n\nOutputs archived in Shock (node {shock_id})."
                    if shock_url:
                        report_message += f"\nDownload: {shock_url}"
                    report_info = report_client.create_extended_report({
                        'message': report_message,
                        'workspace_name': workspace_name,
                        'report_object_name': f"{analysis_name}_report",
                        'objects_created': []
                    })
                    report_name = report_info.get('name', report_name)
                    report_ref = report_info.get('ref', report_ref)
                except Exception as e:
                    logger.warning(f"Failed to create report: {e}")

            normalized = {
                'job_id': final.get('job_id') or wf_result.run_id,
                'analysis_result_ref': '',
                'summary': final.get('summary') or 'Protein query analysis completed',
                'input_parameters': params,
                'start_time': start_time,
                'output_directory': '',  # Don't show internal paths to users
                'general_info_dir': '',  # Don't show internal paths to users
                'network_analysis_dir': '',  # Don't show internal paths to users
                'sequence_analysis_dir': '',  # Don't show internal paths to users
                'embeddings_file_path': '',  # Don't show internal file paths
                'top_proteins_csv_path': '',  # Don't show internal file paths
                'html_report_path': '',  # No HTML reports as requested
                'protein_count': final.get('protein_count') or 0,
                'stages_completed': stages_completed,
                'report_name': report_name,
                'report_ref': report_ref,
                'shock_id': shock_id,
                'shock_url': shock_url
            }
            # Final safety: ensure non-empty placeholders
            if not normalized['report_name']:
                normalized['report_name'] = f"{analysis_name}_report"
            if not normalized['report_ref']:
                normalized['report_ref'] = f"report_{analysis_name}"
            return [normalized]
        except Exception as e:
            try:
                if hasattr(self, '_log_with_kbutillib'):
                    self._log_with_kbutillib('ERROR', f"run_protein_query_analysis failed: {e}")
                else:
                    logger.error(f"run_protein_query_analysis failed: {e}")
            except Exception:
                logger.error(f"run_protein_query_analysis failed: {e}")
            # Create error report without HTML
            try:
                report_client = KBaseReport(self.callback_url)
                error_report = report_client.create_extended_report({
                    'message': f'Analysis failed: {str(e)}',
                    'objects_created': [],
                    'workspace_name': workspace_name,
                    'report_object_name': f"{analysis_name}_error_report"
                })
            except Exception as report_error:
                logger.error(f"Failed to create error report: {report_error}")
                error_report = {
                    'name': f"{analysis_name}_error_report",
                    'ref': f'report_error_{int(time.time())}'
                }
            
            return [{
                'job_id': f'error_{int(time.time())}',
                'analysis_result_ref': 'error',
                'summary': f'Analysis failed: {str(e)}',
                'input_parameters': params,
                'start_time': start_time,
                'protein_count': 0,
                'stages_completed': [],
                'report_name': error_report['name'],
                'report_ref': error_report['ref']
            }]
        except Exception as e:
            logger.error(f"Protein query analysis failed: {e}")
            # Create error report without HTML
            try:
                error_report_client = KBaseReport(self.callback_url)
                error_report = error_report_client.create_extended_report({
                    'message': f'Analysis failed: {str(e)}',
                    'objects_created': [],
                    'workspace_name': workspace_name,
                    'report_object_name': 'error_report'
                })
            except Exception as report_error:
                logger.error(f"Failed to create error report: {report_error}")
                error_report = {
                    'name': 'error_report',
                    'ref': 'error_report_ref'
                }
            
            output = {
                'report_name': error_report['name'],
                'report_ref': error_report['ref'],
                'analysis_result_ref': 'error',
                'summary': f'Analysis failed: {str(e)}',
                'input_parameters': params,
                'start_time': start_time,
                'html_report_path': '',  # No HTML error reports
                'protein_count': 0,
                'stages_completed': [],
                'output_directory': self.shared_folder,
                'general_info_dir': self.shared_folder,
                'network_analysis_dir': self.shared_folder,
                'sequence_analysis_dir': self.shared_folder,
                'job_id': f'error_job_{int(time.time())}',
                'embeddings_file_path': '',
                'top_proteins_csv_path': ''
            }
        #END run_protein_query_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_protein_query_analysis return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def run_kbase_protein_query_module(self, ctx, params):
        """
        Legacy method name for backward compatibility.
        Alias for run_protein_query_analysis.
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ProteinQueryAnalysisResults" (Unified
           Protein Query Analysis Pipeline This method provides a single
           entry point for comprehensive protein analysis, supporting
           multiple input types and configurable analysis stages.) ->
           structure: parameter "report_name" of String, parameter
           "report_ref" of String, parameter "analysis_result_ref" of String,
           parameter "summary" of String, parameter "input_parameters" of
           mapping from String to unspecified object, parameter "start_time"
           of Double, parameter "html_report_path" of String, parameter
           "protein_count" of Long, parameter "stages_completed" of list of
           String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kbase_protein_query_module
        self._log_with_kbutillib('WARNING', 'run_kbase_protein_query_module is deprecated; using run_protein_query_analysis')
        return self.run_protein_query_analysis(ctx, params)
        #END run_kbase_protein_query_module

    def get_available_analyses(self, ctx):
        """
        Get list of available analyses for the UI.
        
        This method returns information about all available analyses that can be
        selected by users in the KBase Narrative interface.
        
        :returns: instance of type "GetAvailableAnalysesResults" ->
           structure: parameter "available_analyses" of mapping from String to
           unspecified object, parameter "summary" of String
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
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
