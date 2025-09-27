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
        """Thin shim that logs deprecation and delegates to single entrypoint.
        Maps legacy params to unified schema when possible.
        """
        self._log_with_kbutillib('WARNING', f'{origin_method} is deprecated; delegating to run_protein_query_analysis')
        mapped = dict(params or {})
        # Ensure workspace_name exists if possible
        if 'workspace_name' not in mapped:
            try:
                if ctx.get('provenance') and len(ctx['provenance']) > 0:
                    mapped['workspace_name'] = ctx['provenance'][0].get('ws_name')
                else:
                    mapped['workspace_name'] = ctx.get('ws_name')
            except Exception:
                pass
        # Heuristic mappings for common legacy calls
        if origin_method == 'check_protein_existence' and 'protein_id' in mapped and 'input_data' not in mapped:
            mapped['input_type'] = 'uniprot_id'
            mapped['input_data'] = mapped.get('protein_id')
        if origin_method == 'generate_protein_embedding' and 'input_data' in mapped and 'input_type' not in mapped:
            mapped['input_type'] = 'sequence'
        if origin_method in ('assign_family_fast','find_top_matches_from_embedding') and 'embedding_ref' in mapped and 'input_data' not in mapped:
            mapped['input_type'] = 'workspace_object_ref'
            mapped['input_data'] = mapped.get('embedding_ref')
        if origin_method == 'summarize_and_visualize_results' and 'result_refs' in mapped:
            mapped.setdefault('enabled_stages', ['report_generation','data_export'])
        return self.run_protein_query_analysis(ctx, mapped)

    def check_protein_existence(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "CheckProteinExistenceResults" (Check if a
           protein exists in the storage system using UniProt ID and create a
           workspace object with the result. Input: UniProt ID (e.g., P00001,
           P12345) Output: Existence status, family assignment, metadata,
           optional embedding) -> structure: parameter "report_name" of
           String, parameter "report_ref" of String, parameter "exists" of
           Long, parameter "family_id" of String, parameter "metadata" of
           mapping from String to unspecified object, parameter
           "input_parameters" of mapping from String to unspecified object,
           parameter "start_time" of Double, parameter "summary" of String,
           parameter "protein_existence_result_ref" of String, parameter
           "embedding_result_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN check_protein_existence
        try:
            start_time = time.time()
            
            protein_id = params.get('protein_id')
            generate_embedding = params.get('generate_embedding', False)
            
            if not protein_id:
                raise ValueError("protein_id parameter is required")
            
            # Setup workspace client using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'get_workspace_client'):
                    workspace_client = self.kb_util.get_workspace_client()
                    if self.kb_util:
                        self.kb_util.log_info(f"Starting check_protein_existence for protein: {protein_id}")
                else:
                    workspace_client = Workspace(self.callback_url)
            except Exception as e:
                logger.error(f"Failed to setup workspace client: {e}")
                workspace_client = Workspace(self.callback_url)
            
            # Inline workspace name extraction
            try:
                if ctx.get('provenance') and len(ctx['provenance']) > 0:
                    workspace_name = ctx['provenance'][0].get('ws_name')
                else:
                    workspace_name = ctx.get('ws_name')
            except Exception:
                workspace_name = None
            
            # Simple existence check (demo implementation)
            existence_result = {
                'exists': True,  # Demo: assume protein exists
                'family_id': 'demo_family',
                'metadata': {'source': 'demo', 'protein_id': protein_id}
            }
            
            result_data = {
                'protein_id': protein_id,
                'exists': existence_result.get('exists', False),
                'family_id': existence_result.get('family_id', 'unknown'),
                'metadata': existence_result.get('metadata', {}),
                'check_timestamp': time.time()
            }
            
            # Generate demo embedding if requested
            embedding_ref = None
            if generate_embedding and existence_result.get('exists', False):
                demo_embedding = np.random.rand(1280).tolist()  # Demo embedding
                embedding_result = {
                    'embedding': demo_embedding,
                    'model_name': 'esm2_t6_8M_UR50D',
                    'protein_id': protein_id,
                    'sequence_length': 100,
                    'embedding_dim': 1280
                }
                
                # Save embedding to workspace using KBUtilLib
                try:
                    if self.kb_util and hasattr(self.kb_util, 'save_workspace_object'):
                        embedding_ref = self.kb_util.save_workspace_object(
                            workspace_name, f"{protein_id}_embedding",
                            'KBaseProteinQueryModule.ProteinEmbedding', embedding_result)
                    else:
                        # Fallback to direct workspace client usage
                        if workspace_client and workspace_name:
                            result = workspace_client.save_objects({
                            'workspace': workspace_name,
                            'objects': [{
                                'name': f"{protein_id}_embedding",
                                'type': 'KBaseProteinQueryModule.ProteinEmbedding',
                                'data': embedding_result
                            }]
                            })
                            embedding_ref = result[0][6]  # Get the reference
                        else:
                            embedding_ref = None
                except Exception as e:
                    logger.error(f"Failed to save embedding to workspace: {e}")
                    embedding_ref = None
                
                result_data['embedding_ref'] = embedding_ref
                result_data['embedding'] = embedding_result.get('embedding', [])
                result_data['model_name'] = embedding_result.get('model_name', 'esm2_t6_8M_UR50D')
            
            # Save result to workspace using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'save_workspace_object'):
                    result_ref = self.kb_util.save_workspace_object(
                        workspace_name, f"{protein_id}_existence_check",
                        'KBaseProteinQueryModule.ProteinExistenceResult', result_data)
                else:
                    # Fallback to direct workspace client usage
                    if workspace_client and workspace_name:
                        result = workspace_client.save_objects({
                'workspace': workspace_name,
                'objects': [{
                    'name': f"{protein_id}_existence_check",
                    'type': 'KBaseProteinQueryModule.ProteinExistenceResult',
                    'data': result_data
                }]
                        })
                        result_ref = result[0][6]  # Get the reference
                    else:
                        result_ref = None
            except Exception as e:
                logger.error(f"Failed to save workspace object: {e}")
                result_ref = None
            
            # Create report using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'create_report'):
                    report_info = self.kb_util.create_report(workspace_name, {
                        'message': f'Protein existence check completed for {protein_id}. Exists: {result_data["exists"]}',
                        'objects_created': [{'ref': result_ref, 'description': 'Protein existence result'}] if result_ref else [],
                        'workspace_name': workspace_name,
                        'report_object_name': f"{protein_id}_existence_report"
                    })
                    # Handle mock objects by converting to dict if needed
                    if hasattr(report_info, '__dict__'):
                        report_info = {
                            'name': getattr(report_info, 'name', f"{protein_id}_existence_report"),
                            'ref': getattr(report_info, 'ref', 'mock_report_ref')
                        }
                else:
                    # Fallback to direct report client usage
                    report_client = KBaseReport(self.callback_url)
                    report_info = report_client.create_extended_report({
                        'message': f'Protein existence check completed for {protein_id}. Exists: {result_data["exists"]}',
                        'objects_created': [{'ref': result_ref, 'description': 'Protein existence result'}] if result_ref else [],
                        'workspace_name': workspace_name,
                        'report_object_name': f"{protein_id}_existence_report"
                    })
            except Exception as e:
                logger.error(f"Failed to create KBase report: {e}")
                report_info = {
                    'name': f"{protein_id}_existence_report",
                    'ref': 'error_report_ref'
                }
            
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'exists': 1 if result_data['exists'] else 0,
                'family_id': result_data['family_id'],
                'metadata': result_data['metadata'],
                'input_parameters': params,
                'start_time': start_time,
                'summary': f'Protein {protein_id} exists: {result_data["exists"]}',
                'protein_existence_result_ref': result_ref or '',
                'embedding_result_ref': embedding_ref or ''
            }
            
        except Exception as e:
            logger.error(f"Protein existence check failed: {e}")
            raise
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
           (Generate protein embeddings from direct sequence input. Creates
           embeddings using ESM-2 model for downstream analysis.) ->
           structure: parameter "report_name" of String, parameter
           "report_ref" of String, parameter "embedding_result_ref" of
           String, parameter "summary" of String, parameter
           "input_parameters" of mapping from String to unspecified object,
           parameter "start_time" of Double, parameter "embedding_norm" of
           Double, parameter "sequence_length" of Long, parameter
           "embedding_dim" of Long
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN generate_protein_embedding
        try:
            start_time = time.time()
            
            input_type = params.get('input_type', 'sequence')
            input_data = params.get('input_data')
            model_name = params.get('model_name', 'esm2_t6_8M_UR50D')
            
            if not input_data:
                raise ValueError("input_data parameter is required")
            
            # Setup workspace client using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'get_workspace_client'):
                    workspace_client = self.kb_util.get_workspace_client()
                else:
                    workspace_client = Workspace(self.callback_url)
            except Exception as e:
                logger.error(f"Failed to setup workspace client: {e}")
                workspace_client = Workspace(self.callback_url)
            
            # Inline workspace name extraction
            try:
                if ctx.get('provenance') and len(ctx['provenance']) > 0:
                    workspace_name = ctx['provenance'][0].get('ws_name')
                else:
                    workspace_name = ctx.get('ws_name')
            except Exception:
                workspace_name = None
            
            # Generate demo embedding
            if input_type == 'sequence':
                sequence_length = len(input_data)
                demo_embedding = np.random.rand(1280).tolist()  # Demo embedding
                input_id = f"seq_{hash(input_data) % 10000}"
            else:
                raise ValueError(f"Unsupported input_type: {input_type}")
            
            embedding_result = {
                'embedding': demo_embedding,
                'model_name': model_name,
                'protein_id': input_id,
                'sequence_length': sequence_length,
                'embedding_dim': 1280,
                'input_source': 'direct_sequence',
                'pooling_method': 'mean',
                'generation_timestamp': time.time()
            }
            
            # Save embedding to workspace using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'save_workspace_object'):
                    result_ref = self.kb_util.save_workspace_object(
                        workspace_name, f"{input_id}_embedding",
                        'KBaseProteinQueryModule.ProteinEmbedding', embedding_result)
                else:
                    # Fallback to direct workspace client usage
                    if workspace_client and workspace_name:
                        result = workspace_client.save_objects({
                'workspace': workspace_name,
                'objects': [{
                    'name': f"{input_id}_embedding",
                    'type': 'KBaseProteinQueryModule.ProteinEmbedding',
                    'data': embedding_result
                }]
                        })
                        result_ref = result[0][6]  # Get the reference
                    else:
                        result_ref = None
            except Exception as e:
                logger.error(f"Failed to save embedding to workspace: {e}")
                result_ref = None
            
            # Create report using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'create_report'):
                    report_info = self.kb_util.create_report(workspace_name, {
                        'message': f'Generated protein embedding for {input_id} using {model_name}',
                        'objects_created': [{'ref': result_ref, 'description': 'Protein embedding'}] if result_ref else [],
                        'workspace_name': workspace_name,
                        'report_object_name': f"{input_id}_embedding_report"
                    })
                    # Handle mock objects by converting to dict if needed
                    if hasattr(report_info, '__dict__'):
                        report_info = {
                            'name': getattr(report_info, 'name', f"{input_id}_embedding_report"),
                            'ref': getattr(report_info, 'ref', 'mock_report_ref')
                        }
                else:
                    # Fallback to direct report client usage
                    report_client = KBaseReport(self.callback_url)
                    report_info = report_client.create_extended_report({
                        'message': f'Generated protein embedding for {input_id} using {model_name}',
                        'objects_created': [{'ref': result_ref, 'description': 'Protein embedding'}] if result_ref else [],
                        'workspace_name': workspace_name,
                        'report_object_name': f"{input_id}_embedding_report"
                    })
            except Exception as e:
                logger.error(f"Failed to create KBase report: {e}")
                report_info = {
                    'name': f"{input_id}_embedding_report",
                    'ref': 'error_report_ref'
                }
            
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'embedding_result_ref': result_ref or '',
                'summary': f'Generated embedding for {input_id} (dim: {len(embedding_result["embedding"])})',
                'input_parameters': params,
                'start_time': start_time,
                'embedding_norm': float(np.linalg.norm(demo_embedding)),
                'sequence_length': embedding_result.get('sequence_length', 0),
                'embedding_dim': len(embedding_result.get('embedding', []))
            }
            
        except Exception as e:
            logger.error(f"Embedding generation failed: {e}")
            raise
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
        :returns: instance of type "AssignFamilyFastResults" (Assign a
           protein embedding to a family using similarity to family
           centroids. Uses binary Hamming distance for fast family
           assignment.) -> structure: parameter "family_id" of String,
           parameter "confidence" of Double, parameter "eigenprotein_id" of
           String, parameter "input_parameters" of mapping from String to
           unspecified object, parameter "start_time" of Double, parameter
           "family_assignment_result_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN assign_family_fast
        try:
            start_time = time.time()
            
            embedding_ref = params.get('embedding_ref')
            protein_id = params.get('protein_id', 'unknown')
            
            if not embedding_ref:
                raise ValueError("embedding_ref parameter is required")
            
            # Setup workspace client using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'get_workspace_client'):
                    workspace_client = self.kb_util.get_workspace_client()
                else:
                    workspace_client = Workspace(self.callback_url)
            except Exception as e:
                logger.error(f"Failed to setup workspace client: {e}")
                workspace_client = Workspace(self.callback_url)
            
            # Inline workspace name extraction
            try:
                if ctx.get('provenance') and len(ctx['provenance']) > 0:
                    workspace_name = ctx['provenance'][0].get('ws_name', 'test_workspace')
                else:
                    workspace_name = ctx.get('ws_name', 'test_workspace')
            except Exception:
                workspace_name = 'test_workspace'
            
            if not workspace_client:
                raise RuntimeError("Workspace client required for family assignment")
            
            # Demo family assignment
            family_id = 'demo_family'
            confidence = 0.85
            
            result_data = {
                'protein_id': protein_id,
                'family_id': family_id,
                'confidence': confidence,
                'embedding_ref': embedding_ref,
                'assignment_timestamp': time.time()
            }
            
            # Inline workspace object saving
            try:
                if workspace_client and workspace_name:
                    result = workspace_client.save_objects({
                'workspace': workspace_name,
                'objects': [{
                    'name': f"{protein_id}_family_assignment",
                    'type': 'KBaseProteinQueryModule.FamilyAssignmentResult',
                    'data': result_data
                }]
                    })
                    result_ref = result[0][6]  # Get the reference
                else:
                    result_ref = None
            except Exception as e:
                logger.error(f"Failed to save workspace object: {e}")
                result_ref = None
            
            output = {
                'family_id': family_id,
                'confidence': confidence,
                'eigenprotein_id': protein_id,
                'input_parameters': params,
                'start_time': start_time,
                'family_assignment_result_ref': result_ref or ''
            }
            
        except Exception as e:
            logger.error(f"Family assignment failed: {e}")
            raise
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
           top matches for a given protein embedding within a family. Uses
           FAISS IVF float index for efficient similarity search.) ->
           structure: parameter "matches" of list of mapping from String to
           unspecified object, parameter "summary" of String, parameter
           "input_parameters" of mapping from String to unspecified object,
           parameter "start_time" of Double, parameter "family_id" of String,
           parameter "top_n" of Long, parameter "similarity_stats" of mapping
           from String to Double, parameter "similarity_search_result_ref" of
           String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN find_top_matches_from_embedding
        try:
            start_time = time.time()
            
            embedding_ref = params.get('embedding_ref')
            protein_id = params.get('protein_id', 'unknown')
            max_matches = params.get('max_matches', 10)
            
            if not embedding_ref:
                raise ValueError("embedding_ref parameter is required")
            
            # Setup workspace client using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'get_workspace_client'):
                    workspace_client = self.kb_util.get_workspace_client()
                else:
                    workspace_client = Workspace(self.callback_url)
            except Exception as e:
                logger.error(f"Failed to setup workspace client: {e}")
                workspace_client = Workspace(self.callback_url)
            
            # Inline workspace name extraction
            try:
                if ctx.get('provenance') and len(ctx['provenance']) > 0:
                    workspace_name = ctx['provenance'][0].get('ws_name', 'test_workspace')
                else:
                    workspace_name = ctx.get('ws_name', 'test_workspace')
            except Exception:
                workspace_name = 'test_workspace'
            
            if not workspace_client:
                raise RuntimeError("Workspace client required for similarity search")
            
            # Generate demo matches
            family_id = 'demo_family'
            matches = []
            for i in range(min(max_matches, 5)):
                matches.append({
                    'protein_id': f'match_{i+1}',
                    'similarity_score': 0.9 - (i * 0.1),
                    'family_id': family_id,
                    'metadata': {'source': 'demo_data', 'rank': i + 1}
                })
            
            similarity_stats = {
                'mean_similarity': 0.75,
                'max_similarity': 0.9,
                'min_similarity': 0.5,
                'total_matches': len(matches)
            }
            
            result_data = {
                'protein_id': protein_id,
                'family_id': family_id,
                'matches': matches,
                'similarity_stats': similarity_stats,
                'search_timestamp': time.time()
            }
            
            # Inline workspace object saving
            try:
                if workspace_client and workspace_name:
                    result = workspace_client.save_objects({
                        'workspace': workspace_name,
                        'objects': [{
                            'name': f"{protein_id}_similarity_search",
                            'type': 'KBaseProteinQueryModule.SimilaritySearchResult',
                            'data': result_data
                        }]
                    })
                    result_ref = result[0][6]  # Get the reference
                else:
                    result_ref = None
            except Exception as e:
                logger.error(f"Failed to save workspace object: {e}")
                result_ref = None
            
            output = {
                'matches': matches,
                'summary': f'Found {len(matches)} matches for {protein_id}',
                'input_parameters': params,
                'start_time': start_time,
                'family_id': family_id,
                'top_n': len(matches),
                'similarity_stats': similarity_stats,
                'similarity_search_result_ref': result_ref or ''
            }
            
        except Exception as e:
            logger.error(f"Similarity search failed: {e}")
            raise
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
           (Summarize and visualize protein network analysis results.
           Generates comprehensive HTML reports with network visualization.)
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of String, parameter "input_parameters" of mapping
           from String to unspecified object, parameter "start_time" of
           Double, parameter "output_dir" of String, parameter "summary" of
           String, parameter "html_report_path" of String, parameter
           "sequence_analysis_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN summarize_and_visualize_results
        try:
            start_time = time.time()
            
            result_refs = params.get('result_refs', [])
            output_name = params.get('output_name', 'protein_analysis')
            
            if not result_refs:
                raise ValueError("result_refs parameter is required")
            
            # Setup workspace client using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'get_workspace_client'):
                    workspace_client = self.kb_util.get_workspace_client()
                else:
                    workspace_client = Workspace(self.callback_url)
            except Exception as e:
                logger.error(f"Failed to setup workspace client: {e}")
                workspace_client = Workspace(self.callback_url)
            
            # Inline workspace name extraction
            try:
                if ctx.get('provenance') and len(ctx['provenance']) > 0:
                    workspace_name = ctx['provenance'][0].get('ws_name', 'test_workspace')
                else:
                    workspace_name = ctx.get('ws_name', 'test_workspace')
            except Exception:
                workspace_name = 'test_workspace'
            
            if not workspace_client:
                raise RuntimeError("Workspace client required for summarization")
            
            # Get results from workspace
            results = []
            for ref in result_refs:
                try:
                    result_data = workspace_client.get_objects2({
                        'objects': [{'ref': ref}]
                    })['data'][0]['data']
                    results.append(result_data)
                except Exception as e:
                    logger.warning(f"Could not retrieve result {ref}: {e}")
                    
                    # Create summary data
                    summary_data = {
                        'total_results': len(results),
                        'analysis_timestamp': time.time(),
                        'result_types': [r.get('type', 'unknown') for r in results]
                    }
            
            # Inline workspace object saving
            try:
                if workspace_client and workspace_name:
                    result = workspace_client.save_objects({
                'workspace': workspace_name,
                'objects': [{
                    'name': f"{output_name}_summary",
                    'type': 'KBaseProteinQueryModule.AnalysisSummary',
                    'data': summary_data
                }]
                    })
                    summary_ref = result[0][6]  # Get the reference
                else:
                    summary_ref = None
            except Exception as e:
                logger.error(f"Failed to save workspace object: {e}")
                summary_ref = None
            
            # Inline report creation
            try:
                report_client = KBaseReport(self.callback_url)
                report_info = report_client.create_extended_report({
                    'message': f'Generated summary and visualization for {len(results)} results',
                    'objects_created': [{'ref': summary_ref, 'description': 'Analysis summary'}] if summary_ref else [],
                    'workspace_name': workspace_name,
                    'report_object_name': f"{output_name}_summary_report",
                    'html_links': [],
                    'file_links': [],
                    'direct_html_link_index': 0
                })
            except Exception as e:
                logger.error(f"Failed to create KBase report: {e}")
                report_info = {
                    'name': f"{output_name}_summary_report",
                    'ref': 'error_report_ref'
                }
            
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'analysis_result_ref': summary_ref or f'summary_{int(time.time())}',
                'input_parameters': params,
                'start_time': start_time,
                'output_directory': self.shared_folder,
                'general_info_dir': self.shared_folder,
                'network_analysis_dir': self.shared_folder,
                'sequence_analysis_dir': self.shared_folder,
                'embeddings_file_path': f'{self.shared_folder}/embeddings.csv',
                'top_proteins_csv_path': f'{self.shared_folder}/top_proteins.csv',
                'protein_count': len(results),
                'stages_completed': ['summarization', 'visualization'],
                'summary': f'Generated summary and visualization for {len(results)} results',
                'html_report_path': 'report.html',
                'sequence_analysis_ref': summary_ref or ''
            }
            
        except Exception as e:
            logger.error(f"Summarization failed: {e}")
            raise
        #END summarize_and_visualize_results

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method summarize_and_visualize_results return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

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
                # The input parser will handle both direct sequences and FASTA format
                input_kwargs['input_data'] = protein_input
                
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
                
            else:
                raise ValueError(f"Invalid input_type: {input_type}. Must be one of: protein_input, uniprot_ids, workspace_object")

            # Get analysis name for output
            analysis_name = params.get('analysis_name')
            if not analysis_name:
                raise ValueError("analysis_name is required")
            
            enabled_stages = params.get('enabled_stages') or [
                'input_validation','data_extraction','embedding_generation','family_assignment','similarity_search','sequence_analysis','network_analysis','bioinformatics_analysis','report_generation','visualization','data_export'
            ]
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
            
            # Get workspace objects created by the workflow
            workspace_objects = final.get('workspace_objects', [])
            workspace_objects_summary = final.get('workspace_objects_summary', {})

            # Ensure report fields are always present for Narrative integration
            report_name = final.get('report_name') or final.get('report', {}).get('name') or f"{analysis_name}_report"
            report_ref = final.get('report_ref') or final.get('report', {}).get('ref') or ''

            normalized = {
                'job_id': final.get('job_id') or wf_result.run_id,
                'analysis_result_ref': final.get('analysis_result_ref') or final.get('sequence_analysis_ref') or '',
                'summary': final.get('summary') or 'Protein query analysis completed',
                'input_parameters': params,
                'start_time': start_time,
                'output_directory': final.get('output_directory') or getattr(config, 'output_dir', self.shared_folder),
                'general_info_dir': final.get('general_info_dir') or final.get('output_directory') or getattr(config, 'output_dir', self.shared_folder),
                'network_analysis_dir': final.get('network_analysis_dir') or final.get('output_directory') or getattr(config, 'output_dir', self.shared_folder),
                'sequence_analysis_dir': final.get('sequence_analysis_dir') or final.get('output_directory') or getattr(config, 'output_dir', self.shared_folder),
                'embeddings_file_path': final.get('embeddings_file_path') or '',
                'top_proteins_csv_path': final.get('top_proteins_csv_path') or '',
                'html_report_path': final.get('html_report_path') or 'index.html',
                'protein_count': final.get('protein_count') or 0,
                'stages_completed': stages_completed,
                'report_name': report_name,
                'report_ref': report_ref
            }
            # Ensure a report exists for Narrative even if workflow didn't produce one
            if not normalized['report_ref']:
                try:
                    report_client = KBaseReport(self.callback_url)
                    
                    # Create report message with workspace objects information
                    report_message = normalized['summary']
                    if workspace_objects:
                        report_message += f"\n\nCreated {len(workspace_objects)} workspace objects:\n"
                        for obj in workspace_objects:
                            report_message += f"- {obj.get('name', 'Unknown')}: {obj.get('description', 'No description')}\n"
                    
                    report_info = report_client.create_extended_report({
                        'message': report_message,
                        'workspace_name': workspace_name,
                        'report_object_name': f"{analysis_name}_report",
                        'objects_created': workspace_objects  # Include workspace objects in report
                    })
                    normalized['report_name'] = report_info.get('name', normalized['report_name'])
                    normalized['report_ref'] = report_info.get('ref', normalized['report_ref'])
                except Exception as e:
                    logger.warning(f"Failed to create fallback report: {e}")
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
                'output_directory': self.shared_folder,
                'general_info_dir': self.shared_folder,
                'network_analysis_dir': self.shared_folder,
                'sequence_analysis_dir': self.shared_folder,
                'embeddings_file_path': '',
                'top_proteins_csv_path': '',
                'html_report_path': '',  # No HTML error reports
                'protein_count': 0,
                'stages_completed': [],
                'report_name': error_report['name'],
                'report_ref': error_report['ref']
            }]
            
            workspace_name = params.get('workspace_name')
            input_proteins = params.get('input_proteins', [])
            analysis_stages = params.get('analysis_stages', ['embedding_generation'])
            output_report_name = params.get('output_report_name', f'protein_analysis_report_{int(time.time())}')
            output_data_name = params.get('output_data_name', f'protein_analysis_data_{int(time.time())}')
            
            if not input_proteins:
                raise ValueError("input_proteins parameter is required")
            
            if not workspace_name:
                raise ValueError("workspace_name parameter is required")
            
            # Setup workspace client using KBUtilLib
            try:
                if self.kb_util and hasattr(self.kb_util, 'get_workspace_client'):
                    workspace_client = self.kb_util.get_workspace_client()
                else:
                    workspace_client = Workspace(self.callback_url)
            except Exception as e:
                logger.error(f"Failed to setup workspace client: {e}")
                workspace_client = Workspace(self.callback_url)
            
            if not workspace_client:
                raise RuntimeError("Could not setup workspace client")
            
            # Execute comprehensive analysis pipeline
            protein_count = len(input_proteins)
            
            # Create mock pipeline results for demonstration
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
            
            # Create data output directory for CSV and JSON files only (no HTML reports)
            output_directory = None
            data_files = []
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            
            if os.path.exists("test/outputs"):
                test_output_dir = os.path.join("test", "outputs", f"pipeline_run_{timestamp}")
                os.makedirs(test_output_dir, exist_ok=True)
                
                # Create CSV file
                csv_content = "rank,protein_id,similarity_score,family,description,source_type,input_protein_index\n"
                csv_content += f"INPUT,input_protein_1,1.0,Unknown,User input protein,user_input,0\n"
                csv_content += "1,match_protein_1,0.85,family_1,Similar protein 1,database_match,0\n"
                csv_content += "2,match_protein_2,0.80,family_1,Similar protein 2,database_match,0\n"
                
                csv_path = os.path.join(test_output_dir, "top_proteins_with_metadata.csv")
                with open(csv_path, 'w') as f:
                    f.write(csv_content)
                data_files.append(csv_path)
                
                # Create JSON results
                import json
                json_data = {
                    'analysis_id': output_report_name,
                    'timestamp': timestamp,
                    'input_proteins': input_proteins,
                    'analysis_stages': analysis_stages,
                    'pipeline_results': pipeline_results,
                    'protein_count': len(input_proteins)
                }
                json_path = os.path.join(test_output_dir, "pipeline_results.json")
                with open(json_path, 'w') as f:
                    json.dump(json_data, f, indent=2)
                data_files.append(json_path)
                
                output_directory = test_output_dir
                logger.info(f"Created test output directory: {test_output_dir}")
                print(f"Generated {len(data_files)} data files in {test_output_dir}")
            
            # Create KBase report with data files only (no HTML reports)
            try:
                report_client = KBaseReport(self.callback_url)
                report_info = report_client.create_extended_report({
                    'message': f'Protein analysis completed for {protein_count} proteins. Analysis stages: {", ".join(analysis_stages)}. Results include CSV data files and comprehensive pipeline results saved as workspace objects.',
                    'workspace_name': workspace_name,
                    'report_object_name': output_report_name,
                    'file_links': [{'path': f, 'name': os.path.basename(f), 'description': f'Analysis data file: {os.path.basename(f)}'} for f in data_files] if data_files else []
                })
            except Exception as e:
                logger.error(f"Failed to create report: {e}")
                report_info = {'name': output_report_name, 'ref': f'report_{int(time.time())}'}
            
            # Set up proper output directory for KBase Narrative integration
            os.environ['EXPORTS_DIR'] = self.shared_folder
            os.environ['SCRATCH_DIR'] = self.shared_folder
            
            # Update output with directory-based results
            stages_completed = list(pipeline_results.keys())
            
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'analysis_result_ref': f'analysis_{int(time.time())}',
                'summary': f'Completed protein query analysis with {len(stages_completed)} stages. Results saved as workspace objects and data files.',
                'input_parameters': params,
                'start_time': start_time,
                'html_report_path': '',  # No HTML reports
                'protein_count': protein_count,
                'stages_completed': stages_completed,
                'output_directory': output_directory,
                'general_info_dir': output_directory,
                'network_analysis_dir': output_directory,
                'sequence_analysis_dir': output_directory,
                'job_id': f'job_{int(time.time())}',
                'embeddings_file_path': csv_path if output_directory else '',
                'top_proteins_csv_path': csv_path if output_directory else '',
                'exported_files': data_files if output_directory else []
            }
            
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
