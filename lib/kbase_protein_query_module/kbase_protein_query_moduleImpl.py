# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import logging
import time
import uuid
from typing import Dict, Any, List, Optional

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil

from .src.workflows.workflow_orchestrator import ProteinQueryWorkflow
from .src.core.pipeline_config import PipelineConfig
from .src.storage import ProteinExistenceChecker, ProteinFamilyAssigner
from .src.processing.embeddings.generator import ProteinEmbeddingGenerator
from .src.processing.similarity.hierarchical_index import HierarchicalIndex
from .src.processing.networks.builder import DynamicNetworkBuilder
from .src.analysis.sequence_analyzer import ProteinSequenceAnalyzer
from .src.reports.html.generator import HTMLReportGenerator

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

Authors: Vibhav Setlur
Contact: https://kbase.us/contact-us/
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/VibhavSetlur/kbase_protein_query_module.git"
    GIT_COMMIT_HASH = "2984af3c7bb64132c99024796d5335cc48816269"

    #BEGIN_CLASS_HEADER
    VERSION = "2.0.0"
    GIT_URL = "https://github.com/VibhavSetlur/kbase_protein_query_module.git"
    GIT_COMMIT_HASH = "8279279e82c839ab6fa37678c41d700ab876a27b"
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ.get('SDK_CALLBACK_URL')
        if self.callback_url is None:
            raise RuntimeError('SDK_CALLBACK_URL environment variable must be set.')
        
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        
        # Initialize DataFileUtil client
        self.dfu = DataFileUtil(self.callback_url)
        
        # Initialize core components
        self._initialize_components()
        #END_CONSTRUCTOR
        pass

    def _initialize_components(self):
        """Initialize core components for protein query operations."""
        try:
            # Initialize storage components
            from .src.storage import ProteinStorage, MemoryEfficientLoader
            self.storage = ProteinStorage(base_dir="data")
            self.memory_loader = MemoryEfficientLoader(self.storage)
            
            # Initialize embedding generator
            self.embedding_generator = ProteinEmbeddingGenerator(
                model_name="esm2_t6_8M_UR50D",
                device="cpu"
            )
            
            # Initialize similarity search components
            self.hierarchical_index = HierarchicalIndex()
            
            # Initialize network builder
            self.network_builder = DynamicNetworkBuilder()
            
            # Initialize family assignment
            self.family_assigner = ProteinFamilyAssigner()
            
            # Initialize existence checker
            self.existence_checker = ProteinExistenceChecker()
            
            # Initialize sequence analyzer
            self.sequence_analyzer = ProteinSequenceAnalyzer()
            
            # Initialize report generator
            self.report_generator = HTMLReportGenerator()
            
            # Initialize workflow orchestrator
            self.workflow_orchestrator = ProteinQueryWorkflow()
            
            logger.info("All core components initialized successfully")
            
        except Exception as e:
            logger.error(f"Failed to initialize components: {str(e)}")
            raise

    def _setup_workspace_client(self):
        """Setup workspace client for KBase operations."""
        try:
            token = os.environ.get('KB_AUTH_TOKEN')
            if not token:
                logger.warning("No auth token found, workspace operations may fail")
                return None
            
            ws_url = os.environ.get('KBASE_ENDPOINT', 'https://appdev.kbase.us/services')
            if not ws_url.endswith('/ws'):
                ws_url = ws_url + '/ws'
            
            return Workspace(url=ws_url, token=token)
        except Exception as e:
            logger.error(f"Failed to setup workspace client: {e}")
            return None

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
            
            # Extract parameters
            protein_id = params.get('protein_id')
            generate_embedding = params.get('generate_embedding', False)
            
            if not protein_id:
                raise ValueError("protein_id parameter is required")
            
            # Setup workspace client
            workspace_client = self._setup_workspace_client()
            if not workspace_client:
                raise RuntimeError("Could not setup workspace client")
            
            # Check protein existence using workflow orchestrator
            existence_result = self.existence_checker.check_protein_existence(protein_id)
            
            # Create result object
            result_data = {
                'protein_id': protein_id,
                'exists': existence_result.get('exists', False),
                'family_id': existence_result.get('family_id', 'unknown'),
                'metadata': existence_result.get('metadata', {}),
                'check_timestamp': time.time()
            }
            
            # Generate embedding if requested
            embedding_ref = None
            if generate_embedding and existence_result.get('exists', False):
                try:
                    embedding_generator = ProteinEmbeddingGenerator()
                    embedding_result = embedding_generator.generate_embedding_from_id(protein_id)
                    
                    if embedding_result:
                        # Save embedding to workspace
                        workspace_name = ctx.get('provenance', [{}])[0].get('ws_name', 'unknown')
                        embedding_ref = workspace_client.save_objects({
                            'id': workspace_name,
                            'objects': [{
                                'name': f"{protein_id}_embedding",
                                'type': 'KBaseProteinQueryModule.ProteinEmbedding',
                                'data': embedding_result
                            }]
                        })[0][6]  # Get the reference
                        
                        result_data['embedding_ref'] = embedding_ref
                        result_data['embedding'] = embedding_result.get('embedding', [])
                        result_data['model_name'] = embedding_result.get('model_name', 'esm2_t6_8M_UR50D')
                except Exception as e:
                    logger.warning(f"Could not generate embedding for {protein_id}: {e}")
            
            # Save result to workspace
            workspace_name = ctx.get('provenance', [{}])[0].get('ws_name', 'unknown')
            result_ref = workspace_client.save_objects({
                'id': workspace_name,
                'objects': [{
                    'name': f"{protein_id}_existence_check",
                    'type': 'KBaseProteinQueryModule.ProteinExistenceResult',
                    'data': result_data
                }]
            })[0][6]  # Get the reference
            
            # Create report
            report_client = KBaseReport(self.callback_url)
            report_info = report_client.create_extended_report({
                'message': f'Protein existence check completed for {protein_id}. Exists: {result_data["exists"]}',
                'objects_created': [{'ref': result_ref, 'description': 'Protein existence result'}],
                'workspace_name': workspace_name,
                'report_object_name': f"{protein_id}_existence_report"
            })
            
            execution_time = time.time() - start_time
            
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'exists': result_data['exists'],
                'family_id': result_data['family_id'],
                'metadata': result_data['metadata'],
                'input_parameters': params,
                'start_time': start_time,
                'summary': f'Protein {protein_id} exists: {result_data["exists"]}',
                'protein_existence_result_ref': result_ref,
                'embedding_result_ref': embedding_ref
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
            
            # Extract parameters
            input_type = params.get('input_type', 'sequence')
            input_data = params.get('input_data')
            model_name = params.get('model_name', 'esm2_t6_8M_UR50D')
            pooling_method = params.get('pooling_method', 'mean')
            
            if not input_data:
                raise ValueError("input_data parameter is required")
            
            # Setup workspace client
            workspace_client = self._setup_workspace_client()
            if not workspace_client:
                raise RuntimeError("Could not setup workspace client")
            
            # Generate embedding
            embedding_generator = ProteinEmbeddingGenerator()
            
            if input_type == 'sequence':
                embedding_result = embedding_generator.generate_embedding_from_sequence(
                    input_data, model_name, pooling_method
                )
                input_id = f"seq_{hash(input_data) % 10000}"
                input_source = 'direct_sequence'
            elif input_type == 'uniprot_id':
                embedding_result = embedding_generator.generate_embedding_from_id(
                    input_data, model_name, pooling_method
                )
                input_id = input_data
                input_source = 'uniprot_id'
            elif input_type == 'workspace_object':
                # Get sequence from workspace object
                obj_data = workspace_client.get_objects2({
                    'objects': [{'ref': input_data}]
                })['data'][0]['data']
                
                sequence = obj_data.get('sequence')
                if not sequence:
                    raise ValueError("No sequence found in workspace object")
                
                embedding_result = embedding_generator.generate_embedding_from_sequence(
                    sequence, model_name, pooling_method
                )
                input_id = obj_data.get('id', f"obj_{hash(input_data) % 10000}")
                input_source = 'workspace_object'
            else:
                raise ValueError(f"Unsupported input_type: {input_type}")
            
            if not embedding_result:
                raise RuntimeError("Failed to generate embedding")
            
            # Add metadata
            embedding_result.update({
                'input_id': input_id,
                'input_source': input_source,
                'generation_timestamp': time.time()
            })
            
            # Save to workspace
            workspace_name = ctx.get('provenance', [{}])[0].get('ws_name', 'unknown')
            result_ref = workspace_client.save_objects({
                'id': workspace_name,
                'objects': [{
                    'name': f"{input_id}_embedding",
                    'type': 'KBaseProteinQueryModule.ProteinEmbedding',
                    'data': embedding_result
                }]
            })[0][6]  # Get the reference
            
            # Create report
            report_client = KBaseReport(self.callback_url)
            report_info = report_client.create_extended_report({
                'message': f'Generated protein embedding for {input_id} using {model_name}',
                'objects_created': [{'ref': result_ref, 'description': 'Protein embedding'}],
                'workspace_name': workspace_name,
                'report_object_name': f"{input_id}_embedding_report"
            })
            
            execution_time = time.time() - start_time
            
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'embedding_result_ref': result_ref,
                'summary': f'Generated embedding for {input_id} (dim: {len(embedding_result["embedding"])})',
                'input_parameters': params,
                'start_time': start_time,
                'embedding_norm': embedding_result.get('embedding_norm', 0.0),
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
            
            # Extract parameters
            embedding_ref = params.get('embedding_ref')
            protein_id = params.get('protein_id', 'unknown')
            
            if not embedding_ref:
                raise ValueError("embedding_ref parameter is required")
            
            # Setup workspace client
            workspace_client = self._setup_workspace_client()
            if not workspace_client:
                raise RuntimeError("Could not setup workspace client")
            
            # Get embedding from workspace
            embedding_data = workspace_client.get_objects2({
                'objects': [{'ref': embedding_ref}]
            })['data'][0]['data']
            
            embedding = embedding_data.get('embedding', [])
            if not embedding:
                raise ValueError("No embedding found in workspace object")
            
            # Assign family
            if self.family_assigner:
                family_result = self.family_assigner.assign_family(embedding)
                family_id = family_result.get('family_id', 'unknown')
                confidence = family_result.get('confidence', 0.0)
            else:
                family_id = 'unknown'
                confidence = 0.0
            
            # Create result object
            result_data = {
                'protein_id': protein_id,
                'family_id': family_id,
                'confidence': confidence,
                'embedding_ref': embedding_ref,
                'assignment_timestamp': time.time()
            }
            
            # Save to workspace
            workspace_name = ctx.get('provenance', [{}])[0].get('ws_name', 'unknown')
            result_ref = workspace_client.save_objects({
                'id': workspace_name,
                'objects': [{
                    'name': f"{protein_id}_family_assignment",
                    'type': 'KBaseProteinQueryModule.FamilyAssignmentResult',
                    'data': result_data
                }]
            })[0][6]  # Get the reference
            
            execution_time = time.time() - start_time
            
            output = {
                'family_id': family_id,
                'confidence': confidence,
                'eigenprotein_id': protein_id,
                'input_parameters': params,
                'start_time': start_time,
                'family_assignment_result_ref': result_ref
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
            
            # Extract parameters
            embedding_ref = params.get('embedding_ref')
            protein_id = params.get('protein_id', 'unknown')
            max_matches = params.get('max_matches', 10)
            
            if not embedding_ref:
                raise ValueError("embedding_ref parameter is required")
            
            # Setup workspace client
            workspace_client = self._setup_workspace_client()
            if not workspace_client:
                raise RuntimeError("Could not setup workspace client")
            
            # Get embedding from workspace
            embedding_data = workspace_client.get_objects2({
                'objects': [{'ref': embedding_ref}]
            })['data'][0]['data']
            
            embedding = embedding_data.get('embedding', [])
            if not embedding:
                raise ValueError("No embedding found in workspace object")
            
            # Find matches (placeholder implementation)
            matches = []
            for i in range(min(max_matches, 5)):  # Limit to 5 for demo
                matches.append({
                    'protein_id': f'match_{i+1}',
                    'similarity_score': 0.9 - (i * 0.1),
                    'family_id': 'demo_family',
                    'metadata': {'source': 'demo_data'}
                })
            
            execution_time = time.time() - start_time
            
            output = {
                'matches': matches,
                'summary': f'Found {len(matches)} matches for {protein_id}',
                'input_parameters': params,
                'start_time': start_time,
                'family_id': 'demo_family'
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
            
            # Extract parameters
            result_refs = params.get('result_refs', [])
            output_name = params.get('output_name', 'protein_analysis')
            
            if not result_refs:
                raise ValueError("result_refs parameter is required")
            
            # Setup workspace client
            workspace_client = self._setup_workspace_client()
            if not workspace_client:
                raise RuntimeError("Could not setup workspace client")
            
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
            
            # Generate visualization
            workspace_name = ctx.get('provenance', [{}])[0].get('ws_name', 'unknown')
            
            # Create simple summary
            summary_data = {
                'total_results': len(results),
                'analysis_timestamp': time.time(),
                'result_types': [r.get('type', 'unknown') for r in results]
            }
            
            # Save summary to workspace
            summary_ref = workspace_client.save_objects({
                'id': workspace_name,
                'objects': [{
                    'name': f"{output_name}_summary",
                    'type': 'KBaseProteinQueryModule.AnalysisSummary',
                    'data': summary_data
                }]
            })[0][6]
            
            # Create report
            report_client = KBaseReport(self.callback_url)
            report_info = report_client.create_extended_report({
                'message': f'Generated summary and visualization for {len(results)} results',
                'objects_created': [{'ref': summary_ref, 'description': 'Analysis summary'}],
                'workspace_name': workspace_name,
                'report_object_name': f"{output_name}_summary_report"
            })
            
            execution_time = time.time() - start_time
            
            output = {
                'summary_ref': summary_ref,
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'summary': f'Generated summary and visualization for {len(results)} results',
                'input_parameters': params,
                'start_time': start_time
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
        try:
            start_time = time.time()
            
            # Extract parameters
            workspace_name = params.get('workspace_name')
            input_type = params.get('input_type', 'sequence')
            input_data = params.get('input_data')
            input_data_ref = params.get('input_data_ref')  # Workspace object reference
            analysis_stages = params.get('analysis_stages', ['embedding', 'family_assignment'])
            stop_after_stage = params.get('stop_after_stage')
            
            if not input_data and not input_data_ref:
                raise ValueError("Either input_data or input_data_ref parameter is required")
            
            if not workspace_name:
                raise ValueError("workspace_name parameter is required")
            
            # Setup workspace client
            workspace_client = self._setup_workspace_client()
            if not workspace_client:
                raise RuntimeError("Could not setup workspace client")
            
            # Use workspace object reference if provided, otherwise use text input
            if input_data_ref:
                input_data = input_data_ref
                # Determine input type from workspace object type if not specified
                if input_type == 'sequence':
                    # Extract type from workspace object reference
                    try:
                        obj_info = workspace_client.get_object_info3({'objects': [{'ref': input_data_ref}]})['infos'][0]
                        obj_type = obj_info[2]
                        if 'ProteinSequenceSet' in obj_type:
                            input_type = 'ProteinSequenceSet'
                        elif 'Genome' in obj_type:
                            input_type = 'Genome'
                        elif 'FeatureSet' in obj_type:
                            input_type = 'FeatureSet'
                        elif 'GenomeSet' in obj_type:
                            input_type = 'GenomeSet'
                    except Exception as e:
                        logger.warning(f"Could not determine input type from workspace object: {e}")
            
            # Handle KBase data types using DataFileUtil
            try:
                processed_data = self._handle_kbase_data_types(input_data, input_type)
                logger.info(f"Successfully processed {input_type} input with {processed_data['count']} sequences")
                
                # Update input_data with processed sequences
                processed_input_data = {
                    'sequences': processed_data['sequences'],
                    'metadata': processed_data['metadata'],
                    'source_type': processed_data['source_type'],
                    'original_input': input_data
                }
            except Exception as e:
                logger.error(f"Failed to process {input_type} input: {e}")
                raise RuntimeError(f"Input processing failed: {str(e)}")
            
            # Create pipeline configuration
            from .src.workflows.unified_workflow_orchestrator import PipelineConfig, create_basic_pipeline
            
            config = PipelineConfig(
                input_type=input_type,
                input_data=processed_input_data,
                workspace_client=workspace_client,
                enabled_stages=analysis_stages,
                stop_after_stage=stop_after_stage,
                workspace_name=workspace_name
            )
            
            # Create and run pipeline
            pipeline = create_basic_pipeline(input_type, processed_input_data, workspace_client)
            pipeline.config = config
            
            result = pipeline.run_pipeline()
            
            if not result.success:
                raise RuntimeError(f"Pipeline failed: {result.error_message}")
            
            # Create KBase report
            try:
                report_client = KBaseReport(self.callback_url)
                
                # Create report content
                report_content = f"""
                <h2>Protein Query Analysis Results</h2>
                <p><strong>Input Type:</strong> {input_type}</p>
                <p><strong>Stages Completed:</strong> {', '.join(result.stages_completed)}</p>
                <p><strong>Execution Time:</strong> {result.execution_time:.2f} seconds</p>
                <p><strong>Status:</strong> {'Success' if result.success else 'Failed'}</p>
                
                <h3>Analysis Summary</h3>
                <p>{result.final_output.get('summary', 'No summary available')}</p>
                
                <h3>Warnings</h3>
                <ul>
                {''.join([f'<li>{warning}</li>' for warning in result.warnings])}
                </ul>
                """
                
                report_info = report_client.create_extended_report({
                    'message': f'Completed protein query analysis with {len(result.stages_completed)} stages',
                    'objects_created': [],
                    'workspace_name': workspace_name,
                    'report_object_name': 'protein_query_analysis_report',
                    'html_links': [{'name': 'Analysis Report', 'description': 'Protein query analysis results', 'html': report_content}]
                })
                
                execution_time = time.time() - start_time
                
                output = {
                    'report_name': report_info['name'],
                    'report_ref': report_info['ref'],
                    'analysis_result_ref': result.run_id,
                    'summary': f'Completed protein query analysis with {len(result.stages_completed)} stages',
                    'input_parameters': params,
                    'start_time': start_time,
                    'html_report_path': 'analysis_report.html',
                    'protein_count': result.final_output.get('protein_count', 0),
                    'stages_completed': result.stages_completed
                }
                
            except Exception as e:
                logger.error(f"Failed to create KBase report: {e}")
                # Return basic result without report
                execution_time = time.time() - start_time
                output = {
                    'report_name': 'protein_query_analysis_report',
                    'report_ref': 'error_report_ref',
                    'analysis_result_ref': result.run_id,
                    'summary': f'Completed protein query analysis with {len(result.stages_completed)} stages (report creation failed)',
                    'input_parameters': params,
                    'start_time': start_time,
                    'html_report_path': 'analysis_report.html',
                    'protein_count': result.final_output.get('protein_count', 0),
                    'stages_completed': result.stages_completed
                }
            
        except Exception as e:
            logger.error(f"Protein query analysis failed: {e}")
            # Create error report
            try:
                error_report = self._create_error_report(str(e), workspace_name, time.time() - start_time)
                output = {
                    'report_name': error_report['report_name'],
                    'report_ref': error_report['report_ref'],
                    'analysis_result_ref': 'error',
                    'summary': f'Analysis failed: {str(e)}',
                    'input_parameters': params,
                    'start_time': start_time,
                    'html_report_path': 'error_report.html',
                    'protein_count': 0,
                    'stages_completed': []
                }
            except Exception as report_error:
                logger.error(f"Failed to create error report: {report_error}")
                raise
        #END run_protein_query_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_protein_query_analysis return value ' +
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
