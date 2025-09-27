"""
Workflow Orchestrator for KBase Protein Query Module

This module provides a comprehensive, scalable pipeline that integrates all analysis stages
with proper dependency management and KBase integration using the new modular structure.
"""

import os
import logging
import time
import uuid
import gc
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field

# Import managers from new structure
from ..input import DataExtractionStage, InputValidationStage, WorkspaceObjectStage
from ..analysis import AnalysisManager, get_enabled_analyses
from ..outputs import OutputManager
from ..util import (
    ProteinEmbeddingGenerator,
    HierarchicalIndex, StreamingIndex,
    ProteinStorage, MemoryEfficientLoader,
    ProteinFamilyAssigner, ProteinExistenceChecker,
    IndexingStrategy
)
from .pipeline_config import PipelineConfig

logger = logging.getLogger(__name__)

@dataclass
class WorkflowResult:
    """Result container for the complete workflow execution."""
    
    success: bool
    run_id: str
    analyses_completed: List[str]
    analysis_results: Dict[str, Any]
    final_output: Dict[str, Any]
    execution_time: float
    output_directory: str
    stages_completed: List[str] = field(default_factory=list)
    error_message: Optional[str] = None

class WorkflowOrchestrator:
    """
    Orchestrates the complete protein query workflow using the new modular structure.
    
    This class coordinates input handling, analysis execution, and output generation
    through dedicated managers for each component.
    """
    
    def __init__(self, config: Optional[Union[Dict[str, Any], PipelineConfig]] = None, kb_util=None):
        """
        Initialize the Workflow Orchestrator.
        
        Args:
            config: Configuration dictionary or PipelineConfig object for the workflow
            kb_util: KBase utility library instance
        """
        # Handle both dict and PipelineConfig objects
        if isinstance(config, PipelineConfig):
            # Convert PipelineConfig to dictionary
            self.config = self._pipeline_config_to_dict(config)
            self.pipeline_config = config
        else:
            self.config = config or {}
            self.pipeline_config = None
            
        self.kb_util = kb_util
        self.run_id = str(uuid.uuid4())[:8]
        
        # Initialize managers
        self.output_manager = None
        self.analysis_manager = None
        
        # Initialize utility components
        self.embedding_generator = None
        self.similarity_index = None
        self.protein_storage = None
        
        # Results tracking
        self.results = {}
        self.execution_start_time = None
        
        logger.info(f"WorkflowOrchestrator initialized with run_id: {self.run_id}")
    
    def _pipeline_config_to_dict(self, config: PipelineConfig) -> Dict[str, Any]:
        """Convert PipelineConfig object to dictionary."""
        import dataclasses
        
        if hasattr(config, '__dict__'):
            return {k: v for k, v in config.__dict__.items() if not k.startswith('_')}
        elif dataclasses.is_dataclass(config):
            return dataclasses.asdict(config)
        else:
            # Fallback: convert to dict manually
            return {
                'input_proteins': getattr(config, 'input_proteins', []),
                'enabled_stages': getattr(config, 'enabled_stages', []),
                'stage_configs': getattr(config, 'stage_configs', {}),
                'storage_config': getattr(config, 'storage_config', {}),
                'similarity_config': getattr(config, 'similarity_config', {}),
                'output_dir': getattr(config, 'output_dir', '/tmp'),
                'workspace_name': getattr(config, 'workspace_name', None),
                'workspace_client': getattr(config, 'workspace_client', None),
                'workspace_url': getattr(config, 'workspace_url', None),
                'auth_token': getattr(config, 'auth_token', None),
                'generate_html_report': getattr(config, 'generate_html_report', True),
                'generate_network_visualization': getattr(config, 'generate_network_visualization', True)
            }
    
    def initialize_components(self, output_dir: str, workspace_name: Optional[str] = None):
        """
        Initialize all workflow components.
        
        Args:
            output_dir: Base directory for outputs
            workspace_name: KBase workspace name if applicable
        """
        try:
            # Use workspace_name from config if not provided
            if not workspace_name:
                workspace_name = self.config.get('workspace_name')
            
            # Initialize output manager with KBUtilLib
            self.output_manager = OutputManager(
                base_output_dir=output_dir,
                run_id=self.run_id,
                workspace_name=workspace_name,
                kb_util=self.kb_util
            )
            
            # Initialize analysis manager with output manager
            self.analysis_manager = AnalysisManager(output_manager=self.output_manager)
            
            # Initialize utility components
            self._initialize_utilities()
            
            logger.info("All workflow components initialized successfully")
            
        except Exception as e:
            logger.error(f"Error initializing workflow components: {e}")
            raise
    
    def _initialize_utilities(self):
        """Initialize utility components."""
        try:
            # Initialize embedding generator
            self.embedding_generator = ProteinEmbeddingGenerator()
            
            # Initialize protein storage
            self.protein_storage = ProteinStorage()
            
            logger.info("Utility components initialized")
            
        except Exception as e:
            logger.error(f"Error initializing utilities: {e}")
            raise
    
    def run_workflow(self, input_data: Dict[str, Any], 
                    output_dir: str,
                    workspace_name: Optional[str] = None,
                    selected_analyses: Optional[List[str]] = None,
                    **kwargs) -> WorkflowResult:
        """
        Run the complete protein query workflow.
        
        Args:
            input_data: Input data for the workflow
            output_dir: Directory to save outputs
            workspace_name: KBase workspace name if applicable
            selected_analyses: List of analyses to run (None for all enabled)
            **kwargs: Additional parameters
            
        Returns:
            WorkflowResult containing execution results
        """
        self.execution_start_time = time.time()
        
        try:
            logger.info(f"Starting workflow execution with run_id: {self.run_id}")
            
            # Initialize components
            self.initialize_components(output_dir, workspace_name)
            
            # Process input data
            processed_data = self._process_input_data(input_data)
            
            # Determine which analyses to run
            analyses_to_run = self._determine_analyses_to_run(selected_analyses)
            
            # Run analyses
            analysis_results = self._run_analyses(processed_data, analyses_to_run, **kwargs)
            
            # Generate final outputs
            final_output = self._generate_final_outputs(analysis_results)
            
            # Calculate execution time
            execution_time = time.time() - self.execution_start_time
            
            # Create workflow result
            result = WorkflowResult(
                success=True,
                run_id=self.run_id,
                analyses_completed=analyses_to_run,
                analysis_results=analysis_results,
                final_output=final_output,
                execution_time=execution_time,
                output_directory=self.output_manager.get_root_dir(),
                stages_completed=analyses_to_run  # Use analyses as stages for now
            )
            
            # Save final metadata
            self._save_final_metadata(result)
            
            logger.info(f"Workflow completed successfully in {execution_time:.2f} seconds")
            return result
            
        except Exception as e:
            execution_time = time.time() - self.execution_start_time if self.execution_start_time else 0
            logger.error(f"Workflow execution failed: {e}")
            
            return WorkflowResult(
                success=False,
                run_id=self.run_id,
                analyses_completed=[],
                analysis_results={},
                final_output={},
                execution_time=execution_time,
                output_directory=output_dir,
                stages_completed=[],
                error_message=str(e)
            )
    
    def _process_input_data(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process and validate input data.
        
        Args:
            input_data: Raw input data
            
        Returns:
            Processed and validated data
        """
        try:
            logger.info("Processing input data")
            
            # Initialize input processing stages
            validation_stage = InputValidationStage()
            extraction_stage = DataExtractionStage()
            workspace_stage = WorkspaceObjectStage()
            
            # Validate input
            validation_result = validation_stage.run(input_data)
            if not validation_result.success:
                raise ValueError(f"Input validation failed: {validation_result.error_message}")
            
            # Extract data
            extraction_result = extraction_stage.run(input_data)
            if not extraction_result.success:
                raise ValueError(f"Data extraction failed: {extraction_result.error_message}")
            
            # Process workspace objects
            workspace_result = workspace_stage.run(input_data)
            if not workspace_result.success:
                raise ValueError(f"Workspace processing failed: {workspace_result.error_message}")
            
            # Combine results
            processed_data = {
                "validation": validation_result.data,
                "extraction": extraction_result.data,
                "workspace": workspace_result.data,
                "original_input": input_data
            }
            
            logger.info("Input data processing completed")
            return processed_data
            
        except Exception as e:
            logger.error(f"Error processing input data: {e}")
            raise
    
    def _determine_analyses_to_run(self, selected_analyses: Optional[List[str]] = None) -> List[str]:
        """
        Determine which analyses to run based on selection and availability.
        
        Args:
            selected_analyses: User-selected analyses (None for all enabled)
            
        Returns:
            List of analysis names to run
        """
        enabled_analyses = get_enabled_analyses()
        
        if selected_analyses is None:
            # Run all enabled analyses
            return list(enabled_analyses.keys())
        
        # Filter selected analyses to only include enabled ones
        available_analyses = []
        for analysis in selected_analyses:
            if analysis in enabled_analyses:
                available_analyses.append(analysis)
            else:
                logger.warning(f"Analysis '{analysis}' is not enabled or available")
        
        return available_analyses
    
    def _run_analyses(self, processed_data: Dict[str, Any], 
                     analyses_to_run: List[str], **kwargs) -> Dict[str, Any]:
        """
        Run the selected analyses.
        
        Args:
            processed_data: Processed input data
            analyses_to_run: List of analyses to run
            **kwargs: Additional parameters
            
        Returns:
            Dictionary of analysis results
        """
        try:
            logger.info(f"Running analyses: {analyses_to_run}")
            
            # Prepare data for analyses
            analysis_data = self._prepare_analysis_data(processed_data)
            
            # Run analyses through analysis manager
            results = self.analysis_manager.run_multiple_analyses(
                analysis_names=analyses_to_run,
                input_data=analysis_data,
                output_dir=self.output_manager.get_root_dir(),
                **kwargs
            )
            
            logger.info(f"Completed {len(results)} analyses")
            return results
            
        except Exception as e:
            logger.error(f"Error running analyses: {e}")
            raise
    
    def _prepare_analysis_data(self, processed_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Prepare data for analysis execution.
        
        Args:
            processed_data: Processed input data
            
        Returns:
            Data prepared for analysis
        """
        # Extract relevant data for analyses
        analysis_data = {
            "proteins": processed_data.get("extraction", {}).get("proteins", []),
            "workspace_info": processed_data.get("workspace", {}),
            "config": self.config
        }
        
        # Add embeddings if available
        if self.embedding_generator:
            try:
                proteins = analysis_data["proteins"]
                if proteins:
                    embeddings = self.embedding_generator.generate_embeddings(proteins)
                    analysis_data["embeddings"] = embeddings
            except Exception as e:
                logger.warning(f"Could not generate embeddings: {e}")
        
        return analysis_data
    
    def _generate_final_outputs(self, analysis_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate final consolidated outputs and create workspace objects.
        
        Args:
            analysis_results: Results from all analyses
            
        Returns:
            Consolidated final output
        """
        try:
            logger.info("Generating final outputs")
            
            final_output = {
                "run_id": self.run_id,
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "analyses_completed": list(analysis_results.keys()),
                "summary": self._generate_summary(analysis_results),
                "analysis_results": analysis_results
            }
            
            # Save final output
            self.output_manager.write_json(
                "",
                "final_output.json",
                final_output,
                description="Final consolidated output from all analyses"
            )
            
            # Create workspace objects for all outputs
            logger.info("Creating workspace objects from analysis results")
            workspace_objects = self.output_manager.create_workspace_objects(analysis_results)
            
            # Add workspace objects info to final output
            final_output["workspace_objects"] = workspace_objects
            final_output["workspace_objects_summary"] = self.output_manager.get_workspace_objects_summary()
            
            logger.info(f"Created {len(workspace_objects)} workspace objects")
            
            return final_output
            
        except Exception as e:
            logger.error(f"Error generating final outputs: {e}")
            raise
    
    def _generate_summary(self, analysis_results: Dict[str, Any]) -> str:
        """
        Generate a summary of the analysis results.
        
        Args:
            analysis_results: Results from all analyses
            
        Returns:
            Summary string
        """
        summary_parts = [
            f"Protein Query Analysis Run {self.run_id}",
            f"Completed {len(analysis_results)} analyses:",
        ]
        
        for analysis_name, result in analysis_results.items():
            if "error" in result:
                summary_parts.append(f"  - {analysis_name}: FAILED ({result['error']})")
            else:
                summary_parts.append(f"  - {analysis_name}: SUCCESS")
        
        return "\n".join(summary_parts)
    
    def _save_final_metadata(self, result: WorkflowResult):
        """
        Save final metadata about the workflow execution.
        
        Args:
            result: Final workflow result
        """
        try:
            # Save metadata
            self.output_manager.save_metadata(
                config=self.config,
                analyses_run=result.analyses_completed,
                summary=result.final_output.get("summary", "")
            )
            
            # Save process info
            self.output_manager.save_process_info(
                stages_completed=result.analyses_completed,
                execution_time=result.execution_time
            )
            
            # Create final manifest
            self.output_manager.finalize_manifest()
            
            logger.info("Final metadata saved")
            
        except Exception as e:
            logger.error(f"Error saving final metadata: {e}")
    
    def execute(self, input_data: Optional[Union[Dict[str, Any], PipelineConfig]] = None) -> WorkflowResult:
        """Execute the workflow with the given input data.
        
        This is the main entry point for the workflow orchestrator.
        If no input_data is provided, it will use the config data.
        """
        if input_data is None:
            input_data = self.config or {}
        
        # Handle both dict and PipelineConfig objects
        if isinstance(input_data, PipelineConfig):
            # It's a PipelineConfig object
            output_dir = getattr(input_data, 'output_dir', None) or '/tmp/protein_query_output'
            workspace_name = getattr(input_data, 'workspace_name', None)
            # Convert to dict for run_workflow
            input_dict = self._pipeline_config_to_dict(input_data)
        else:
            # It's a dictionary
            output_dir = input_data.get('output_dir', '/tmp/protein_query_output')
            workspace_name = input_data.get('workspace_name')
            input_dict = input_data
        
        return self.run_workflow(input_dict, output_dir, workspace_name)
    
    def get_available_analyses(self) -> Dict[str, Dict[str, Any]]:
        """
        Get list of available analyses.
        
        Returns:
            Dictionary of available analyses
        """
        if self.analysis_manager:
            return self.analysis_manager.get_available_analyses()
        return get_enabled_analyses()
    
    def cleanup(self):
        """Clean up resources and temporary files."""
        try:
            # Clear results
            self.results.clear()
            
            # Force garbage collection
            gc.collect()
            
            logger.info("Workflow cleanup completed")
            
        except Exception as e:
            logger.error(f"Error during cleanup: {e}")
