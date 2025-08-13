"""
Workflow Orchestrator for KBase Protein Query Module

This module provides a comprehensive, scalable pipeline that integrates all analysis stages
with proper dependency management and KBase integration. It combines the best features
from both the unified and legacy orchestrators.
"""

import os
import logging
import time
import uuid
import gc
import yaml
import numpy as np
import pandas as pd
import networkx as nx
from typing import Dict, Any, List, Optional, Union, Tuple
from dataclasses import dataclass, field
from pathlib import Path

from ..core import BaseStage, StageResult, PipelineConfig
from ..stages import (
    STAGE_REGISTRY, STAGE_DEPENDENCIES, get_stage_class, get_stage_dependencies,
    InputValidationStage, DataExtractionStage, WorkspaceObjectStage,
    EmbeddingGenerationStage, FamilyAssignmentStage, SimilaritySearchStage,
    SequenceAnalysisStage, NetworkAnalysisStage, BioinformaticsAnalysisStage,
    ReportGenerationStage, VisualizationStage, DataExportStage
)
from ..utils import input_parser
from ..processing.embeddings.generator import ProteinEmbeddingGenerator
from ..storage import ProteinFamilyAssigner
from ..processing.similarity.hierarchical_index import HierarchicalIndex, StreamingIndex
from ..processing.networks.builder import DynamicNetworkBuilder
from ..storage import ProteinStorage, MemoryEfficientLoader

logger = logging.getLogger(__name__)


@dataclass
class WorkflowResult:
    """Result container for the complete workflow execution."""
    
    success: bool
    run_id: str
    stages_completed: List[str]
    stage_results: Dict[str, Any]
    final_output: Dict[str, Any]
    execution_time: float
    error_message: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


class ProteinQueryWorkflow:
    """
    Comprehensive workflow orchestrator for protein query analysis.
    
    This orchestrator provides:
    - Modular stage-based architecture
    - Automatic dependency resolution
    - KBase workspace integration
    - Performance optimization for large datasets
    - Comprehensive error handling and logging
    - Support for both small and massive datasets
    """
    
    def __init__(self, config: PipelineConfig = None):
        """
        Initialize the workflow orchestrator.
        
        Args:
            config: Pipeline configuration
        """
        self.config = config or PipelineConfig()
        self.run_id = str(uuid.uuid4())
        self.logger = logging.getLogger(f"{__name__}.{self.run_id}")
        
        # Initialize components
        self._initialize_components()
        
        # Stage execution tracking
        self.stages_completed = []
        self.stage_results = {}
        self.performance_metrics = {}
        
        # Load pre-computed data if available
        self._load_precomputed_data()
    
    def _initialize_components(self):
        """Initialize all workflow components."""
        try:
            # Initialize storage components
            self.storage = ProteinStorage(**self.config.storage_config)
            self.memory_loader = MemoryEfficientLoader(self.storage)
            
            # Initialize embedding generator
            self.embedding_generator = ProteinEmbeddingGenerator(
                model_name=self.config.embedding_model,
                device=self.config.embedding_device
            )
            
            # Initialize similarity search components
            self.hierarchical_index = HierarchicalIndex()
            self.streaming_index = StreamingIndex()
            
            # Initialize network builder
            self.network_builder = DynamicNetworkBuilder()
            
            # Initialize family assignment
            self.family_assigner = ProteinFamilyAssigner()
            
            self.logger.info("All workflow components initialized successfully")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize components: {str(e)}")
            raise
    
    def _load_precomputed_data(self):
        """Load pre-computed data for performance optimization."""
        try:
            # Load family centroids
            centroids_path = self._find_centroids_file()
            if centroids_path:
                self.logger.info(f"Loading family centroids from: {centroids_path}")
                self.family_assigner.load_family_centroids(centroids_path)
            else:
                self.logger.warning("Family centroids file not found. Family assignment will not be available.")
            
            # Load pre-computed indexes
            self._load_precomputed_indexes()
            
        except Exception as e:
            self.logger.warning(f"Failed to load pre-computed data: {str(e)}")
    
    def _find_centroids_file(self) -> Optional[str]:
        """Find the family centroids file."""
        possible_paths = [
            str(self.storage.base_dir / "family_centroids" / "family_centroids_binary.npz") if self.storage else None,
            str(self.storage.base_dir / "family_centroids_binary.npz") if self.storage else None,
            "data/family_centroids/family_centroids_binary.npz",
            "data/family_centroids_binary.npz",
            "/kb/module/data/family_centroids/family_centroids_binary.npz",
            "/kb/module/data/family_centroids_binary.npz"
        ]
        
        for path in possible_paths:
            if path and os.path.exists(path):
                return path
        
        # Search recursively
        import glob
        search_paths = ["data", "/kb/module/data", "."]
        for search_path in search_paths:
            if os.path.exists(search_path):
                centroids_files = glob.glob(os.path.join(search_path, "**/*family_centroids*.npz"), recursive=True)
                if centroids_files:
                    return centroids_files[0]
        
        return None
    
    def _load_precomputed_indexes(self):
        """Load pre-computed similarity indexes."""
        try:
            # Load hierarchical index if available
            index_path = self.config.similarity_config.get('index_path')
            if index_path and os.path.exists(index_path):
                self.hierarchical_index.load_index(index_path)
                self.logger.info(f"Loaded hierarchical index from: {index_path}")
            
            # Load streaming index if available
            streaming_path = self.config.similarity_config.get('streaming_index_path')
            if streaming_path and os.path.exists(streaming_path):
                self.streaming_index.load_index(streaming_path)
                self.logger.info(f"Loaded streaming index from: {streaming_path}")
                
        except Exception as e:
            self.logger.warning(f"Failed to load pre-computed indexes: {str(e)}")
    
    def execute(self, input_data: Dict[str, Any] = None) -> WorkflowResult:
        """
        Execute the complete protein query analysis workflow.
        
        Args:
            input_data: Input data for the workflow
            
        Returns:
            WorkflowResult containing execution results
        """
        start_time = time.time()
        
        try:
            self.logger.info(f"Starting protein query analysis workflow: {self.run_id}")
            
            # Prepare input data
            if input_data is None:
                input_data = self._prepare_input_data()
            
            # Execute stages in dependency order
            result = self._execute_stages(input_data)
            
            # Calculate execution time
            execution_time = time.time() - start_time
            
            # Create workflow result
            workflow_result = WorkflowResult(
                success=result.success,
                run_id=self.run_id,
                stages_completed=self.stages_completed,
                stage_results=self.stage_results,
                final_output=result.output_data,
                execution_time=execution_time,
                error_message=result.error_message,
                warnings=result.warnings,
                metadata={
                    'performance_metrics': self.performance_metrics,
                    'config': self.config.to_dict()
                }
            )
            
            self.logger.info(f"Workflow completed successfully in {execution_time:.2f} seconds")
            return workflow_result
            
        except Exception as e:
            execution_time = time.time() - start_time
            self.logger.error(f"Workflow failed: {str(e)}")
            
            return WorkflowResult(
                success=False,
                run_id=self.run_id,
                stages_completed=self.stages_completed,
                stage_results=self.stage_results,
                final_output={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _prepare_input_data(self) -> Dict[str, Any]:
        """Prepare input data from configuration."""
        input_data = {
            'config': self.config,
            'run_id': self.run_id
        }
        
        # Add input proteins
        if self.config.input_proteins:
            input_data['proteins'] = self.config.input_proteins
        elif self.config.input_file_path:
            input_data['input_file'] = self.config.input_file_path
        elif self.config.workspace_object_ref:
            input_data['workspace_ref'] = self.config.workspace_object_ref
        
        return input_data
    
    def _execute_stages(self, input_data: Dict[str, Any]) -> StageResult:
        """Execute all stages in dependency order."""
        current_data = input_data.copy()
        
        # Get stage execution order
        execution_order = self._get_stage_execution_order()
        
        for stage_name in execution_order:
            if stage_name not in self.config.enabled_stages:
                self.logger.info(f"Skipping disabled stage: {stage_name}")
                continue
            
            try:
                self.logger.info(f"Executing stage: {stage_name}")
                stage_start_time = time.time()
                
                # Get stage class and execute
                stage_class = get_stage_class(stage_name)
                stage = stage_class(self.config.stage_configs.get(stage_name, {}))
                
                # Execute stage
                result = stage.execute(current_data, self.config.workspace_client)
                
                # Record performance metrics
                stage_time = time.time() - stage_start_time
                self.performance_metrics[stage_name] = {
                    'execution_time': stage_time,
                    'success': result.success,
                    'output_size': len(str(result.output_data))
                }
                
                # Store result
                self.stage_results[stage_name] = result
                
                if result.success:
                    self.stages_completed.append(stage_name)
                    current_data.update(result.output_data)
                    self.logger.info(f"Stage {stage_name} completed successfully in {stage_time:.2f}s")
                else:
                    self.logger.error(f"Stage {stage_name} failed: {result.error_message}")
                    return result
                
                # Check if we should stop after this stage
                if self.config.stop_after_stage == stage_name:
                    self.logger.info(f"Stopping after stage: {stage_name}")
                    break
                
            except Exception as e:
                self.logger.error(f"Failed to execute stage {stage_name}: {str(e)}")
                return StageResult(
                    success=False,
                    output_data={},
                    metadata={},
                    execution_time=time.time() - time.time(),
                    error_message=f"Stage {stage_name} failed: {str(e)}"
                )
        
        return StageResult(
            success=True,
            output_data=current_data,
            metadata={'stages_completed': self.stages_completed},
            execution_time=0.0
        )
    
    def _get_stage_execution_order(self) -> List[str]:
        """Get the order of stage execution based on dependencies."""
        # Define stage dependencies
        dependencies = {
            'input_validation': [],
            'data_extraction': ['input_validation'],
            'workspace_object': ['input_validation'],
            'embedding_generation': ['data_extraction', 'workspace_object'],
            'family_assignment': ['embedding_generation'],
            'similarity_search': ['embedding_generation'],
            'sequence_analysis': ['data_extraction', 'workspace_object'],
            'network_analysis': ['similarity_search', 'family_assignment'],
            'bioinformatics_analysis': ['sequence_analysis'],
            'report_generation': ['network_analysis', 'bioinformatics_analysis'],
            'visualization': ['network_analysis'],
            'data_export': ['report_generation', 'visualization']
        }
        
        # Topological sort
        execution_order = []
        visited = set()
        temp_visited = set()
        
        def visit(stage):
            if stage in temp_visited:
                raise ValueError(f"Circular dependency detected: {stage}")
            if stage in visited:
                return
            
            temp_visited.add(stage)
            
            for dep in dependencies.get(stage, []):
                visit(dep)
            
            temp_visited.remove(stage)
            visited.add(stage)
            execution_order.append(stage)
        
        for stage in dependencies.keys():
            if stage not in visited:
                visit(stage)
        
        return execution_order
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get performance summary of the workflow execution."""
        total_time = sum(metrics['execution_time'] for metrics in self.performance_metrics.values())
        
        return {
            'total_execution_time': total_time,
            'stages_completed': len(self.stages_completed),
            'stage_performance': self.performance_metrics,
            'success_rate': len([m for m in self.performance_metrics.values() if m['success']]) / len(self.performance_metrics) if self.performance_metrics else 0
        }
    
    def cleanup(self):
        """Clean up resources after workflow execution."""
        try:
            # Clear memory
            gc.collect()
            
            # Clear stage results
            self.stage_results.clear()
            
            # Clear performance metrics
            self.performance_metrics.clear()
            
            self.logger.info("Workflow cleanup completed")
            
        except Exception as e:
            self.logger.warning(f"Cleanup failed: {str(e)}")