"""
Base Stage for Protein Query Analysis Pipeline

This module defines the base class for all pipeline stages in the Protein Query Analysis workflow.
Each stage must inherit from this class and implement the required methods.
"""

import logging
import time
from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class StageResult:
    """Result container for pipeline stage execution."""
    success: bool
    output_data: Dict[str, Any]
    metadata: Dict[str, Any]
    execution_time: float
    error_message: Optional[str] = None
    warnings: List[str] = None
    
    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


class BaseStage(ABC):
    """
    Base class for all pipeline stages.
    
    Each stage must implement:
    - run(): Execute the stage logic
    - validate_input(): Validate input data
    - get_output_schema(): Define output data structure
    - get_stage_name(): Return stage identifier
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        Initialize the stage with configuration.
        
        Args:
            config: Stage-specific configuration parameters
        """
        self.config = config or {}
        self.stage_name = self.get_stage_name()
        self.logger = logging.getLogger(f"{__name__}.{self.stage_name}")
    
    @abstractmethod
    def get_stage_name(self) -> str:
        """Return the stage name/identifier."""
        pass
    
    @abstractmethod
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """
        Validate input data for this stage.
        
        Args:
            input_data: Input data dictionary
            
        Returns:
            True if input is valid, False otherwise
        """
        pass
    
    @abstractmethod
    def get_output_schema(self) -> Dict[str, Any]:
        """
        Define the output data schema for this stage.
        
        Returns:
            Dictionary describing the expected output structure
        """
        pass
    
    @abstractmethod
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """
        Execute the stage logic.
        
        Args:
            input_data: Input data from previous stages
            workspace_client: KBase workspace client for object storage
            
        Returns:
            StageResult containing execution results
        """
        pass
    
    def pre_run_hook(self, input_data: Dict[str, Any]) -> None:
        """
        Hook called before stage execution.
        Override to add pre-execution logic.
        
        Args:
            input_data: Input data for the stage
        """
        self.logger.info(f"Starting {self.stage_name} stage")
    
    def post_run_hook(self, result: StageResult) -> None:
        """
        Hook called after stage execution.
        Override to add post-execution logic.
        
        Args:
            result: Stage execution result
        """
        if result.success:
            self.logger.info(f"Completed {self.stage_name} stage in {result.execution_time:.2f}s")
        else:
            self.logger.error(f"Failed {self.stage_name} stage: {result.error_message}")
    
    def execute(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """
        Execute the stage with timing and error handling.
        
        Args:
            input_data: Input data from previous stages
            workspace_client: KBase workspace client
            
        Returns:
            StageResult with execution results
        """
        start_time = time.time()
        
        try:
            # Pre-execution hook
            self.pre_run_hook(input_data)
            
            # Validate input
            if not self.validate_input(input_data):
                return StageResult(
                    success=False,
                    output_data={},
                    metadata={},
                    execution_time=time.time() - start_time,
                    error_message="Input validation failed"
                )
            
            # Execute stage
            result = self.run(input_data, workspace_client)
            result.execution_time = time.time() - start_time
            
            # Post-execution hook
            self.post_run_hook(result)
            
            return result
            
        except Exception as e:
            execution_time = time.time() - start_time
            self.logger.error(f"Stage {self.stage_name} failed with exception: {e}")
            
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def get_required_inputs(self) -> List[str]:
        """
        Get list of required input keys for this stage.
        
        Returns:
            List of required input keys
        """
        return []
    
    def get_optional_inputs(self) -> List[str]:
        """
        Get list of optional input keys for this stage.
        
        Returns:
            List of optional input keys
        """
        return []
    
    def get_config_schema(self) -> Dict[str, Any]:
        """
        Define the configuration schema for this stage.
        
        Returns:
            Dictionary describing expected configuration parameters
        """
        return {}
    
    def validate_config(self, config: Dict[str, Any]) -> bool:
        """
        Validate stage configuration.
        
        Args:
            config: Configuration dictionary
            
        Returns:
            True if configuration is valid
        """
        return True
    
    def get_stage_description(self) -> str:
        """
        Get human-readable description of this stage.
        
        Returns:
            Stage description
        """
        return f"Stage: {self.stage_name}"
    
    def get_stage_dependencies(self) -> List[str]:
        """
        Get list of stage names this stage depends on.
        
        Returns:
            List of stage names that must run before this stage
        """
        return []
