"""
Workflow Orchestrator for KBase Protein Query Module

Orchestrates the complete protein analysis pipeline with modular components.
"""

import os
import logging
import time
import uuid
from typing import Dict, Any, List, Optional

# Import managers from modular structure
from ..input import InputManager
from ..analysis import AnalysisManager, get_enabled_analyses
from ..output import OutputManager

logger = logging.getLogger(__name__)


class WorkflowOrchestrator:
    """
    Orchestrates the complete protein query workflow using modular components.
    
    Coordinates input handling, analysis execution, and output generation.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None, kb_util=None):
        """
        Initialize the Workflow Orchestrator.
        
        Args:
            config: Configuration dictionary
            kb_util: KBase utility library instance
        """
        self.config = config or {}
        self.kb_util = kb_util
        self.run_id = str(uuid.uuid4())[:8]
        logger.debug(f"Initializing WorkflowOrchestrator with run_id: {self.run_id}")

        # Get available analyses directly
        self.analyses = get_enabled_analyses()
        logger.debug(f"Available analyses: {list(self.analyses.keys())}")

        # Initialize component managers (lazy load)
        self.input_manager: Optional[InputManager] = None
        self.analysis_manager: Optional[AnalysisManager] = None
        self.output_manager: Optional[OutputManager] = None
        
        # Results tracking
        self.results: Dict[str, Any] = {}
        self.execution_start_time: Optional[float] = None
        
        logger.info(f"WorkflowOrchestrator initialized with run_id: {self.run_id}")
    
    def initialize_components(self, output_dir: str, workspace_name: Optional[str] = None):
        """
        Initialize all workflow components.
        
        Args:
            output_dir: Base directory for outputs
            workspace_name: KBase workspace name if applicable
        """
        try:
            logger.debug(f"Initializing components with output_dir={output_dir}, workspace_name={workspace_name}")
            # Use workspace_name from config if not provided
            if not workspace_name:
                workspace_name = self.config.get('workspace_name')
                logger.debug(f"Using workspace_name from config: {workspace_name}")
            
            # Initialize input manager
            logger.debug("Initializing InputManager")
            self.input_manager = InputManager(
                config=self.config,
                kb_util=self.kb_util
            )
            
            # Initialize output manager with KBUtilLib
            logger.debug("Initializing OutputManager")
            self.output_manager = OutputManager(
                base_output_dir=output_dir,
                run_id=self.run_id,
                workspace_name=workspace_name,
                kb_util=self.kb_util
            )
            
            # Initialize analysis manager with output manager and config
            logger.debug("Initializing AnalysisManager")
            self.analysis_manager = AnalysisManager(output_manager=self.output_manager, config=self.config)
            
            logger.info("All workflow components initialized successfully")
            
        except Exception as e:
            logger.error(f"Error initializing workflow components: {e}", exc_info=True)
            raise
    
    
    def run_workflow(
        self,
        input_data: Dict[str, Any],
        output_dir: Optional[str] = None,
        workspace_name: Optional[str] = None,
        selected_analyses: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        """
        Run the complete protein query workflow: input -> analyses -> outputs.

        Returns a concise result dictionary suitable for KBase report generation.
        """
        self.execution_start_time = time.time()

        # Resolve runtime configuration
        resolved_output_dir = output_dir or self.config.get('output_dir', '/tmp/protein_query_output')
        resolved_workspace = workspace_name or self.config.get('workspace_name')

        try:
            logger.info(f"Starting workflow execution run_id={self.run_id}")
            logger.debug(f"Input data keys: {list(input_data.keys()) if isinstance(input_data, dict) else 'non-dict'}")
            logger.debug(f"Resolved output_dir: {resolved_output_dir}, workspace: {resolved_workspace}")

            # Initialize components once
            if not (self.input_manager and self.analysis_manager and self.output_manager):
                logger.debug("Components not initialized, initializing now")
                self.initialize_components(resolved_output_dir, resolved_workspace)
            else:
                logger.debug("Components already initialized")

            # 1) Process input
            logger.debug("Processing input data")
            processed_data = self.input_manager.process_input(input_data)
            logger.debug(f"Input processing completed, success: {processed_data.get('success') if isinstance(processed_data, dict) else 'unknown'}")
            if isinstance(processed_data, dict) and processed_data.get("success") is False:
                exec_time = float(processed_data.get("processing_time", 0.0) or 0.0)
                result = {
                    "success": False,
                    "run_id": self.run_id,
                    "analyses_completed": [],
                    "analysis_results": {},
                    "final_output": {},
                    "execution_time": exec_time,
                    "output_directory": resolved_output_dir,
                    "stages_completed": [],
                    "error_message": processed_data.get("error_message", "Input processing failed"),
                }
                # Persist minimal metadata
                try:
                    if self.output_manager:
                        self.output_manager.save_process_info([], exec_time)
                except Exception:
                    pass
                return result

            # 2) Determine analyses to run
            if selected_analyses is None:
                # Handle both 'analyses' and 'analysis_stages' parameter names
                selected_analyses = self.config.get('selected_analyses') or self.config.get('analysis_stages') or list(self.analyses.keys())
            logger.debug(f"Selected analyses to run: {selected_analyses}")

            analyses_to_run = list(selected_analyses or [])

            # 3) Run analyses sequentially and save outputs
            logger.debug(f"Running {len(analyses_to_run)} analyses")
            analysis_results: Dict[str, Any] = {}
            failed_analyses: List[str] = []
            for analysis_name in analyses_to_run:
                logger.debug(f"Running analysis: {analysis_name}")
                
                # Get analysis-specific output directory
                analysis_output_dir = self.output_manager.get_analysis_dir(analysis_name)
                
                # Ensure processed_data has output_dir set to analysis directory
                analysis_input_data = processed_data.copy() if processed_data else {}
                analysis_input_data['output_dir'] = analysis_output_dir
                
                # Run analysis with analysis-specific output directory
                result = self.analysis_manager.run_analyses(analysis_name, analysis_input_data, **kwargs)
                if result is None:
                    logger.error(f"Analysis '{analysis_name}' returned None - analysis did not run")
                    failed_analyses.append(analysis_name)
                    # Create error result for failed analysis
                    analysis_results[analysis_name] = {
                        'success': False,
                        'error_message': f"Analysis '{analysis_name}' is not available (likely missing dependencies)",
                        'processing_time': 0.0
                    }
                else:
                    analysis_results[analysis_name] = result
                    logger.debug(f"Analysis {analysis_name} completed")
                    if isinstance(result, dict) and result.get('success') is not False:
                        try:
                            logger.debug(f"Saving output for analysis {analysis_name}")
                            # Save results.json to analysis directory
                            self.output_manager.save_analysis_output(
                                analysis_name, result, self.output_manager.get_root_dir()
                            )
                            
                            # Copy output files to analysis directory if they're not already there
                            output_files = result.get('output_files', [])
                            if output_files:
                                import shutil
                                for file_path in output_files:
                                    if isinstance(file_path, str) and os.path.exists(file_path):
                                        # If file is not in analysis directory, copy it
                                        file_basename = os.path.basename(file_path)
                                        target_path = os.path.join(analysis_output_dir, file_basename)
                                        if file_path != target_path:
                                            shutil.copy2(file_path, target_path)
                                            logger.debug(f"Copied output file to analysis directory: {file_basename}")
                                            # Update file path in result
                                            result['output_files'] = [
                                                os.path.join(analysis_output_dir, os.path.basename(f)) 
                                                if isinstance(f, str) and os.path.exists(f) else f
                                                for f in output_files
                                            ]
                        except Exception as e:
                            logger.error(f"Error saving analysis output for {analysis_name}: {e}", exc_info=True)
                            raise
                    else:
                        logger.warning(f"Analysis {analysis_name} returned unsuccessful result: {result}")
            
            # If any analyses failed to run, log warning
            if failed_analyses:
                logger.warning(f"The following analyses failed to run (likely missing dependencies): {failed_analyses}")
                logger.warning("Required dependencies: transformers (for embeddings), networkx and sklearn (for network analysis)")

            # 4) Build final outputs and persist run metadata
            protein_count = int(len((processed_data or {}).get("proteins", []) or []))
            final_output = {
                "summary": f"Completed {len(analyses_to_run)} analysis step(s)",
                "protein_count": protein_count,
            }

            # Compute execution time (prefer component timings if present)
            exec_time_wall = time.time() - self.execution_start_time
            input_time = float((processed_data or {}).get("processing_time", 0.0) or 0.0)
            analyses_time = 0.0
            for v in (analysis_results or {}).values():
                if isinstance(v, dict):
                    analyses_time += float(v.get("processing_time", 0.0) or 0.0)
            execution_time = max(exec_time_wall, input_time + analyses_time)

            # Persist metadata
            try:
                if self.output_manager:
                    self.output_manager.save_metadata(
                        config=self.config,
                        analyses_run=analyses_to_run,
                        summary=final_output.get("summary", ""),
                    )
                    self.output_manager.save_process_info(analyses_to_run, execution_time)
            except Exception as e:
                logger.warning(f"Failed to save run metadata: {e}")

            # Decide which directory to expose for external zipping/use
            output_directory_value = self.output_manager.get_root_dir()

            result = {
                "success": True,
                "run_id": self.run_id,
                "analyses_completed": analyses_to_run,
                "analysis_results": analysis_results,
                "final_output": final_output,
                "execution_time": execution_time,
                "output_directory": output_directory_value,
                "stages_completed": analyses_to_run,
            }

            logger.info(
                f"Workflow completed successfully in {execution_time:.2f}s; analyses={analyses_to_run}"
            )
            return result

        except Exception as e:
            execution_time = time.time() - self.execution_start_time if self.execution_start_time else 0.0
            logger.error(f"Workflow execution failed: {e}")
            return {
                "success": False,
                "run_id": self.run_id,
                "analyses_completed": [],
                "analysis_results": {},
                "final_output": {},
                "execution_time": execution_time,
                "output_directory": resolved_output_dir,
                "stages_completed": [],
                "error_message": str(e),
            }
    


    def save_results(self, results: Dict[str, Any]) -> Optional[str]:
        """
        Persist a final results snapshot under metadata; returns path if written.
        """
        try:
            if not self.output_manager:
                return None
            return self.output_manager.write_json("metadata", "final_results.json", results)
        except Exception as e:
            logger.warning(f"Failed to save final results snapshot: {e}")
            return None
    

def main() -> int:
    """
    Test WorkflowOrchestrator with a simple workflow.
    """
    ok = True
    try:
        import tempfile
        import shutil
        
        test_dir = tempfile.mkdtemp(prefix='workflow_test_')
        try:
            config = {
                'output_dir': test_dir,
                'workspace_name': 'test_workspace'
            }
            orchestrator = WorkflowOrchestrator(config=config)
            
            # Test with minimal input
            input_data = {
                'input_type': 'protein_sequence',
                'protein_sequence': 'ACDEFGHIKLMNPQRSTVWY'
            }
            
            result = orchestrator.run_workflow(
                input_data=input_data,
                output_dir=test_dir,
                selected_analyses=[]  # Empty to avoid requiring data files
            )
            
            if not isinstance(result, dict):
                raise RuntimeError("Result is not a dict")
            
            if result.get('success') is not True:
                raise RuntimeError(f"Workflow failed: {result.get('error_message', 'Unknown error')}")
            
            print("WorkflowOrchestrator test: SUCCESS")
        finally:
            try:
                shutil.rmtree(test_dir, ignore_errors=True)
            except Exception:
                pass
    except Exception as e:
        ok = False
        print(f"WorkflowOrchestrator test: FAILED - {e}")
        import traceback
        traceback.print_exc()
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())