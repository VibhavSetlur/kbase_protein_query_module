"""
Base stage abstractions for the Protein Query workflow.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass
class StageResult:
    success: bool
    output_data: Dict[str, Any]
    metadata: Dict[str, Any]
    execution_time: float
    error_message: Optional[str] = None


class BaseStage:
    """
    Minimal abstract base for pipeline stages. Concrete stages should subclass
    and implement run() and any schema/validation as needed.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        self.config = config or {}

    def get_stage_name(self) -> str:
        return self.__class__.__name__.lower()

    def get_required_inputs(self) -> List[str]:
        return []

    def get_optional_inputs(self) -> List[str]:
        return []

    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        return True

    def get_output_schema(self) -> Dict[str, Any]:
        return {}

    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        raise NotImplementedError("Stages must implement run()")


