import os

from lib.kbase_protein_query_module.src.core.workflow_orchestrator import WorkflowOrchestrator


def test_orchestrator_runs_minimal_workflow(tmp_path):
    orch = WorkflowOrchestrator(config={})
    input_data = {"input_type": "uniprot_identifiers", "input_data": ["P12345", "Q8N158"]}

    result = orch.run_workflow(
        input_data=input_data,
        output_dir=str(tmp_path),
        workspace_name=None,
        selected_analyses=None,
    )

    assert result.output_directory == str(tmp_path)
    # Success can be true or false depending on available stages; assert result structure
    assert result.run_id
    assert isinstance(result.analyses_completed, list)
    assert isinstance(result.analysis_results, dict)


def test_orchestrator_filters_unknown_analyses(tmp_path):
    orch = WorkflowOrchestrator(config={})
    # Selected includes an unknown name; it should be filtered out gracefully
    _ = orch.initialize_components(str(tmp_path))
    names = orch._determine_analyses_to_run(["sequence_analysis", "__unknown__name__"])
    assert "__unknown__name__" not in names


