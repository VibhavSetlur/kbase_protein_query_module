from lib.kbase_protein_query_module.src.core.workflow_orchestrator import WorkflowOrchestrator


def test_stage_selection_filters_unknown_and_returns_summary(tmp_path):
    orch = WorkflowOrchestrator(config={})
    _ = orch.initialize_components(str(tmp_path))
    selected = orch._determine_analyses_to_run(["sequence_analysis", "__bad__"])
    assert "__bad__" not in selected

    # Generate a minimal final output to exercise summary aggregation
    summary = orch._generate_summary({"sequence_analysis": {"ok": True}})
    assert "Completed 1 analyses" in summary


