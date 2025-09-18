import pytest

from lib.kbase_protein_query_module.src.core.pipeline_config import PipelineConfig


def test_pipeline_config_minimal_inputs_accepts_one_of_three():
    # input_proteins path
    cfg = PipelineConfig(input_proteins=["AAA"])
    assert cfg.input_proteins == ["AAA"]

    # file path
    cfg2 = PipelineConfig(input_file_path="/tmp/fake.faa")
    assert cfg2.input_file_path.endswith("fake.faa")

    # workspace ref
    cfg3 = PipelineConfig(workspace_object_ref="1/1/1")
    assert cfg3.workspace_object_ref == "1/1/1"


def test_pipeline_config_validates_thresholds():
    with pytest.raises(ValueError):
        PipelineConfig(input_proteins=["AAA"], similarity_threshold=1.5)

    with pytest.raises(ValueError):
        PipelineConfig(input_proteins=["AAA"], max_similar_proteins=0)


def test_pipeline_config_auto_scaling_fields_present():
    cfg = PipelineConfig(input_proteins=["AAA"])
    # Ensure auto-config filled sane minimums
    assert cfg.max_workers >= 1
    assert cfg.batch_size_proteins >= 50
    assert cfg.cache_size_mb >= 128


def test_pipeline_config_to_from_dict_roundtrip():
    cfg = PipelineConfig(input_proteins=["AAA"], stage_configs={"sequence_analysis": {"x": 1}})
    d = cfg.to_dict()
    restored = PipelineConfig.from_dict(d)
    assert restored.stage_configs["sequence_analysis"]["x"] == 1


