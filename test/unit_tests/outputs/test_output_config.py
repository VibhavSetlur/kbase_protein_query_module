from lib.kbase_protein_query_module.src.outputs.config import (
    get_output_config,
    get_analysis_output_config,
    get_enabled_output_analyses,
    is_output_enabled_for_analysis,
    get_output_formats_for_analysis,
    get_file_naming_config,
    get_directory_config,
    validate_output_config,
)


def test_output_config_has_expected_keys():
    cfg = get_output_config()
    assert cfg["base_output_dir"]
    assert isinstance(cfg["formats"], list)


def test_analysis_output_toggles_and_formats():
    enabled = get_enabled_output_analyses()
    assert isinstance(enabled, dict)

    # Known analysis names from config
    assert is_output_enabled_for_analysis("network_analysis") is True
    formats = get_output_formats_for_analysis("sequence_analysis")
    assert "json" in formats


def test_file_and_dir_configs_present():
    files = get_file_naming_config()
    assert "final_output_file" in files

    dirs = get_directory_config()
    assert "analysis_subdir" in dirs


def test_validate_output_config_true():
    assert validate_output_config() is True


