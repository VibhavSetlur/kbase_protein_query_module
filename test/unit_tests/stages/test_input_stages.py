import types

from lib.kbase_protein_query_module.src.input.input_validation import InputValidationStage
from lib.kbase_protein_query_module.src.input.data_extraction import DataExtractionStage


def test_validate_uniprot_identifiers_and_extracts_minimally():
    validator = InputValidationStage()
    data = {"input_type": "uniprot_identifiers", "input_data": ["P12345", "BAD"]}
    vres = validator.run(data)
    assert vres.success
    assert vres.output_data["validation_stats"]["total_inputs"] == 2
    assert vres.output_data["validation_stats"]["valid_inputs"] >= 1

    extractor = DataExtractionStage(config={"batch_size": 10})
    eres = extractor.run({"validated_input": vres.output_data["validated_input"]})
    assert eres.success
    assert "extracted_data" in eres.output_data


def test_validate_fasta_string_and_extracts():
    fasta = ">seq1\nMTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRK\n>seq2\nGPPPGKKKSKKAVKKK\n"
    validator = InputValidationStage()
    vres = validator.run({"input_type": "fasta", "input_data": fasta})
    assert vres.success
    vdict = vres.output_data.get("validated_data")
    assert isinstance(vdict, dict)

    extractor = DataExtractionStage()
    eres = extractor.run({"validated_data": vdict})
    assert eres.success
    recs = eres.output_data["extracted_data"]["protein_records"]
    assert len(recs) >= 1


def test_workspace_ref_validation_requires_client():
    validator = InputValidationStage()
    # Missing client should fail validation at stage runtime
    vres = validator.run({"input_type": "workspace_object", "input_data": "1/1"})
    assert vres.success is False or "error" in (vres.error_message or "")


