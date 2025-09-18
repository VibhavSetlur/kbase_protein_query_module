import pytest

from lib.kbase_protein_query_module.src.util import input_parser as ip


def test_detect_input_type_variants():
    p = ip.InputParser()
    assert p.detect_input_type("P12345") == 'Uniprot'
    assert p.detect_input_type("MTEYKLVVVGAGGVGKSALTIQLI") == 'SingleProtein'
    assert p.detect_input_type({"ref": "1/1"}) == 'WorkspaceObject'
    # Mixed list
    t = p.detect_input_type(["P12345", "MTEYKLVVVGAGG"])
    assert t == 'MixedInput'


def test_parse_fasta_string_helper():
    fasta = ">a\nMTEYKLVVVGAGGVGKSALTI\n>b\nACDEFGHIKLMNPQRSTVWY\n"
    recs = ip.parse_fasta_string(fasta)
    ids = {r.protein_id for r in recs}
    assert ids == {"a", "b"}
    assert all(len(r.sequence) > 0 for r in recs)


def test_single_protein_validation_and_cleanup():
    p = ip.InputParser()
    # Use public API to validate equivalent behavior
    assert p._is_protein_sequence("BADSEQ") is False
    recs = p.parse_input('SingleProtein', "MTEYKLVVVGAGGVGKSALTIQLIQ")
    assert recs[0].protein_id == 'single_protein'


