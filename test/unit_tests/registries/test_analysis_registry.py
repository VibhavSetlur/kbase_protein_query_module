import typing as _t
import os

from lib.kbase_protein_query_module.src.core.analysis_registry import (
    BaseAnalysis,
    AnalysisMetadata,
    AnalysisRegistry,
    register_analysis,
    create_analysis,
    list_available_analyses,
    get_registry,
)


class _DummyAnalysis(BaseAnalysis):
    def analyze(self, proteins: _t.List[_t.Any], **kwargs) -> dict:
        return {"num_proteins": len(proteins)}

    def get_output_files(self, output_dir: str) -> _t.List[str]:
        return [f"{output_dir}/dummy.json"]

    def get_metadata(self) -> AnalysisMetadata:
        return AnalysisMetadata(
            name="dummy",
            description="Dummy analysis for tests",
            version="0.0.1",
            author="tests",
            output_files=["dummy.json"],
            dependencies=[],
            category="test",
            computational_complexity="low",
            memory_requirements="low",
            supports_batch=True,
            supports_streaming=False,
        )


def test_registry_manual_registration_is_retrievable():
    registry = AnalysisRegistry()
    registry.register_analysis("dummy_manual", _DummyAnalysis)

    instance = registry.get_analysis("dummy_manual", config={"x": 1})
    assert isinstance(instance, _DummyAnalysis)
    files = instance.get_output_files("/tmp")
    assert any(os.path.basename(p) == "dummy.json" for p in files)

    meta = registry.list_analyses()["dummy_manual"]
    assert meta.name == "dummy"


def test_register_analysis_decorator_registers_globally():
    @register_analysis("dummy_decorated")
    class _Decorated(_DummyAnalysis):
        pass

    # Global registry should now know about it
    reg = get_registry()
    names = reg.list_analyses().keys()
    assert "dummy_decorated" in names

    # Convenience helpers work
    created = create_analysis("dummy_decorated", config={})
    assert isinstance(created, _Decorated)


def test_list_available_analyses_returns_metadata_dict():
    analyses = list_available_analyses()
    assert isinstance(analyses, dict)
    # At least contains our decorated test analysis
    assert any(name.startswith("dummy_") for name in analyses)


