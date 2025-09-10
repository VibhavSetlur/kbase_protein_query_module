# Add a New Analysis (Quick Steps)

Goal: plug a new analysis into the pipeline and outputs.

1) Create analysis class
- Location: `lib/.../src/stages/analysis/`
- Inherit from `BaseAnalysis` (see `core/analysis_registry.py`)
- Implement: `analyze(proteins, **kwargs)`, `get_output_files(output_dir)`, `get_metadata()`

2) Register analysis
- Use `@register_analysis("my_analysis")` decorator
- Provide meaningful `AnalysisMetadata`

3) Produce outputs
- Write files into a stage/output directory
- Return data structures consumable by `ReportGenerationStage`/`FileOrganizer`

4) Wire into workflow
- Add config flags or include in `enabled_stages`
- Integrate in `ProteinQueryWorkflow` if needed

5) Test
- Unit tests for analyze() on small inputs (TEST_MODE)
- Verify files appear in `output_directory` and summary has entries
