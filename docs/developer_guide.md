# Developer Guide (Concise)

Code layout (`lib/kbase_protein_query_module/src/`):
- `core/`: pipeline config, registry, resource/parallel/performance utils
- `stages/`: input, processing, analysis, output stages
- `processing/`: embeddings, similarity, networks
- `storage/`: family/embeddings storage, metadata, indexing
- `reports/`: HTML and file organization helpers
- `utils/`: input parsing, documentation helpers
- `workflows/`: `ProteinQueryWorkflow` orchestrates stages

Entry points
- KIDL: `kbase_protein_query_module.spec` (single method)
- Implementation: `lib/kbase_protein_query_module/kbase_protein_query_moduleImpl.py`
- UI: `ui/narrative/methods/ProteinQueryAnalysis/spec.json`

Extend safely
- Add new analysis: see `analysis_extensions.md`
- Add new index strategy: see `indexing_and_similarity.md`
- Update outputs: use `FileOrganizer` and Report generation stages
- Keep changes inside `#BEGIN/#END` in Impl to survive `kb-sdk compile`

Testing
- Unit: `pytest` with `TEST_MODE=true`
- Integration: `kb-sdk test` (requires KBase runtime)

Conventions
- Prefer single config map into workflow
- Validate inputs early (`InputValidationStage`)
- Organize outputs under one directory for Narrative
