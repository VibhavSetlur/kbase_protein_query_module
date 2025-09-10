# kbase_protein_query_module – Overview

Purpose: Unified protein query analysis for KBase Narrative.
- Inputs: UniProt IDs, FASTA file, Workspace objects, or direct sequences
- Pipeline: input parsing → local ESM embeddings → family assignment → similarity search → sequence analysis → network analysis → report + organized outputs
- Output: KBase report + output directory (CSV/TSV, HTML, networks)

Key locations:
- Method: `run_protein_query_analysis` (see `kbase_protein_query_module.spec`)
- UI: `ui/narrative/methods/ProteinQueryAnalysis/spec.json`
- Code: `lib/kbase_protein_query_module/src/`

Next:
- See `usage.md` (how to run in Narrative)
- See `developer_guide.md` (code layout)
- See `analysis_extensions.md` (add new analyses)
- See `embedding_and_models.md` (local model)
- See `indexing_and_similarity.md` (search/index)
- See `ui_and_spec.md` (UI + KIDL)
- See `testing.md` (pytest + kb-sdk)
