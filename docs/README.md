# Documentation

This directory contains documentation for the kbase_protein_query_module.

## Getting Started

- **[Quick Start](QUICK_START.md)** - Get started in 5 minutes with a simple example
- **[Creating Analyses](CREATING_ANALYSES.md)** - Complete guide to create new analyses
- **[Creating Utilities](CREATING_UTILITIES.md)** - Guide to create reusable utility modules
- **[Config Guide](CONFIG_GUIDE.md)** - How to configure analyses and system settings
- **[Workflow Guide](WORKFLOW_GUIDE.md)** - Understanding the workflow system

## Development Guides

- **[Analysis Extensions](analysis_extensions.md)** - Advanced analysis development
- **[Testing](testing.md)** - Testing guidelines and practices

## Quick Reference

### Creating a New Analysis

1. Create directory: `src/analysis/your_analysis/`
2. Create analysis class with `run_network_analysis()` method
3. Add self-test with `main()` function
4. Register in `src/analysis/config.py`
5. Test with self-test: `python your_analysis.py`

### Creating a New Utility

1. Create directory: `src/util/your_utility/`
2. Create utility class
3. Add `__init__.py` with exports
4. Import and use in analyses

### Adding Configuration

1. Add analysis to `_ANALYSIS_BASE` in `config.py`
2. Set `module_path` and `class_name`
3. Add optional `requires_deps` for dependencies
4. Access config in analysis `__init__` method

### Understanding Workflow

1. InputManager processes user input
2. WorkflowOrchestrator coordinates execution
3. AnalysisManager runs selected analyses
4. OutputManager saves results
5. Results returned to KBase

## Testing

All analyses must include self-tests. Run:

```bash
python lib/kbase_protein_query_module/src/analysis/your_analysis/your_analysis.py
```

Expected output: `ANALYSIS_OK` or `ANALYSIS_FAIL: <error>`

## Module Overview

Purpose: Unified protein query analysis for KBase Narrative.

- **Inputs**: UniProt IDs, FASTA files, Workspace objects, or direct sequences
- **Pipeline**: input parsing → embeddings → similarity search → network analysis → report
- **Output**: KBase report + organized output files (TSV, HTML, networks)

Key locations:
- Method: `run_protein_query_analysis` (see `kbase_protein_query_module.spec`)
- UI: `ui/narrative/methods/ProteinQueryAnalysis/spec.json`
- Code: `lib/kbase_protein_query_module/src/`
