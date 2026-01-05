# KBase Protein Query Module Source

This directory contains the source code for the KBase Protein Query Module.

## Directory Structure

*   **`analysis/`**: Contains logic for running various analyses (e.g., network analysis).
*   **`core/`**: Core workflow orchestration and pipeline configuration.
*   **`input/`**: Handlers for different input types (Protein Sequence, UniProt ID, Workspace Object).
*   **`output/`**: Manages output generation (HTML reports, JSON).
*   **`util/`**: Utility modules for embeddings, FAISS, similarity search, storage, and UniProt API interactions.

## Usage

This package is designed to be used as part of the KBase SDK module. The main entry point is typically through the `WorkflowOrchestrator` in `core/`.
