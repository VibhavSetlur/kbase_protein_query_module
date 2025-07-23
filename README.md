# kbase_protein_network_analysis_toolkit

This repository contains the KBase Protein Network Analysis Toolkit, a modular suite for large-scale protein network analysis, embedding generation, family assignment, and efficient storage/indexing. The toolkit is designed for scalability and integration with the KBase platform.

## Repository Structure

- `lib/` - Source code for all modules and installed clients. See `lib/README.md` for details.
- `test/` - Unit and integration tests. See `test/README.md` for test structure and instructions.
- `scripts/` - Utility and deployment scripts.
- `data/` - Example and generated data for development/testing.
- `ui/` - KBase app UI specifications. See `ui/README.md` for method and UI details.

## Setup and Testing

1. Install the [KBase SDK](https://github.com/kbase/kb_sdk) and its dependencies.
2. Add your KBase developer token to `test_local/test.cfg`.
3. Run:
   ```bash
   make
   kb-sdk test
   ```
4. After making changes, rerun `kb-sdk test` to verify functionality.

## Usage

- For module usage in other SDK modules:
  ```bash
  kb-sdk install kbase_protein_network_analysis_toolkit
  ```
- For app usage, see the [KBase catalog page](https://narrative.kbase.us/#catalog/modules/kbase_protein_network_analysis_toolkit).

## Documentation

- All modules, scripts, and tests are documented with clear docstrings and comments.
- See subdirectory README files for detailed usage and contribution guidelines.
- For UI and method specifications, see `ui/README.md`.

## Contribution Guidelines

- Follow modular, clear commenting practices.
- Update `RELEASE_NOTES.md` for all changes.
- Ensure all new code is documented and tested.

## Help

- [KBase SDK Docs](https://kbase.github.io/kb_sdk_docs/)
- [FAQ](https://kbase.github.io/kb_sdk_docs/references/questions_and_answers.html)
- [Troubleshooting Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html)
