# Usage (KBase Narrative)

1) Add app: ProteinQueryAnalysis
2) Choose input_type:
   - `uniprot_identifiers`: comma/newline-separated IDs
   - `fasta_file`: upload or select a FASTA file
   - `workspace_object`: select supported WS objects
   - `sequence_strings`: paste sequences (ACDEFGHIKLMNPQRSTVWY)
3) Optional: set `output_directory` and `output_report`.
4) Run. The app returns:
   - `report_name`, `report_ref`: KBase Report
   - Organized files in `output_directory` (CSV/TSV, HTML, network files)

Notes
- Model is local; no model selection parameter.
- Outputs are designed for downstream Narrative inspection and downloads.
