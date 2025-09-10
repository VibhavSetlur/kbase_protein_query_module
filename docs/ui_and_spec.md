# UI and Spec (Mapping)

KIDL
- Single funcdef: `run_protein_query_analysis(params)` → `ProteinQueryAnalysisResults`
- Results: `report_name`, `report_ref`, `analysis_result_ref`, optional metadata

UI (`spec.json`)
- No output widget: `widgets.output = "no-display"`
- Inputs mapped via `behavior.service-mapping.input_mapping`:
  `input_type`, `input_proteins`, `uniprot_ids`, `fasta_file`, `workspace_object`, `sequence_strings`, `analysis_type`, `enabled_stages`, `stop_after_stage`, `output_directory`, `output_report`, `output_data`, `input_data`
- Outputs mapped: `report_name`, `report_ref`, `analysis_result_ref`

Keep UI minimal. Route everything to one backend params map.
