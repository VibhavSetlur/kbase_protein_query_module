/*
A KBase module: kbase_protein_query_module

ProteinQueryAnalysis app.

*/

module kbase_protein_query_module {
    
    typedef structure {
        string job_id;
        string analysis_result_ref;
        string summary;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        int protein_count;
        list<string> stages_completed;
        string output_directory;
        string general_info_dir;
        string network_analysis_dir;
        string sequence_analysis_dir;
        string embeddings_file_path;
        string top_proteins_csv_path;
    } ProteinQueryAnalysisResults;
    
    funcdef run_protein_query_analysis(mapping<string, UnspecifiedObject> params) returns (ProteinQueryAnalysisResults output) authentication required;
};
