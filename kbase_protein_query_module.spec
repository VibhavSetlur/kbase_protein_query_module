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
    
    typedef structure {
        string report_name;
        string report_ref;
        int exists;
        string family_id;
        mapping<string, UnspecifiedObject> metadata;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        string summary;
    } CheckProteinExistenceResults;
    
    typedef structure {
        string embedding_result_ref;
        string summary;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        float embedding_norm;
        int sequence_length;
        int embedding_dim;
    } GenerateProteinEmbeddingResults;
    
    typedef structure {
        string family_id;
        float confidence;
        string eigenprotein_id;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        string family_assignment_result_ref;
    } AssignFamilyFastResults;
    
    typedef structure {
        list<UnspecifiedObject> matches;
        string summary;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        string family_id;
        int top_n;
        mapping<string, UnspecifiedObject> similarity_stats;
        string similarity_search_result_ref;
    } FindTopMatchesResults;
    
    typedef structure {
        string analysis_result_ref;
        string output_directory;
        string general_info_dir;
        string network_analysis_dir;
        string sequence_analysis_dir;
        string embeddings_file_path;
        string top_proteins_csv_path;
        string summary;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        int protein_count;
        list<string> stages_completed;
    } SummarizeAndVisualizeResults;
    
    typedef structure {
        mapping<string, UnspecifiedObject> available_analyses;
        string summary;
    } GetAvailableAnalysesResults;
    
    funcdef run_protein_query_analysis(mapping<string, UnspecifiedObject> params) returns (ProteinQueryAnalysisResults output) authentication required;
    funcdef check_protein_existence(mapping<string, UnspecifiedObject> params) returns (CheckProteinExistenceResults output) authentication required;
    funcdef generate_protein_embedding(mapping<string, UnspecifiedObject> params) returns (GenerateProteinEmbeddingResults output) authentication required;
    funcdef assign_family_fast(mapping<string, UnspecifiedObject> params) returns (AssignFamilyFastResults output) authentication required;
    funcdef find_top_matches_from_embedding(mapping<string, UnspecifiedObject> params) returns (FindTopMatchesResults output) authentication required;
    funcdef summarize_and_visualize_results(mapping<string, UnspecifiedObject> params) returns (SummarizeAndVisualizeResults output) authentication required;
    funcdef get_available_analyses() returns (GetAvailableAnalysesResults output) authentication required;
};
