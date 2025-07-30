/*
A KBase module: kbase_protein_query_module
*/

module kbase_protein_query_module {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_kbase_protein_query_module(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

    /*
        Check if a protein exists in the storage system.
    */
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
        string input_id;
        string input_type;
        int exists;
        string matched_id;
        string matched_db;
        mapping<string, string> metadata;
    } ProteinExistenceResult;
    funcdef check_protein_existence(mapping<string, UnspecifiedObject> params) returns (CheckProteinExistenceResults output) authentication required;

    /*
        Generate a protein embedding from a sequence.
    */
    typedef structure {
        string report_name;
        string report_ref;
        string embedding_result_ref;
        string summary;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        float embedding_norm;
        int sequence_length;
        int embedding_dim;
    } GenerateProteinEmbeddingResults;
    funcdef generate_protein_embedding(mapping<string, UnspecifiedObject> params) returns (GenerateProteinEmbeddingResults output) authentication required;

    /*
        Quickly assign a protein embedding to a family by similarity to the medoid.
    */
    typedef structure {
        string family_id;
        float confidence;
        string eigenprotein_id;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
    } AssignFamilyFastResults;
    funcdef assign_family_fast(mapping<string, UnspecifiedObject> params) returns (AssignFamilyFastResults output) authentication required;

    /*
        Find top matches for a given protein embedding.
    */
    typedef structure {
        list<mapping<string, UnspecifiedObject>> matches;
        string summary;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        string family_id;
        int top_n;
        mapping<string, float> similarity_stats;
    } FindTopMatchesFromEmbeddingResults;
    funcdef find_top_matches_from_embedding(mapping<string, UnspecifiedObject> params) returns (FindTopMatchesFromEmbeddingResults output) authentication required;

    /*
        Summarize and visualize protein network analysis results.
    */
    typedef structure {
        string report_name;
        string report_ref;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        string output_dir;
        string summary;
    } SummarizeAndVisualizeResultsResults;
    funcdef summarize_and_visualize_results(mapping<string, UnspecifiedObject> params) returns (SummarizeAndVisualizeResultsResults output) authentication required;

    typedef structure {
        string input_id;
        string input_type;
        string embedding_ref;
        list<float> embedding;
        string model_name;
        string pooling_method;
        mapping<string, string> metadata;
    } ProteinEmbeddingResult;

    typedef structure {
        string input_id;
        string input_type;
        string embedding_ref;
        string assigned_family_id;
        float similarity_score;
        mapping<string, string> metadata;
    } ProteinFamilyAssignmentResult;

    typedef structure {
        string input_id;
        string input_type;
        string top_matches_result_ref;
        string summary_html;
        mapping<string, string> metadata;
    } SummarizeVisualizeResult;
};
