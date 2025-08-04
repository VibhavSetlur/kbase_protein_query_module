/*
A KBase module: kbase_protein_query_module

This module provides comprehensive protein network analysis capabilities including:
- Protein existence checking in database
- Protein embedding generation using ESM-2 models
- Family assignment using similarity to centroids
- Similarity search and network building
- Sequence analysis with bioinformatics integration
- Comprehensive HTML report generation

Authors: Vibhav Setlur
Contact: https://kbase.us/contact-us/
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
        Check if a protein exists in the storage system and create a workspace object with the result.
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
        string protein_existence_result_ref;
    } CheckProteinExistenceResults;
    
    typedef structure {
        string protein_id;
        int exists;
        string family_id;
        mapping<string, UnspecifiedObject> metadata;
        string embedding_ref;
        list<float> embedding;
        string model_name;
        float search_timestamp;
        string summary;
    } ProteinExistenceResult;
    
    funcdef check_protein_existence(mapping<string, UnspecifiedObject> params) returns (CheckProteinExistenceResults output) authentication required;

    /*
        Generate a protein embedding from a sequence or workspace object.
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
    
    typedef structure {
        string input_id;
        string input_type;
        string embedding_ref;
        list<float> embedding;
        string model_name;
        string pooling_method;
        mapping<string, string> metadata;
        int sequence_length;
        float embedding_norm;
        int embedding_dim;
    } ProteinEmbeddingResult;
    
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
        string family_assignment_result_ref;
    } AssignFamilyFastResults;
    
    typedef structure {
        string input_id;
        string input_type;
        string embedding_ref;
        string assigned_family_id;
        float similarity_score;
        mapping<string, string> metadata;
        string eigenprotein_id;
        float confidence;
    } ProteinFamilyAssignmentResult;
    
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
        string similarity_search_result_ref;
    } FindTopMatchesFromEmbeddingResults;
    
    typedef structure {
        string input_id;
        string input_type;
        string embedding_ref;
        string family_id;
        int top_n;
        list<mapping<string, UnspecifiedObject>> matches;
        mapping<string, float> similarity_stats;
        mapping<string, string> metadata;
    } ProteinSimilaritySearchResult;
    
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
        string html_report_path;
        string sequence_analysis_ref;
    } SummarizeAndVisualizeResultsResults;
    
    typedef structure {
        string input_id;
        string input_type;
        string top_matches_result_ref;
        string summary_html;
        mapping<string, string> metadata;
        string embedding_ref;
        string family_assignment_ref;
        string protein_existence_ref;
        list<mapping<string, UnspecifiedObject>> matches;
        string family_id;
        int top_n;
    } SummarizeVisualizeResult;
    
    funcdef summarize_and_visualize_results(mapping<string, UnspecifiedObject> params) returns (SummarizeAndVisualizeResultsResults output) authentication required;
};
