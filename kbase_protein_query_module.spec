/*
A KBase module: kbase_protein_query_module

This module provides comprehensive protein query analysis capabilities using UniProt IDs as the canonical identifier:

COMPREHENSIVE ANALYSIS WORKFLOW:
1. CheckProteinExistence: Verify protein exists using UniProt ID, optionally generate embedding
2. GenerateProteinEmbeddings: Create embeddings from sequence input or protein check results
3. AssignProteinFamily: Assign proteins to families using similarity to centroids
4. FindTopMatches: Perform similarity search within families
5. SummarizeAndVisualize: Generate comprehensive HTML reports with network analysis
6. RunProteinQueryAnalysis: Unified pipeline for comprehensive protein analysis

ADVANCED CAPABILITIES:
- UniProt ID canonical identifier system (exact match only)
- ESM-2 protein language model for embedding generation
- Efficient FAISS-based similarity search and clustering
- Family assignment using binary centroid similarity
- Comprehensive metadata storage and retrieval
- HTML report generation with network visualization
- Workspace object management for downstream analysis
- Bioinformatics integration with protein databases
- Network analysis and protein relationship mapping
- Advanced similarity metrics and statistical analysis
- Modular pipeline architecture with configurable stages
- Real-time performance monitoring and error handling

Authors: Vibhav Setlur
Contact: https://kbase.us/contact-us/
*/

module kbase_protein_query_module {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        Check if a protein exists in the storage system using UniProt ID and create a workspace object with the result.
        Input: UniProt ID (e.g., P00001, P12345)
        Output: Existence status, family assignment, metadata, optional embedding
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
        string embedding_result_ref;
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
        Generate protein embeddings from direct sequence input.
        Creates embeddings using ESM-2 model for downstream analysis.
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
        string input_source;
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
        Assign a protein embedding to a family using similarity to family centroids.
        Uses binary Hamming distance for fast family assignment.
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
        Find top matches for a given protein embedding within a family.
        Uses FAISS IVF float index for efficient similarity search.
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
        Generates comprehensive HTML reports with network visualization.
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

    /*
        Unified Protein Query Analysis Pipeline
        
        This method provides a single entry point for comprehensive protein analysis,
        supporting multiple input types and configurable analysis stages.
    */
    typedef structure {
        string report_name;
        string report_ref;
        string analysis_result_ref;
        string summary;
        mapping<string, UnspecifiedObject> input_parameters;
        float start_time;
        string html_report_path;
        int protein_count;
        list<string> stages_completed;
    } ProteinQueryAnalysisResults;
    
    funcdef run_protein_query_analysis(mapping<string, UnspecifiedObject> params) returns (ProteinQueryAnalysisResults output) authentication required;
};
