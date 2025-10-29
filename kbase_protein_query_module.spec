/*
A KBase module: kbase_protein_query_module

Protein query and analysis module with comprehensive network analysis capabilities.
*/

module kbase_protein_query_module {
    
    typedef structure {
        string report_name;
        string report_ref;
    } ProteinQueryAnalysisResults;
    
    typedef structure {
        mapping<string, UnspecifiedObject> available_analyses;
        string summary;
    } GetAvailableAnalysesResults;
    
    typedef structure {
        string state;
        string message;
        string version;
        string git_url;
        string git_commit_hash;
    } StatusResults;
    
    funcdef run_protein_query_analysis(mapping<string, UnspecifiedObject> params) returns (ProteinQueryAnalysisResults output) authentication required;
    funcdef get_available_analyses() returns (GetAvailableAnalysesResults output) authentication required;
    funcdef status() returns (StatusResults output) authentication required;
};
