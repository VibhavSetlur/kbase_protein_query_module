/*
A KBase module: kbase_protein_query_module

ProteinQueryAnalysis app.

*/

module kbase_protein_query_module {
    
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
