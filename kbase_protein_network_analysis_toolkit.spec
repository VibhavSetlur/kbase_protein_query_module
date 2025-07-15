/*
A KBase module: kbase_protein_network_analysis_toolkit
*/

module kbase_protein_network_analysis_toolkit {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_kbase_protein_network_analysis_toolkit(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
