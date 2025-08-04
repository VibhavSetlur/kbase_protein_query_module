package us.kbase.kbaseproteinquerymodule;

import java.util.List;
import java.util.Map;

/**
 * <p>Original spec-file type: FindTopMatchesFromEmbeddingResults</p>
 * <pre>
 * Find top matches for a given protein embedding.
 * FindTopMatchesFromEmbeddingResults is a reference to a hash where the following keys are defined:
 * report_name has a value which is a string
 * report_ref has a value which is a string
 * top_matches has a value which is a list of ProteinSimilaritySearchResult
 * summary has a value which is a string
 * input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
 * start_time has a value which is a float
 * similarity_search_result_ref has a value which is a string
 * </pre>
 */
public class FindTopMatchesFromEmbeddingResults {
    private String reportName;
    private String reportRef;
    private List<ProteinSimilaritySearchResult> topMatches;
    private String summary;
    private Map<String, Object> inputParameters;
    private Double startTime;
    private String similaritySearchResultRef;

    public FindTopMatchesFromEmbeddingResults() {
    }

    public String getReportName() {
        return reportName;
    }

    public void setReportName(String reportName) {
        this.reportName = reportName;
    }

    public String getReportRef() {
        return reportRef;
    }

    public void setReportRef(String reportRef) {
        this.reportRef = reportRef;
    }

    public List<ProteinSimilaritySearchResult> getTopMatches() {
        return topMatches;
    }

    public void setTopMatches(List<ProteinSimilaritySearchResult> topMatches) {
        this.topMatches = topMatches;
    }

    public String getSummary() {
        return summary;
    }

    public void setSummary(String summary) {
        this.summary = summary;
    }

    public Map<String, Object> getInputParameters() {
        return inputParameters;
    }

    public void setInputParameters(Map<String, Object> inputParameters) {
        this.inputParameters = inputParameters;
    }

    public Double getStartTime() {
        return startTime;
    }

    public void setStartTime(Double startTime) {
        this.startTime = startTime;
    }

    public String getSimilaritySearchResultRef() {
        return similaritySearchResultRef;
    }

    public void setSimilaritySearchResultRef(String similaritySearchResultRef) {
        this.similaritySearchResultRef = similaritySearchResultRef;
    }
} 