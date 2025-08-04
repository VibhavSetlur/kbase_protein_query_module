package us.kbase.kbaseproteinquerymodule;

import java.util.List;
import java.util.Map;

/**
 * <p>Original spec-file type: SummarizeAndVisualizeResultsResults</p>
 * <pre>
 * Summarize and visualize protein network analysis results.
 * SummarizeAndVisualizeResultsResults is a reference to a hash where the following keys are defined:
 * report_name has a value which is a string
 * report_ref has a value which is a string
 * summary has a value which is a string
 * input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
 * start_time has a value which is a float
 * visualization_result_ref has a value which is a string
 * network_html has a value which is a string
 * </pre>
 */
public class SummarizeAndVisualizeResultsResults {
    private String reportName;
    private String reportRef;
    private String summary;
    private Map<String, Object> inputParameters;
    private Double startTime;
    private String visualizationResultRef;
    private String networkHtml;

    public SummarizeAndVisualizeResultsResults() {
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

    public String getVisualizationResultRef() {
        return visualizationResultRef;
    }

    public void setVisualizationResultRef(String visualizationResultRef) {
        this.visualizationResultRef = visualizationResultRef;
    }

    public String getNetworkHtml() {
        return networkHtml;
    }

    public void setNetworkHtml(String networkHtml) {
        this.networkHtml = networkHtml;
    }
} 