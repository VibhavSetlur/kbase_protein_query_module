
package us.kbase.kbaseproteinquerymodule;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import us.kbase.common.service.UObject;


/**
 * <p>Original spec-file type: SummarizeAndVisualizeResultsResults</p>
 * <pre>
 * Summarize and visualize protein network analysis results.
 * Generates comprehensive HTML reports with network visualization.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "report_name",
    "report_ref",
    "input_parameters",
    "start_time",
    "output_dir",
    "summary",
    "html_report_path",
    "sequence_analysis_ref"
})
public class SummarizeAndVisualizeResultsResults {

    @JsonProperty("report_name")
    private java.lang.String reportName;
    @JsonProperty("report_ref")
    private java.lang.String reportRef;
    @JsonProperty("input_parameters")
    private Map<String, UObject> inputParameters;
    @JsonProperty("start_time")
    private Double startTime;
    @JsonProperty("output_dir")
    private java.lang.String outputDir;
    @JsonProperty("summary")
    private java.lang.String summary;
    @JsonProperty("html_report_path")
    private java.lang.String htmlReportPath;
    @JsonProperty("sequence_analysis_ref")
    private java.lang.String sequenceAnalysisRef;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("report_name")
    public java.lang.String getReportName() {
        return reportName;
    }

    @JsonProperty("report_name")
    public void setReportName(java.lang.String reportName) {
        this.reportName = reportName;
    }

    public SummarizeAndVisualizeResultsResults withReportName(java.lang.String reportName) {
        this.reportName = reportName;
        return this;
    }

    @JsonProperty("report_ref")
    public java.lang.String getReportRef() {
        return reportRef;
    }

    @JsonProperty("report_ref")
    public void setReportRef(java.lang.String reportRef) {
        this.reportRef = reportRef;
    }

    public SummarizeAndVisualizeResultsResults withReportRef(java.lang.String reportRef) {
        this.reportRef = reportRef;
        return this;
    }

    @JsonProperty("input_parameters")
    public Map<String, UObject> getInputParameters() {
        return inputParameters;
    }

    @JsonProperty("input_parameters")
    public void setInputParameters(Map<String, UObject> inputParameters) {
        this.inputParameters = inputParameters;
    }

    public SummarizeAndVisualizeResultsResults withInputParameters(Map<String, UObject> inputParameters) {
        this.inputParameters = inputParameters;
        return this;
    }

    @JsonProperty("start_time")
    public Double getStartTime() {
        return startTime;
    }

    @JsonProperty("start_time")
    public void setStartTime(Double startTime) {
        this.startTime = startTime;
    }

    public SummarizeAndVisualizeResultsResults withStartTime(Double startTime) {
        this.startTime = startTime;
        return this;
    }

    @JsonProperty("output_dir")
    public java.lang.String getOutputDir() {
        return outputDir;
    }

    @JsonProperty("output_dir")
    public void setOutputDir(java.lang.String outputDir) {
        this.outputDir = outputDir;
    }

    public SummarizeAndVisualizeResultsResults withOutputDir(java.lang.String outputDir) {
        this.outputDir = outputDir;
        return this;
    }

    @JsonProperty("summary")
    public java.lang.String getSummary() {
        return summary;
    }

    @JsonProperty("summary")
    public void setSummary(java.lang.String summary) {
        this.summary = summary;
    }

    public SummarizeAndVisualizeResultsResults withSummary(java.lang.String summary) {
        this.summary = summary;
        return this;
    }

    @JsonProperty("html_report_path")
    public java.lang.String getHtmlReportPath() {
        return htmlReportPath;
    }

    @JsonProperty("html_report_path")
    public void setHtmlReportPath(java.lang.String htmlReportPath) {
        this.htmlReportPath = htmlReportPath;
    }

    public SummarizeAndVisualizeResultsResults withHtmlReportPath(java.lang.String htmlReportPath) {
        this.htmlReportPath = htmlReportPath;
        return this;
    }

    @JsonProperty("sequence_analysis_ref")
    public java.lang.String getSequenceAnalysisRef() {
        return sequenceAnalysisRef;
    }

    @JsonProperty("sequence_analysis_ref")
    public void setSequenceAnalysisRef(java.lang.String sequenceAnalysisRef) {
        this.sequenceAnalysisRef = sequenceAnalysisRef;
    }

    public SummarizeAndVisualizeResultsResults withSequenceAnalysisRef(java.lang.String sequenceAnalysisRef) {
        this.sequenceAnalysisRef = sequenceAnalysisRef;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((("SummarizeAndVisualizeResultsResults"+" [reportName=")+ reportName)+", reportRef=")+ reportRef)+", inputParameters=")+ inputParameters)+", startTime=")+ startTime)+", outputDir=")+ outputDir)+", summary=")+ summary)+", htmlReportPath=")+ htmlReportPath)+", sequenceAnalysisRef=")+ sequenceAnalysisRef)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
