package us.kbase.kbaseproteinquerymodule;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import us.kbase.common.service.UObject;

/**
 * <p>Original spec-file type: ProteinQueryAnalysisResults</p>
 * <pre>
 * Unified Protein Query Analysis Pipeline
 * 
 * This method provides a single entry point for comprehensive protein analysis,
 * supporting multiple input types and configurable analysis stages.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "report_name",
    "report_ref",
    "analysis_result_ref",
    "summary",
    "input_parameters",
    "start_time",
    "html_report_path",
    "protein_count",
    "stages_completed"
})
public class ProteinQueryAnalysisResults {

    @JsonProperty("report_name")
    private String reportName;
    @JsonProperty("report_ref")
    private String reportRef;
    @JsonProperty("analysis_result_ref")
    private String analysisResultRef;
    @JsonProperty("summary")
    private String summary;
    @JsonProperty("input_parameters")
    private Map<String, UObject> inputParameters;
    @JsonProperty("start_time")
    private Double startTime;
    @JsonProperty("html_report_path")
    private String htmlReportPath;
    @JsonProperty("protein_count")
    private Long proteinCount;
    @JsonProperty("stages_completed")
    private List<String> stagesCompleted;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("report_name")
    public String getReportName() {
        return reportName;
    }

    @JsonProperty("report_name")
    public void setReportName(String reportName) {
        this.reportName = reportName;
    }

    public ProteinQueryAnalysisResults withReportName(String reportName) {
        this.reportName = reportName;
        return this;
    }

    @JsonProperty("report_ref")
    public String getReportRef() {
        return reportRef;
    }

    @JsonProperty("report_ref")
    public void setReportRef(String reportRef) {
        this.reportRef = reportRef;
    }

    public ProteinQueryAnalysisResults withReportRef(String reportRef) {
        this.reportRef = reportRef;
        return this;
    }

    @JsonProperty("analysis_result_ref")
    public String getAnalysisResultRef() {
        return analysisResultRef;
    }

    @JsonProperty("analysis_result_ref")
    public void setAnalysisResultRef(String analysisResultRef) {
        this.analysisResultRef = analysisResultRef;
    }

    public ProteinQueryAnalysisResults withAnalysisResultRef(String analysisResultRef) {
        this.analysisResultRef = analysisResultRef;
        return this;
    }

    @JsonProperty("summary")
    public String getSummary() {
        return summary;
    }

    @JsonProperty("summary")
    public void setSummary(String summary) {
        this.summary = summary;
    }

    public ProteinQueryAnalysisResults withSummary(String summary) {
        this.summary = summary;
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

    public ProteinQueryAnalysisResults withInputParameters(Map<String, UObject> inputParameters) {
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

    public ProteinQueryAnalysisResults withStartTime(Double startTime) {
        this.startTime = startTime;
        return this;
    }

    @JsonProperty("html_report_path")
    public String getHtmlReportPath() {
        return htmlReportPath;
    }

    @JsonProperty("html_report_path")
    public void setHtmlReportPath(String htmlReportPath) {
        this.htmlReportPath = htmlReportPath;
    }

    public ProteinQueryAnalysisResults withHtmlReportPath(String htmlReportPath) {
        this.htmlReportPath = htmlReportPath;
        return this;
    }

    @JsonProperty("protein_count")
    public Long getProteinCount() {
        return proteinCount;
    }

    @JsonProperty("protein_count")
    public void setProteinCount(Long proteinCount) {
        this.proteinCount = proteinCount;
    }

    public ProteinQueryAnalysisResults withProteinCount(Long proteinCount) {
        this.proteinCount = proteinCount;
        return this;
    }

    @JsonProperty("stages_completed")
    public List<String> getStagesCompleted() {
        return stagesCompleted;
    }

    @JsonProperty("stages_completed")
    public void setStagesCompleted(List<String> stagesCompleted) {
        this.stagesCompleted = stagesCompleted;
    }

    public ProteinQueryAnalysisResults withStagesCompleted(List<String> stagesCompleted) {
        this.stagesCompleted = stagesCompleted;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((("ProteinQueryAnalysisResults"+" [reportName=")+ reportName)+", reportRef=")+ reportRef)+", analysisResultRef=")+ analysisResultRef)+", summary=")+ summary)+", inputParameters=")+ inputParameters)+", startTime=")+ startTime)+", htmlReportPath=")+ htmlReportPath)+", proteinCount=")+ proteinCount)+", stagesCompleted=")+ stagesCompleted)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
