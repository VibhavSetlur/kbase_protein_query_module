
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
 * <p>Original spec-file type: GenerateProteinEmbeddingResults</p>
 * <pre>
 * Generate protein embeddings from direct sequence input.
 * Creates embeddings using ESM-2 model for downstream analysis.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "report_name",
    "report_ref",
    "embedding_result_ref",
    "summary",
    "input_parameters",
    "start_time",
    "embedding_norm",
    "sequence_length",
    "embedding_dim"
})
public class GenerateProteinEmbeddingResults {

    @JsonProperty("report_name")
    private java.lang.String reportName;
    @JsonProperty("report_ref")
    private java.lang.String reportRef;
    @JsonProperty("embedding_result_ref")
    private java.lang.String embeddingResultRef;
    @JsonProperty("summary")
    private java.lang.String summary;
    @JsonProperty("input_parameters")
    private Map<String, UObject> inputParameters;
    @JsonProperty("start_time")
    private Double startTime;
    @JsonProperty("embedding_norm")
    private Double embeddingNorm;
    @JsonProperty("sequence_length")
    private Long sequenceLength;
    @JsonProperty("embedding_dim")
    private Long embeddingDim;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("report_name")
    public java.lang.String getReportName() {
        return reportName;
    }

    @JsonProperty("report_name")
    public void setReportName(java.lang.String reportName) {
        this.reportName = reportName;
    }

    public GenerateProteinEmbeddingResults withReportName(java.lang.String reportName) {
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

    public GenerateProteinEmbeddingResults withReportRef(java.lang.String reportRef) {
        this.reportRef = reportRef;
        return this;
    }

    @JsonProperty("embedding_result_ref")
    public java.lang.String getEmbeddingResultRef() {
        return embeddingResultRef;
    }

    @JsonProperty("embedding_result_ref")
    public void setEmbeddingResultRef(java.lang.String embeddingResultRef) {
        this.embeddingResultRef = embeddingResultRef;
    }

    public GenerateProteinEmbeddingResults withEmbeddingResultRef(java.lang.String embeddingResultRef) {
        this.embeddingResultRef = embeddingResultRef;
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

    public GenerateProteinEmbeddingResults withSummary(java.lang.String summary) {
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

    public GenerateProteinEmbeddingResults withInputParameters(Map<String, UObject> inputParameters) {
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

    public GenerateProteinEmbeddingResults withStartTime(Double startTime) {
        this.startTime = startTime;
        return this;
    }

    @JsonProperty("embedding_norm")
    public Double getEmbeddingNorm() {
        return embeddingNorm;
    }

    @JsonProperty("embedding_norm")
    public void setEmbeddingNorm(Double embeddingNorm) {
        this.embeddingNorm = embeddingNorm;
    }

    public GenerateProteinEmbeddingResults withEmbeddingNorm(Double embeddingNorm) {
        this.embeddingNorm = embeddingNorm;
        return this;
    }

    @JsonProperty("sequence_length")
    public Long getSequenceLength() {
        return sequenceLength;
    }

    @JsonProperty("sequence_length")
    public void setSequenceLength(Long sequenceLength) {
        this.sequenceLength = sequenceLength;
    }

    public GenerateProteinEmbeddingResults withSequenceLength(Long sequenceLength) {
        this.sequenceLength = sequenceLength;
        return this;
    }

    @JsonProperty("embedding_dim")
    public Long getEmbeddingDim() {
        return embeddingDim;
    }

    @JsonProperty("embedding_dim")
    public void setEmbeddingDim(Long embeddingDim) {
        this.embeddingDim = embeddingDim;
    }

    public GenerateProteinEmbeddingResults withEmbeddingDim(Long embeddingDim) {
        this.embeddingDim = embeddingDim;
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
        return ((((((((((((((((((((("GenerateProteinEmbeddingResults"+" [reportName=")+ reportName)+", reportRef=")+ reportRef)+", embeddingResultRef=")+ embeddingResultRef)+", summary=")+ summary)+", inputParameters=")+ inputParameters)+", startTime=")+ startTime)+", embeddingNorm=")+ embeddingNorm)+", sequenceLength=")+ sequenceLength)+", embeddingDim=")+ embeddingDim)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
