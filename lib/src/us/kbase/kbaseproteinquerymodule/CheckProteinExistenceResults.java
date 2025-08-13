
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
 * <p>Original spec-file type: CheckProteinExistenceResults</p>
 * <pre>
 * Check if a protein exists in the storage system using UniProt ID and create a workspace object with the result.
 * Input: UniProt ID (e.g., P00001, P12345)
 * Output: Existence status, family assignment, metadata, optional embedding
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "report_name",
    "report_ref",
    "exists",
    "family_id",
    "metadata",
    "input_parameters",
    "start_time",
    "summary",
    "protein_existence_result_ref",
    "embedding_result_ref"
})
public class CheckProteinExistenceResults {

    @JsonProperty("report_name")
    private java.lang.String reportName;
    @JsonProperty("report_ref")
    private java.lang.String reportRef;
    @JsonProperty("exists")
    private Long exists;
    @JsonProperty("family_id")
    private java.lang.String familyId;
    @JsonProperty("metadata")
    private Map<String, UObject> metadata;
    @JsonProperty("input_parameters")
    private Map<String, UObject> inputParameters;
    @JsonProperty("start_time")
    private Double startTime;
    @JsonProperty("summary")
    private java.lang.String summary;
    @JsonProperty("protein_existence_result_ref")
    private java.lang.String proteinExistenceResultRef;
    @JsonProperty("embedding_result_ref")
    private java.lang.String embeddingResultRef;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("report_name")
    public java.lang.String getReportName() {
        return reportName;
    }

    @JsonProperty("report_name")
    public void setReportName(java.lang.String reportName) {
        this.reportName = reportName;
    }

    public CheckProteinExistenceResults withReportName(java.lang.String reportName) {
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

    public CheckProteinExistenceResults withReportRef(java.lang.String reportRef) {
        this.reportRef = reportRef;
        return this;
    }

    @JsonProperty("exists")
    public Long getExists() {
        return exists;
    }

    @JsonProperty("exists")
    public void setExists(Long exists) {
        this.exists = exists;
    }

    public CheckProteinExistenceResults withExists(Long exists) {
        this.exists = exists;
        return this;
    }

    @JsonProperty("family_id")
    public java.lang.String getFamilyId() {
        return familyId;
    }

    @JsonProperty("family_id")
    public void setFamilyId(java.lang.String familyId) {
        this.familyId = familyId;
    }

    public CheckProteinExistenceResults withFamilyId(java.lang.String familyId) {
        this.familyId = familyId;
        return this;
    }

    @JsonProperty("metadata")
    public Map<String, UObject> getMetadata() {
        return metadata;
    }

    @JsonProperty("metadata")
    public void setMetadata(Map<String, UObject> metadata) {
        this.metadata = metadata;
    }

    public CheckProteinExistenceResults withMetadata(Map<String, UObject> metadata) {
        this.metadata = metadata;
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

    public CheckProteinExistenceResults withInputParameters(Map<String, UObject> inputParameters) {
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

    public CheckProteinExistenceResults withStartTime(Double startTime) {
        this.startTime = startTime;
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

    public CheckProteinExistenceResults withSummary(java.lang.String summary) {
        this.summary = summary;
        return this;
    }

    @JsonProperty("protein_existence_result_ref")
    public java.lang.String getProteinExistenceResultRef() {
        return proteinExistenceResultRef;
    }

    @JsonProperty("protein_existence_result_ref")
    public void setProteinExistenceResultRef(java.lang.String proteinExistenceResultRef) {
        this.proteinExistenceResultRef = proteinExistenceResultRef;
    }

    public CheckProteinExistenceResults withProteinExistenceResultRef(java.lang.String proteinExistenceResultRef) {
        this.proteinExistenceResultRef = proteinExistenceResultRef;
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

    public CheckProteinExistenceResults withEmbeddingResultRef(java.lang.String embeddingResultRef) {
        this.embeddingResultRef = embeddingResultRef;
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
        return ((((((((((((((((((((((("CheckProteinExistenceResults"+" [reportName=")+ reportName)+", reportRef=")+ reportRef)+", exists=")+ exists)+", familyId=")+ familyId)+", metadata=")+ metadata)+", inputParameters=")+ inputParameters)+", startTime=")+ startTime)+", summary=")+ summary)+", proteinExistenceResultRef=")+ proteinExistenceResultRef)+", embeddingResultRef=")+ embeddingResultRef)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
