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
    private String reportName;
    @JsonProperty("report_ref")
    private String reportRef;
    @JsonProperty("exists")
    private Long exists;
    @JsonProperty("family_id")
    private String familyId;
    @JsonProperty("metadata")
    private Map<String, UObject> metadata;
    @JsonProperty("input_parameters")
    private Map<String, UObject> inputParameters;
    @JsonProperty("start_time")
    private Double startTime;
    @JsonProperty("summary")
    private String summary;
    @JsonProperty("protein_existence_result_ref")
    private String proteinExistenceResultRef;
    @JsonProperty("embedding_result_ref")
    private String embeddingResultRef;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("report_name")
    public String getReportName() {
        return reportName;
    }

    @JsonProperty("report_name")
    public void setReportName(String reportName) {
        this.reportName = reportName;
    }

    public CheckProteinExistenceResults withReportName(String reportName) {
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

    public CheckProteinExistenceResults withReportRef(String reportRef) {
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
    public String getFamilyId() {
        return familyId;
    }

    @JsonProperty("family_id")
    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    public CheckProteinExistenceResults withFamilyId(String familyId) {
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
    public String getSummary() {
        return summary;
    }

    @JsonProperty("summary")
    public void setSummary(String summary) {
        this.summary = summary;
    }

    public CheckProteinExistenceResults withSummary(String summary) {
        this.summary = summary;
        return this;
    }

    @JsonProperty("protein_existence_result_ref")
    public String getProteinExistenceResultRef() {
        return proteinExistenceResultRef;
    }

    @JsonProperty("protein_existence_result_ref")
    public void setProteinExistenceResultRef(String proteinExistenceResultRef) {
        this.proteinExistenceResultRef = proteinExistenceResultRef;
    }

    public CheckProteinExistenceResults withProteinExistenceResultRef(String proteinExistenceResultRef) {
        this.proteinExistenceResultRef = proteinExistenceResultRef;
        return this;
    }

    @JsonProperty("embedding_result_ref")
    public String getEmbeddingResultRef() {
        return embeddingResultRef;
    }

    @JsonProperty("embedding_result_ref")
    public void setEmbeddingResultRef(String embeddingResultRef) {
        this.embeddingResultRef = embeddingResultRef;
    }

    public CheckProteinExistenceResults withEmbeddingResultRef(String embeddingResultRef) {
        this.embeddingResultRef = embeddingResultRef;
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
        return ((((((((((((((((((((((("CheckProteinExistenceResults"+" [reportName=")+ reportName)+", reportRef=")+ reportRef)+", exists=")+ exists)+", familyId=")+ familyId)+", metadata=")+ metadata)+", inputParameters=")+ inputParameters)+", startTime=")+ startTime)+", summary=")+ summary)+", proteinExistenceResultRef=")+ proteinExistenceResultRef)+", embeddingResultRef=")+ embeddingResultRef)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
