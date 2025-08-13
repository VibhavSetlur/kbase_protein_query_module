
package us.kbase.kbaseproteinquerymodule;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: ProteinFamilyAssignmentResult</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_id",
    "input_type",
    "embedding_ref",
    "assigned_family_id",
    "similarity_score",
    "metadata",
    "eigenprotein_id",
    "confidence"
})
public class ProteinFamilyAssignmentResult {

    @JsonProperty("input_id")
    private java.lang.String inputId;
    @JsonProperty("input_type")
    private java.lang.String inputType;
    @JsonProperty("embedding_ref")
    private java.lang.String embeddingRef;
    @JsonProperty("assigned_family_id")
    private java.lang.String assignedFamilyId;
    @JsonProperty("similarity_score")
    private Double similarityScore;
    @JsonProperty("metadata")
    private Map<String, String> metadata;
    @JsonProperty("eigenprotein_id")
    private java.lang.String eigenproteinId;
    @JsonProperty("confidence")
    private Double confidence;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("input_id")
    public java.lang.String getInputId() {
        return inputId;
    }

    @JsonProperty("input_id")
    public void setInputId(java.lang.String inputId) {
        this.inputId = inputId;
    }

    public ProteinFamilyAssignmentResult withInputId(java.lang.String inputId) {
        this.inputId = inputId;
        return this;
    }

    @JsonProperty("input_type")
    public java.lang.String getInputType() {
        return inputType;
    }

    @JsonProperty("input_type")
    public void setInputType(java.lang.String inputType) {
        this.inputType = inputType;
    }

    public ProteinFamilyAssignmentResult withInputType(java.lang.String inputType) {
        this.inputType = inputType;
        return this;
    }

    @JsonProperty("embedding_ref")
    public java.lang.String getEmbeddingRef() {
        return embeddingRef;
    }

    @JsonProperty("embedding_ref")
    public void setEmbeddingRef(java.lang.String embeddingRef) {
        this.embeddingRef = embeddingRef;
    }

    public ProteinFamilyAssignmentResult withEmbeddingRef(java.lang.String embeddingRef) {
        this.embeddingRef = embeddingRef;
        return this;
    }

    @JsonProperty("assigned_family_id")
    public java.lang.String getAssignedFamilyId() {
        return assignedFamilyId;
    }

    @JsonProperty("assigned_family_id")
    public void setAssignedFamilyId(java.lang.String assignedFamilyId) {
        this.assignedFamilyId = assignedFamilyId;
    }

    public ProteinFamilyAssignmentResult withAssignedFamilyId(java.lang.String assignedFamilyId) {
        this.assignedFamilyId = assignedFamilyId;
        return this;
    }

    @JsonProperty("similarity_score")
    public Double getSimilarityScore() {
        return similarityScore;
    }

    @JsonProperty("similarity_score")
    public void setSimilarityScore(Double similarityScore) {
        this.similarityScore = similarityScore;
    }

    public ProteinFamilyAssignmentResult withSimilarityScore(Double similarityScore) {
        this.similarityScore = similarityScore;
        return this;
    }

    @JsonProperty("metadata")
    public Map<String, String> getMetadata() {
        return metadata;
    }

    @JsonProperty("metadata")
    public void setMetadata(Map<String, String> metadata) {
        this.metadata = metadata;
    }

    public ProteinFamilyAssignmentResult withMetadata(Map<String, String> metadata) {
        this.metadata = metadata;
        return this;
    }

    @JsonProperty("eigenprotein_id")
    public java.lang.String getEigenproteinId() {
        return eigenproteinId;
    }

    @JsonProperty("eigenprotein_id")
    public void setEigenproteinId(java.lang.String eigenproteinId) {
        this.eigenproteinId = eigenproteinId;
    }

    public ProteinFamilyAssignmentResult withEigenproteinId(java.lang.String eigenproteinId) {
        this.eigenproteinId = eigenproteinId;
        return this;
    }

    @JsonProperty("confidence")
    public Double getConfidence() {
        return confidence;
    }

    @JsonProperty("confidence")
    public void setConfidence(Double confidence) {
        this.confidence = confidence;
    }

    public ProteinFamilyAssignmentResult withConfidence(Double confidence) {
        this.confidence = confidence;
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
        return ((((((((((((((((((("ProteinFamilyAssignmentResult"+" [inputId=")+ inputId)+", inputType=")+ inputType)+", embeddingRef=")+ embeddingRef)+", assignedFamilyId=")+ assignedFamilyId)+", similarityScore=")+ similarityScore)+", metadata=")+ metadata)+", eigenproteinId=")+ eigenproteinId)+", confidence=")+ confidence)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
