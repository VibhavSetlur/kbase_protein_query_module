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
    private String inputId;
    @JsonProperty("input_type")
    private String inputType;
    @JsonProperty("embedding_ref")
    private String embeddingRef;
    @JsonProperty("assigned_family_id")
    private String assignedFamilyId;
    @JsonProperty("similarity_score")
    private Double similarityScore;
    @JsonProperty("metadata")
    private Map<String, String> metadata;
    @JsonProperty("eigenprotein_id")
    private String eigenproteinId;
    @JsonProperty("confidence")
    private Double confidence;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("input_id")
    public String getInputId() {
        return inputId;
    }

    @JsonProperty("input_id")
    public void setInputId(String inputId) {
        this.inputId = inputId;
    }

    public ProteinFamilyAssignmentResult withInputId(String inputId) {
        this.inputId = inputId;
        return this;
    }

    @JsonProperty("input_type")
    public String getInputType() {
        return inputType;
    }

    @JsonProperty("input_type")
    public void setInputType(String inputType) {
        this.inputType = inputType;
    }

    public ProteinFamilyAssignmentResult withInputType(String inputType) {
        this.inputType = inputType;
        return this;
    }

    @JsonProperty("embedding_ref")
    public String getEmbeddingRef() {
        return embeddingRef;
    }

    @JsonProperty("embedding_ref")
    public void setEmbeddingRef(String embeddingRef) {
        this.embeddingRef = embeddingRef;
    }

    public ProteinFamilyAssignmentResult withEmbeddingRef(String embeddingRef) {
        this.embeddingRef = embeddingRef;
        return this;
    }

    @JsonProperty("assigned_family_id")
    public String getAssignedFamilyId() {
        return assignedFamilyId;
    }

    @JsonProperty("assigned_family_id")
    public void setAssignedFamilyId(String assignedFamilyId) {
        this.assignedFamilyId = assignedFamilyId;
    }

    public ProteinFamilyAssignmentResult withAssignedFamilyId(String assignedFamilyId) {
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
    public String getEigenproteinId() {
        return eigenproteinId;
    }

    @JsonProperty("eigenprotein_id")
    public void setEigenproteinId(String eigenproteinId) {
        this.eigenproteinId = eigenproteinId;
    }

    public ProteinFamilyAssignmentResult withEigenproteinId(String eigenproteinId) {
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
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((("ProteinFamilyAssignmentResult"+" [inputId=")+ inputId)+", inputType=")+ inputType)+", embeddingRef=")+ embeddingRef)+", assignedFamilyId=")+ assignedFamilyId)+", similarityScore=")+ similarityScore)+", metadata=")+ metadata)+", eigenproteinId=")+ eigenproteinId)+", confidence=")+ confidence)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
