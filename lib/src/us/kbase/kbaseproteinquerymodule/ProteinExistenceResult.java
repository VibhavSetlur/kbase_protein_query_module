
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
 * <p>Original spec-file type: ProteinExistenceResult</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "protein_id",
    "exists",
    "family_id",
    "metadata",
    "embedding_ref",
    "embedding",
    "model_name",
    "search_timestamp",
    "summary"
})
public class ProteinExistenceResult {

    @JsonProperty("protein_id")
    private java.lang.String proteinId;
    @JsonProperty("exists")
    private Long exists;
    @JsonProperty("family_id")
    private java.lang.String familyId;
    @JsonProperty("metadata")
    private Map<String, UObject> metadata;
    @JsonProperty("embedding_ref")
    private java.lang.String embeddingRef;
    @JsonProperty("embedding")
    private List<Double> embedding;
    @JsonProperty("model_name")
    private java.lang.String modelName;
    @JsonProperty("search_timestamp")
    private java.lang.Double searchTimestamp;
    @JsonProperty("summary")
    private java.lang.String summary;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("protein_id")
    public java.lang.String getProteinId() {
        return proteinId;
    }

    @JsonProperty("protein_id")
    public void setProteinId(java.lang.String proteinId) {
        this.proteinId = proteinId;
    }

    public ProteinExistenceResult withProteinId(java.lang.String proteinId) {
        this.proteinId = proteinId;
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

    public ProteinExistenceResult withExists(Long exists) {
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

    public ProteinExistenceResult withFamilyId(java.lang.String familyId) {
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

    public ProteinExistenceResult withMetadata(Map<String, UObject> metadata) {
        this.metadata = metadata;
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

    public ProteinExistenceResult withEmbeddingRef(java.lang.String embeddingRef) {
        this.embeddingRef = embeddingRef;
        return this;
    }

    @JsonProperty("embedding")
    public List<Double> getEmbedding() {
        return embedding;
    }

    @JsonProperty("embedding")
    public void setEmbedding(List<Double> embedding) {
        this.embedding = embedding;
    }

    public ProteinExistenceResult withEmbedding(List<Double> embedding) {
        this.embedding = embedding;
        return this;
    }

    @JsonProperty("model_name")
    public java.lang.String getModelName() {
        return modelName;
    }

    @JsonProperty("model_name")
    public void setModelName(java.lang.String modelName) {
        this.modelName = modelName;
    }

    public ProteinExistenceResult withModelName(java.lang.String modelName) {
        this.modelName = modelName;
        return this;
    }

    @JsonProperty("search_timestamp")
    public java.lang.Double getSearchTimestamp() {
        return searchTimestamp;
    }

    @JsonProperty("search_timestamp")
    public void setSearchTimestamp(java.lang.Double searchTimestamp) {
        this.searchTimestamp = searchTimestamp;
    }

    public ProteinExistenceResult withSearchTimestamp(java.lang.Double searchTimestamp) {
        this.searchTimestamp = searchTimestamp;
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

    public ProteinExistenceResult withSummary(java.lang.String summary) {
        this.summary = summary;
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
        return ((((((((((((((((((((("ProteinExistenceResult"+" [proteinId=")+ proteinId)+", exists=")+ exists)+", familyId=")+ familyId)+", metadata=")+ metadata)+", embeddingRef=")+ embeddingRef)+", embedding=")+ embedding)+", modelName=")+ modelName)+", searchTimestamp=")+ searchTimestamp)+", summary=")+ summary)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
