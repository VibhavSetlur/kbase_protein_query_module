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
    private String proteinId;
    @JsonProperty("exists")
    private Long exists;
    @JsonProperty("family_id")
    private String familyId;
    @JsonProperty("metadata")
    private Map<String, UObject> metadata;
    @JsonProperty("embedding_ref")
    private String embeddingRef;
    @JsonProperty("embedding")
    private List<Double> embedding;
    @JsonProperty("model_name")
    private String modelName;
    @JsonProperty("search_timestamp")
    private Double searchTimestamp;
    @JsonProperty("summary")
    private String summary;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("protein_id")
    public String getProteinId() {
        return proteinId;
    }

    @JsonProperty("protein_id")
    public void setProteinId(String proteinId) {
        this.proteinId = proteinId;
    }

    public ProteinExistenceResult withProteinId(String proteinId) {
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
    public String getFamilyId() {
        return familyId;
    }

    @JsonProperty("family_id")
    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    public ProteinExistenceResult withFamilyId(String familyId) {
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
    public String getEmbeddingRef() {
        return embeddingRef;
    }

    @JsonProperty("embedding_ref")
    public void setEmbeddingRef(String embeddingRef) {
        this.embeddingRef = embeddingRef;
    }

    public ProteinExistenceResult withEmbeddingRef(String embeddingRef) {
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
    public String getModelName() {
        return modelName;
    }

    @JsonProperty("model_name")
    public void setModelName(String modelName) {
        this.modelName = modelName;
    }

    public ProteinExistenceResult withModelName(String modelName) {
        this.modelName = modelName;
        return this;
    }

    @JsonProperty("search_timestamp")
    public Double getSearchTimestamp() {
        return searchTimestamp;
    }

    @JsonProperty("search_timestamp")
    public void setSearchTimestamp(Double searchTimestamp) {
        this.searchTimestamp = searchTimestamp;
    }

    public ProteinExistenceResult withSearchTimestamp(Double searchTimestamp) {
        this.searchTimestamp = searchTimestamp;
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

    public ProteinExistenceResult withSummary(String summary) {
        this.summary = summary;
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
        return ((((((((((((((((((((("ProteinExistenceResult"+" [proteinId=")+ proteinId)+", exists=")+ exists)+", familyId=")+ familyId)+", metadata=")+ metadata)+", embeddingRef=")+ embeddingRef)+", embedding=")+ embedding)+", modelName=")+ modelName)+", searchTimestamp=")+ searchTimestamp)+", summary=")+ summary)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
