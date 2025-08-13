
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
 * <p>Original spec-file type: ProteinSimilaritySearchResult</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_id",
    "input_type",
    "embedding_ref",
    "family_id",
    "top_n",
    "matches",
    "similarity_stats",
    "metadata"
})
public class ProteinSimilaritySearchResult {

    @JsonProperty("input_id")
    private java.lang.String inputId;
    @JsonProperty("input_type")
    private java.lang.String inputType;
    @JsonProperty("embedding_ref")
    private java.lang.String embeddingRef;
    @JsonProperty("family_id")
    private java.lang.String familyId;
    @JsonProperty("top_n")
    private Long topN;
    @JsonProperty("matches")
    private List<Map<String, UObject>> matches;
    @JsonProperty("similarity_stats")
    private Map<String, Double> similarityStats;
    @JsonProperty("metadata")
    private Map<String, String> metadata;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("input_id")
    public java.lang.String getInputId() {
        return inputId;
    }

    @JsonProperty("input_id")
    public void setInputId(java.lang.String inputId) {
        this.inputId = inputId;
    }

    public ProteinSimilaritySearchResult withInputId(java.lang.String inputId) {
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

    public ProteinSimilaritySearchResult withInputType(java.lang.String inputType) {
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

    public ProteinSimilaritySearchResult withEmbeddingRef(java.lang.String embeddingRef) {
        this.embeddingRef = embeddingRef;
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

    public ProteinSimilaritySearchResult withFamilyId(java.lang.String familyId) {
        this.familyId = familyId;
        return this;
    }

    @JsonProperty("top_n")
    public Long getTopN() {
        return topN;
    }

    @JsonProperty("top_n")
    public void setTopN(Long topN) {
        this.topN = topN;
    }

    public ProteinSimilaritySearchResult withTopN(Long topN) {
        this.topN = topN;
        return this;
    }

    @JsonProperty("matches")
    public List<Map<String, UObject>> getMatches() {
        return matches;
    }

    @JsonProperty("matches")
    public void setMatches(List<Map<String, UObject>> matches) {
        this.matches = matches;
    }

    public ProteinSimilaritySearchResult withMatches(List<Map<String, UObject>> matches) {
        this.matches = matches;
        return this;
    }

    @JsonProperty("similarity_stats")
    public Map<String, Double> getSimilarityStats() {
        return similarityStats;
    }

    @JsonProperty("similarity_stats")
    public void setSimilarityStats(Map<String, Double> similarityStats) {
        this.similarityStats = similarityStats;
    }

    public ProteinSimilaritySearchResult withSimilarityStats(Map<String, Double> similarityStats) {
        this.similarityStats = similarityStats;
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

    public ProteinSimilaritySearchResult withMetadata(Map<String, String> metadata) {
        this.metadata = metadata;
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
        return ((((((((((((((((((("ProteinSimilaritySearchResult"+" [inputId=")+ inputId)+", inputType=")+ inputType)+", embeddingRef=")+ embeddingRef)+", familyId=")+ familyId)+", topN=")+ topN)+", matches=")+ matches)+", similarityStats=")+ similarityStats)+", metadata=")+ metadata)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
