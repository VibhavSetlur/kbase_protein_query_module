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
 * <p>Original spec-file type: FindTopMatchesFromEmbeddingResults</p>
 * <pre>
 * Find top matches for a given protein embedding within a family.
 * Uses FAISS IVF float index for efficient similarity search.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "matches",
    "summary",
    "input_parameters",
    "start_time",
    "family_id",
    "top_n",
    "similarity_stats",
    "similarity_search_result_ref"
})
public class FindTopMatchesFromEmbeddingResults {

    @JsonProperty("matches")
    private List<Map<String, UObject>> matches;
    @JsonProperty("summary")
    private String summary;
    @JsonProperty("input_parameters")
    private Map<String, UObject> inputParameters;
    @JsonProperty("start_time")
    private Double startTime;
    @JsonProperty("family_id")
    private String familyId;
    @JsonProperty("top_n")
    private Long topN;
    @JsonProperty("similarity_stats")
    private Map<String, Double> similarityStats;
    @JsonProperty("similarity_search_result_ref")
    private String similaritySearchResultRef;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("matches")
    public List<Map<String, UObject>> getMatches() {
        return matches;
    }

    @JsonProperty("matches")
    public void setMatches(List<Map<String, UObject>> matches) {
        this.matches = matches;
    }

    public FindTopMatchesFromEmbeddingResults withMatches(List<Map<String, UObject>> matches) {
        this.matches = matches;
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

    public FindTopMatchesFromEmbeddingResults withSummary(String summary) {
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

    public FindTopMatchesFromEmbeddingResults withInputParameters(Map<String, UObject> inputParameters) {
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

    public FindTopMatchesFromEmbeddingResults withStartTime(Double startTime) {
        this.startTime = startTime;
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

    public FindTopMatchesFromEmbeddingResults withFamilyId(String familyId) {
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

    public FindTopMatchesFromEmbeddingResults withTopN(Long topN) {
        this.topN = topN;
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

    public FindTopMatchesFromEmbeddingResults withSimilarityStats(Map<String, Double> similarityStats) {
        this.similarityStats = similarityStats;
        return this;
    }

    @JsonProperty("similarity_search_result_ref")
    public String getSimilaritySearchResultRef() {
        return similaritySearchResultRef;
    }

    @JsonProperty("similarity_search_result_ref")
    public void setSimilaritySearchResultRef(String similaritySearchResultRef) {
        this.similaritySearchResultRef = similaritySearchResultRef;
    }

    public FindTopMatchesFromEmbeddingResults withSimilaritySearchResultRef(String similaritySearchResultRef) {
        this.similaritySearchResultRef = similaritySearchResultRef;
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
        return ((((((((((((((((((("FindTopMatchesFromEmbeddingResults"+" [matches=")+ matches)+", summary=")+ summary)+", inputParameters=")+ inputParameters)+", startTime=")+ startTime)+", familyId=")+ familyId)+", topN=")+ topN)+", similarityStats=")+ similarityStats)+", similaritySearchResultRef=")+ similaritySearchResultRef)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
