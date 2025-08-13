
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
 * <p>Original spec-file type: SummarizeVisualizeResult</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_id",
    "input_type",
    "top_matches_result_ref",
    "summary_html",
    "metadata",
    "embedding_ref",
    "family_assignment_ref",
    "protein_existence_ref",
    "matches",
    "family_id",
    "top_n"
})
public class SummarizeVisualizeResult {

    @JsonProperty("input_id")
    private java.lang.String inputId;
    @JsonProperty("input_type")
    private java.lang.String inputType;
    @JsonProperty("top_matches_result_ref")
    private java.lang.String topMatchesResultRef;
    @JsonProperty("summary_html")
    private java.lang.String summaryHtml;
    @JsonProperty("metadata")
    private Map<String, String> metadata;
    @JsonProperty("embedding_ref")
    private java.lang.String embeddingRef;
    @JsonProperty("family_assignment_ref")
    private java.lang.String familyAssignmentRef;
    @JsonProperty("protein_existence_ref")
    private java.lang.String proteinExistenceRef;
    @JsonProperty("matches")
    private List<Map<String, UObject>> matches;
    @JsonProperty("family_id")
    private java.lang.String familyId;
    @JsonProperty("top_n")
    private Long topN;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("input_id")
    public java.lang.String getInputId() {
        return inputId;
    }

    @JsonProperty("input_id")
    public void setInputId(java.lang.String inputId) {
        this.inputId = inputId;
    }

    public SummarizeVisualizeResult withInputId(java.lang.String inputId) {
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

    public SummarizeVisualizeResult withInputType(java.lang.String inputType) {
        this.inputType = inputType;
        return this;
    }

    @JsonProperty("top_matches_result_ref")
    public java.lang.String getTopMatchesResultRef() {
        return topMatchesResultRef;
    }

    @JsonProperty("top_matches_result_ref")
    public void setTopMatchesResultRef(java.lang.String topMatchesResultRef) {
        this.topMatchesResultRef = topMatchesResultRef;
    }

    public SummarizeVisualizeResult withTopMatchesResultRef(java.lang.String topMatchesResultRef) {
        this.topMatchesResultRef = topMatchesResultRef;
        return this;
    }

    @JsonProperty("summary_html")
    public java.lang.String getSummaryHtml() {
        return summaryHtml;
    }

    @JsonProperty("summary_html")
    public void setSummaryHtml(java.lang.String summaryHtml) {
        this.summaryHtml = summaryHtml;
    }

    public SummarizeVisualizeResult withSummaryHtml(java.lang.String summaryHtml) {
        this.summaryHtml = summaryHtml;
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

    public SummarizeVisualizeResult withMetadata(Map<String, String> metadata) {
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

    public SummarizeVisualizeResult withEmbeddingRef(java.lang.String embeddingRef) {
        this.embeddingRef = embeddingRef;
        return this;
    }

    @JsonProperty("family_assignment_ref")
    public java.lang.String getFamilyAssignmentRef() {
        return familyAssignmentRef;
    }

    @JsonProperty("family_assignment_ref")
    public void setFamilyAssignmentRef(java.lang.String familyAssignmentRef) {
        this.familyAssignmentRef = familyAssignmentRef;
    }

    public SummarizeVisualizeResult withFamilyAssignmentRef(java.lang.String familyAssignmentRef) {
        this.familyAssignmentRef = familyAssignmentRef;
        return this;
    }

    @JsonProperty("protein_existence_ref")
    public java.lang.String getProteinExistenceRef() {
        return proteinExistenceRef;
    }

    @JsonProperty("protein_existence_ref")
    public void setProteinExistenceRef(java.lang.String proteinExistenceRef) {
        this.proteinExistenceRef = proteinExistenceRef;
    }

    public SummarizeVisualizeResult withProteinExistenceRef(java.lang.String proteinExistenceRef) {
        this.proteinExistenceRef = proteinExistenceRef;
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

    public SummarizeVisualizeResult withMatches(List<Map<String, UObject>> matches) {
        this.matches = matches;
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

    public SummarizeVisualizeResult withFamilyId(java.lang.String familyId) {
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

    public SummarizeVisualizeResult withTopN(Long topN) {
        this.topN = topN;
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
        return ((((((((((((((((((((((((("SummarizeVisualizeResult"+" [inputId=")+ inputId)+", inputType=")+ inputType)+", topMatchesResultRef=")+ topMatchesResultRef)+", summaryHtml=")+ summaryHtml)+", metadata=")+ metadata)+", embeddingRef=")+ embeddingRef)+", familyAssignmentRef=")+ familyAssignmentRef)+", proteinExistenceRef=")+ proteinExistenceRef)+", matches=")+ matches)+", familyId=")+ familyId)+", topN=")+ topN)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
