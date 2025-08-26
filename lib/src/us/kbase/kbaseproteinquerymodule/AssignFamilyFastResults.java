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
 * <p>Original spec-file type: AssignFamilyFastResults</p>
 * <pre>
 * Assign a protein embedding to a family using similarity to family centroids.
 * Uses binary Hamming distance for fast family assignment.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "family_id",
    "confidence",
    "eigenprotein_id",
    "input_parameters",
    "start_time",
    "family_assignment_result_ref"
})
public class AssignFamilyFastResults {

    @JsonProperty("family_id")
    private String familyId;
    @JsonProperty("confidence")
    private Double confidence;
    @JsonProperty("eigenprotein_id")
    private String eigenproteinId;
    @JsonProperty("input_parameters")
    private Map<String, UObject> inputParameters;
    @JsonProperty("start_time")
    private Double startTime;
    @JsonProperty("family_assignment_result_ref")
    private String familyAssignmentResultRef;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("family_id")
    public String getFamilyId() {
        return familyId;
    }

    @JsonProperty("family_id")
    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    public AssignFamilyFastResults withFamilyId(String familyId) {
        this.familyId = familyId;
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

    public AssignFamilyFastResults withConfidence(Double confidence) {
        this.confidence = confidence;
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

    public AssignFamilyFastResults withEigenproteinId(String eigenproteinId) {
        this.eigenproteinId = eigenproteinId;
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

    public AssignFamilyFastResults withInputParameters(Map<String, UObject> inputParameters) {
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

    public AssignFamilyFastResults withStartTime(Double startTime) {
        this.startTime = startTime;
        return this;
    }

    @JsonProperty("family_assignment_result_ref")
    public String getFamilyAssignmentResultRef() {
        return familyAssignmentResultRef;
    }

    @JsonProperty("family_assignment_result_ref")
    public void setFamilyAssignmentResultRef(String familyAssignmentResultRef) {
        this.familyAssignmentResultRef = familyAssignmentResultRef;
    }

    public AssignFamilyFastResults withFamilyAssignmentResultRef(String familyAssignmentResultRef) {
        this.familyAssignmentResultRef = familyAssignmentResultRef;
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
        return ((((((((((((((((("AssignFamilyFastResults"+" [familyId=")+ familyId)+", confidence=")+ confidence)+", eigenproteinId=")+ eigenproteinId)+", inputParameters=")+ inputParameters)+", startTime=")+ startTime)+", familyAssignmentResultRef=")+ familyAssignmentResultRef)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
