package us.kbase.kbaseproteinquerymodule;

import java.util.List;
import java.util.Map;

/**
 * <p>Original spec-file type: AssignFamilyFastResults</p>
 * <pre>
 * Quickly assign a protein embedding to a family by similarity to the medoid.
 * AssignFamilyFastResults is a reference to a hash where the following keys are defined:
 * report_name has a value which is a string
 * report_ref has a value which is a string
 * family_id has a value which is a string
 * similarity_score has a value which is a float
 * summary has a value which is a string
 * input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
 * start_time has a value which is a float
 * family_assignment_result_ref has a value which is a string
 * </pre>
 */
public class AssignFamilyFastResults {
    private String reportName;
    private String reportRef;
    private String familyId;
    private Double similarityScore;
    private String summary;
    private Map<String, Object> inputParameters;
    private Double startTime;
    private String familyAssignmentResultRef;

    public AssignFamilyFastResults() {
    }

    public String getReportName() {
        return reportName;
    }

    public void setReportName(String reportName) {
        this.reportName = reportName;
    }

    public String getReportRef() {
        return reportRef;
    }

    public void setReportRef(String reportRef) {
        this.reportRef = reportRef;
    }

    public String getFamilyId() {
        return familyId;
    }

    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    public Double getSimilarityScore() {
        return similarityScore;
    }

    public void setSimilarityScore(Double similarityScore) {
        this.similarityScore = similarityScore;
    }

    public String getSummary() {
        return summary;
    }

    public void setSummary(String summary) {
        this.summary = summary;
    }

    public Map<String, Object> getInputParameters() {
        return inputParameters;
    }

    public void setInputParameters(Map<String, Object> inputParameters) {
        this.inputParameters = inputParameters;
    }

    public Double getStartTime() {
        return startTime;
    }

    public void setStartTime(Double startTime) {
        this.startTime = startTime;
    }

    public String getFamilyAssignmentResultRef() {
        return familyAssignmentResultRef;
    }

    public void setFamilyAssignmentResultRef(String familyAssignmentResultRef) {
        this.familyAssignmentResultRef = familyAssignmentResultRef;
    }
} 