package us.kbase.kbaseproteinquerymodule;

import java.util.List;
import java.util.Map;

/**
 * <p>Original spec-file type: CheckProteinExistenceResults</p>
 * <pre>
 * Check if a protein exists in the storage system and create a workspace object with the result.
 * CheckProteinExistenceResults is a reference to a hash where the following keys are defined:
 * report_name has a value which is a string
 * report_ref has a value which is a string
 * exists has a value which is an int
 * family_id has a value which is a string
 * metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
 * input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
 * start_time has a value which is a float
 * summary has a value which is a string
 * protein_existence_result_ref has a value which is a string
 * </pre>
 */
public class CheckProteinExistenceResults {
    private String reportName;
    private String reportRef;
    private Long exists;
    private String familyId;
    private Map<String, Object> metadata;
    private Map<String, Object> inputParameters;
    private Double startTime;
    private String summary;
    private String proteinExistenceResultRef;

    public CheckProteinExistenceResults() {
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

    public Long getExists() {
        return exists;
    }

    public void setExists(Long exists) {
        this.exists = exists;
    }

    public String getFamilyId() {
        return familyId;
    }

    public void setFamilyId(String familyId) {
        this.familyId = familyId;
    }

    public Map<String, Object> getMetadata() {
        return metadata;
    }

    public void setMetadata(Map<String, Object> metadata) {
        this.metadata = metadata;
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

    public String getSummary() {
        return summary;
    }

    public void setSummary(String summary) {
        this.summary = summary;
    }

    public String getProteinExistenceResultRef() {
        return proteinExistenceResultRef;
    }

    public void setProteinExistenceResultRef(String proteinExistenceResultRef) {
        this.proteinExistenceResultRef = proteinExistenceResultRef;
    }
} 