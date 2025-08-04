package us.kbase.kbaseproteinquerymodule;

import java.util.List;
import java.util.Map;

/**
 * <p>Original spec-file type: GenerateProteinEmbeddingResults</p>
 * <pre>
 * Generate a protein embedding from a sequence or workspace object.
 * GenerateProteinEmbeddingResults is a reference to a hash where the following keys are defined:
 * report_name has a value which is a string
 * report_ref has a value which is a string
 * embedding_result_ref has a value which is a string
 * summary has a value which is a string
 * input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
 * start_time has a value which is a float
 * embedding_norm has a value which is a float
 * sequence_length has a value which is an int
 * embedding_dim has a value which is an int
 * </pre>
 */
public class GenerateProteinEmbeddingResults {
    private String reportName;
    private String reportRef;
    private String embeddingResultRef;
    private String summary;
    private Map<String, Object> inputParameters;
    private Double startTime;
    private Double embeddingNorm;
    private Long sequenceLength;
    private Long embeddingDim;

    public GenerateProteinEmbeddingResults() {
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

    public String getEmbeddingResultRef() {
        return embeddingResultRef;
    }

    public void setEmbeddingResultRef(String embeddingResultRef) {
        this.embeddingResultRef = embeddingResultRef;
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

    public Double getEmbeddingNorm() {
        return embeddingNorm;
    }

    public void setEmbeddingNorm(Double embeddingNorm) {
        this.embeddingNorm = embeddingNorm;
    }

    public Long getSequenceLength() {
        return sequenceLength;
    }

    public void setSequenceLength(Long sequenceLength) {
        this.sequenceLength = sequenceLength;
    }

    public Long getEmbeddingDim() {
        return embeddingDim;
    }

    public void setEmbeddingDim(Long embeddingDim) {
        this.embeddingDim = embeddingDim;
    }
} 