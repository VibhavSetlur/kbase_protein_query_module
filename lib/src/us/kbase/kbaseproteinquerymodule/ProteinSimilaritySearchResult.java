package us.kbase.kbaseproteinquerymodule;

import java.util.List;
import java.util.Map;

/**
 * <p>Original spec-file type: ProteinSimilaritySearchResult</p>
 * <pre>
 * ProteinSimilaritySearchResult is a reference to a hash where the following keys are defined:
 * protein_id has a value which is a string
 * family_id has a value which is a string
 * similarity_score has a value which is a float
 * metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
 * </pre>
 */
public class ProteinSimilaritySearchResult {
    private String proteinId;
    private String familyId;
    private Double similarityScore;
    private Map<String, Object> metadata;

    public ProteinSimilaritySearchResult() {
    }

    public ProteinSimilaritySearchResult(String proteinId, String familyId, Double similarityScore, Map<String, Object> metadata) {
        this.proteinId = proteinId;
        this.familyId = familyId;
        this.similarityScore = similarityScore;
        this.metadata = metadata;
    }

    public String getProteinId() {
        return proteinId;
    }

    public void setProteinId(String proteinId) {
        this.proteinId = proteinId;
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

    public Map<String, Object> getMetadata() {
        return metadata;
    }

    public void setMetadata(Map<String, Object> metadata) {
        this.metadata = metadata;
    }
} 