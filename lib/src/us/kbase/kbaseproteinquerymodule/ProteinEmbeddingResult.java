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

/**
 * <p>Original spec-file type: ProteinEmbeddingResult</p>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_id",
    "input_source",
    "embedding_ref",
    "embedding",
    "model_name",
    "pooling_method",
    "metadata",
    "sequence_length",
    "embedding_norm",
    "embedding_dim"
})
public class ProteinEmbeddingResult {

    @JsonProperty("input_id")
    private String inputId;
    @JsonProperty("input_source")
    private String inputSource;
    @JsonProperty("embedding_ref")
    private String embeddingRef;
    @JsonProperty("embedding")
    private List<Double> embedding;
    @JsonProperty("model_name")
    private String modelName;
    @JsonProperty("pooling_method")
    private String poolingMethod;
    @JsonProperty("metadata")
    private Map<String, String> metadata;
    @JsonProperty("sequence_length")
    private Long sequenceLength;
    @JsonProperty("embedding_norm")
    private Double embeddingNorm;
    @JsonProperty("embedding_dim")
    private Long embeddingDim;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("input_id")
    public String getInputId() {
        return inputId;
    }

    @JsonProperty("input_id")
    public void setInputId(String inputId) {
        this.inputId = inputId;
    }

    public ProteinEmbeddingResult withInputId(String inputId) {
        this.inputId = inputId;
        return this;
    }

    @JsonProperty("input_source")
    public String getInputSource() {
        return inputSource;
    }

    @JsonProperty("input_source")
    public void setInputSource(String inputSource) {
        this.inputSource = inputSource;
    }

    public ProteinEmbeddingResult withInputSource(String inputSource) {
        this.inputSource = inputSource;
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

    public ProteinEmbeddingResult withEmbeddingRef(String embeddingRef) {
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

    public ProteinEmbeddingResult withEmbedding(List<Double> embedding) {
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

    public ProteinEmbeddingResult withModelName(String modelName) {
        this.modelName = modelName;
        return this;
    }

    @JsonProperty("pooling_method")
    public String getPoolingMethod() {
        return poolingMethod;
    }

    @JsonProperty("pooling_method")
    public void setPoolingMethod(String poolingMethod) {
        this.poolingMethod = poolingMethod;
    }

    public ProteinEmbeddingResult withPoolingMethod(String poolingMethod) {
        this.poolingMethod = poolingMethod;
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

    public ProteinEmbeddingResult withMetadata(Map<String, String> metadata) {
        this.metadata = metadata;
        return this;
    }

    @JsonProperty("sequence_length")
    public Long getSequenceLength() {
        return sequenceLength;
    }

    @JsonProperty("sequence_length")
    public void setSequenceLength(Long sequenceLength) {
        this.sequenceLength = sequenceLength;
    }

    public ProteinEmbeddingResult withSequenceLength(Long sequenceLength) {
        this.sequenceLength = sequenceLength;
        return this;
    }

    @JsonProperty("embedding_norm")
    public Double getEmbeddingNorm() {
        return embeddingNorm;
    }

    @JsonProperty("embedding_norm")
    public void setEmbeddingNorm(Double embeddingNorm) {
        this.embeddingNorm = embeddingNorm;
    }

    public ProteinEmbeddingResult withEmbeddingNorm(Double embeddingNorm) {
        this.embeddingNorm = embeddingNorm;
        return this;
    }

    @JsonProperty("embedding_dim")
    public Long getEmbeddingDim() {
        return embeddingDim;
    }

    @JsonProperty("embedding_dim")
    public void setEmbeddingDim(Long embeddingDim) {
        this.embeddingDim = embeddingDim;
    }

    public ProteinEmbeddingResult withEmbeddingDim(Long embeddingDim) {
        this.embeddingDim = embeddingDim;
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
        return ((((((((((((((((((((((("ProteinEmbeddingResult"+" [inputId=")+ inputId)+", inputSource=")+ inputSource)+", embeddingRef=")+ embeddingRef)+", embedding=")+ embedding)+", modelName=")+ modelName)+", poolingMethod=")+ poolingMethod)+", metadata=")+ metadata)+", sequenceLength=")+ sequenceLength)+", embeddingNorm=")+ embeddingNorm)+", embeddingDim=")+ embeddingDim)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
