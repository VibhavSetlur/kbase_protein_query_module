package us.kbase.kbaseproteinquerymodule;

import com.fasterxml.jackson.core.type.TypeReference;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import us.kbase.auth.AuthToken;
import us.kbase.common.service.JsonClientCaller;
import us.kbase.common.service.JsonClientException;
import us.kbase.common.service.RpcContext;
import us.kbase.common.service.UObject;
import us.kbase.common.service.UnauthorizedException;

/**
 * <p>Original spec-file module name: kbase_protein_query_module</p>
 * <pre>
 * A KBase module: kbase_protein_query_module
 * This module provides comprehensive protein query analysis capabilities using UniProt IDs as the canonical identifier:
 * COMPREHENSIVE ANALYSIS WORKFLOW:
 * 1. CheckProteinExistence: Verify protein exists using UniProt ID, optionally generate embedding
 * 2. GenerateProteinEmbeddings: Create embeddings from sequence input or protein check results
 * 3. AssignProteinFamily: Assign proteins to families using similarity to centroids
 * 4. FindTopMatches: Perform similarity search within families
 * 5. SummarizeAndVisualize: Generate comprehensive HTML reports with network analysis
 * ADVANCED CAPABILITIES:
 * - UniProt ID canonical identifier system (exact match only)
 * - ESM-2 protein language model for embedding generation
 * - Efficient FAISS-based similarity search and clustering
 * - Family assignment using binary centroid similarity
 * - Comprehensive metadata storage and retrieval
 * - HTML report generation with network visualization
 * - Workspace object management for downstream analysis
 * - Bioinformatics integration with protein databases
 * - Network analysis and protein relationship mapping
 * - Advanced similarity metrics and statistical analysis
 * Authors: Vibhav Setlur
 * Contact: https://kbase.us/contact-us/
 * </pre>
 */
public class KbaseProteinQueryModuleClient {
    private JsonClientCaller caller;
    private String serviceVersion = null;


    /** Constructs a client with a custom URL and no user credentials.
     * @param url the URL of the service.
     */
    public KbaseProteinQueryModuleClient(URL url) {
        caller = new JsonClientCaller(url);
    }
    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param token the user's authorization token.
     * @throws UnauthorizedException if the token is not valid.
     * @throws IOException if an IOException occurs when checking the token's
     * validity.
     */
    public KbaseProteinQueryModuleClient(URL url, AuthToken token) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, token);
    }

    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KbaseProteinQueryModuleClient(URL url, String user, String password) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password);
    }

    /** Constructs a client with a custom URL
     * and a custom authorization service URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @param auth the URL of the authorization server.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KbaseProteinQueryModuleClient(URL url, String user, String password, URL auth) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password, auth);
    }

    /** Get the token this client uses to communicate with the server.
     * @return the authorization token.
     */
    public AuthToken getToken() {
        return caller.getToken();
    }

    /** Get the URL of the service with which this client communicates.
     * @return the service URL.
     */
    public URL getURL() {
        return caller.getURL();
    }

    /** Set the timeout between establishing a connection to a server and
     * receiving a response. A value of zero or null implies no timeout.
     * @param milliseconds the milliseconds to wait before timing out when
     * attempting to read from a server.
     */
    public void setConnectionReadTimeOut(Integer milliseconds) {
        this.caller.setConnectionReadTimeOut(milliseconds);
    }

    /** Check if this client allows insecure http (vs https) connections.
     * @return true if insecure connections are allowed.
     */
    public boolean isInsecureHttpConnectionAllowed() {
        return caller.isInsecureHttpConnectionAllowed();
    }

    /** Deprecated. Use isInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public boolean isAuthAllowedForHttp() {
        return caller.isAuthAllowedForHttp();
    }

    /** Set whether insecure http (vs https) connections should be allowed by
     * this client.
     * @param allowed true to allow insecure connections. Default false
     */
    public void setIsInsecureHttpConnectionAllowed(boolean allowed) {
        caller.setInsecureHttpConnectionAllowed(allowed);
    }

    /** Deprecated. Use setIsInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public void setAuthAllowedForHttp(boolean isAuthAllowedForHttp) {
        caller.setAuthAllowedForHttp(isAuthAllowedForHttp);
    }

    /** Set whether all SSL certificates, including self-signed certificates,
     * should be trusted.
     * @param trustAll true to trust all certificates. Default false.
     */
    public void setAllSSLCertificatesTrusted(final boolean trustAll) {
        caller.setAllSSLCertificatesTrusted(trustAll);
    }
    
    /** Check if this client trusts all SSL certificates, including
     * self-signed certificates.
     * @return true if all certificates are trusted.
     */
    public boolean isAllSSLCertificatesTrusted() {
        return caller.isAllSSLCertificatesTrusted();
    }
    /** Sets streaming mode on. In this case, the data will be streamed to
     * the server in chunks as it is read from disk rather than buffered in
     * memory. Many servers are not compatible with this feature.
     * @param streamRequest true to set streaming mode on, false otherwise.
     */
    public void setStreamingModeOn(boolean streamRequest) {
        caller.setStreamingModeOn(streamRequest);
    }

    /** Returns true if streaming mode is on.
     * @return true if streaming mode is on.
     */
    public boolean isStreamingModeOn() {
        return caller.isStreamingModeOn();
    }

    public void _setFileForNextRpcResponse(File f) {
        caller.setFileForNextRpcResponse(f);
    }

    public String getServiceVersion() {
        return this.serviceVersion;
    }

    public void setServiceVersion(String newValue) {
        this.serviceVersion = newValue;
    }

    /**
     * <p>Original spec-file function name: check_protein_existence</p>
     * <pre>
     * </pre>
     * @param   params   instance of mapping from String to unspecified object
     * @return   parameter "output" of type {@link us.kbase.kbaseproteinquerymodule.CheckProteinExistenceResults CheckProteinExistenceResults}
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public CheckProteinExistenceResults checkProteinExistence(Map<String,UObject> params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<CheckProteinExistenceResults>> retType = new TypeReference<List<CheckProteinExistenceResults>>() {};
        List<CheckProteinExistenceResults> res = caller.jsonrpcCall("kbase_protein_query_module.check_protein_existence", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: generate_protein_embedding</p>
     * <pre>
     * </pre>
     * @param   params   instance of mapping from String to unspecified object
     * @return   parameter "output" of type {@link us.kbase.kbaseproteinquerymodule.GenerateProteinEmbeddingResults GenerateProteinEmbeddingResults}
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public GenerateProteinEmbeddingResults generateProteinEmbedding(Map<String,UObject> params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<GenerateProteinEmbeddingResults>> retType = new TypeReference<List<GenerateProteinEmbeddingResults>>() {};
        List<GenerateProteinEmbeddingResults> res = caller.jsonrpcCall("kbase_protein_query_module.generate_protein_embedding", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: assign_family_fast</p>
     * <pre>
     * </pre>
     * @param   params   instance of mapping from String to unspecified object
     * @return   parameter "output" of type {@link us.kbase.kbaseproteinquerymodule.AssignFamilyFastResults AssignFamilyFastResults}
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public AssignFamilyFastResults assignFamilyFast(Map<String,UObject> params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<AssignFamilyFastResults>> retType = new TypeReference<List<AssignFamilyFastResults>>() {};
        List<AssignFamilyFastResults> res = caller.jsonrpcCall("kbase_protein_query_module.assign_family_fast", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: find_top_matches_from_embedding</p>
     * <pre>
     * </pre>
     * @param   params   instance of mapping from String to unspecified object
     * @return   parameter "output" of type {@link us.kbase.kbaseproteinquerymodule.FindTopMatchesFromEmbeddingResults FindTopMatchesFromEmbeddingResults}
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public FindTopMatchesFromEmbeddingResults findTopMatchesFromEmbedding(Map<String,UObject> params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<FindTopMatchesFromEmbeddingResults>> retType = new TypeReference<List<FindTopMatchesFromEmbeddingResults>>() {};
        List<FindTopMatchesFromEmbeddingResults> res = caller.jsonrpcCall("kbase_protein_query_module.find_top_matches_from_embedding", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: summarize_and_visualize_results</p>
     * <pre>
     * </pre>
     * @param   params   instance of mapping from String to unspecified object
     * @return   parameter "output" of type {@link us.kbase.kbaseproteinquerymodule.SummarizeAndVisualizeResultsResults SummarizeAndVisualizeResultsResults}
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public SummarizeAndVisualizeResultsResults summarizeAndVisualizeResults(Map<String,UObject> params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<SummarizeAndVisualizeResultsResults>> retType = new TypeReference<List<SummarizeAndVisualizeResultsResults>>() {};
        List<SummarizeAndVisualizeResultsResults> res = caller.jsonrpcCall("kbase_protein_query_module.summarize_and_visualize_results", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: run_protein_query_analysis</p>
     * <pre>
     * </pre>
     * @param   params   instance of mapping from String to unspecified object
     * @return   parameter "output" of type {@link us.kbase.kbaseproteinquerymodule.ProteinQueryAnalysisResults ProteinQueryAnalysisResults}
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public ProteinQueryAnalysisResults runProteinQueryAnalysis(Map<String,UObject> params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<ProteinQueryAnalysisResults>> retType = new TypeReference<List<ProteinQueryAnalysisResults>>() {};
        List<ProteinQueryAnalysisResults> res = caller.jsonrpcCall("kbase_protein_query_module.run_protein_query_analysis", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    public Map<String, Object> status(RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        TypeReference<List<Map<String, Object>>> retType = new TypeReference<List<Map<String, Object>>>() {};
        List<Map<String, Object>> res = caller.jsonrpcCall("kbase_protein_query_module.status", args, retType, true, false, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }
}
