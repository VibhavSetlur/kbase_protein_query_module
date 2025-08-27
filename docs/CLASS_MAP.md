# Class Map and Dependencies


## Core Framework

### `AnalysisMetadata`
- **Location**: `lib/kbase_protein_query_module/src/core/analysis_registry.py:48`
- **Description**: Metadata for registered analyses.

### `AnalysisRegistry`
- **Location**: `lib/kbase_protein_query_module/src/core/analysis_registry.py:153`
- **Description**: 
    Registry for managing available analysis types.
    
    This registry allows dynamic discovery and registration of analysis types,
    making it easy to extend the module with new analyses.
    ...
- **Key Methods**: register_analysis, get_analysis, list_analyses, get_analyses_by_category, generate_analysis_documentation

### `BaseAnalysis`
- **Location**: `lib/kbase_protein_query_module/src/core/analysis_registry.py:62`
- **Extends**: ABC
- **Description**: 
    Base class for all protein analyses.
    
    All new analysis types MUST inherit from this class and implement
    the required abstract methods.
    
    CLASS LOCATION: lib/kbase_protein_query...
- **Key Methods**: analyze, get_output_files, get_metadata, validate_input, estimate_resources

### `BatchQueue`
- **Location**: `lib/kbase_protein_query_module/src/core/parallel_processor.py:371`
- **Description**: 
    Intelligent batch queue for managing large-scale processing jobs.
    
- **Key Methods**: add_task, get_next_batch, is_empty, get_queue_status

### `MemoryOptimizer`
- **Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:257`
- **Description**: Memory optimization utilities for large-scale processing.
- **Key Methods**: optimize_numpy_arrays, batch_iterator

### `ParallelProcessor`
- **Location**: `lib/kbase_protein_query_module/src/core/parallel_processor.py:36`
- **Used in**: 2 files
- **Description**: 
    High-performance parallel processor for protein analysis workflows.
    
    Features:
    - Adaptive thread/process pool management
    - Resource-aware task scheduling
    - Fault tolerance wit...
- **Key Methods**: process_batch, process_proteins_parallel, get_performance_metrics

### `PerformanceMetric`
- **Location**: `lib/kbase_protein_query_module/src/core/performance_monitor.py:27`
- **Description**: A single performance metric measurement.

### `PerformanceProfiler`
- **Location**: `lib/kbase_protein_query_module/src/core/performance_monitor.py:38`
- **Used in**: 1 files
- **Description**: 
    Comprehensive performance profiler for protein analysis operations.
    
    Features:
    - Real-time performance monitoring
    - Memory usage tracking
    - CPU utilization analysis
    - I/O ...
- **Key Methods**: profile_operation, get_performance_summary, export_metrics, generate_performance_report, set_baseline

### `PipelineConfig`
- **Location**: `lib/kbase_protein_query_module/src/core/pipeline_config.py:18`
- **Used in**: 4 files
- **Description**: 
    Configuration for the protein query analysis pipeline.
    
    This dataclass contains all configuration parameters needed to run
    the protein query analysis workflow.
    
- **Key Methods**: get_resource_limits, to_dict, from_dict

### `ProcessingTask`
- **Location**: `lib/kbase_protein_query_module/src/core/parallel_processor.py:26`
- **Description**: A task for parallel processing.

### `ResourceLimits`
- **Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:22`
- **Used in**: 2 files
- **Description**: 
    Server-aware resource limits for KBase DOE servers.
    
    These limits are designed to be respectful of shared server resources
    and ensure efficient operation without impacting other proce...

### `ResourceManager`
- **Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:62`
- **Used in**: 3 files
- **Description**: 
    Server-aware resource manager for KBase DOE server environments.
    
    This manager is designed for shared server environments and uses percentage-based
    limits to ensure respectful resourc...
- **Key Methods**: get_current_metrics, check_resource_availability, optimize_batch_size, resource_context, get_performance_summary

### `ResourceMetrics`
- **Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:52`
- **Description**: Current resource usage metrics.

### `ScalabilityConfig`
- **Location**: `lib/kbase_protein_query_module/src/core/resource_manager.py:300`
- **Description**: Configuration for scalability parameters.

### `StreamingProcessor`
- **Location**: `lib/kbase_protein_query_module/src/core/parallel_processor.py:248`
- **Description**: 
    Streaming processor for handling large datasets that don't fit in memory.
    
- **Key Methods**: stream_process


## Pipeline Stages

### `BaseStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/base_stage.py:32`
- **Extends**: ABC
- **Extended by**: IntegratedAnalysisStage, ReportGenerationStage, VisualizationStage, DataExportStage, EmbeddingGenerationStage, FamilyAssignmentStage, SimilaritySearchStage, InputValidationStage, DataExtractionStage, WorkspaceObjectStage, DataExportStage, ReportGenerationStage, VisualizationStage, EmbeddingGenerationStage, FamilyAssignmentStage, SimilaritySearchStage, SequenceAnalysisStage, BioinformaticsAnalysisStage, NetworkAnalysisStage, MultiProteinAnalysisStage, WorkspaceObjectStage, DataExtractionStage, InputValidationStage
- **Used in**: 21 files
- **Description**: 
    Base class for all pipeline stages.
    
    Each stage must implement:
    - run(): Execute the stage logic
    - validate_input(): Validate input data
    - get_output_schema(): Define output d...
- **Key Methods**: get_stage_name, validate_input, get_output_schema, run, pre_run_hook

### `BioinformaticsAnalysisStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/analysis/bioinformatics_analysis.py:12`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: Bioinformatics analysis stage leveraging SequenceAnalysisStage links and motifs.
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `DataExportStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/output/data_export.py:15`
- **Extends**: BaseStage
- **Used in**: 5 files
- **Description**: Data export stage that writes selected outputs to disk as JSON/CSV/HTML.
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `DataExtractionStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/input/data_extraction.py:19`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: 
    Extracts protein data from various sources.
    
    Handles:
    - FASTA file/string parsing
    - UniProt data retrieval
    - Workspace object extraction
    - Data format conversion
    
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `EmbeddingGenerationStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/processing/embedding_generation.py:17`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: 
    Generates protein embeddings using deep learning models.
    
    Handles:
    - Protein sequence embedding generation
    - Model loading and caching
    - Batch processing
    - Embedding stora...
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `FamilyAssignmentStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/processing/family_assignment.py:16`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: Family assignment stage using FAISS binary centroid search.
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `InputData`
- **Location**: `lib/kbase_protein_query_module/src/stages/input/input_validation.py:21`
- **Description**: Container for validated input data.

### `InputValidationStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/input/input_validation.py:29`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: 
    Validates and preprocesses input data from various sources.
    
    Handles:
    - Input type detection
    - Basic validation
    - Format standardization
    
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `IntegratedAnalysisStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/analysis_stages.py:21`
- **Extends**: BaseStage
- **Used in**: 2 files
- **Description**: 
    Integrated wrapper stage that orchestrates all analysis sub-stages without performing
    any analysis itself. It delegates to:
      - `SequenceAnalysisStage` (from analysis package)
      - `Ne...
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `MultiProteinAnalysisStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/analysis/multi_protein_analysis.py:31`
- **Extends**: BaseStage
- **Used in**: 1 files
- **Description**: 
    Multi-protein analysis stage that performs MSA and phylogenetic analysis.
    
    This stage analyzes relationships between multiple proteins using:
    - MAFFT for multiple sequence alignment
 ...
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `NetworkAnalysisStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/analysis/network_analysis.py:41`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: Pipeline stage for network analysis of similarity search results.
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `ReportGenerationStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/output/report_generation.py:16`
- **Extends**: BaseStage
- **Used in**: 5 files
- **Description**: Report generation stage that produces comprehensive HTML reports.
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `SequenceAnalysisStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/analysis/sequence_analysis.py:33`
- **Extends**: BaseStage
- **Used in**: 7 files
- **Description**: 
    Pipeline stage for comprehensive protein sequence analysis with external references.
    
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `SimilaritySearchStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/processing/similarity_search.py:16`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: Similarity search stage using FAISS float indices.
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `StageResult`
- **Location**: `lib/kbase_protein_query_module/src/stages/base_stage.py:18`
- **Used in**: 21 files
- **Description**: Result container for pipeline stage execution.

### `VisualizationStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/output/visualization.py:18`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: Visualization stage for generating interactive network visualizations.
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `WorkspaceObjectStage`
- **Location**: `lib/kbase_protein_query_module/src/stages/input/workspace_object.py:16`
- **Extends**: BaseStage
- **Used in**: 4 files
- **Description**: 
    Handles workspace object retrieval and processing.
    
    Handles:
    - Workspace object validation
    - Object data retrieval
    - Data format conversion
    
- **Key Methods**: get_stage_name, get_required_inputs, get_optional_inputs, validate_input, get_output_schema

### `_VD`
- **Location**: `lib/kbase_protein_query_module/src/stages/input/data_extraction.py:84`


## Storage & Indexing

### `CompressedMetadataStorage`
- **Location**: `lib/kbase_protein_query_module/src/storage/protein_storage.py:1137`
- **Used in**: 2 files
- **Description**: 
    Efficient metadata storage with compression and indexing.
    
- **Key Methods**: store_metadata, load_metadata

### `FAISSIndexingStrategy`
- **Location**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:322`
- **Extends**: IndexingStrategy
- **Extended by**: RegisteredFAISSStrategy
- **Description**: 
    FAISS-based indexing strategy with advanced features.
    
    Supports various FAISS index types including IVF, PQ, and HNSW for different
    use cases and performance requirements.
    
    CL...
- **Key Methods**: build_index, search, get_index_info, save_index, load_index

### `IndexingConfig`
- **Location**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:58`
- **Description**: 
    Configuration for indexing strategies.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/storage/indexing_strategy.py:50
    USED BY: All indexing strategy implementations
    

### `IndexingRegistry`
- **Location**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:567`
- **Description**: Registry implementation for indexing strategies.
- **Key Methods**: register_strategy, get_strategy, list_strategies

### `IndexingStrategy`
- **Location**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:75`
- **Extends**: ABC
- **Extended by**: FAISSIndexingStrategy
- **Description**: 
    Abstract base class for all indexing strategies.
    
    This class defines the interface that all indexing implementations must follow.
    New indexing systems should inherit from this class a...
- **Key Methods**: build_index, search, get_index_info, save_index, load_index

### `MemoryEfficientLoader`
- **Location**: `lib/kbase_protein_query_module/src/storage/protein_storage.py:1268`
- **Used in**: 2 files
- **Description**: 
    Memory-efficient loading for large datasets with mapped families only.
    
- **Key Methods**: load_families_batch

### `ProteinExistenceChecker`
- **Location**: `lib/kbase_protein_query_module/src/storage/protein_existence_checker.py:16`
- **Used in**: 5 files
- **Description**: 
    Checks if a protein exists in the storage and returns its family and metadata.
    Uses efficient protein IDs index for fast searching (exact UniProt ID match).
    
- **Key Methods**: check_protein_existence, check_protein_exists, check_protein_with_metadata, load_database, check_multiple_proteins

### `ProteinFamilyAssigner`
- **Location**: `lib/kbase_protein_query_module/src/storage/protein_family_assigner.py:18`
- **Used in**: 5 files
- **Description**: 
    Assigns protein embeddings to precomputed protein families using centroid similarity.
    All protein identifiers are UniProt IDs (exact match only).
    
- **Key Methods**: load_family_centroids, assign_family, predict_family, assign_families_batch, assign_family_with_threshold

### `ProteinIDsIndex`
- **Location**: `lib/kbase_protein_query_module/src/storage/protein_storage.py:1310`
- **Used in**: 2 files
- **Description**: 
    Efficient protein IDs index for fast searching (exact UniProt ID match only).
    Optimized for mapped families only.
    
- **Key Methods**: search_protein, get_protein_family, get_all_proteins, get_proteins_by_family

### `ProteinStorage`
- **Location**: `lib/kbase_protein_query_module/src/storage/protein_storage.py:24`
- **Used in**: 6 files
- **Description**: 
    Advanced storage system for massive protein datasets with hybrid FAISS indexing.
    
    Features:
    - Selective loading of mapped families only (memory efficient)
    - Hybrid FAISS indexing:...
- **Key Methods**: create_hybrid_family_index, search_within_family_hybrid, create_family_centroid_binary_advanced, assign_family_advanced, save_html_file

### `RegisteredFAISSStrategy`
- **Location**: `lib/kbase_protein_query_module/src/storage/indexing_strategy.py:593`
- **Extends**: FAISSIndexingStrategy
- **Description**: Registered FAISS strategy for the registry.


## Utilities

### `ClassTracker`
- **Location**: `lib/kbase_protein_query_module/src/utils/documentation_generator.py:27`
- **Description**: 
    Tracks class definitions, usage, and dependencies across the module.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/utils/documentation_generator.py:25
    USED BY: DocumentationGene...
- **Key Methods**: scan_codebase, get_class_info, find_class_usage, generate_class_map

### `DataFormatter`
- **Location**: `lib/kbase_protein_query_module/src/reports/utils/data_formatter.py:14`
- **Used in**: 1 files
- **Description**: Formats data for HTML report generation.
- **Key Methods**: format_protein_data, format_similarity_scores, format_family_data, create_summary_table

### `DocumentationGenerator`
- **Location**: `lib/kbase_protein_query_module/src/utils/documentation_generator.py:197`
- **Description**: 
    Generates comprehensive documentation for module extensions.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/utils/documentation_generator.py:120
    USED BY: Developers for understan...
- **Key Methods**: generate_extension_guide

### `InputParser`
- **Location**: `lib/kbase_protein_query_module/src/utils/input_parser.py:33`
- **Used in**: 1 files
- **Description**: 
    Parser for various protein input formats.
    
    Supports:
    - FASTA files (local and remote)
    - UniProt identifiers
    - ProteinSequenceSet objects
    - Genome references
    - Single p...
- **Key Methods**: parse_input, validate_records

### `ProteinRecord`
- **Location**: `lib/kbase_protein_query_module/src/utils/input_parser.py:19`
- **Used in**: 1 files
- **Description**: Data class representing a protein record.


## Other

### `Application`
- **Location**: `lib/kbase_protein_query_module/kbase_protein_query_moduleServer.py:316`
- **Extends**: object
- **Key Methods**: logcallback, log, process_error, now_in_utc

### `BaseClient`
- **Location**: `lib/kbase_protein_query_module/baseclient.py:100`
- **Extends**: object
- **Used in**: 1 files
- **Description**: 
    The KBase base client.
    Required initialization arguments (positional):
    url - the url of the the service to contact:
        For SDK methods: either the url of the callback service or the
...
- **Key Methods**: run_job, call_method

### `ChartGenerator`
- **Location**: `lib/kbase_protein_query_module/src/reports/html/components/charts.py:12`
- **Used in**: 1 files
- **Description**: Generates chart components for HTML reports.
- **Key Methods**: generate_similarity_chart, generate_family_chart

### `DashboardGenerator`
- **Location**: `lib/kbase_protein_query_module/src/reports/html/components/dashboard.py:12`
- **Used in**: 1 files
- **Description**: Generates dashboard components for HTML reports.
- **Key Methods**: generate_dashboard, generate_summary_stats

### `DynamicNetworkBuilder`
- **Location**: `lib/kbase_protein_query_module/src/processing/networks/builder.py:483`
- **Extended by**: NetworkBuilder
- **Used in**: 7 files
- **Description**: 
    Builds dynamic, localized protein networks from similarity search results.
    
    This class constructs targeted networks around query proteins using
    various network construction methods in...
- **Key Methods**: build_mutual_knn_network, build_threshold_network, build_hybrid_network, build_network_from_similar_proteins, create_interactive_visualization

### `HTMLReportGenerator`
- **Location**: `lib/kbase_protein_query_module/src/reports/html/generator.py:26`
- **Used in**: 4 files
- **Description**: 
    Generates comprehensive HTML reports with multi-layer design.
    
    Features:
    - Dashboard layer with aggregated visualizations
    - Detail layer for individual protein analysis
    - Mult...
- **Key Methods**: generate_comprehensive_report

### `HierarchicalIndex`
- **Location**: `lib/kbase_protein_query_module/src/processing/similarity/hierarchical_index.py:24`
- **Extended by**: SimilarityIndex
- **Used in**: 6 files
- **Description**: 
    Hierarchical index structure for massive protein datasets.
    All search and indexing is performed using UniProt IDs as the canonical identifier (exact match only).
    
- **Key Methods**: create_family_index, create_family_index_float, search_family, search_family_float, search_all_families

### `JSONObjectEncoder`
- **Location**: `lib/kbase_protein_query_module/kbase_protein_query_moduleServer.py:58`
- **Extends**: json.JSONEncoder
- **Key Methods**: default

### `JSONRPCServiceCustom`
- **Location**: `lib/kbase_protein_query_module/kbase_protein_query_moduleServer.py:70`
- **Extends**: JSONRPCService
- **Key Methods**: call, call_py

### `KBaseAuth`
- **Location**: `lib/kbase_protein_query_module/authclient.py:58`
- **Extends**: object
- **Used in**: 1 files
- **Description**: 
    A very basic KBase auth client for the Python server.
    
- **Key Methods**: get_user

### `MethodContext`
- **Location**: `lib/kbase_protein_query_module/kbase_protein_query_moduleServer.py:200`
- **Extends**: dict
- **Key Methods**: log_err, log_info, log_debug, set_log_level, get_log_level

### `NetworkBuilder`
- **Location**: `lib/kbase_protein_query_module/src/processing/networks/__init__.py:10`
- **Extends**: DynamicNetworkBuilder
- **Used in**: 7 files

### `NetworkVisualizationGenerator`
- **Location**: `lib/kbase_protein_query_module/src/reports/html/components/network_viz.py:12`
- **Used in**: 1 files
- **Description**: Generates network visualization components for HTML reports.
- **Key Methods**: generate_network_visualizations, generate_interactive_network

### `ProteinDetailsGenerator`
- **Location**: `lib/kbase_protein_query_module/src/reports/html/components/protein_details.py:12`
- **Used in**: 1 files
- **Description**: Generates protein details components for HTML reports.
- **Key Methods**: generate_protein_details, generate_protein_table

### `ProteinEmbeddingGenerator`
- **Location**: `lib/kbase_protein_query_module/src/processing/embeddings/generator.py:32`
- **Used in**: 6 files
- **Description**: 
    Generates protein embeddings using ESM-2 models.
    
    This class handles the conversion of protein sequences to high-dimensional
    embeddings that capture structural and functional informat...
- **Key Methods**: generate_embedding, generate_embeddings_batch, generate_embeddings, save_embeddings, load_embeddings

### `ProteinQueryWorkflow`
- **Location**: `lib/kbase_protein_query_module/src/workflows/workflow_orchestrator.py:54`
- **Used in**: 3 files
- **Description**: 
    Comprehensive workflow orchestrator for protein query analysis.
    
    This orchestrator provides:
    - Modular stage-based architecture
    - Automatic dependency resolution
    - KBase works...
- **Key Methods**: load_family_subset, perform_optimized_similarity_search, build_optimized_network, generate_query_embedding, classify_query_family

### `ReferenceDataLoader`
- **Location**: `lib/kbase_protein_query_module/src/data/reference_loader.py:18`
- **Used in**: 2 files
- **Description**: 
    Loader for reference data used in protein sequence analysis.
    
    This class loads all reference data from external JSON files, making the system
    more maintainable and scalable. It follow...
- **Key Methods**: amino_acids, motif_patterns, bioinformatics_databases, physicochemical_constants, get_amino_acid_properties

### `ServerError`
- **Location**: `lib/kbase_protein_query_module/kbase_protein_query_moduleServer.py:280`
- **Extends**: Exception
- **Used in**: 1 files
- **Description**: 
    The call returned an error. Fields:
    name - the name of the error.
    code - the error code.
    message - a human readable error message.
    data - the server side stacktrace.
    

### `SimilarityIndex`
- **Location**: `lib/kbase_protein_query_module/src/processing/similarity/__init__.py:10`
- **Extends**: HierarchicalIndex

### `StreamingIndex`
- **Location**: `lib/kbase_protein_query_module/src/processing/similarity/hierarchical_index.py:505`
- **Used in**: 3 files
- **Description**: 
    Streaming index for memory-efficient processing of massive datasets.
    
- **Key Methods**: create_streaming_index, stream_search

### `TokenCache`
- **Location**: `lib/kbase_protein_query_module/authclient.py:14`
- **Extends**: object
- **Description**:  A basic cache for tokens. 
- **Key Methods**: get_user, add_valid_token

### `WorkflowResult`
- **Location**: `lib/kbase_protein_query_module/src/workflows/workflow_orchestrator.py:40`
- **Used in**: 3 files
- **Description**: Result container for the complete workflow execution.

### `_JSONObjectEncoder`
- **Location**: `lib/kbase_protein_query_module/baseclient.py:90`
- **Extends**: _json.JSONEncoder
- **Key Methods**: default

### `kbase_protein_query_module`
- **Location**: `lib/kbase_protein_query_module/kbase_protein_query_moduleImpl.py:20`
- **Used in**: 1 files
- **Description**: 
    Module Name:
    kbase_protein_query_module

    Module Description:
    A KBase module: kbase_protein_query_module

This module provides comprehensive protein query analysis capabilities using U...
- **Key Methods**: check_protein_existence, generate_protein_embedding, assign_family_fast, find_top_matches_from_embedding, summarize_and_visualize_results

