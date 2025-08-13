"""
Analysis Stages for KBase Protein Query Module

This module contains all analysis stages that perform comprehensive protein analysis.
"""

import logging
import time
import numpy as np
import pandas as pd
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass

from .base_stage import BaseStage, StageResult
from ..analysis import sequence_analyzer
from ..processing.networks.builder import DynamicNetworkBuilder

logger = logging.getLogger(__name__)

class SequenceAnalysisStage(BaseStage):
    """
    Performs comprehensive sequence analysis.
    
    Handles:
    - Amino acid composition analysis
    - Physicochemical property calculation
    - Secondary structure prediction
    - Sequence motif identification
    - Bioinformatics database integration
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.sequence_analyzer = sequence_analyzer.ProteinSequenceAnalyzer()
        self.include_bioinformatics = config.get('include_bioinformatics', True) if config else True
        self.include_structure_prediction = config.get('include_structure_prediction', True) if config else True
    
    def get_stage_name(self) -> str:
        return "sequence_analysis"
    
    def get_required_inputs(self) -> List[str]:
        return ['protein_records']
    
    def get_optional_inputs(self) -> List[str]:
        return ['analysis_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['data_extraction']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'protein_records' not in input_data:
            logger.error("Missing protein_records")
            return False
        
        protein_records = input_data['protein_records']
        if not isinstance(protein_records, list) or len(protein_records) == 0:
            logger.error("protein_records must be a non-empty list")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'sequence_analyses': {
                'type': 'object',
                'description': 'Dictionary mapping protein IDs to sequence analysis results'
            },
            'analysis_stats': {
                'type': 'object',
                'properties': {
                    'total_analyses': {'type': 'integer'},
                    'successful_analyses': {'type': 'integer'},
                    'failed_analyses': {'type': 'integer'},
                    'property_distributions': {'type': 'object'},
                    'motif_statistics': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Perform sequence analysis for proteins."""
        start_time = time.time()
        
        try:
            protein_records = input_data['protein_records']
            analysis_config = input_data.get('analysis_config', {})
            
            # Perform sequence analysis
            sequence_analyses = {}
            total_analyses = len(protein_records)
            successful_analyses = 0
            failed_analyses = 0
            
            # Collect property distributions
            property_distributions = {
                'molecular_weights': [],
                'isoelectric_points': [],
                'instability_indices': [],
                'aliphatic_indices': [],
                'extinction_coefficients': []
            }
            
            motif_statistics = {
                'total_motifs': 0,
                'motif_types': {},
                'proteins_with_motifs': 0
            }
            
            for record in protein_records:
                try:
                    analysis_result = self.sequence_analyzer.analyze_sequence(
                        record.sequence, record.protein_id
                    )
                    
                    sequence_analyses[record.protein_id] = analysis_result
                    successful_analyses += 1
                    
                    # Collect property distributions
                    if 'physicochemical_properties' in analysis_result:
                        props = analysis_result['physicochemical_properties']
                        property_distributions['molecular_weights'].append(props.get('molecular_weight', 0))
                        property_distributions['isoelectric_points'].append(props.get('isoelectric_point', 0))
                        property_distributions['instability_indices'].append(props.get('instability_index', 0))
                        property_distributions['aliphatic_indices'].append(props.get('aliphatic_index', 0))
                        property_distributions['extinction_coefficients'].append(props.get('extinction_coefficient', 0))
                    
                    # Collect motif statistics
                    if 'sequence_motifs' in analysis_result:
                        motifs = analysis_result['sequence_motifs']
                        if motifs.get('motifs'):
                            motif_statistics['proteins_with_motifs'] += 1
                            motif_statistics['total_motifs'] += len(motifs['motifs'])
                            
                            for motif in motifs['motifs']:
                                motif_type = motif.get('type', 'unknown')
                                motif_statistics['motif_types'][motif_type] = motif_statistics['motif_types'].get(motif_type, 0) + 1
                    
                except Exception as e:
                    logger.warning(f"Failed to analyze sequence for {record.protein_id}: {e}")
                    sequence_analyses[record.protein_id] = {
                        'status': 'error',
                        'error': str(e),
                        'protein_id': record.protein_id
                    }
                    failed_analyses += 1
            
            # Calculate property statistics
            for prop_name, values in property_distributions.items():
                if values:
                    property_distributions[prop_name] = {
                        'mean': np.mean(values),
                        'median': np.median(values),
                        'min': np.min(values),
                        'max': np.max(values),
                        'std': np.std(values)
                    }
                else:
                    property_distributions[prop_name] = {}
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_analyses > 0,
                output_data={
                    'sequence_analyses': sequence_analyses,
                    'analysis_stats': {
                        'total_analyses': total_analyses,
                        'successful_analyses': successful_analyses,
                        'failed_analyses': failed_analyses,
                        'property_distributions': property_distributions,
                        'motif_statistics': motif_statistics
                    }
                },
                metadata={
                    'include_bioinformatics': self.include_bioinformatics,
                    'include_structure_prediction': self.include_structure_prediction,
                    'analysis_success_rate': successful_analyses / total_analyses if total_analyses > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Sequence analysis failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class StructureAnalysisStage(BaseStage):
    """
    Performs protein structure analysis and prediction.
    
    Handles:
    - Secondary structure prediction
    - Domain identification
    - Structure classification
    - Structural motif analysis
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.include_secondary_structure = config.get('include_secondary_structure', True) if config else True
        self.include_domain_prediction = config.get('include_domain_prediction', True) if config else True
        self.include_structure_classification = config.get('include_structure_classification', True) if config else True
    
    def get_stage_name(self) -> str:
        return "structure_analysis"
    
    def get_required_inputs(self) -> List[str]:
        return ['sequence_analyses']
    
    def get_optional_inputs(self) -> List[str]:
        return ['structure_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['sequence_analysis']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'sequence_analyses' not in input_data:
            logger.error("Missing sequence_analyses")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'structure_analyses': {
                'type': 'object',
                'description': 'Dictionary mapping protein IDs to structure analysis results'
            },
            'structure_stats': {
                'type': 'object',
                'properties': {
                    'total_analyses': {'type': 'integer'},
                    'successful_analyses': {'type': 'integer'},
                    'failed_analyses': {'type': 'integer'},
                    'structure_distributions': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Perform structure analysis for proteins."""
        start_time = time.time()
        
        try:
            sequence_analyses = input_data['sequence_analyses']
            structure_config = input_data.get('structure_config', {})
            
            # Perform structure analysis
            structure_analyses = {}
            total_analyses = len(sequence_analyses)
            successful_analyses = 0
            failed_analyses = 0
            
            # Collect structure distributions
            structure_distributions = {
                'secondary_structure': {'helix': 0, 'sheet': 0, 'coil': 0},
                'structure_classes': {},
                'domain_counts': []
            }
            
            for protein_id, sequence_analysis in sequence_analyses.items():
                try:
                    if sequence_analysis.get('status') == 'error':
                        structure_analyses[protein_id] = {
                            'status': 'error',
                            'error': sequence_analysis.get('error', 'Sequence analysis failed'),
                            'protein_id': protein_id
                        }
                        failed_analyses += 1
                        continue
                    
                    # Extract secondary structure prediction
                    structure_result = {
                        'protein_id': protein_id,
                        'status': 'success',
                        'secondary_structure': {},
                        'structure_classification': {},
                        'domain_analysis': {}
                    }
                    
                    # Secondary structure (from sequence analysis)
                    if 'secondary_structure' in sequence_analysis:
                        structure_result['secondary_structure'] = sequence_analysis['secondary_structure']
                        
                        # Update distributions
                        ss_pred = sequence_analysis['secondary_structure'].get('prediction', {})
                        structure_distributions['secondary_structure']['helix'] += ss_pred.get('helix_percentage', 0)
                        structure_distributions['secondary_structure']['sheet'] += ss_pred.get('sheet_percentage', 0)
                        structure_distributions['secondary_structure']['coil'] += ss_pred.get('coil_percentage', 0)
                    
                    # Structure classification based on sequence properties
                    if 'physicochemical_properties' in sequence_analysis:
                        props = sequence_analysis['physicochemical_properties']
                        
                        # Simple structure classification based on properties
                        structure_class = self._classify_structure(props)
                        structure_result['structure_classification'] = {
                            'predicted_class': structure_class,
                            'confidence': 0.7,  # Placeholder confidence
                            'classification_method': 'property_based'
                        }
                        
                        # Update distributions
                        structure_distributions['structure_classes'][structure_class] = \
                            structure_distributions['structure_classes'].get(structure_class, 0) + 1
                    
                    # Domain analysis (placeholder for future implementation)
                    structure_result['domain_analysis'] = {
                        'predicted_domains': [],
                        'domain_count': 0,
                        'analysis_method': 'placeholder'
                    }
                    
                    structure_analyses[protein_id] = structure_result
                    successful_analyses += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to analyze structure for {protein_id}: {e}")
                    structure_analyses[protein_id] = {
                        'status': 'error',
                        'error': str(e),
                        'protein_id': protein_id
                    }
                    failed_analyses += 1
            
            # Normalize secondary structure distributions
            if successful_analyses > 0:
                for ss_type in structure_distributions['secondary_structure']:
                    structure_distributions['secondary_structure'][ss_type] /= successful_analyses
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_analyses > 0,
                output_data={
                    'structure_analyses': structure_analyses,
                    'structure_stats': {
                        'total_analyses': total_analyses,
                        'successful_analyses': successful_analyses,
                        'failed_analyses': failed_analyses,
                        'structure_distributions': structure_distributions
                    }
                },
                metadata={
                    'include_secondary_structure': self.include_secondary_structure,
                    'include_domain_prediction': self.include_domain_prediction,
                    'analysis_success_rate': successful_analyses / total_analyses if total_analyses > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Structure analysis failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
    
    def _classify_structure(self, properties: Dict[str, Any]) -> str:
        """Classify protein structure based on physicochemical properties."""
        # Simple classification logic
        mw = properties.get('molecular_weight', 0)
        pi = properties.get('isoelectric_point', 7.0)
        instab = properties.get('instability_index', 40.0)
        
        if mw < 10000:
            return 'small_protein'
        elif instab < 40:
            return 'stable_protein'
        elif pi < 6:
            return 'acidic_protein'
        elif pi > 8:
            return 'basic_protein'
        else:
            return 'neutral_protein'

class NetworkAnalysisStage(BaseStage):
    """
    Performs network analysis and visualization.
    
    Handles:
    - Network construction from similarity data
    - Network topology analysis
    - Community detection
    - Network visualization generation
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.network_builder = DynamicNetworkBuilder()
        self.network_method = config.get('network_method', 'mutual_knn') if config else 'mutual_knn'
        self.k_neighbors = config.get('k_neighbors', 8) if config else 8
        self.similarity_threshold = config.get('similarity_threshold', 0.5) if config else 0.5
        self.min_network_size = config.get('min_network_size', 5) if config else 5
        self.max_network_size = config.get('max_network_size', 100) if config else 100
    
    def get_stage_name(self) -> str:
        return "network_analysis"
    
    def get_required_inputs(self) -> List[str]:
        return ['embeddings', 'similarity_results']
    
    def get_optional_inputs(self) -> List[str]:
        return ['network_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['similarity_search']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        required = ['embeddings', 'similarity_results']
        for field in required:
            if field not in input_data:
                logger.error(f"Missing required input: {field}")
                return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'networks': {
                'type': 'object',
                'description': 'Dictionary mapping protein IDs to network analysis results'
            },
            'network_stats': {
                'type': 'object',
                'properties': {
                    'total_networks': {'type': 'integer'},
                    'successful_networks': {'type': 'integer'},
                    'failed_networks': {'type': 'integer'},
                    'network_properties': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Perform network analysis for proteins."""
        start_time = time.time()
        
        try:
            embeddings = input_data['embeddings']
            similarity_results = input_data['similarity_results']
            network_config = input_data.get('network_config', {})
            
            # Perform network analysis
            networks = {}
            total_networks = len(embeddings)
            successful_networks = 0
            failed_networks = 0
            
            # Collect network properties
            network_properties = {
                'node_counts': [],
                'edge_counts': [],
                'densities': [],
                'diameters': [],
                'clustering_coefficients': []
            }
            
            for protein_id, embedding in embeddings.items():
                try:
                    similarity_result = similarity_results.get(protein_id, {})
                    
                    if similarity_result.get('status') != 'success':
                        networks[protein_id] = {
                            'status': 'no_similarity_data',
                            'protein_id': protein_id,
                            'network': None,
                            'network_properties': {}
                        }
                        failed_networks += 1
                        continue
                    
                    # Build network from similar proteins
                    similar_proteins = similarity_result.get('similar_proteins', [])
                    
                    if len(similar_proteins) < self.min_network_size:
                        networks[protein_id] = {
                            'status': 'insufficient_similar_proteins',
                            'protein_id': protein_id,
                            'network': None,
                            'network_properties': {
                                'node_count': 1,
                                'edge_count': 0,
                                'reason': f'Only {len(similar_proteins)} similar proteins found'
                            }
                        }
                        failed_networks += 1
                        continue
                    
                    # Create network
                    network, properties = self.network_builder.build_network_from_similar_proteins(
                        similar_proteins=similar_proteins,
                        embeddings=embeddings,
                        protein_ids=list(embeddings.keys()),
                        query_embedding=embedding,
                        query_protein_id=protein_id,
                        method=self.network_method
                    )
                    
                    # Analyze network properties
                    if network and len(network.nodes()) > 0:
                        network_analysis = self.network_builder.analyze_network_properties(network)
                        
                        # Update property distributions
                        network_properties['node_counts'].append(len(network.nodes()))
                        network_properties['edge_counts'].append(len(network.edges()))
                        network_properties['densities'].append(network_analysis.get('density', 0))
                        network_properties['diameters'].append(network_analysis.get('diameter', 0))
                        network_properties['clustering_coefficients'].append(network_analysis.get('clustering_coefficient', 0))
                        
                        networks[protein_id] = {
                            'status': 'success',
                            'protein_id': protein_id,
                            'network': network,
                            'network_properties': network_analysis,
                            'similar_proteins_count': len(similar_proteins)
                        }
                        successful_networks += 1
                    else:
                        networks[protein_id] = {
                            'status': 'network_construction_failed',
                            'protein_id': protein_id,
                            'network': None,
                            'network_properties': {}
                        }
                        failed_networks += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to analyze network for {protein_id}: {e}")
                    networks[protein_id] = {
                        'status': 'error',
                        'error': str(e),
                        'protein_id': protein_id,
                        'network': None,
                        'network_properties': {}
                    }
                    failed_networks += 1
            
            # Calculate network property statistics
            for prop_name, values in network_properties.items():
                if values:
                    network_properties[prop_name] = {
                        'mean': np.mean(values),
                        'median': np.median(values),
                        'min': np.min(values),
                        'max': np.max(values),
                        'std': np.std(values)
                    }
                else:
                    network_properties[prop_name] = {}
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_networks > 0,
                output_data={
                    'networks': networks,
                    'network_stats': {
                        'total_networks': total_networks,
                        'successful_networks': successful_networks,
                        'failed_networks': failed_networks,
                        'network_properties': network_properties
                    }
                },
                metadata={
                    'network_method': self.network_method,
                    'k_neighbors': self.k_neighbors,
                    'similarity_threshold': self.similarity_threshold,
                    'network_success_rate': successful_networks / total_networks if total_networks > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Network analysis failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )

class BioinformaticsAnalysisStage(BaseStage):
    """
    Performs bioinformatics analysis and database integration.
    
    Handles:
    - Database searches and annotations
    - Functional analysis
    - Pathway analysis
    - Evolutionary analysis
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__(config)
        self.include_functional_analysis = config.get('include_functional_analysis', True) if config else True
        self.include_pathway_analysis = config.get('include_pathway_analysis', True) if config else True
        self.include_evolutionary_analysis = config.get('include_evolutionary_analysis', True) if config else True
    
    def get_stage_name(self) -> str:
        return "bioinformatics_analysis"
    
    def get_required_inputs(self) -> List[str]:
        return ['sequence_analyses']
    
    def get_optional_inputs(self) -> List[str]:
        return ['bioinformatics_config', 'workspace_client']
    
    def get_stage_dependencies(self) -> List[str]:
        return ['sequence_analysis']
    
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """Validate input data."""
        if 'sequence_analyses' not in input_data:
            logger.error("Missing sequence_analyses")
            return False
        
        return True
    
    def get_output_schema(self) -> Dict[str, Any]:
        return {
            'bioinformatics_analyses': {
                'type': 'object',
                'description': 'Dictionary mapping protein IDs to bioinformatics analysis results'
            },
            'bioinformatics_stats': {
                'type': 'object',
                'properties': {
                    'total_analyses': {'type': 'integer'},
                    'successful_analyses': {'type': 'integer'},
                    'failed_analyses': {'type': 'integer'},
                    'database_coverage': {'type': 'object'}
                }
            }
        }
    
    def run(self, input_data: Dict[str, Any], workspace_client=None) -> StageResult:
        """Perform bioinformatics analysis for proteins."""
        start_time = time.time()
        
        try:
            sequence_analyses = input_data['sequence_analyses']
            bioinformatics_config = input_data.get('bioinformatics_config', {})
            
            # Perform bioinformatics analysis
            bioinformatics_analyses = {}
            total_analyses = len(sequence_analyses)
            successful_analyses = 0
            failed_analyses = 0
            
            # Collect database coverage
            database_coverage = {
                'uniprot_annotations': 0,
                'pdb_structures': 0,
                'pfam_domains': 0,
                'go_terms': 0,
                'kegg_pathways': 0
            }
            
            for protein_id, sequence_analysis in sequence_analyses.items():
                try:
                    if sequence_analysis.get('status') == 'error':
                        bioinformatics_analyses[protein_id] = {
                            'status': 'error',
                            'error': sequence_analysis.get('error', 'Sequence analysis failed'),
                            'protein_id': protein_id
                        }
                        failed_analyses += 1
                        continue
                    
                    # Extract bioinformatics links from sequence analysis
                    bioinformatics_result = {
                        'protein_id': protein_id,
                        'status': 'success',
                        'database_links': {},
                        'functional_analysis': {},
                        'pathway_analysis': {},
                        'evolutionary_analysis': {}
                    }
                    
                    # Database links (from sequence analysis)
                    if 'bioinformatics_links' in sequence_analysis:
                        bioinformatics_result['database_links'] = sequence_analysis['bioinformatics_links']
                        
                        # Update database coverage
                        for db_name in database_coverage.keys():
                            if db_name in sequence_analysis['bioinformatics_links']:
                                database_coverage[db_name] += 1
                    
                    # Functional analysis (placeholder)
                    if self.include_functional_analysis:
                        bioinformatics_result['functional_analysis'] = {
                            'predicted_functions': [],
                            'functional_domains': [],
                            'analysis_method': 'placeholder'
                        }
                    
                    # Pathway analysis (placeholder)
                    if self.include_pathway_analysis:
                        bioinformatics_result['pathway_analysis'] = {
                            'predicted_pathways': [],
                            'metabolic_roles': [],
                            'analysis_method': 'placeholder'
                        }
                    
                    # Evolutionary analysis (placeholder)
                    if self.include_evolutionary_analysis:
                        bioinformatics_result['evolutionary_analysis'] = {
                            'evolutionary_conservation': 0.0,
                            'orthologs': [],
                            'paralogs': [],
                            'analysis_method': 'placeholder'
                        }
                    
                    bioinformatics_analyses[protein_id] = bioinformatics_result
                    successful_analyses += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to perform bioinformatics analysis for {protein_id}: {e}")
                    bioinformatics_analyses[protein_id] = {
                        'status': 'error',
                        'error': str(e),
                        'protein_id': protein_id
                    }
                    failed_analyses += 1
            
            # Normalize database coverage
            if successful_analyses > 0:
                for db_name in database_coverage:
                    database_coverage[db_name] = database_coverage[db_name] / successful_analyses
            
            execution_time = time.time() - start_time
            
            return StageResult(
                success=successful_analyses > 0,
                output_data={
                    'bioinformatics_analyses': bioinformatics_analyses,
                    'bioinformatics_stats': {
                        'total_analyses': total_analyses,
                        'successful_analyses': successful_analyses,
                        'failed_analyses': failed_analyses,
                        'database_coverage': database_coverage
                    }
                },
                metadata={
                    'include_functional_analysis': self.include_functional_analysis,
                    'include_pathway_analysis': self.include_pathway_analysis,
                    'include_evolutionary_analysis': self.include_evolutionary_analysis,
                    'analysis_success_rate': successful_analyses / total_analyses if total_analyses > 0 else 0
                },
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Bioinformatics analysis failed: {e}")
            return StageResult(
                success=False,
                output_data={},
                metadata={},
                execution_time=execution_time,
                error_message=str(e)
            )
