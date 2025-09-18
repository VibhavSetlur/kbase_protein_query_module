"""
Network Analysis Output Handler

Handles output generation for protein network analysis results.
"""

import json
from typing import Dict, Any, List
from ...base_analysis_output import BaseAnalysisOutput, AnalysisOutputResult
from ...analysis_output_registry import register_analysis_output


@register_analysis_output("network_analysis")
class NetworkAnalysisOutput(BaseAnalysisOutput):
    """Output handler for network analysis."""
    
    def generate_outputs(self, analysis_data: Dict[str, Any], stage_name: str) -> AnalysisOutputResult:
        """Generate network analysis outputs."""
        output_files = []
        metadata = {}
        
        try:
            # Extract network analysis results
            network_graph = analysis_data.get('network_graph', {})
            node_centrality = analysis_data.get('node_centrality', {})
            community_structure = analysis_data.get('community_structure', {})
            pathway_enrichment = analysis_data.get('pathway_enrichment', [])
            
            # Generate network summary
            network_summary = {
                'network_graph': network_graph,
                'node_centrality': node_centrality,
                'community_structure': community_structure,
                'pathway_enrichment': pathway_enrichment,
                'analysis_timestamp': analysis_data.get('timestamp'),
                'stage_name': stage_name
            }
            
            summary_path = self.output_manager.write_json(
                'network_analysis',
                'network_summary.json',
                network_summary,
                stage=stage_name,
                description='Network analysis summary with graph structure and centrality'
            )
            output_files.append(summary_path)
            
            # Generate centrality analysis
            if node_centrality:
                centrality_data = {
                    'node_centrality': node_centrality,
                    'centrality_metrics': list(node_centrality.keys()) if node_centrality else [],
                    'analysis_method': analysis_data.get('centrality_method', 'unknown'),
                    'stage_name': stage_name
                }
                centrality_path = self.output_manager.write_json(
                    'network_analysis',
                    'centrality_analysis.json',
                    centrality_data,
                    stage=stage_name,
                    description='Node centrality analysis results'
                )
                output_files.append(centrality_path)
            
            # Generate community structure
            if community_structure:
                community_data = {
                    'community_structure': community_structure,
                    'community_count': len(community_structure.get('communities', [])),
                    'modularity': community_structure.get('modularity', 0.0),
                    'stage_name': stage_name
                }
                community_path = self.output_manager.write_json(
                    'network_analysis',
                    'community_structure.json',
                    community_data,
                    stage=stage_name,
                    description='Network community structure analysis'
                )
                output_files.append(community_path)
            
            # Generate pathway enrichment
            if pathway_enrichment:
                pathway_data = {
                    'pathway_enrichment': pathway_enrichment,
                    'enriched_pathways': len(pathway_enrichment),
                    'enrichment_method': analysis_data.get('enrichment_method', 'unknown'),
                    'stage_name': stage_name
                }
                pathway_path = self.output_manager.write_json(
                    'network_analysis',
                    'pathway_enrichment.json',
                    pathway_data,
                    stage=stage_name,
                    description='Pathway enrichment analysis results'
                )
                output_files.append(pathway_path)
            
            # Generate network statistics
            network_stats = self._generate_network_statistics(network_graph, node_centrality, community_structure)
            stats_path = self.output_manager.write_json(
                'network_analysis',
                'network_statistics.json',
                network_stats,
                stage=stage_name,
                description='Network topology and analysis statistics'
            )
            output_files.append(stats_path)
            
            metadata = {
                'nodes': len(network_graph.get('nodes', [])),
                'edges': len(network_graph.get('edges', [])),
                'communities': len(community_structure.get('communities', [])),
                'enriched_pathways': len(pathway_enrichment),
                'output_files': len(output_files),
                'analysis_type': 'network_analysis'
            }
            
            summary = f"Network analysis completed: {len(network_graph.get('nodes', []))} nodes, {len(network_graph.get('edges', []))} edges, {len(community_structure.get('communities', []))} communities"
            
            return AnalysisOutputResult(
                success=True,
                output_files=output_files,
                metadata=metadata,
                summary=summary
            )
            
        except Exception as e:
            return AnalysisOutputResult(
                success=False,
                output_files=[],
                metadata={},
                summary=f"Error generating network analysis outputs: {str(e)}",
                error_message=str(e)
            )
    
    def _generate_network_statistics(self, network_graph: Dict[str, Any], node_centrality: Dict[str, Any], community_structure: Dict[str, Any]) -> Dict[str, Any]:
        """Generate network statistics."""
        nodes = network_graph.get('nodes', [])
        edges = network_graph.get('edges', [])
        
        stats = {
            'network_size': {
                'nodes': len(nodes),
                'edges': len(edges),
                'density': len(edges) / (len(nodes) * (len(nodes) - 1) / 2) if len(nodes) > 1 else 0.0
            },
            'centrality_summary': {},
            'community_summary': {
                'communities': len(community_structure.get('communities', [])),
                'modularity': community_structure.get('modularity', 0.0)
            }
        }
        
        # Calculate centrality statistics
        if node_centrality:
            for metric, values in node_centrality.items():
                if isinstance(values, dict) and values:
                    metric_values = list(values.values())
                    stats['centrality_summary'][metric] = {
                        'mean': sum(metric_values) / len(metric_values),
                        'max': max(metric_values),
                        'min': min(metric_values)
                    }
        
        return stats
    
    def get_output_description(self) -> str:
        """Get description of network analysis outputs."""
        return "Network topology, centrality analysis, community structure, and pathway enrichment"
    
    def get_supported_formats(self) -> List[str]:
        """Get supported output formats."""
        return ['json']
    
    def validate_analysis_data(self, analysis_data: Dict[str, Any]) -> bool:
        """Validate network analysis data."""
        # At least network graph should be present
        return 'network_graph' in analysis_data
