# Pipeline HTML Test Run Summary

**Timestamp:** 20250827_160759
**Run Directory:** /home/vibhav/Downloads/Work/ANL/Research/kbase_protein_query_module/test/outputs/pipeline_run_20250827_160759

## Generated Reports

### Single Protein Analysis
- `demo_single_protein_report.html` - Demo single protein report with basic tabs
- `single_protein_report.html` - Generated from actual pipeline (if successful)
- `single_protein_workflow.html` - Generated from workflow orchestrator (if successful)

### Multi-Protein Analysis  
- `demo_multi_protein_report.html` - Demo multi-protein report with all tabs including:
  - Multi-Protein Analysis tab (MSA, phylogeny, clustering)
  - Network Analysis tab (centrality measures, network properties)
  - Enhanced dashboard with multi-protein statistics
- `multi_protein_report.html` - Generated from actual pipeline (if successful)  
- `multi_protein_workflow.html` - Generated from workflow orchestrator (if successful)

## Key Features Tested

### Single Protein Mode
- Basic sequence analysis
- Family assignment
- Similarity search
- Standard dashboard view
- Bioinformatics links

### Multi-Protein Mode
- All single protein features PLUS:
- Multiple Sequence Alignment (MSA)
- Phylogenetic tree construction
- Conservation analysis
- Protein clustering
- Network analysis with centrality measures
- Enhanced dashboard with comparative statistics

## Validation Steps

1. Open each HTML file in a web browser
2. Verify all tabs are present and functional
3. Check that multi-protein reports have the additional "Multi-Protein Analysis" tab
4. Confirm responsive design and Bootstrap styling
5. Validate that data displays correctly in tables and metrics

## Expected Behavior

- **Single protein input**: Should show Dashboard, Protein Details, and Bioinformatics tabs
- **Multi-protein input (≥2 proteins)**: Should show Dashboard, Multi-Protein Analysis, Network Analysis, and Bioinformatics tabs
- All reports should be responsive and researcher-friendly
