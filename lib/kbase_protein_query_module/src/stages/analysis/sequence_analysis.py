"""
Refactored Protein Sequence Analyzer Stage

This stage provides comprehensive analysis of protein sequences using external reference data
for maintainability and scalability. All reference data is loaded from JSON files and the
stage follows the BaseStage/StageResult interface used across the pipeline.

Features:
- Amino acid composition and statistics
- Physicochemical properties (verified calculations)
- Secondary structure predictions
- Functional domain analysis
- Bioinformatics database links
- Sequence motifs and patterns
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
import logging
from collections import Counter
import re
import requests
from urllib.parse import quote
import json

from ..base_stage import BaseStage, StageResult
from ...data import ReferenceDataLoader

logger = logging.getLogger(__name__)


class SequenceAnalysisStage(BaseStage):
    """
    Pipeline stage for comprehensive protein sequence analysis with external references.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None, data_loader: Optional[ReferenceDataLoader] = None):
        super().__init__(config or {})
        if data_loader is None:
            self.data_loader = ReferenceDataLoader()
        else:
            self.data_loader = data_loader
        if not self.data_loader.validate_data_integrity():
            raise RuntimeError("Reference data integrity validation failed")
        self.include_bioinformatics = self.config.get('include_bioinformatics', True)
        self.include_structure_prediction = self.config.get('include_structure_prediction', True)
        logger.info("SequenceAnalysisStage initialized with external reference data")

    def get_stage_name(self) -> str:
        return "sequence_analysis"

    def get_required_inputs(self) -> List[str]:
        return ['protein_records']

    def get_optional_inputs(self) -> List[str]:
        return ['analysis_config']

    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        # Accept either direct 'protein_records' or via 'extracted_data'
        records = input_data.get('protein_records')
        if records is None and isinstance(input_data.get('extracted_data'), dict):
            records = input_data['extracted_data'].get('protein_records')
        return isinstance(records, list) and len(records) > 0

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
        results: Dict[str, Any] = {}
        total = 0
        success = 0
        failed = 0
        property_distributions = {
            'molecular_weights': [],
            'isoelectric_points': [],
            'instability_indices': [],
            'aliphatic_indices': [],
            'extinction_coefficients_280nm': []
        }
        motif_statistics = {
            'total_motifs': 0,
            'motif_types': {},
            'proteins_with_motifs': 0
        }

        # Get protein records from direct field or extracted_data alias
        records = input_data.get('protein_records')
        if records is None and isinstance(input_data.get('extracted_data'), dict):
            records = input_data['extracted_data'].get('protein_records', [])

        for record in records:
            total += 1
            try:
                analysis = self.analyze_sequence(record.sequence, record.protein_id)
                results[record.protein_id] = analysis
                props = analysis.get('physicochemical_properties', {})
                if props:
                    property_distributions['molecular_weights'].append(props.get('molecular_weight', 0))
                    property_distributions['isoelectric_points'].append(props.get('isoelectric_point', 0))
                    property_distributions['instability_indices'].append(props.get('instability_index', 0))
                    property_distributions['aliphatic_indices'].append(props.get('aliphatic_index', 0))
                    property_distributions['extinction_coefficients_280nm'].append(props.get('extinction_coefficient_280nm', 0))
                motifs = analysis.get('sequence_motifs', {})
                if motifs:
                    motif_statistics['proteins_with_motifs'] += 1
                    counts = sum(len(v or []) for v in motifs.values())
                    motif_statistics['total_motifs'] += counts
                    for mtype, lst in motifs.items():
                        motif_statistics['motif_types'][mtype] = motif_statistics['motif_types'].get(mtype, 0) + len(lst or [])
                success += 1
            except Exception as e:
                results[getattr(record, 'protein_id', f'protein_{total}')] = {
                    'status': 'error',
                    'error': str(e)
                }
                failed += 1

        # summarize distributions
        import numpy as _np
        for k, values in property_distributions.items():
            if values:
                property_distributions[k] = {
                    'mean': float(_np.mean(values)),
                    'median': float(_np.median(values)),
                    'min': float(_np.min(values)),
                    'max': float(_np.max(values)),
                    'std': float(_np.std(values))
                }
            else:
                property_distributions[k] = {}

        return StageResult(
            success=success > 0,
            output_data={
                'sequence_analyses': results,
                'analysis_stats': {
                    'total_analyses': total,
                    'successful_analyses': success,
                    'failed_analyses': failed,
                    'property_distributions': property_distributions,
                    'motif_statistics': motif_statistics
                }
            },
            metadata={
                'include_bioinformatics': self.include_bioinformatics,
                'include_structure_prediction': self.include_structure_prediction
            },
            execution_time=0.0
        )

    def analyze_sequence(self, sequence: str, protein_id: str = "UNKNOWN") -> Dict[str, Any]:
        """
        Perform comprehensive sequence analysis.
        
        Args:
            sequence: Protein sequence string
            protein_id: Protein identifier
            
        Returns:
            Dictionary containing all analysis results
        """
        if not sequence or not isinstance(sequence, str):
            raise ValueError("Sequence must be a non-empty string")
        
        sequence = sequence.upper()
        
        ss_pred = self._predict_secondary_structure(sequence)

        results = {
            'protein_id': protein_id,
            'sequence': sequence,
            'length': len(sequence),
            'amino_acid_composition': self._analyze_amino_acid_composition(sequence),
            'physicochemical_properties': self._analyze_physicochemical_properties(sequence),
            'secondary_structure_prediction': ss_pred,
            'sequence_motifs': self._find_sequence_motifs(sequence),
            'bioinformatics_links': self._generate_bioinformatics_links(sequence, protein_id),
            'statistics': self._calculate_statistics(sequence)
        }
        
        return results

    def _analyze_amino_acid_composition(self, sequence: str) -> Dict[str, Any]:
        """Analyze amino acid composition using external reference data."""
        aa_counts = Counter(sequence)
        total_aa = len(sequence)
        
        # Get amino acid data from reference files
        amino_acids_data = self.data_loader.amino_acids.get('amino_acids', {})
        amino_acid_groups = self.data_loader.get_amino_acid_groups()
        
        composition = {}
        for aa_code in amino_acids_data:
            count = aa_counts.get(aa_code, 0)
            percentage = (count / total_aa) * 100 if total_aa > 0 else 0
            aa_info = amino_acids_data[aa_code]
            composition[aa_code] = {
                'count': count,
                'percentage': round(percentage, 2),
                'name': aa_info.get('name', 'Unknown'),
                'three_letter': aa_info.get('three_letter', ''),
                'symbol': aa_info.get('symbol', aa_code)
            }
        
        # Group by properties using reference data
        group_composition = {}
        for group_name, group_aa in amino_acid_groups.items():
            group_count = sum(aa_counts.get(aa, 0) for aa in group_aa)
            group_percentage = (group_count / total_aa) * 100 if total_aa > 0 else 0
            group_composition[group_name] = {
                'count': group_count,
                'percentage': round(group_percentage, 2),
                'amino_acids': group_aa
            }
        
        return {
            'individual': composition,
            'groups': group_composition,
            'total_amino_acids': total_aa
        }

    def _analyze_physicochemical_properties(self, sequence: str) -> Dict[str, Any]:
        """Analyze physicochemical properties using external reference data."""
        total_mw = 0
        total_charge = 0
        hydrophobicity_scores = []
        
        # Get amino acid properties from reference data
        amino_acids_data = self.data_loader.amino_acids.get('amino_acids', {})
        physicochemical_constants = self.data_loader.physicochemical_constants.get('physicochemical_constants', {})
        
        # Count amino acids for charge calculation
        aa_counts = Counter(sequence)
        
        for aa in sequence:
            if aa in amino_acids_data:
                aa_info = amino_acids_data[aa]
                total_mw += aa_info.get('molecular_weight', 0)
                hydrophobicity_scores.append(aa_info.get('hydrophobicity', 0))
        
        # Calculate net charge at physiological pH
        ph_physiological = physicochemical_constants.get('ph_physiological', 7.0)
        net_charge = 0
        for aa, count in aa_counts.items():
            if aa in amino_acids_data:
                pka = amino_acids_data[aa].get('pka', 7.0)
                # Calculate charge contribution at physiological pH
                if pka < ph_physiological:  # Acidic side chains
                    charge_contribution = -count * (1 / (1 + 10**(ph_physiological - pka)))
                elif pka > ph_physiological:  # Basic side chains
                    charge_contribution = count * (1 / (1 + 10**(pka - ph_physiological)))
                else:
                    charge_contribution = 0
                net_charge += charge_contribution
        
        # Calculate averages
        avg_hydrophobicity = np.mean(hydrophobicity_scores) if hydrophobicity_scores else 0
        
        # Calculate isoelectric point using iterative method
        pi_estimate = self._calculate_isoelectric_point(sequence)
        
        # Calculate extinction coefficient at 280nm (for Trp, Tyr, Cys)
        extinction_coefficient = self._calculate_extinction_coefficient(sequence)
        
        # Calculate instability index
        instability_index = self._calculate_instability_index(sequence)
        
        # Calculate aliphatic index
        aliphatic_index = self._calculate_aliphatic_index(sequence)
        
        # Calculate gravy (Grand average of hydropathy)
        gravy = avg_hydrophobicity
        
        return {
            'molecular_weight': round(total_mw, 2),
            'net_charge_ph7': round(net_charge, 2),
            'average_hydrophobicity': round(avg_hydrophobicity, 3),
            'isoelectric_point': round(pi_estimate, 2),
            'extinction_coefficient_280nm': round(extinction_coefficient, 0),
            'instability_index': round(instability_index, 2),
            'aliphatic_index': round(aliphatic_index, 2),
            'gravy': round(gravy, 3),
            'amino_acid_properties': {aa: amino_acids_data.get(aa, {}) for aa in sequence}
        }

    def _predict_secondary_structure(self, sequence: str) -> Dict[str, Any]:
        """Predict secondary structure using Chou-Fasman parameters from reference data."""
        if len(sequence) < 3:
            return {
                'helix': 0,
                'sheet': 0,
                'turn': 0,
                'coil': 100,
                'dominant_structure': 'coil',
                'dominant_score': 100
            }
        
        # Get secondary structure parameters from reference data
        amino_acids_data = self.data_loader.amino_acids.get('amino_acids', {})
        physicochemical_constants = self.data_loader.physicochemical_constants.get('physicochemical_constants', {})
        ss_parameters = physicochemical_constants.get('secondary_structure_parameters', {})
        chou_fasman = ss_parameters.get('chou_fasman', {})
        thresholds = ss_parameters.get('thresholds', {'helix': 1.0, 'sheet': 1.0, 'turn': 1.0})
        
        # Calculate average preferences
        helix_pref = []
        sheet_pref = []
        turn_pref = []
        
        for aa in sequence:
            if aa in amino_acids_data and aa in chou_fasman:
                ss_prefs = chou_fasman[aa]
                helix_pref.append(ss_prefs.get('helix', 1.0))
                sheet_pref.append(ss_prefs.get('sheet', 1.0))
                turn_pref.append(ss_prefs.get('turn', 1.0))
        
        avg_helix = np.mean(helix_pref) if helix_pref else 1.0
        avg_sheet = np.mean(sheet_pref) if sheet_pref else 1.0
        avg_turn = np.mean(turn_pref) if turn_pref else 1.0
        
        # Determine dominant structure
        structures = {
            'helix': avg_helix,
            'sheet': avg_sheet,
            'turn': avg_turn
        }
        
        dominant_structure = max(structures, key=structures.get)
        dominant_score = structures[dominant_structure]
        
        # Calculate percentages (simplified)
        total_pref = avg_helix + avg_sheet + avg_turn
        if total_pref > 0:
            helix_pct = (avg_helix / total_pref) * 100
            sheet_pct = (avg_sheet / total_pref) * 100
            turn_pct = (avg_turn / total_pref) * 100
            coil_pct = max(0, 100 - helix_pct - sheet_pct - turn_pct)
        else:
            helix_pct = sheet_pct = turn_pct = 0
            coil_pct = 100
        
        return {
            'helix': round(helix_pct, 1),
            'sheet': round(sheet_pct, 1),
            'turn': round(turn_pct, 1),
            'coil': round(coil_pct, 1),
            'dominant_structure': dominant_structure,
            'dominant_score': round(dominant_score, 3)
        }

    def _find_sequence_motifs(self, sequence: str) -> Dict[str, Any]:
        """Find sequence motifs and patterns using external reference data."""
        motifs = {}
        
        # Get motif patterns from reference data
        motif_patterns_data = self.data_loader.motif_patterns.get('motif_patterns', {})
        
        for motif_name, motif_info in motif_patterns_data.items():
            pattern = motif_info.get('pattern', '')
            if pattern:
                matches = []
                for match in re.finditer(pattern, sequence):
                    matches.append({
                        'start': match.start(),
                        'end': match.end(),
                        'sequence': match.group(),
                        'position': f"{match.start() + 1}-{match.end()}",
                        'description': motif_info.get('description', ''),
                        'source': motif_info.get('source', ''),
                        'biological_significance': motif_info.get('biological_significance', '')
                    })
                motifs[motif_name] = matches
        
        return motifs

    def _generate_bioinformatics_links(self, sequence: str, protein_id: str) -> Dict[str, str]:
        """Generate links to various bioinformatics databases using external reference data."""
        encoded_sequence = quote(sequence)
        encoded_protein_id = quote(protein_id)
        
        # Get database information from reference data
        databases_data = self.data_loader.bioinformatics_databases.get('databases', {})
        
        links = {}
        for db_key, db_info in databases_data.items():
            base_url = db_info.get('base_url', '')
            search_type = db_info.get('search_type', 'protein_id')
            
            if search_type == 'protein_id':
                links[db_key] = f"{base_url}{encoded_protein_id}"
            elif search_type == 'sequence_search':
                # For sequence search tools, we'll use a generic approach
                links[db_key] = base_url
            else:
                links[db_key] = base_url
        
        return links

    def _calculate_statistics(self, sequence: str) -> Dict[str, Any]:
        """Calculate general sequence statistics."""
        return {
            'length': len(sequence),
            'unique_amino_acids': len(set(sequence)),
            'most_common_aa': Counter(sequence).most_common(1)[0] if sequence else None,
            'gc_content': 0,  # Not applicable for proteins
            'n_content': sequence.count('N') / len(sequence) * 100 if sequence else 0,
            'c_content': sequence.count('C') / len(sequence) * 100 if sequence else 0
        }

    def _calculate_isoelectric_point(self, sequence: str) -> float:
        """Calculate isoelectric point using iterative method."""
        # Simplified pI calculation based on amino acid composition
        aa_counts = Counter(sequence)
        
        # Get pKa values from reference data
        amino_acids_data = self.data_loader.amino_acids.get('amino_acids', {})
        
        # Calculate charge at different pH values
        ph_values = np.arange(3.0, 12.0, 0.1)
        charges = []
        
        for ph in ph_values:
            charge = 0
            # N-terminus
            charge += 1 / (1 + 10**(ph - 8.0))
            # C-terminus  
            charge -= 1 / (1 + 10**(3.1 - ph))
            # Side chains
            for aa, count in aa_counts.items():
                if aa in amino_acids_data:
                    pka = amino_acids_data[aa].get('pka', 7.0)
                    if aa in ['D', 'E']:  # Acidic
                        charge += count * (1 / (1 + 10**(ph - pka)))
                    elif aa in ['K', 'R', 'H']:  # Basic
                        charge -= count * (1 / (1 + 10**(pka - ph)))
            charges.append(abs(charge))
        
        # Find pH where charge is closest to zero
        min_charge_idx = np.argmin(charges)
        pi = ph_values[min_charge_idx]
        
        return max(3.0, min(11.0, pi))  # Clamp between 3-11

    def _calculate_extinction_coefficient(self, sequence: str) -> float:
        """Calculate extinction coefficient at 280nm."""
        # Get extinction coefficient parameters from reference data
        physicochemical_constants = self.data_loader.physicochemical_constants.get('physicochemical_constants', {})
        extinction_params = physicochemical_constants.get('extinction_coefficient_parameters', {})
        
        trp_count = sequence.count('W')
        tyr_count = sequence.count('Y')
        cys_count = sequence.count('C')
        
        # Use parameters from reference data
        trp_coeff = extinction_params.get('tryptophan', 5500)
        tyr_coeff = extinction_params.get('tyrosine', 1490)
        cys_coeff = extinction_params.get('cysteine', 125)
        
        # Assume half of Cys are in disulfide bonds
        cys_contribution = cys_count * (cys_coeff / 2)
        
        extinction = trp_count * trp_coeff + tyr_count * tyr_coeff + cys_contribution
        return extinction

    def _calculate_instability_index(self, sequence: str) -> float:
        """Calculate instability index using external reference data."""
        if len(sequence) < 2:
            return 0.0
        
        # Get instability index parameters from reference data
        physicochemical_constants = self.data_loader.physicochemical_constants.get('physicochemical_constants', {})
        instability_params = physicochemical_constants.get('instability_index_parameters', {})
        dipeptide_values = instability_params.get('dipeptide_values', {})
        
        # Calculate instability index
        total_score = 0
        for i in range(len(sequence) - 1):
            dipeptide = sequence[i:i+2]
            score = dipeptide_values.get(dipeptide, 0.0)
            total_score += score
        
        # Normalize by sequence length
        instability_index = (10 / len(sequence)) * total_score
        return instability_index

    def _calculate_aliphatic_index(self, sequence: str) -> float:
        """Calculate aliphatic index using external reference data."""
        if len(sequence) == 0:
            return 0.0
        
        # Get aliphatic index parameters from reference data
        physicochemical_constants = self.data_loader.physicochemical_constants.get('physicochemical_constants', {})
        aliphatic_params = physicochemical_constants.get('aliphatic_index_parameters', {})
        aliphatic_aa = aliphatic_params.get('aliphatic_amino_acids', ['A', 'I', 'L', 'M', 'V'])
        weight_factors = aliphatic_params.get('weight_factors', {})
        
        # Calculate aliphatic index
        aliphatic_count = sum(sequence.count(aa) for aa in aliphatic_aa)
        aliphatic_index = (aliphatic_count / len(sequence)) * 100
        
        return aliphatic_index

    def get_analysis_summary(self, analysis_results: Dict[str, Any]) -> str:
        """Generate a human-readable summary of analysis results."""
        report = []
        
        # Basic info
        report.append(f"Protein ID: {analysis_results['protein_id']}")
        report.append(f"Sequence Length: {analysis_results['length']} amino acids")
        report.append("")
        
        # Amino acid composition
        comp = analysis_results['amino_acid_composition']
        report.append("Amino Acid Composition:")
        report.append(f"  Total: {comp['total_amino_acids']} amino acids")
        report.append(f"  Unique: {len(comp['individual'])} different types")
        
        # Most abundant amino acids
        sorted_aa = sorted(comp['individual'].items(), key=lambda x: x[1]['count'], reverse=True)
        report.append("  Most abundant:")
        for aa, info in sorted_aa[:5]:
            report.append(f"    {aa} ({info['name']}): {info['count']} ({info['percentage']}%)")
        report.append("")
        
        # Physicochemical properties
        props = analysis_results['physicochemical_properties']
        report.append("Physicochemical Properties:")
        report.append(f"  Molecular Weight: {props['molecular_weight']} Da")
        report.append(f"  Isoelectric Point: {props['isoelectric_point']}")
        report.append(f"  Net Charge (pH 7): {props['net_charge_ph7']}")
        report.append(f"  Average Hydrophobicity: {props['average_hydrophobicity']}")
        report.append(f"  Instability Index: {props['instability_index']}")
        report.append(f"  Aliphatic Index: {props['aliphatic_index']}")
        report.append("")
        
        # Secondary structure
        ss = analysis_results['secondary_structure_prediction']
        report.append("Secondary Structure Prediction:")
        report.append(f"  Helix: {ss['helix']}%")
        report.append(f"  Sheet: {ss['sheet']}%")
        report.append(f"  Turn: {ss['turn']}%")
        report.append(f"  Coil: {ss['coil']}%")
        report.append(f"  Dominant Structure: {ss['dominant_structure']} ({ss['dominant_score']:.3f})")
        report.append("")
        
        # Motifs
        motifs = analysis_results['sequence_motifs']
        report.append("Sequence Motifs:")
        for motif_type, motif_list in motifs.items():
            if motif_list:
                report.append(f"  {motif_type.replace('_', ' ').title()}: {len(motif_list)} found")
        
        return "\n".join(report)


