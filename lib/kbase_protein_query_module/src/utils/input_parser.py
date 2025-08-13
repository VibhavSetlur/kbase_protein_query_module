"""
Input Parser for KBase Protein Query Module

This module handles parsing various input formats for protein analysis,
including FASTA files, UniProt identifiers, and KBase workspace objects.
"""

import os
import requests
import tempfile
from typing import List, Dict, Any, Union, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class ProteinRecord:
    """Data class representing a protein record."""
    protein_id: str
    source: str
    sequence: str
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}
        # Calculate sequence length
        self.metadata['length'] = str(len(self.sequence))


class InputParser:
    """
    Parser for various protein input formats.
    
    Supports:
    - FASTA files (local and remote)
    - UniProt identifiers
    - ProteinSequenceSet objects
    - Genome references
    - Single protein sequences
    """
    
    def __init__(self, workspace_client=None):
        """Initialize the parser with optional workspace client."""
        self.workspace_client = workspace_client
        self.supported_formats = {
            'FASTA', 'Uniprot', 'ProteinSequenceSet', 
            'Genome', 'FeatureSet', 'GenomeSet', 'SingleProtein'
        }
    
    def parse_input(self, input_type: str, input_data: Union[str, List[str]]) -> List[ProteinRecord]:
        """
        Parse input data based on the specified type.
        
        Args:
            input_type: Type of input ('FASTA', 'Uniprot', etc.)
            input_data: Input data (file path, URL, identifier, etc.)
            
        Returns:
            List of ProteinRecord objects
        """
        if input_type not in self.supported_formats:
            raise ValueError(f"Unsupported input type: {input_type}")
        
        if input_type == 'FASTA':
            return self._parse_fasta(input_data)
        elif input_type == 'Uniprot':
            return self._parse_uniprot_identifiers(input_data)
        elif input_type == 'ProteinSequenceSet':
            return self._parse_protein_sequence_set(input_data)
        elif input_type == 'Genome':
            return self._parse_genome_reference(input_data)
        elif input_type == 'FeatureSet':
            return self._parse_feature_set(input_data)
        elif input_type == 'GenomeSet':
            return self._parse_genome_set(input_data)
        elif input_type == 'SingleProtein':
            return self._parse_single_protein(input_data)
        else:
            raise ValueError(f"Input type {input_type} not implemented")
    
    def _parse_fasta(self, fasta_input: str) -> List[ProteinRecord]:
        """Parse FASTA input (file path or URL)."""
        if os.path.exists(fasta_input):
            return self._parse_fasta_file(fasta_input)
        elif fasta_input.startswith(('http://', 'https://')):
            return self._parse_fasta_url(fasta_input)
        else:
            raise FileNotFoundError(f"FASTA file not found: {fasta_input}")
    
    def _parse_fasta_file(self, file_path: str) -> List[ProteinRecord]:
        """Parse a local FASTA file."""
        records = []
        current_id = None
        current_sequence = []
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous record
                        if current_id and current_sequence:
                            sequence = ''.join(current_sequence)
                            records.append(ProteinRecord(
                                protein_id=current_id,
                                source='FASTA',
                                sequence=sequence
                            ))
                        
                        # Start new record
                        current_id = line[1:].split()[0]  # Take first word after >
                        current_sequence = []
                    elif line and current_id:
                        current_sequence.append(line)
                
                # Save last record
                if current_id and current_sequence:
                    sequence = ''.join(current_sequence)
                    records.append(ProteinRecord(
                        protein_id=current_id,
                        source='FASTA',
                        sequence=sequence
                    ))
        except Exception as e:
            logger.error(f"Error parsing FASTA file {file_path}: {e}")
            raise
        
        return records
    
    def _parse_fasta_url(self, url: str) -> List[ProteinRecord]:
        """Parse FASTA from URL."""
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Create temporary file and parse
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(response.text)
                temp_path = f.name
            
            try:
                records = self._parse_fasta_file(temp_path)
                return records
            finally:
                os.unlink(temp_path)
        except Exception as e:
            logger.error(f"Error fetching FASTA from URL {url}: {e}")
            raise
    
    def _parse_uniprot_identifiers(self, identifiers: str) -> List[ProteinRecord]:
        """Parse UniProt identifiers."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for UniProt parsing")
        
        id_list = [id.strip() for id in identifiers.split(',')]
        records = []
        
        for uniprot_id in id_list:
            try:
                # Fetch from UniProt API
                url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                
                # Parse the single FASTA entry
                lines = response.text.strip().split('\n')
                if len(lines) >= 2:
                    header = lines[0]
                    sequence = ''.join(lines[1:])
                    
                    # Extract protein ID from header
                    protein_id = header.split('|')[1] if '|' in header else uniprot_id
                    
                    records.append(ProteinRecord(
                        protein_id=protein_id,
                        source='Uniprot',
                        sequence=sequence
                    ))
            except Exception as e:
                logger.warning(f"Failed to fetch UniProt ID {uniprot_id}: {e}")
                continue
        
        return records
    
    def _parse_protein_sequence_set(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse ProteinSequenceSet from workspace."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for ProteinSequenceSet parsing")
        
        try:
            # Get object info
            obj_info = self.workspace_client.get_object_info3({
                'objects': [{'ref': ws_ref}]
            })
            
            # Get object data
            obj_data = self.workspace_client.get_objects2({
                'objects': [{'ref': ws_ref}]
            })
            
            records = []
            proteins = obj_data['data'][0]['data'].get('proteins', [])
            
            for protein in proteins:
                records.append(ProteinRecord(
                    protein_id=protein.get('id', 'unknown'),
                    source='ProteinSequenceSet',
                    sequence=protein.get('sequence', ''),
                    metadata={
                        'description': protein.get('description', ''),
                        'length': str(len(protein.get('sequence', '')))
                    }
                ))
            
            return records
        except Exception as e:
            logger.error(f"Error parsing ProteinSequenceSet {ws_ref}: {e}")
            raise
    
    def _parse_genome_reference(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse genome reference to extract proteins."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for genome parsing")
        
        try:
            # Get object info
            obj_info = self.workspace_client.get_object_info3({
                'objects': [{'ref': ws_ref}]
            })
            
            # Get object data
            obj_data = self.workspace_client.get_objects2({
                'objects': [{'ref': ws_ref}]
            })
            
            genome_data = obj_data['data'][0]['data']
            features = genome_data.get('features', [])
            
            records = []
            for feature in features:
                # Extract CDS features (coding sequences)
                if feature.get('type') == 'CDS':
                    protein_id = feature.get('id', 'unknown')
                    
                    # Extract protein sequence
                    protein_translation = feature.get('protein_translation', '')
                    if protein_translation:
                        records.append(ProteinRecord(
                            protein_id=protein_id,
                            source='Genome',
                            sequence=protein_translation,
                            metadata={
                                'feature_type': feature.get('type', ''),
                                'location': feature.get('location', []),
                                'function': feature.get('function', ''),
                                'length': str(len(protein_translation))
                            }
                        ))
            
            logger.info(f"Extracted {len(records)} protein sequences from genome")
            return records
            
        except Exception as e:
            logger.error(f"Error parsing genome {ws_ref}: {e}")
            raise
    
    def _parse_single_protein(self, sequence: str) -> List[ProteinRecord]:
        """Parse a single protein sequence."""
        return [ProteinRecord(
            protein_id='single_protein',
            source='SingleProtein',
            sequence=sequence
        )]
    
    def _parse_feature_set(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse feature set to extract proteins."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for feature set parsing")
        
        try:
            # Get object data
            obj_data = self.workspace_client.get_objects2({
                'objects': [{'ref': ws_ref}]
            })
            
            feature_set_data = obj_data['data'][0]['data']
            features = feature_set_data.get('elements', [])
            
            records = []
            for feature in features:
                # Extract protein sequences from features
                protein_translation = feature.get('protein_translation', '')
                if protein_translation:
                    protein_id = feature.get('id', 'unknown')
                    records.append(ProteinRecord(
                        protein_id=protein_id,
                        source='FeatureSet',
                        sequence=protein_translation,
                        metadata={
                            'feature_type': feature.get('type', ''),
                            'function': feature.get('function', ''),
                            'length': str(len(protein_translation))
                        }
                    ))
            
            logger.info(f"Extracted {len(records)} protein sequences from feature set")
            return records
            
        except Exception as e:
            logger.error(f"Error parsing feature set {ws_ref}: {e}")
            raise
    
    def _parse_genome_set(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse genome set to extract proteins from multiple genomes."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for genome set parsing")
        
        try:
            # Get object data
            obj_data = self.workspace_client.get_objects2({
                'objects': [{'ref': ws_ref}]
            })
            
            genome_set_data = obj_data['data'][0]['data']
            genome_refs = genome_set_data.get('elements', [])
            
            all_records = []
            for genome_ref in genome_refs:
                try:
                    # Parse each genome in the set
                    genome_records = self._parse_genome_reference(genome_ref)
                    all_records.extend(genome_records)
                except Exception as e:
                    logger.warning(f"Failed to parse genome {genome_ref}: {e}")
                    continue
            
            logger.info(f"Extracted {len(all_records)} protein sequences from genome set")
            return all_records
            
        except Exception as e:
            logger.error(f"Error parsing genome set {ws_ref}: {e}")
            raise
    
    def validate_records(self, records: List[ProteinRecord]) -> Dict[str, Any]:
        """
        Validate protein records.
        
        Returns:
            Dictionary with validation results
        """
        validation_results = {
            'total_records': len(records),
            'valid_records': 0,
            'invalid_records': 0,
            'errors': []
        }
        
        for record in records:
            if self._is_valid_sequence(record.sequence):
                validation_results['valid_records'] += 1
            else:
                validation_results['invalid_records'] += 1
                validation_results['errors'].append({
                    'protein_id': record.protein_id,
                    'error': 'Invalid sequence'
                })
        
        return validation_results
    
    def _is_valid_sequence(self, sequence: str) -> bool:
        """Check if a protein sequence is valid."""
        if not sequence:
            return False
        
        # Check for valid amino acid characters
        valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
        sequence_chars = set(sequence.upper())
        
        if not sequence_chars.issubset(valid_chars):
            return False
        
        # Check reasonable length (1-50000 amino acids)
        if len(sequence) < 1 or len(sequence) > 50000:
            return False
        
        return True
