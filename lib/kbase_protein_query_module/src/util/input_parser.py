"""
Enhanced Input Parser for KBase Protein Query Module

This module handles parsing and standardizing various input formats for protein analysis,
including FASTA files, UniProt identifiers, KBase workspace objects, genome sets, and more.
All inputs are standardized to ProteinRecord format for consistent pipeline processing.
"""

import os
import requests
import tempfile
import re
from typing import List, Dict, Any, Union, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class ProteinRecord:
    """Standardized data class representing a protein record for pipeline processing."""
    protein_id: str
    source: str
    sequence: str
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}
        # Calculate sequence length
        self.metadata['length'] = str(len(self.sequence))
        # Ensure protein_id is standardized
        if not self.protein_id:
            self.protein_id = f"protein_{len(self.sequence)}"
        # Clean sequence (remove whitespace, validate amino acids)
        self.sequence = self._clean_sequence(self.sequence)
    
    def _clean_sequence(self, sequence: str) -> str:
        """Clean and validate protein sequence."""
        # Remove whitespace and convert to uppercase
        cleaned = ''.join(sequence.split()).upper()
        
        # Validate amino acid characters
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(c in valid_aa for c in cleaned):
            logger.warning(f"Invalid amino acid characters found in sequence for {self.protein_id}")
            # Remove invalid characters
            cleaned = ''.join(c for c in cleaned if c in valid_aa)
        
        return cleaned


class InputParser:
    """
    Enhanced parser for various protein input formats with standardization.
    
    Supports:
    - FASTA files (local and remote)
    - UniProt identifiers (single or multiple)
    - ProteinSequenceSet objects
    - Genome references (single or genome sets)
    - FeatureSet objects
    - Single protein sequences
    - Workspace object references
    
    All inputs are standardized to ProteinRecord format for consistent pipeline processing.
    """
    
    def __init__(self, workspace_client=None):
        """Initialize the parser with optional workspace client."""
        self.workspace_client = workspace_client
        self.supported_formats = {
            'FASTA', 'Uniprot', 'ProteinSequenceSet', 
            'Genome', 'FeatureSet', 'GenomeSet', 'SingleProtein',
            'WorkspaceObject', 'MixedInput'
        }
        
        # UniProt ID patterns for validation
        self.uniprot_patterns = [
            r'^[A-Z][0-9A-Z]{5}$',  # UniProtKB/Swiss-Prot format (P12345)
            r'^[A-Z][0-9A-Z]{9}$',  # UniProtKB/TrEMBL format (A0A0A0A0A0)
            r'^[A-Z][0-9A-Z]{4}$',  # Short format (P1234)
        ]
    
    def parse_input(self, input_type: str, input_data: Union[str, List[str], Dict[str, Any]]) -> List[ProteinRecord]:
        """
        Parse and standardize input data based on the specified type.
        
        Args:
            input_type: Type of input ('FASTA', 'Uniprot', 'Genome', etc.)
            input_data: Input data (file path, URL, identifier, workspace ref, etc.)
            
        Returns:
            List of standardized ProteinRecord objects
        """
        if input_type not in self.supported_formats:
            raise ValueError(f"Unsupported input type: {input_type}")
        
        logger.info(f"Parsing input type: {input_type}")
        
        try:
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
                # Prefer instance method; fallback to module helper if unavailable
                if hasattr(self, '_parse_single_protein'):
                    return self._parse_single_protein(input_data)
                else:
                    return parse_single_protein(input_data)
            elif input_type == 'WorkspaceObject':
                return self._parse_workspace_object(input_data)
            elif input_type == 'MixedInput':
                return self._parse_mixed_input(input_data)
            else:
                raise ValueError(f"Input type {input_type} not implemented")
        except Exception as e:
            logger.error(f"Error parsing input type {input_type}: {e}")
            raise
    
    def detect_input_type(self, input_data: Union[str, List[str], Dict[str, Any]]) -> str:
        """
        Automatically detect the input type based on the data format.
        
        Args:
            input_data: Input data to analyze
            
        Returns:
            Detected input type string
        """
        if isinstance(input_data, str):
            # Check if it's a file path
            if os.path.exists(input_data) or input_data.startswith(('http://', 'https://')):
                if input_data.endswith('.fasta') or input_data.endswith('.fa'):
                    return 'FASTA'
            
            # Check if it's a UniProt ID
            if self._is_uniprot_id(input_data):
                return 'Uniprot'
            
            # Check if it's a workspace reference
            if '/' in input_data and self._is_workspace_ref(input_data):
                return 'WorkspaceObject'
            
            # Check if it's a protein sequence
            if self._is_protein_sequence(input_data):
                return 'SingleProtein'
        
        elif isinstance(input_data, list):
            # Check if all items are UniProt IDs
            if all(self._is_uniprot_id(item) for item in input_data):
                return 'Uniprot'
            # Check if all items are protein sequences
            elif all(self._is_protein_sequence(item) for item in input_data):
                return 'SingleProtein'
            # Mixed input
            else:
                return 'MixedInput'
        
        elif isinstance(input_data, dict):
            # Check for workspace object indicators
            if 'ref' in input_data or 'workspace_ref' in input_data:
                return 'WorkspaceObject'
        
        # Default to mixed input for complex cases
        return 'MixedInput'
    
    def _is_uniprot_id(self, identifier: str) -> bool:
        """Check if identifier is a valid UniProt ID."""
        if not identifier or len(identifier) < 3:
            return False
        
        for pattern in self.uniprot_patterns:
            if re.match(pattern, identifier):
                return True
        
        return False
    
    def _is_workspace_ref(self, ref: str) -> bool:
        """Check if reference is a valid workspace reference."""
        if not ref or '/' not in ref:
            return False
        
        parts = ref.split('/')
        if len(parts) != 2:
            return False
        
        try:
            int(parts[0])  # Workspace ID should be numeric
            int(parts[1])  # Object ID should be numeric
            return True
        except ValueError:
            return False
    
    def _is_protein_sequence(self, sequence: str) -> bool:
        """Check if string is a valid protein sequence."""
        if not sequence or len(sequence) < 10:  # Minimum reasonable protein length
            return False
        
        # Check if it contains mostly amino acid characters
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        sequence_upper = sequence.upper()
        aa_count = sum(1 for c in sequence_upper if c in valid_aa)
        aa_ratio = aa_count / len(sequence_upper)
        
        return aa_ratio > 0.8  # At least 80% should be valid amino acids
    
    def _parse_fasta(self, fasta_input: str) -> List[ProteinRecord]:
        """Parse FASTA input (file path, URL, or raw FASTA string)."""
        if os.path.exists(fasta_input):
            return self._parse_fasta_file(fasta_input)
        elif fasta_input.startswith(('http://', 'https://')):
            return self._parse_fasta_url(fasta_input)
        else:
            # Treat as raw FASTA content
            return self._parse_fasta_string(fasta_input)

    def _parse_fasta_string(self, fasta_text: str) -> List[ProteinRecord]:
        """Parse FASTA content from a string."""
        records = []
        current_id = None
        current_sequence = []
        for line in fasta_text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id and current_sequence:
                    sequence = ''.join(current_sequence)
                    records.append(ProteinRecord(
                        protein_id=current_id,
                        source='FASTA_STRING',
                        sequence=sequence,
                        metadata={}
                    ))
                current_id = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_id and current_sequence:
            sequence = ''.join(current_sequence)
            records.append(ProteinRecord(
                protein_id=current_id,
                source='FASTA_STRING',
                sequence=sequence,
                metadata={}
            ))
        return records
    
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
                                sequence=sequence,
                                metadata={'file_path': file_path}
                            ))
                        
                        # Start new record
                        current_id = line[1:].split()[0]  # Take first word after >
                        current_sequence = []
                    else:
                        current_sequence.append(line)
                
                # Save last record
                if current_id and current_sequence:
                    sequence = ''.join(current_sequence)
                    records.append(ProteinRecord(
                        protein_id=current_id,
                        source='FASTA',
                        sequence=sequence,
                        metadata={'file_path': file_path}
                    ))
            
            logger.info(f"Parsed {len(records)} proteins from FASTA file: {file_path}")
            return records
            
        except Exception as e:
            logger.error(f"Error parsing FASTA file {file_path}: {e}")
            raise
    
    def _parse_fasta_url(self, url: str) -> List[ProteinRecord]:
        """Parse a remote FASTA file."""
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Create temporary file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(response.text)
                temp_path = f.name
            
            try:
                records = self._parse_fasta_file(temp_path)
                # Update source to indicate URL
                for record in records:
                    record.source = 'FASTA_URL'
                    record.metadata['url'] = url
                return records
            finally:
                os.unlink(temp_path)
                
        except Exception as e:
            logger.error(f"Error parsing FASTA URL {url}: {e}")
            raise

# --- Module-level convenience functions expected by some stages/tests ---
def parse_fasta_string(fasta_text: str) -> List[ProteinRecord]:
    parser = InputParser()
    return parser._parse_fasta_string(fasta_text)


def parse_fasta_file(file_path: str) -> List[ProteinRecord]:
    parser = InputParser()
    return parser._parse_fasta_file(file_path)


def parse_single_protein(sequence: str) -> List[ProteinRecord]:
    parser = InputParser()
    if not parser._is_protein_sequence(sequence):
        raise ValueError(f"Invalid protein sequence: {sequence[:50]}...")
    rec = ProteinRecord(
        protein_id='single_protein',
        source='SingleProtein',
        sequence=sequence,
        metadata={'input_type': 'direct_sequence'}
    )
    return [rec]
    
    def _parse_uniprot_identifiers(self, identifiers: Union[str, List[str]]) -> List[ProteinRecord]:
        """Parse UniProt identifiers and fetch sequences."""
        if isinstance(identifiers, str):
            # Handle comma-separated string
            if ',' in identifiers:
                identifiers = [id.strip() for id in identifiers.split(',')]
            else:
                identifiers = [identifiers]
        
        records = []
        for uniprot_id in identifiers:
            if not self._is_uniprot_id(uniprot_id):
                logger.warning(f"Invalid UniProt ID format: {uniprot_id}")
                continue
            
            try:
                # Fetch sequence from UniProt
                sequence = self._fetch_uniprot_sequence(uniprot_id)
                if sequence:
                    records.append(ProteinRecord(
                        protein_id=uniprot_id,
                        source='UniProt',
                        sequence=sequence,
                        metadata={'uniprot_id': uniprot_id}
                    ))
                else:
                    logger.warning(f"Could not fetch sequence for UniProt ID: {uniprot_id}")
            except Exception as e:
                logger.error(f"Error fetching UniProt ID {uniprot_id}: {e}")
        
        logger.info(f"Successfully parsed {len(records)} UniProt identifiers")
        return records
    
    def _fetch_uniprot_sequence(self, uniprot_id: str) -> Optional[str]:
        """Fetch protein sequence from UniProt."""
        try:
            # Use UniProt REST API
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Parse FASTA response
            lines = response.text.strip().split('\n')
            if len(lines) < 2:
                return None
            
            # Extract sequence (all lines except header)
            sequence = ''.join(lines[1:])
            return sequence
            
        except Exception as e:
            logger.error(f"Error fetching UniProt sequence for {uniprot_id}: {e}")
            return None
    
    def _parse_protein_sequence_set(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse ProteinSequenceSet workspace object."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for ProteinSequenceSet parsing")
        
        try:
            # Get object data
            obj_data = self.workspace_client.get_objects2({
                'objects': [{'ref': ws_ref}]
            })
            
            data = obj_data['data'][0]['data']
            # Handle both 'sequences' and 'proteins' fields for compatibility
            sequences = data.get('sequences', data.get('proteins', []))
            
            records = []
            for seq_data in sequences:
                records.append(ProteinRecord(
                    protein_id=seq_data.get('id', 'unknown'),
                    source='ProteinSequenceSet',
                    sequence=seq_data.get('sequence', ''),
                    metadata={
                        'workspace_ref': ws_ref,
                        'object_type': 'ProteinSequenceSet'
                    }
                ))
            
            logger.info(f"Parsed {len(records)} proteins from ProteinSequenceSet: {ws_ref}")
            return records
            
        except Exception as e:
            logger.error(f"Error parsing ProteinSequenceSet {ws_ref}: {e}")
            raise
    
    def _parse_genome_reference(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse genome reference to extract proteins."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for genome parsing")
        
        try:
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
                                'workspace_ref': ws_ref,
                                'object_type': 'Genome',
                                'feature_type': feature.get('type', ''),
                                'location': feature.get('location', []),
                                'function': feature.get('function', '')
                            }
                        ))
            
            logger.info(f"Extracted {len(records)} protein sequences from genome: {ws_ref}")
            return records
            
        except Exception as e:
            logger.error(f"Error parsing genome {ws_ref}: {e}")
            raise
    
    def _parse_genome_set(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse GenomeSet workspace object."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for GenomeSet parsing")
        
        try:
            # Get object data
            obj_data = self.workspace_client.get_objects2({
                'objects': [{'ref': ws_ref}]
            })
            
            data = obj_data['data'][0]['data']
            genome_refs = data.get('genome_refs', [])
            
            records = []
            for genome_ref in genome_refs:
                genome_records = self._parse_genome_reference(genome_ref)
                records.extend(genome_records)
            
            logger.info(f"Extracted {len(records)} protein sequences from GenomeSet: {ws_ref}")
            return records
            
        except Exception as e:
            logger.error(f"Error parsing GenomeSet {ws_ref}: {e}")
            raise
    
    def _parse_feature_set(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse FeatureSet workspace object."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for FeatureSet parsing")
        
        try:
            # Get object data
            obj_data = self.workspace_client.get_objects2({
                'objects': [{'ref': ws_ref}]
            })
            
            data = obj_data['data'][0]['data']
            features = data.get('features', [])
            
            records = []
            for feature in features:
                if feature.get('type') == 'CDS':
                    protein_id = feature.get('id', 'unknown')
                    protein_translation = feature.get('protein_translation', '')
                    
                    if protein_translation:
                        records.append(ProteinRecord(
                            protein_id=protein_id,
                            source='FeatureSet',
                            sequence=protein_translation,
                            metadata={
                                'workspace_ref': ws_ref,
                                'object_type': 'FeatureSet',
                                'feature_type': feature.get('type', '')
                            }
                        ))
            
            logger.info(f"Extracted {len(records)} protein sequences from FeatureSet: {ws_ref}")
            return records
            
        except Exception as e:
            logger.error(f"Error parsing FeatureSet {ws_ref}: {e}")
            raise
    
    def _parse_single_protein(self, sequence: str) -> List[ProteinRecord]:
        """Parse a single protein sequence."""
        if not self._is_protein_sequence(sequence):
            raise ValueError(f"Invalid protein sequence: {sequence[:50]}...")
        
        record = ProteinRecord(
            protein_id='single_protein',
            source='SingleProtein',
            sequence=sequence,
            metadata={'input_type': 'direct_sequence'}
        )
        
        logger.info(f"Parsed single protein sequence (length: {len(sequence)})")
        return [record]
    
    def _parse_workspace_object(self, ws_ref: str) -> List[ProteinRecord]:
        """Parse workspace object by detecting its type."""
        if not self.workspace_client:
            raise ValueError("Workspace client required for workspace object parsing")
        
        try:
            # Get object info to determine type
            obj_info = self.workspace_client.get_object_info3({
                'objects': [{'ref': ws_ref}]
            })
            
            object_type = obj_info['infos'][0][2]  # Type name
            
            # Route to appropriate parser based on type
            if 'ProteinSequenceSet' in object_type:
                return self._parse_protein_sequence_set(ws_ref)
            elif 'Genome' in object_type:
                return self._parse_genome_reference(ws_ref)
            elif 'GenomeSet' in object_type:
                return self._parse_genome_set(ws_ref)
            elif 'FeatureSet' in object_type:
                return self._parse_feature_set(ws_ref)
            else:
                logger.warning(f"Unknown workspace object type: {object_type}")
                return []
                
        except Exception as e:
            logger.error(f"Error parsing workspace object {ws_ref}: {e}")
            raise
    
    def _parse_mixed_input(self, input_data: Union[List[str], Dict[str, Any]]) -> List[ProteinRecord]:
        """Parse mixed input types."""
        records = []
        
        if isinstance(input_data, list):
            for item in input_data:
                if self._is_uniprot_id(item):
                    uniprot_records = self._parse_uniprot_identifiers([item])
                    records.extend(uniprot_records)
                elif self._is_protein_sequence(item):
                    seq_records = self._parse_single_protein(item)
                    records.extend(seq_records)
                elif self._is_workspace_ref(item):
                    ws_records = self._parse_workspace_object(item)
                    records.extend(ws_records)
                else:
                    logger.warning(f"Unknown input item type: {item}")
        
        elif isinstance(input_data, dict):
            # Handle dictionary input
            if 'proteins' in input_data:
                for protein in input_data['proteins']:
                    if isinstance(protein, str):
                        if self._is_uniprot_id(protein):
                            uniprot_records = self._parse_uniprot_identifiers([protein])
                            records.extend(uniprot_records)
                        elif self._is_protein_sequence(protein):
                            seq_records = self._parse_single_protein(protein)
                            records.extend(seq_records)
                    elif isinstance(protein, dict):
                        # Handle protein dictionary
                        protein_id = protein.get('id', 'unknown')
                        sequence = protein.get('sequence', '')
                        if sequence:
                            records.append(ProteinRecord(
                                protein_id=protein_id,
                                source='MixedInput',
                                sequence=sequence,
                                metadata=protein.get('metadata', {})
                            ))
        
        logger.info(f"Parsed {len(records)} proteins from mixed input")
        return records
    
    def validate_records(self, records: List[ProteinRecord]) -> Dict[str, Any]:
        """
        Validate a list of protein records.
        
        Args:
            records: List of ProteinRecord objects
            
        Returns:
            Validation summary dictionary
        """
        validation_results = {
            'total_records': len(records),
            'valid_records': 0,
            'invalid_records': 0,
            'errors': [],
            'warnings': []
        }
        
        for record in records:
            try:
                # Check sequence length
                if len(record.sequence) < 10:
                    validation_results['warnings'].append(
                        f"Short sequence for {record.protein_id}: {len(record.sequence)} aa"
                    )
                
                # Check for invalid characters
                valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
                invalid_chars = set(record.sequence) - valid_aa
                if invalid_chars:
                    validation_results['warnings'].append(
                        f"Invalid characters in {record.protein_id}: {invalid_chars}"
                    )
                
                # Check protein ID
                if not record.protein_id or record.protein_id.strip() == '':
                    validation_results['errors'].append(
                        f"Empty protein ID for sequence: {record.sequence[:20]}..."
                    )
                    validation_results['invalid_records'] += 1
                    continue
                
                validation_results['valid_records'] += 1
                
            except Exception as e:
                validation_results['errors'].append(
                    f"Error validating record {record.protein_id}: {e}"
                )
                validation_results['invalid_records'] += 1
        
        return validation_results
