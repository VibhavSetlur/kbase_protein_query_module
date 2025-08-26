"""
Reference Data Loader Module

This module provides functionality to load and manage reference data for protein analysis,
including amino acid properties, sequence motifs, bioinformatics databases, and
physicochemical constants. All data is loaded from external JSON files for maintainability
and scalability.
"""

import json
import os
from typing import Dict, Any, Optional
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class ReferenceDataLoader:
    """
    Loader for reference data used in protein sequence analysis.
    
    This class loads all reference data from external JSON files, making the system
    more maintainable and scalable. It follows KBase SDK guidelines for data management.
    """
    
    def __init__(self, data_dir: Optional[str] = None):
        """
        Initialize the reference data loader.
        
        Args:
            data_dir: Directory containing reference data files. If None, uses default location.
        """
        if data_dir is None:
            # Default to the reference directory within this module
            current_dir = Path(__file__).parent
            self.data_dir = current_dir / "reference"
        else:
            self.data_dir = Path(data_dir)
        
        # Initialize data storage
        self._amino_acids_data = None
        self._motif_patterns_data = None
        self._bioinformatics_databases_data = None
        self._physicochemical_constants_data = None
        
        # Load all reference data
        self._load_all_reference_data()
    
    def _load_all_reference_data(self) -> None:
        """Load all reference data files."""
        try:
            self._load_amino_acids_data()
            self._load_motif_patterns_data()
            self._load_bioinformatics_databases_data()
            self._load_physicochemical_constants_data()
            logger.info("All reference data loaded successfully")
        except Exception as e:
            logger.error(f"Failed to load reference data: {e}")
            raise
    
    def _load_amino_acids_data(self) -> None:
        """Load amino acid reference data."""
        file_path = self.data_dir / "amino_acids.json"
        try:
            with open(file_path, 'r') as f:
                self._amino_acids_data = json.load(f)
            logger.debug(f"Loaded amino acids data from {file_path}")
        except Exception as e:
            logger.error(f"Failed to load amino acids data from {file_path}: {e}")
            raise
    
    def _load_motif_patterns_data(self) -> None:
        """Load sequence motif patterns data."""
        file_path = self.data_dir / "sequence_motifs.json"
        try:
            with open(file_path, 'r') as f:
                self._motif_patterns_data = json.load(f)
            logger.debug(f"Loaded motif patterns data from {file_path}")
        except Exception as e:
            logger.error(f"Failed to load motif patterns data from {file_path}: {e}")
            raise
    
    def _load_bioinformatics_databases_data(self) -> None:
        """Load bioinformatics databases data."""
        file_path = self.data_dir / "bioinformatics_databases.json"
        try:
            with open(file_path, 'r') as f:
                self._bioinformatics_databases_data = json.load(f)
            logger.debug(f"Loaded bioinformatics databases data from {file_path}")
        except Exception as e:
            logger.error(f"Failed to load bioinformatics databases data from {file_path}: {e}")
            raise
    
    def _load_physicochemical_constants_data(self) -> None:
        """Load physicochemical constants data."""
        file_path = self.data_dir / "physicochemical_constants.json"
        try:
            with open(file_path, 'r') as f:
                self._physicochemical_constants_data = json.load(f)
            logger.debug(f"Loaded physicochemical constants data from {file_path}")
        except Exception as e:
            logger.error(f"Failed to load physicochemical constants data from {file_path}: {e}")
            raise
    
    @property
    def amino_acids(self) -> Dict[str, Any]:
        """Get amino acid data."""
        if self._amino_acids_data is None:
            raise RuntimeError("Amino acids data not loaded")
        return self._amino_acids_data
    
    @property
    def motif_patterns(self) -> Dict[str, Any]:
        """Get motif patterns data."""
        if self._motif_patterns_data is None:
            raise RuntimeError("Motif patterns data not loaded")
        return self._motif_patterns_data
    
    @property
    def bioinformatics_databases(self) -> Dict[str, Any]:
        """Get bioinformatics databases data."""
        if self._bioinformatics_databases_data is None:
            raise RuntimeError("Bioinformatics databases data not loaded")
        return self._bioinformatics_databases_data
    
    @property
    def physicochemical_constants(self) -> Dict[str, Any]:
        """Get physicochemical constants data."""
        if self._physicochemical_constants_data is None:
            raise RuntimeError("Physicochemical constants data not loaded")
        return self._physicochemical_constants_data
    
    def get_amino_acid_properties(self, amino_acid: str) -> Dict[str, Any]:
        """
        Get properties for a specific amino acid.
        
        Args:
            amino_acid: Single letter amino acid code
            
        Returns:
            Dictionary containing amino acid properties
        """
        amino_acids = self.amino_acids.get('amino_acids', {})
        return amino_acids.get(amino_acid.upper(), {})
    
    def get_amino_acid_groups(self) -> Dict[str, list]:
        """
        Get amino acid groups.
        
        Returns:
            Dictionary mapping group names to lists of amino acid codes
        """
        return self.amino_acids.get('amino_acid_groups', {})
    
    def get_motif_pattern(self, motif_name: str) -> Dict[str, Any]:
        """
        Get a specific motif pattern.
        
        Args:
            motif_name: Name of the motif pattern
            
        Returns:
            Dictionary containing motif pattern information
        """
        patterns = self.motif_patterns.get('motif_patterns', {})
        return patterns.get(motif_name, {})
    
    def get_motif_categories(self) -> Dict[str, list]:
        """
        Get motif categories.
        
        Returns:
            Dictionary mapping category names to lists of motif names
        """
        return self.motif_patterns.get('motif_categories', {})
    
    def get_database_info(self, database_name: str) -> Dict[str, Any]:
        """
        Get information for a specific database.
        
        Args:
            database_name: Name of the database
            
        Returns:
            Dictionary containing database information
        """
        databases = self.bioinformatics_databases.get('databases', {})
        return databases.get(database_name, {})
    
    def get_database_categories(self) -> Dict[str, list]:
        """
        Get database categories.
        
        Returns:
            Dictionary mapping category names to lists of database names
        """
        return self.bioinformatics_databases.get('database_categories', {})
    
    def get_physicochemical_constant(self, constant_name: str) -> Any:
        """
        Get a specific physicochemical constant.
        
        Args:
            constant_name: Name of the constant (can be nested using dot notation)
            
        Returns:
            Value of the constant
        """
        constants = self.physicochemical_constants.get('physicochemical_constants', {})
        
        # Handle nested access using dot notation
        keys = constant_name.split('.')
        value = constants
        for key in keys:
            if isinstance(value, dict):
                value = value.get(key)
            else:
                return None
        
        return value
    
    def get_calculation_methods(self) -> Dict[str, str]:
        """
        Get calculation methods.
        
        Returns:
            Dictionary mapping calculation names to method descriptions
        """
        return self.physicochemical_constants.get('calculation_methods', {})
    
    def get_references(self) -> Dict[str, str]:
        """
        Get scientific references.
        
        Returns:
            Dictionary mapping reference names to citation strings
        """
        return self.physicochemical_constants.get('references', {})
    
    def reload_data(self) -> None:
        """Reload all reference data from files."""
        logger.info("Reloading all reference data")
        self._load_all_reference_data()
    
    def validate_data_integrity(self) -> bool:
        """
        Validate the integrity of loaded reference data.
        
        Returns:
            True if all data is valid, False otherwise
        """
        try:
            # Check that all required data is loaded
            required_data = [
                self._amino_acids_data,
                self._motif_patterns_data,
                self._bioinformatics_databases_data,
                self._physicochemical_constants_data
            ]
            
            if any(data is None for data in required_data):
                logger.error("Some reference data is not loaded")
                return False
            
            # Validate amino acids data structure
            amino_acids = self.amino_acids.get('amino_acids', {})
            if not amino_acids:
                logger.error("Amino acids data is empty")
                return False
            
            # Validate motif patterns data structure
            motif_patterns = self.motif_patterns.get('motif_patterns', {})
            if not motif_patterns:
                logger.error("Motif patterns data is empty")
                return False
            
            # Validate bioinformatics databases data structure
            databases = self.bioinformatics_databases.get('databases', {})
            if not databases:
                logger.error("Bioinformatics databases data is empty")
                return False
            
            # Validate physicochemical constants data structure
            constants = self.physicochemical_constants.get('physicochemical_constants', {})
            if not constants:
                logger.error("Physicochemical constants data is empty")
                return False
            
            logger.info("All reference data validation passed")
            return True
            
        except Exception as e:
            logger.error(f"Data integrity validation failed: {e}")
            return False
    
    def get_data_summary(self) -> Dict[str, Any]:
        """
        Get a summary of loaded reference data.
        
        Returns:
            Dictionary containing summary information about loaded data
        """
        summary = {
            'amino_acids': {
                'count': len(self.amino_acids.get('amino_acids', {})),
                'groups_count': len(self.amino_acids.get('amino_acid_groups', {}))
            },
            'motif_patterns': {
                'count': len(self.motif_patterns.get('motif_patterns', {})),
                'categories_count': len(self.motif_patterns.get('motif_categories', {}))
            },
            'bioinformatics_databases': {
                'count': len(self.bioinformatics_databases.get('databases', {})),
                'categories_count': len(self.bioinformatics_databases.get('database_categories', {}))
            },
            'physicochemical_constants': {
                'calculation_methods_count': len(self.physicochemical_constants.get('calculation_methods', {})),
                'references_count': len(self.physicochemical_constants.get('references', {}))
            }
        }
        
        return summary
