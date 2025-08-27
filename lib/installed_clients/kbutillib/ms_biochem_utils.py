"""ModelSEED Biochemistry Database utilities for compound and reaction searches."""

import os
import re
import subprocess
import string
from typing import Any, Optional, Dict, List
from collections import defaultdict

from .shared_env_utils import SharedEnvUtils


class MSBiochemUtils(SharedEnvUtils):
    """Utilities for searching the ModelSEED Biochemistry Database.

    Provides methods for searching compounds and reactions by name, equation,
    formula, and other properties. Automatically manages the ModelSEED Database
    git repository and provides convenient search interfaces.
    """

    def __init__(
        self,
        modelseed_db_path: Optional[str] = None,
        auto_download: bool = True,
        **kwargs: Any,
    ) -> None:
        """Initialize ModelSEED Biochemistry utilities.

        Args:
            modelseed_db_path: Path to ModelSEED database directory
            auto_download: Whether to automatically download database if not found
            **kwargs: Additional keyword arguments passed to SharedEnvironment
        """
        super().__init__(**kwargs)

        # Set default database path to dependencies directory
        if modelseed_db_path is None:
            # First check if user has configured a path
            config_path = self.get_config(
                "ModelSEEDBiochem", "modelseedbiochem_directory"
            )
            if config_path:
                modelseed_db_path = config_path
            else:
                # Use dependencies directory by default
                from pathlib import Path

                dependencies_dir = Path(__file__).parent / "dependencies"
                modelseed_db_path = str(dependencies_dir / "ModelSEEDDatabase")

        self.modelseed_db_path = modelseed_db_path
        self.auto_download = auto_download
        self._biochem_db = None
        self._identifier_hash = None
        self._structure_hash = None
        self._element_hashes = None

        # Initialize database
        self._ensure_database_available()

    @property
    def biochem_db(self) -> Any:
        """Get the ModelSEED biochemistry database instance."""
        if self._biochem_db is None:
            self._ensure_database_available()
        return self._biochem_db
    
    @property
    def identifier_hash(self) -> Dict:
        """Index the biochemistry compounds by their identifiers."""
        if self._identifier_hash is None:
            self._identifier_hash = {}
            for cpd in self._biochem_db.compounds:
                if cpd.is_obsolete == False:
                    item = self._standardize_string(cpd.id)
                    self._identifier_hash.setdefault(item, {"type": "msid", "ids": []})
                    self._identifier_hash[item]["ids"].append(cpd.id)
                    item = self._standardize_string(cpd.name)
                    self._identifier_hash.setdefault(item, {"type": "name", "ids": []})
                    self._identifier_hash[item]["ids"].append(cpd.id)   
                    for name in cpd.names:
                        name = self._standardize_string(str(name))
                        self._identifier_hash.setdefault(name, {"type": "synonym", "ids": []})
                        self._identifier_hash[name]["ids"].append(cpd.id)
                    for anno_type in cpd.annotation:
                        if isinstance(cpd.annotation[anno_type], set):
                            for item in cpd.annotation[anno_type]:
                                item = self._standardize_string(item)
                                self._identifier_hash.setdefault(item, {"type": anno_type, "ids": []})
                                self._identifier_hash[item]["ids"].append(cpd.id)
        return self._identifier_hash

    @property
    def structure_hash(self) -> Dict:
        """Index the biochemistry compounds by their structures."""
        if self._structure_hash is None:
            self._structure_hash = {}
            for cpd in self._biochem_db.compounds:
                if cpd.is_obsolete == False:
                    #All structures are unified in a single hash, and we store the type and ids for each struct in the hash values
                    if "InChI" in cpd.annotation:
                        self._structure_hash.setdefault(cpd.annotation["InChI"], {"type": "InChI", "ids": []})
                        self._structure_hash[cpd.annotation["InChI"]]["ids"].append(cpd.id)
                    if "SMILE" in cpd.annotation:
                        self._structure_hash.setdefault(cpd.annotation["SMILE"], {"type": "SMILE", "ids": []})
                        self._structure_hash[cpd.annotation["SMILE"]]["ids"].append(cpd.id)
                    if "InChIKey" in cpd.annotation:
                        key = cpd.annotation["InChIKey"]
                        self._structure_hash.setdefault(key, {"type": "InChIKey", "ids": []})
                        self._structure_hash[key]["ids"].append(cpd.id)
                        key_components = key.split("-")
                        self._structure_hash.setdefault(key_components[0], {"type": "InChIKeyBaseOne", "ids": []})
                        self._structure_hash[key_components[0]]["ids"].append(cpd.id)
                        self._structure_hash.setdefault(key_components[0]+"-"+key_components[1], {"type": "InChIKeyBaseTwo", "ids": []})
                        self._structure_hash[key_components[0]+"-"+key_components[1]]["ids"].append(cpd.id)
        return self._structure_hash

    @property
    def element_hashes(self) -> Dict:
        """Index the biochemistry compounds by the elements they contain."""
        if self._element_hashes is None:
            self._element_hashes = {"element_hash":{},"element_count":{},"cpd_elements":{}}
            for cpd in self._biochem_db.compounds:
                if cpd.is_obsolete == False and cpd.formula is not None and len(cpd.formula) > 0:
                    elements = self._parse_formula(cpd.formula)
                    self._element_hashes["cpd_elements"][cpd.id] = elements
                    self._element_hashes["element_count"][cpd.id] = len(elements)
                    if "H" in elements:
                        self._element_hashes["element_count"][cpd.id] += -1
                    for element, count in elements.items():
                        self._element_hashes["element_hash"].setdefault(element, {})
                        self._element_hashes["element_hash"][element].setdefault(count, [])
                        self._element_hashes["element_hash"][element][count].append(cpd.id)
        return self._element_hashes

    def _standardize_string(self, input_string: str) -> str:
        """Standardize a string by lowercasing and stripping whitespace."""
        input_string = input_string.translate(str.maketrans('', '', string.punctuation + string.whitespace))
        return input_string.lower().strip()

    def _parse_formula(self,formula: str) -> dict:
        # Match elements (capital letter + optional lowercase letters) followed by optional digits
        tokens = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        
        counts = defaultdict(int)
        for element, num in tokens:
            counts[element] += int(num) if num else 1
        return dict(counts)

    def _ensure_database_available(self) -> None:
        """Ensure ModelSEED database is available, using dependency management."""
        # First ensure ModelSEEDpy is available
        if not self.ensure_modelseed_py():
            raise ImportError("Failed to obtain ModelSEEDpy dependency")

        # Ensure ModelSEEDDatabase is available
        if not self.ensure_modelseed_database():
            if self.auto_download:
                self.log_warning(
                    "Failed to automatically obtain ModelSEEDDatabase, falling back to manual download"
                )
                self._download_database_manual()
            else:
                raise FileNotFoundError(
                    "ModelSEED database not found and auto_download is disabled"
                )

        # Load the database
        try:
            from modelseedpy.biochem.modelseed_biochem import ModelSEEDBiochem

            self._biochem_db = ModelSEEDBiochem.get(path=self.modelseed_db_path)
            self.log_info(f"ModelSEED database loaded from {self.modelseed_db_path}")
        except Exception as e:
            self.log_error(f"Failed to load ModelSEED database: {e}")
            raise

    def _download_database_manual(self) -> None:
        """Download the ModelSEED database from GitHub as fallback."""
        git_url = "https://github.com/ModelSEED/ModelSEEDDatabase.git"

        # Remove existing directory if it exists and is problematic
        if os.path.exists(self.modelseed_db_path):
            import shutil

            self.log_info(f"Removing existing directory {self.modelseed_db_path}")
            shutil.rmtree(self.modelseed_db_path)

        try:
            # Create parent directory if it doesn't exist
            parent_dir = os.path.dirname(self.modelseed_db_path)
            os.makedirs(parent_dir, exist_ok=True)

            # Clone the repository
            self.log_info(f"Cloning ModelSEED database from {git_url}")
            subprocess.run(
                ["git", "clone", git_url, self.modelseed_db_path],
                capture_output=True,
                text=True,
                check=True,
                timeout=180,  # Longer timeout for this large repository
            )

            self.log_info("ModelSEED database downloaded successfully")

        except subprocess.TimeoutExpired:
            error_msg = f"Timeout while cloning ModelSEED database from {git_url}"
            self.log_error(error_msg)
            raise RuntimeError(error_msg)
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to clone ModelSEED database: {e.stderr}"
            self.log_error(error_msg)
            raise RuntimeError(error_msg)
        except Exception as e:
            error_msg = f"Error downloading ModelSEED database: {e}"
            self.log_error(error_msg)
            raise

    def update_database(self) -> None:
        """Update the ModelSEED database to the latest version."""
        if not os.path.exists(self.modelseed_db_path):
            self._download_database()
            return

        try:
            self.log_info("Updating ModelSEED database...")
            subprocess.run(
                ["git", "pull"],
                cwd=self.modelseed_db_path,
                capture_output=True,
                text=True,
                check=True,
            )

            self.log_info("ModelSEED database updated successfully")

            # Reload the database
            self._biochem_db = None
            self._ensure_database_available()

        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to update ModelSEED database: {e.stderr}"
            self.log_error(error_msg)
            raise RuntimeError(error_msg)          

    def search_compounds(
        self,
        query_identifiers: List[str] = [],
        query_structures: List[str] = [],
        query_formula: str = None
    ) -> List[Dict[str, Any]]:
        """Search for compounds in the ModelSEED database."""
        identifier_hash = self.identifier_hash
        element_hashes = self.element_hashes
        structure_hash = self.structure_hash
        matches = {}
        for item in query_identifiers:
            item = self._standardize_string(item)
            if item in identifier_hash:
                for hit in identifier_hash[item]["ids"]:
                    matches.setdefault(hit, {"score": 0, "identifier_hits": {}, "formula_hits": {},"structure_hits": {}})
                    matches[hit]["identifier_hits"][item] = 1 * identifier_hash[item]["type"]
                    matches[hit]["score"] += 10
        for item in query_structures:
            if item in structure_hash:
                for hit in structure_hash[item]["ids"]:
                    matches.setdefault(hit, {"score": 0, "identifier_hits": {}, "formula_hits": {},"structure_hits": {}})
                    matches[hit]["structure_hits"][item] = 1 * structure_hash[item]["type"]
                    matches[hit]["score"] += 8
        if query_formula is not None and len(query_formula) > 0:
            elements = self._parse_formula(query_formula)
            accumulated_hits = {}
            element_count = 0
            for element in elements:
                if element != "H":
                    element_count += 1
                    if element in element_hashes["element_hash"] and elements[element] in element_hashes["element_hash"][element]:
                        for hit in element_hashes["element_hash"][element][elements[element]]:
                            accumulated_hits.setdefault(hit, 0)
                            accumulated_hits[hit] += 1
            #Checking if the hits have the same number of elements and the same number of instances of each element
            for hit in accumulated_hits:
                if accumulated_hits[hit] == element_count and element_hashes["element_count"][hit] == element_count:
                    matches.setdefault(hit, {"score": 0, "identifier_hits": {}, "formula_hits": {}, "structure_hits": {}})
                    hittype = "No H"
                    matches[hit]["score"] += 2
                    if "H" in elements and "H" in element_hashes["cpd_elements"][hit] and elements["H"] == element_hashes["cpd_elements"][hit]["H"]:
                        hittype = "H match"
                        matches[hit]["score"] += 2
                    matches[hit]["formula_hits"][query_formula+":"+self.biochem_db.compounds.get_by_id(hit).formula] = hittype
        return matches

    def search_reactions(
        self,
        query: str,
        search_fields: Optional[List[str]] = None,
        exact_match: bool = False,
        case_sensitive: bool = False,
        max_results: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """Search for reactions in the ModelSEED database.

        Args:
            query: Search query string
            search_fields: Fields to search in ['name', 'equation', 'ec_numbers', 'aliases']
            exact_match: Whether to require exact matches
            case_sensitive: Whether search should be case sensitive
            max_results: Maximum number of results to return

        Returns:
            List of reaction dictionaries with match information
        """
        if search_fields is None:
            search_fields = ["name", "equation", "aliases"]

        results = []
        query_lower = query.lower() if not case_sensitive else query

        for reaction in self.biochem_db.reactions:
            match_score = 0
            matched_fields = []

            # Search in reaction name
            if "name" in search_fields and reaction.name:
                search_text = reaction.name if case_sensitive else reaction.name.lower()
                if exact_match:
                    if search_text == query_lower:
                        match_score += 10
                        matched_fields.append("name")
                else:
                    if query_lower in search_text:
                        match_score += 10
                        matched_fields.append("name")

            # Search in equation
            if "equation" in search_fields and hasattr(reaction, "equation"):
                search_text = (
                    reaction.equation if case_sensitive else reaction.equation.lower()
                )
                if exact_match:
                    if search_text == query_lower:
                        match_score += 15
                        matched_fields.append("equation")
                else:
                    if query_lower in search_text:
                        match_score += 8
                        matched_fields.append("equation")

            # Search in aliases
            if "aliases" in search_fields and hasattr(reaction, "aliases"):
                for alias in reaction.aliases:
                    search_text = alias if case_sensitive else alias.lower()
                    if exact_match:
                        if search_text == query_lower:
                            match_score += 8
                            matched_fields.append("aliases")
                    else:
                        if query_lower in search_text:
                            match_score += 8
                            matched_fields.append("aliases")

            # Search in EC numbers
            if "ec_numbers" in search_fields and hasattr(reaction, "ec_numbers"):
                for ec in reaction.ec_numbers:
                    if query_lower in ec.lower():
                        match_score += 12
                        matched_fields.append("ec_numbers")

            if match_score > 0:
                result = {
                    "reaction_id": reaction.id,
                    "name": reaction.name,
                    "equation": getattr(reaction, "equation", ""),
                    "direction": getattr(reaction, "direction", ""),
                    "reversibility": getattr(reaction, "reversibility", ""),
                    "match_score": match_score,
                    "matched_fields": list(set(matched_fields)),
                    "reaction_object": reaction,
                }

                # Add additional attributes if available
                if hasattr(reaction, "aliases"):
                    result["aliases"] = reaction.aliases
                if hasattr(reaction, "ec_numbers"):
                    result["ec_numbers"] = reaction.ec_numbers
                if hasattr(reaction, "pathways"):
                    result["pathways"] = reaction.pathways

                results.append(result)

        # Sort by match score (highest first)
        results.sort(key=lambda x: x["match_score"], reverse=True)

        # Limit results if specified
        if max_results:
            results = results[:max_results]

        self.log_info(f"Found {len(results)} reactions matching '{query}'")
        return results

    def get_compound_by_id(self, compound_id: str) -> Optional[Any]:
        """Get a compound by its ModelSEED ID.

        Args:
            compound_id: ModelSEED compound ID

        Returns:
            Compound object or None if not found
        """
        try:
            return self.biochem_db.compounds.get_by_id(compound_id)
        except:
            return None

    def get_reaction_by_id(self, reaction_id: str) -> Optional[Any]:
        """Get a reaction by its ModelSEED ID.

        Args:
            reaction_id: ModelSEED reaction ID

        Returns:
            Reaction object or None if not found
        """
        try:
            return self.biochem_db.reactions.get_by_id(reaction_id)
        except:
            return None

    def get_database_statistics(self) -> Dict[str, Any]:
        """Get statistics about the ModelSEED database.

        Returns:
            Dictionary with database statistics
        """
        stats = {
            "database_path": self.modelseed_db_path,
            "total_compounds": len(self.biochem_db.compounds)
            if self.biochem_db.compounds
            else 0,
            "total_reactions": len(self.biochem_db.reactions)
            if self.biochem_db.reactions
            else 0,
        }

        # Count compounds by charge
        charge_distribution = {}
        formula_patterns = {}

        for compound in self.biochem_db.compounds:
            charge = getattr(compound, "charge", 0)
            charge_distribution[charge] = charge_distribution.get(charge, 0) + 1

            # Count common formula patterns
            if compound.formula:
                # Extract major elements
                elements = re.findall(r"[A-Z][a-z]?", compound.formula)
                for element in elements:
                    formula_patterns[element] = formula_patterns.get(element, 0) + 1

        stats["charge_distribution"] = charge_distribution
        stats["common_elements"] = dict(
            sorted(formula_patterns.items(), key=lambda x: x[1], reverse=True)[:10]
        )

        # Count reactions by reversibility
        reversibility_stats = {}
        for reaction in self.biochem_db.reactions:
            reversibility = getattr(reaction, "reversibility", "unknown")
            reversibility_stats[reversibility] = (
                reversibility_stats.get(reversibility, 0) + 1
            )

        stats["reversibility_distribution"] = reversibility_stats

        return stats
