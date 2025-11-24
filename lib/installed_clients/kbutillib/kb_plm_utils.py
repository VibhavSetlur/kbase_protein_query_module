"""KBase Protein Language Model (PLM) utilities for protein homology search and annotation.

This module provides utilities for interfacing with the KBase Protein Language Model API
to find protein homologs, construct BLAST databases, and perform sequence similarity searches.
"""

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests

from kb_genome_utils import KBGenomeUtils


class KBPLMUtils(KBGenomeUtils):
    """Utilities for protein homology search using the KBase Protein Language Model API.

    This class provides methods to:
    - Query the PLM API for protein homologs
    - Retrieve Uniprot sequences for hits
    - Create custom BLAST databases
    - Perform BLAST searches
    - Match features to their best homologs
    """

    def __init__(
        self,
        plm_api_url: str = "https://kbase.us/services/llm_homology_api",
        **kwargs: Any
    ) -> None:
        """Initialize KBase PLM utilities.

        Args:
            plm_api_url: Base URL for the PLM API (default: KBase PLM API)
            **kwargs: Additional keyword arguments passed to KBGenomeUtils
        """
        super().__init__(**kwargs)
        self.plm_api_url = plm_api_url.rstrip("/")
        self.plm_search_endpoint = f"{self.plm_api_url}/search"
        self.plm_result_endpoint = f"{self.plm_api_url}/result"
        self.plm_sequence_endpoint = f"{self.plm_api_url}/sequences"

        # Check if BLAST is available
        self._check_blast_availability()

    def _check_blast_availability(self) -> bool:
        """Check if BLAST tools are available in the system.

        Returns:
            bool: True if BLAST is available, False otherwise
        """
        try:
            subprocess.run(
                ["blastp", "-version"],
                capture_output=True,
                check=True,
                timeout=5
            )
            self.blast_available = True
            self.log_info("BLAST tools are available")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            self.blast_available = False
            self.log_warning(
                "BLAST tools not found. Install NCBI BLAST+ to use BLAST functionality. "
                "On Ubuntu/Debian: sudo apt-get install ncbi-blast+, "
                "On MacOS: brew install blast"
            )
            return False

    def query_plm_api(
        self,
        query_sequences: List[Dict[str, str]],
        max_hits: int = 100,
        similarity_threshold: float = 0.0,
        return_embeddings: bool = False,
        poll_interval: float = 2.0,
        max_wait_time: float = 300.0
    ) -> Dict[str, Any]:
        """Query the PLM API for protein homologs.

        This function submits a job to the PLM API and polls for results.
        The API returns a job_id immediately, then the function polls the
        result endpoint until the job completes.

        Args:
            query_sequences: List of dicts with 'id' and 'sequence' keys
            max_hits: Maximum number of hits to return per query (1-100)
            similarity_threshold: Minimum similarity score threshold
            return_embeddings: Whether to return embeddings for hits
            poll_interval: Time in seconds between polling attempts (default: 2.0)
            max_wait_time: Maximum time in seconds to wait for results (default: 300.0)

        Returns:
            Dict containing API response with hits for each query

        Raises:
            ValueError: If input validation fails
            requests.RequestException: If API request fails
            TimeoutError: If job doesn't complete within max_wait_time
        """
        # Validate inputs
        if not query_sequences:
            raise ValueError("query_sequences cannot be empty")

        if not 1 <= max_hits <= 100:
            raise ValueError("max_hits must be between 1 and 100")

        # Prepare request payload
        payload = {
            "query_sequences": query_sequences,
            "max_hits": max_hits,
            "similarity_threshold": similarity_threshold,
            "best_hit_only": False,
            "return_query_embeddings": return_embeddings,
            "return_hit_embeddings": return_embeddings
        }

        self.log_info(
            f"Querying PLM API with {len(query_sequences)} sequences, "
            f"requesting up to {max_hits} hits per sequence"
        )

        # Step 1: Submit the job to the search endpoint
        try:
            response = requests.post(
                self.plm_search_endpoint,
                json=payload,
                timeout=30,  # Short timeout for job submission
                verify=False  # Following the pattern in base_utils.py
            )
            response.raise_for_status()

            job_response = response.json()

            # Extract job_id from response
            if "job_id" not in job_response:
                raise ValueError(
                    f"Expected 'job_id' in response, got: {job_response.keys()}"
                )

            job_id = job_response["job_id"]
            self.log_info(f"PLM job submitted successfully with job_id: {job_id}")

        except requests.exceptions.Timeout:
            self.log_error("PLM API job submission timed out")
            raise
        except requests.exceptions.RequestException as e:
            self.log_error(f"PLM API job submission failed: {str(e)}")
            raise

        # Step 2: Poll the result endpoint for job completion
        # Note: The result endpoint requires POST with job_id in the body, not GET
        elapsed_time = 0.0

        self.log_info(
            f"Polling for results at {self.plm_result_endpoint} for job {job_id} "
            f"(interval: {poll_interval}s, max_wait: {max_wait_time}s)"
        )

        while elapsed_time < max_wait_time:
            try:
                time.sleep(poll_interval)
                elapsed_time += poll_interval

                result_response = requests.post(
                    self.plm_result_endpoint,
                    json={"job_id": job_id},
                    timeout=30,
                    verify=False
                )

                # Check if job is complete
                if result_response.status_code == 200:
                    response_json = result_response.json()

                    # API returns: {"status": "done|pending|running|failed", "result": {...}, "error": ...}
                    status = response_json.get("status", "unknown")

                    if status == "done":
                        # Job completed successfully - return the result field
                        result = response_json.get("result")
                        if result is None:
                            raise RuntimeError("PLM job completed but result is None")

                        self.log_info(
                            f"PLM job completed successfully after {elapsed_time:.1f}s"
                        )
                        return result

                    elif status == "failed":
                        error_msg = response_json.get("error", "Unknown error")
                        raise RuntimeError(f"PLM job failed: {error_msg}")

                    elif status in ["pending", "running"]:
                        self.log_debug(
                            f"Job still {status} after {elapsed_time:.1f}s"
                        )
                        continue

                    else:
                        # Unknown status - log and continue polling
                        self.log_warning(
                            f"Unknown job status '{status}' after {elapsed_time:.1f}s"
                        )
                        continue

                elif result_response.status_code == 202:
                    # Job still processing
                    self.log_debug(
                        f"Job still processing (HTTP 202) after {elapsed_time:.1f}s"
                    )
                    continue
                elif result_response.status_code == 404:
                    raise RuntimeError(f"Job {job_id} not found (HTTP 404)")
                else:
                    result_response.raise_for_status()

            except requests.exceptions.Timeout:
                self.log_warning(
                    f"Timeout while polling for results (elapsed: {elapsed_time:.1f}s)"
                )
                continue
            except requests.exceptions.RequestException as e:
                self.log_error(f"Error polling for results: {str(e)}")
                raise

        # If we've exhausted max_wait_time
        raise TimeoutError(
            f"PLM job {job_id} did not complete within {max_wait_time}s"
        )

    def get_uniprot_sequences(
        self,
        uniprot_ids: List[str]
    ) -> Dict[str, str]:
        """Retrieve protein sequences from UniProt.

        Args:
            uniprot_ids: List of UniProt IDs

        Returns:
            Dict mapping UniProt IDs to their sequences
        """
        sequences = {}

        self.log_info(f"Retrieving {len(uniprot_ids)} sequences from UniProt")

        # UniProt REST API endpoint
        base_url = "https://rest.uniprot.org/uniprotkb"

        for uniprot_id in uniprot_ids:
            try:
                # Fetch FASTA format
                url = f"{base_url}/{uniprot_id}.fasta"
                response = requests.get(url, timeout=30)

                if response.status_code == 200:
                    fasta_text = response.text
                    # Parse FASTA (skip header line, join sequence lines)
                    lines = fasta_text.strip().split('\n')
                    if len(lines) > 1:
                        sequence = ''.join(lines[1:])
                        sequences[uniprot_id] = sequence
                else:
                    self.log_warning(
                        f"Could not retrieve sequence for {uniprot_id}: "
                        f"HTTP {response.status_code}"
                    )

            except requests.exceptions.RequestException as e:
                self.log_warning(f"Failed to retrieve {uniprot_id}: {str(e)}")
                continue

        self.log_info(f"Successfully retrieved {len(sequences)} sequences from UniProt")
        return sequences

    def create_blast_database(
        self,
        sequences: Dict[str, str],
        db_path: str
    ) -> bool:
        """Create a BLAST database from protein sequences.

        Args:
            sequences: Dict mapping sequence IDs to sequences
            db_path: Path where the BLAST database should be created

        Returns:
            bool: True if database was created successfully

        Raises:
            RuntimeError: If BLAST is not available or database creation fails
        """
        if not self.blast_available:
            raise RuntimeError(
                "BLAST tools are not available. Please install NCBI BLAST+."
            )

        # Create FASTA file for input
        fasta_path = f"{db_path}.fasta"

        try:
            with open(fasta_path, 'w') as f:
                for seq_id, sequence in sequences.items():
                    f.write(f">{seq_id}\n{sequence}\n")

            self.log_info(f"Created FASTA file with {len(sequences)} sequences")

            # Run makeblastdb
            cmd = [
                "makeblastdb",
                "-in", fasta_path,
                "-dbtype", "prot",
                "-out", db_path,
                "-parse_seqids"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120
            )

            if result.returncode == 0:
                self.log_info(f"Successfully created BLAST database at {db_path}")
                return True
            else:
                self.log_error(f"makeblastdb failed: {result.stderr}")
                return False

        except Exception as e:
            self.log_error(f"Failed to create BLAST database: {str(e)}")
            return False

    def run_blastp(
        self,
        query_sequences: Dict[str, str],
        database_path: str,
        evalue: float = 0.001,
        max_target_seqs: int = 1
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Run BLASTP search against a custom database.

        Args:
            query_sequences: Dict mapping query IDs to sequences
            database_path: Path to the BLAST database
            evalue: E-value threshold
            max_target_seqs: Maximum number of target sequences to return per query

        Returns:
            Dict mapping query IDs to list of hit dictionaries

        Raises:
            RuntimeError: If BLAST is not available or search fails
        """
        if not self.blast_available:
            raise RuntimeError(
                "BLAST tools are not available. Please install NCBI BLAST+."
            )

        # Create temporary FASTA file for queries
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            query_fasta = f.name
            for query_id, sequence in query_sequences.items():
                f.write(f">{query_id}\n{sequence}\n")

        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            output_file = f.name

        try:
            # Run BLASTP with JSON output
            cmd = [
                "blastp",
                "-query", query_fasta,
                "-db", database_path,
                "-evalue", str(evalue),
                "-max_target_seqs", str(max_target_seqs),
                "-outfmt", "15",  # JSON output
                "-out", output_file
            ]

            self.log_info(
                f"Running BLASTP for {len(query_sequences)} queries "
                f"against database {database_path}"
            )

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0:
                self.log_error(f"BLASTP failed: {result.stderr}")
                raise RuntimeError(f"BLASTP failed: {result.stderr}")

            # Parse JSON output
            with open(output_file, 'r') as f:
                blast_results = json.load(f)

            # Process results into simpler format
            hits_by_query = {}

            if "BlastOutput2" in blast_results:
                for result_item in blast_results["BlastOutput2"]:
                    report = result_item.get("report", {})
                    results = report.get("results", {})
                    search = results.get("search", {})
                    query_title = search.get("query_title", "")

                    hits = []
                    for hit in search.get("hits", []):
                        description = hit.get("description", [{}])[0]
                        hit_id = description.get("id", "")
                        hit_title = description.get("title", "")

                        hsps = hit.get("hsps", [])
                        if hsps:
                            hsp = hsps[0]  # Take best HSP
                            hits.append({
                                "hit_id": hit_id,
                                "hit_title": hit_title,
                                "evalue": hsp.get("evalue", 1.0),
                                "bit_score": hsp.get("bit_score", 0.0),
                                "identity": hsp.get("identity", 0),
                                "align_len": hsp.get("align_len", 0),
                                "query_from": hsp.get("query_from", 0),
                                "query_to": hsp.get("query_to", 0),
                                "hit_from": hsp.get("hit_from", 0),
                                "hit_to": hsp.get("hit_to", 0)
                            })

                    hits_by_query[query_title] = hits

            self.log_info(f"BLASTP completed successfully")
            return hits_by_query

        finally:
            # Clean up temporary files
            try:
                os.unlink(query_fasta)
                os.unlink(output_file)
            except:
                pass

    def find_best_hits_for_features(
        self,
        feature_container_name: str,
        max_plm_hits: int = 100,
        similarity_threshold: float = 0.0,
        blast_evalue: float = 0.001,
        temp_dir: Optional[str] = None
    ) -> Dict[str, Dict[str, Any]]:
        """Find the best protein homologs for features in a KBase object.

        This is the main workflow method that:
        1. Extracts protein sequences from the feature container
        2. Queries the PLM API for top hits
        3. Retrieves UniProt sequences for those hits
        4. Creates a BLAST database
        5. Runs BLAST to find best matches
        6. Returns the best UniProt IDs for each feature

        Args:
            feature_container_name: Name of the loaded feature container object
            max_plm_hits: Maximum number of PLM hits to retrieve per feature (1-100)
            similarity_threshold: Minimum PLM similarity score threshold
            blast_evalue: E-value threshold for BLAST
            temp_dir: Directory for temporary files (uses system temp if None)

        Returns:
            Dict mapping feature IDs to their best hit information:
            {
                "feature_id": {
                    "best_uniprot_id": str,
                    "plm_score": float,
                    "blast_evalue": float,
                    "blast_bit_score": float,
                    "blast_identity": int,
                    "all_plm_hits": List[str]  # All UniProt IDs from PLM
                }
            }

        Raises:
            ValueError: If feature container is not loaded or has no features
            RuntimeError: If BLAST is not available
        """
        self.initialize_call(
            "find_best_hits_for_features",
            {
                "feature_container_name": feature_container_name,
                "max_plm_hits": max_plm_hits,
                "similarity_threshold": similarity_threshold,
                "blast_evalue": blast_evalue
            },
            print_params=True
        )

        # Step 1: Extract features and sequences
        self.log_info(f"Extracting features from {feature_container_name}")

        if not self._check_for_object(feature_container_name):
            raise ValueError(
                f"Feature container '{feature_container_name}' not loaded. "
                f"Use load_kbase_gene_container() first."
            )

        features = self.object_to_features(feature_container_name)
        if not features:
            raise ValueError(
                f"No features found in {feature_container_name}"
            )

        # Extract protein sequences
        query_sequences = []
        feature_sequences = {}

        for feature in features:
            ftr_id = feature.get("id", "")
            protein_seq = feature.get("protein_translation", "")

            if protein_seq and len(protein_seq) > 0:
                query_sequences.append({
                    "id": ftr_id,
                    "sequence": protein_seq
                })
                feature_sequences[ftr_id] = protein_seq

        self.log_info(
            f"Found {len(feature_sequences)} features with protein sequences"
        )

        if not query_sequences:
            raise ValueError(
                f"No protein sequences found in features of {feature_container_name}"
            )

        # Step 2: Query PLM API
        self.log_info("Querying PLM API for homologs")
        plm_results = self.query_plm_api(
            query_sequences,
            max_hits=max_plm_hits,
            similarity_threshold=similarity_threshold
        )

        # Step 3: Collect all unique UniProt IDs and build mapping
        all_uniprot_ids = set()
        plm_hits_by_feature = {}

        for hits_data in plm_results.get("hits", []):
            query_id = hits_data.get("query_id", "")
            hits = hits_data.get("hits", [])

            feature_hits = []
            for hit in hits:
                uniprot_id = hit.get("id", "")
                score = hit.get("score", 0.0)

                if uniprot_id:
                    all_uniprot_ids.add(uniprot_id)
                    feature_hits.append({
                        "uniprot_id": uniprot_id,
                        "plm_score": score
                    })

            plm_hits_by_feature[query_id] = feature_hits

        self.log_info(
            f"Collected {len(all_uniprot_ids)} unique UniProt IDs "
            f"from PLM results"
        )

        # Step 4: Retrieve UniProt sequences
        self.log_info("Retrieving sequences from UniProt")
        uniprot_sequences = self.get_uniprot_sequences(list(all_uniprot_ids))

        if not uniprot_sequences:
            self.log_error("Could not retrieve any UniProt sequences")
            raise RuntimeError("Failed to retrieve UniProt sequences")

        # Step 5: Create BLAST database
        if temp_dir is None:
            temp_dir = tempfile.gettempdir()

        db_path = os.path.join(temp_dir, f"plm_hits_db_{os.getpid()}")

        self.log_info("Creating BLAST database from UniProt sequences")
        if not self.create_blast_database(uniprot_sequences, db_path):
            raise RuntimeError("Failed to create BLAST database")

        # Step 6: Run BLAST
        self.log_info("Running BLAST to find best matches")
        blast_results = self.run_blastp(
            feature_sequences,
            db_path,
            evalue=blast_evalue,
            max_target_seqs=1
        )

        # Step 7: Compile final results
        final_results = {}

        for feature_id, feature_seq in feature_sequences.items():
            result_entry = {
                "best_uniprot_id": None,
                "plm_score": None,
                "blast_evalue": None,
                "blast_bit_score": None,
                "blast_identity": None,
                "all_plm_hits": []
            }

            # Get PLM hits for this feature
            plm_hits = plm_hits_by_feature.get(feature_id, [])
            result_entry["all_plm_hits"] = [
                hit["uniprot_id"] for hit in plm_hits
            ]

            # Get BLAST result for this feature
            blast_hits = blast_results.get(feature_id, [])

            if blast_hits:
                best_blast_hit = blast_hits[0]
                best_uniprot_id = best_blast_hit["hit_id"]

                result_entry["best_uniprot_id"] = best_uniprot_id
                result_entry["blast_evalue"] = best_blast_hit["evalue"]
                result_entry["blast_bit_score"] = best_blast_hit["bit_score"]
                result_entry["blast_identity"] = best_blast_hit["identity"]

                # Find the PLM score for this hit
                for plm_hit in plm_hits:
                    if plm_hit["uniprot_id"] == best_uniprot_id:
                        result_entry["plm_score"] = plm_hit["plm_score"]
                        break
            else:
                self.log_warning(f"No BLAST hits found for feature {feature_id}")

            final_results[feature_id] = result_entry

        # Clean up BLAST database files
        try:
            for ext in ['.phr', '.pin', '.psq', '.fasta']:
                file_path = f"{db_path}{ext}"
                if os.path.exists(file_path):
                    os.unlink(file_path)
        except Exception as e:
            self.log_warning(f"Could not clean up temporary files: {str(e)}")

        self.log_info(
            f"Successfully found best hits for {len(final_results)} features"
        )

        return final_results

    def get_best_uniprot_ids(
        self,
        feature_container_name: str,
        **kwargs: Any
    ) -> Dict[str, str]:
        """Simplified method to get just the best UniProt ID for each feature.

        Args:
            feature_container_name: Name of the loaded feature container object
            **kwargs: Additional arguments passed to find_best_hits_for_features

        Returns:
            Dict mapping feature IDs to their best UniProt IDs
        """
        results = self.find_best_hits_for_features(
            feature_container_name,
            **kwargs
        )

        return {
            feature_id: hit_data["best_uniprot_id"]
            for feature_id, hit_data in results.items()
            if hit_data["best_uniprot_id"] is not None
        }
