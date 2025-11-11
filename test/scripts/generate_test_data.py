"""
Generate real UniProt test data for pipeline.

Fetches real proteins from a small PFAM family, generates embeddings using ESM-2,
and creates the data files needed for testing.

Creates:
- data/embeddings/embeddings.tsv (uniprot_id \t comma-separated floats)
- data/indexes/id_to_row.json (uniprot_id -> row index)
- data/indexes/binary_embeddings.npy (uint8 matrix: value>0 -> 1 else 0)

This uses real UniProt data so metadata fetching works correctly in tests.

Recommended PFAM families for testing (small, well-characterized):
- PF00096: Zinc finger, C2H2 type (default) - ~1000+ reviewed proteins, small domain
- PF00004: ATP synthase domain - ~500+ reviewed proteins
- PF00005: ABC transporter - ~400+ reviewed proteins
- PF00106: Short chain dehydrogenase - ~300+ reviewed proteins
- PF00076: RNA recognition motif - ~200+ reviewed proteins

Usage:
    # Generate data with default settings (PF00096, 50 proteins)
    python test/scripts/generate_test_data.py
    
    # Use a different PFAM family
    python test/scripts/generate_test_data.py --pfam_id PF00004 --max_proteins 100
    
    # Skip embedding generation (faster, for testing script only)
    python test/scripts/generate_test_data.py --skip_embeddings

Note: Embedding generation can take 5-10 minutes for 50 proteins depending on hardware.
The script will download the ESM-2 model on first use (~50MB).
"""

import os
import argparse
import json
import time
import sys
from typing import List, Dict, Optional
import logging

import numpy as np
import requests

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Ensure project lib/ is on path for module imports when run directly
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
LIB_SRC_DIR = os.path.join(PROJECT_ROOT, "lib", "kbase_protein_query_module", "src")
if LIB_SRC_DIR not in sys.path:
    sys.path.insert(0, LIB_SRC_DIR)

# Import after path setup
try:
    from util.embeddings.generator import ProteinEmbeddingGenerator
    from util.uniprot.api import fetch_protein_sequence, fetch_metadata
except ImportError as e:
    logger.error(f"Failed to import required modules: {e}")
    logger.error("Make sure you're running from the project root and dependencies are installed")
    raise

UNIPROT_REST = "https://rest.uniprot.org"


def query_pfam_proteins(pfam_id: str, max_proteins: int = 100, reviewed_only: bool = True) -> List[str]:
    """
    Query UniProt for proteins in a PFAM family.
    
    Args:
        pfam_id: PFAM family ID (e.g., "PF00096")
        max_proteins: Maximum number of proteins to return
        reviewed_only: If True, only return reviewed (Swiss-Prot) proteins
    
    Returns:
        List of UniProt accession IDs
    """
    logger.info(f"Querying UniProt for PFAM family {pfam_id}...")
    
    # Build query - UniProt uses xref:pfam-{pfam_id} syntax
    # Note: Use 'database:(type:pfam AND id:{pfam_id})' as alternative if xref doesn't work
    if reviewed_only:
        # For reviewed proteins, we need to be more explicit
        query = f"(xref:pfam-{pfam_id}) AND (reviewed:true)"
    else:
        query = f"xref:pfam-{pfam_id}"
    
    # Query parameters
    params = {
        "query": query,
        "fields": "accession",
        "format": "tsv",
        "size": max_proteins,
        "sort": "accession asc",  # Deterministic ordering (ascending)
    }
    
    try:
        url = f"{UNIPROT_REST}/uniprotkb/search"
        logger.debug(f"UniProt query: {query}")
        resp = requests.get(url, params=params, timeout=30)
        
        if resp.status_code != 200:
            error_msg = resp.text[:500] if resp.text else "No error message"
            logger.error(f"UniProt API error response: {error_msg}")
            raise RuntimeError(f"UniProt query failed: HTTP {resp.status_code}\nQuery: {query}\nError: {error_msg}")
        
        lines = resp.text.strip().splitlines()
        if len(lines) < 2:
            logger.warning(f"No proteins found for PFAM {pfam_id}")
            logger.warning(f"Response: {resp.text[:200]}")
            return []
        
        # Parse accessions (skip header)
        accessions = [line.strip() for line in lines[1:] if line.strip()]
        logger.info(f"Found {len(accessions)} proteins for PFAM {pfam_id}")
        
        if len(accessions) == 0:
            logger.warning("No accessions parsed from response")
            return []
        
        return accessions[:max_proteins]
    
    except requests.exceptions.RequestException as e:
        logger.error(f"Network error querying UniProt for PFAM {pfam_id}: {e}")
        raise
    except Exception as e:
        logger.error(f"Failed to query UniProt for PFAM {pfam_id}: {e}")
        raise


def fetch_protein_sequences(uniprot_ids: List[str], batch_size: int = 10) -> Dict[str, str]:
    """
    Fetch protein sequences for UniProt IDs.
    
    Args:
        uniprot_ids: List of UniProt accession IDs
        batch_size: Number of sequences to fetch in parallel (not used, but kept for future optimization)
    
    Returns:
        Dictionary mapping UniProt ID to sequence
    """
    logger.info(f"Fetching sequences for {len(uniprot_ids)} proteins...")
    sequences = {}
    failed = []
    
    for i, uniprot_id in enumerate(uniprot_ids):
        try:
            sequence = fetch_protein_sequence(uniprot_id)
            if sequence and len(sequence) > 0:
                sequences[uniprot_id] = sequence
                if (i + 1) % 10 == 0:
                    logger.info(f"Fetched {i + 1}/{len(uniprot_ids)} sequences...")
            else:
                logger.warning(f"Empty sequence for {uniprot_id}")
                failed.append(uniprot_id)
            # Small delay to avoid rate limiting
            time.sleep(0.1)
        except Exception as e:
            logger.warning(f"Failed to fetch sequence for {uniprot_id}: {e}")
            failed.append(uniprot_id)
    
    logger.info(f"Successfully fetched {len(sequences)} sequences, {len(failed)} failed")
    if failed:
        logger.warning(f"Failed IDs: {failed[:10]}...")  # Show first 10
    
    return sequences


def generate_embeddings(sequences: Dict[str, str], model_name: str = "esm2_t6_8M_UR50D", batch_size: int = 10) -> Dict[str, np.ndarray]:
    """
    Generate embeddings for protein sequences using ESM-2.
    
    Args:
        sequences: Dictionary mapping UniProt ID to sequence
        model_name: ESM model name to use
        batch_size: Number of sequences to process before logging progress
    
    Returns:
        Dictionary mapping UniProt ID to embedding vector
    """
    logger.info(f"Generating embeddings for {len(sequences)} proteins using {model_name}...")
    
    # Initialize embedding generator
    try:
        generator = ProteinEmbeddingGenerator(model_name=model_name, device="cpu")
        logger.info(f"Loaded model with embedding dimension: {generator.embedding_dim}")
    except Exception as e:
        logger.error(f"Failed to initialize embedding generator: {e}")
        raise
    
    embeddings = {}
    failed = []
    
    # Process sequences with error handling and progress logging
    sequences_list = list(sequences.items())
    total = len(sequences_list)
    
    for i, (uniprot_id, sequence) in enumerate(sequences_list):
        try:
            # Validate sequence before processing
            if not sequence or len(sequence.strip()) < 3:
                logger.warning(f"Skipping {uniprot_id}: sequence too short or empty")
                failed.append(uniprot_id)
                continue
            
            # Generate mean-pooled embedding
            embedding = generator.generate_embedding(sequence, pooling="mean")
            if embedding is not None and embedding.size > 0:
                embeddings[uniprot_id] = embedding
                if (i + 1) % batch_size == 0 or (i + 1) == total:
                    logger.info(f"Generated {i + 1}/{total} embeddings...")
            else:
                logger.warning(f"Empty embedding for {uniprot_id}")
                failed.append(uniprot_id)
        except KeyboardInterrupt:
            logger.warning("Embedding generation interrupted by user")
            raise
        except Exception as e:
            logger.warning(f"Failed to generate embedding for {uniprot_id}: {e}")
            logger.debug(f"Sequence length: {len(sequence) if sequence else 0}")
            failed.append(uniprot_id)
            # Continue with next sequence instead of failing completely
    
    logger.info(f"Successfully generated {len(embeddings)} embeddings, {len(failed)} failed")
    if failed:
        logger.warning(f"Failed IDs (showing first 10): {failed[:10]}")
        if len(failed) > 10:
            logger.warning(f"... and {len(failed) - 10} more")
    
    if len(embeddings) == 0:
        raise RuntimeError("No embeddings were generated successfully")
    
    return embeddings


def write_embeddings_file(embeddings: Dict[str, np.ndarray], output_path: str) -> None:
    """
    Write embeddings to TSV file.
    Overwrites existing file if it exists.
    
    Args:
        embeddings: Dictionary mapping UniProt ID to embedding vector
        output_path: Output file path
    """
    logger.info(f"Writing embeddings to {output_path}...")
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Remove existing file if it exists to ensure clean overwrite
    if os.path.exists(output_path):
        logger.info(f"Removing existing file: {output_path}")
        os.remove(output_path)
    
    # Write embeddings file
    with open(output_path, "w") as f:
        # Write header
        f.write("uniprot_id\tembedding\n")
        
        # Write embeddings (sorted for deterministic output)
        for uniprot_id in sorted(embeddings.keys()):
            embedding = embeddings[uniprot_id]
            # Convert numpy array to comma-separated string
            emb_str = ",".join(f"{x:.6f}" for x in embedding.tolist())
            f.write(f"{uniprot_id}\t{emb_str}\n")
    
    logger.info(f"Written {len(embeddings)} embeddings to {output_path}")


def write_indexes(embeddings: Dict[str, np.ndarray], output_dir: str) -> Dict[str, str]:
    """
    Write index files for embeddings.
    Overwrites existing files if they exist.
    
    Args:
        embeddings: Dictionary mapping UniProt ID to embedding vector
        output_dir: Output directory for index files
    
    Returns:
        Dictionary with paths to created index files
    """
    logger.info(f"Writing indexes to {output_dir}...")
    
    # Ensure directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get sorted list of IDs for deterministic ordering
    uniprot_ids = sorted(embeddings.keys())
    
    # Create id_to_row mapping
    id_to_row = {uniprot_id: i for i, uniprot_id in enumerate(uniprot_ids)}
    
    # Write id_to_row.json (overwrite if exists)
    id_map_path = os.path.join(output_dir, "id_to_row.json")
    if os.path.exists(id_map_path):
        logger.info(f"Removing existing file: {id_map_path}")
        os.remove(id_map_path)
    with open(id_map_path, "w") as f:
        json.dump(id_to_row, f, indent=2)
    logger.info(f"Written id_to_row.json with {len(id_to_row)} entries")
    
    # Create binary embeddings matrix
    # Stack all embeddings into a matrix
    embedding_matrix = np.vstack([embeddings[uid] for uid in uniprot_ids]).astype(np.float32)
    
    # Create binary matrix (value > 0 -> 1, else 0)
    binary_matrix = (embedding_matrix > 0).astype(np.uint8)
    
    # Save binary embeddings (overwrite if exists)
    binary_path = os.path.join(output_dir, "binary_embeddings.npy")
    if os.path.exists(binary_path):
        logger.info(f"Removing existing file: {binary_path}")
        os.remove(binary_path)
    np.save(binary_path, binary_matrix)
    logger.info(f"Written binary_embeddings.npy with shape {binary_matrix.shape}")
    
    return {
        "id_to_row": id_map_path,
        "binary_embeddings": binary_path,
    }


def main() -> int:
    """Main function to generate real UniProt test data."""
    parser = argparse.ArgumentParser(
        description="Generate real UniProt test data from a PFAM family"
    )
    parser.add_argument(
        "--pfam_id",
        type=str,
        default="PF00096",
        help="PFAM family ID (default: PF00096 - Zinc finger, C2H2 type)"
    )
    parser.add_argument(
        "--max_proteins",
        type=int,
        default=50,
        help="Maximum number of proteins to fetch (default: 50)"
    )
    parser.add_argument(
        "--reviewed_only",
        action="store_true",
        help="Only fetch reviewed (Swiss-Prot) proteins (default: True, unless --include_unreviewed is used)"
    )
    parser.add_argument(
        "--include_unreviewed",
        action="store_true",
        help="Include unreviewed (TrEMBL) proteins (default: False)"
    )
    parser.add_argument(
        "--model",
        type=str,
        default="esm2_t6_8M_UR50D",
        help="ESM model name (default: esm2_t6_8M_UR50D)"
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="data",
        help="Output directory (default: data)"
    )
    parser.add_argument(
        "--skip_embeddings",
        action="store_true",
        help="Skip embedding generation (only fetch sequences)"
    )
    
    args = parser.parse_args()
    
    try:
        logger.info("=" * 60)
        logger.info("Generating real UniProt test data")
        logger.info("=" * 60)
        # Determine if we should fetch reviewed only
        # Default is True (reviewed only) unless --include_unreviewed is specified
        reviewed_only = not args.include_unreviewed
        
        logger.info(f"PFAM family: {args.pfam_id}")
        logger.info(f"Max proteins: {args.max_proteins}")
        logger.info(f"Reviewed only: {reviewed_only}")
        logger.info(f"Model: {args.model}")
        logger.info(f"Output directory: {args.out_dir}")
        logger.info("=" * 60)
        
        # Step 1: Query UniProt for proteins in PFAM family
        uniprot_ids = query_pfam_proteins(
            args.pfam_id,
            max_proteins=args.max_proteins,
            reviewed_only=reviewed_only
        )
        
        if not uniprot_ids:
            logger.error("No proteins found. Exiting.")
            return 1
        
        logger.info(f"Found {len(uniprot_ids)} proteins")
        logger.info(f"Sample IDs: {uniprot_ids[:5]}")
        
        # Step 2: Fetch protein sequences
        sequences = fetch_protein_sequences(uniprot_ids)
        
        if not sequences:
            logger.error("No sequences fetched. Exiting.")
            return 1
        
        logger.info(f"Fetched {len(sequences)} sequences")
        
        # Step 3: Generate embeddings (if not skipped)
        if args.skip_embeddings:
            logger.info("Skipping embedding generation")
            embeddings = {}
            # Create dummy embeddings for testing (will be replaced later)
            for uniprot_id in sequences.keys():
                embeddings[uniprot_id] = np.random.normal(0, 1, size=(320,)).astype(np.float32)
        else:
            embeddings = generate_embeddings(sequences, model_name=args.model)
        
        if not embeddings:
            logger.error("No embeddings generated. Exiting.")
            return 1
        
        logger.info(f"Generated {len(embeddings)} embeddings")
        
        # Step 4: Prepare output directories
        embeddings_dir = os.path.join(args.out_dir, "embeddings")
        indexes_dir = os.path.join(args.out_dir, "indexes")
        
        # Ensure directories exist (will be created if needed)
        os.makedirs(embeddings_dir, exist_ok=True)
        os.makedirs(indexes_dir, exist_ok=True)
        
        logger.info(f"Output directories prepared: {embeddings_dir}, {indexes_dir}")
        
        # Step 5: Write embeddings file (overwrites if exists)
        embeddings_file = os.path.join(embeddings_dir, "embeddings.tsv")
        write_embeddings_file(embeddings, embeddings_file)
        
        # Step 6: Write index files (overwrites if exist)
        index_paths = write_indexes(embeddings, indexes_dir)
        
        # Step 7: Verify with metadata fetch
        logger.info("Verifying with metadata fetch...")
        sample_ids = list(embeddings.keys())[:5]
        metadata = fetch_metadata(sample_ids)
        if metadata:
            logger.info(f"Successfully fetched metadata for {len(metadata)} proteins")
            for meta in metadata[:2]:
                logger.info(f"  {meta.get('Entry')}: {meta.get('Protein names', 'N/A')}")
        
        # Summary
        logger.info("=" * 60)
        logger.info("Summary")
        logger.info("=" * 60)
        logger.info(f"Proteins: {len(embeddings)}")
        logger.info(f"Embeddings file: {embeddings_file}")
        logger.info(f"Index files: {index_paths}")
        logger.info("=" * 60)
        
        # Print file sizes
        if os.path.exists(embeddings_file):
            size_mb = os.path.getsize(embeddings_file) / (1024 * 1024)
            logger.info(f"Embeddings file size: {size_mb:.2f} MB")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error generating test data: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
