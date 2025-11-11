#!/usr/bin/env python3
"""
FAISS Indexing Guide for Protein Embeddings Database

This script provides a comprehensive guide and implementation for creating
FAISS indexes from protein embeddings, with support for family-based subsetting
for efficient hierarchical similarity search.

Key Concepts:
1. Flat Index: Simple brute-force search, best for small datasets (<10K proteins)
2. IVF Index: Inverted File Index, efficient for medium-large datasets (10K-1M proteins)
3. Family-based Indexing: Organize proteins into families for hierarchical search
   - Level 1: Map query to protein family using family centroids
   - Level 2: Search within selected family using family-specific index

Usage:
    # Create indexes for all families
    python test/scripts/create_faiss_indexing_guide.py --embeddings data/embeddings/embeddings.tsv --output data/indexes
    
    # Create indexes with family grouping
    python test/scripts/create_faiss_indexing_guide.py --embeddings data/embeddings/embeddings.tsv --output data/indexes --families data/families/
    
    # Create index for specific family
    python test/scripts/create_faiss_indexing_guide.py --embeddings data/embeddings/embeddings.tsv --output data/indexes --family "kinase_enzyme"
"""

import os
import sys
import argparse
import logging
import json
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "lib"))

try:
    import faiss
    HAS_FAISS = True
except ImportError:
    HAS_FAISS = False
    print("Warning: faiss-cpu not installed. Install with: pip install faiss-cpu")

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class FAISSIndexingGuide:
    """
    Guide and implementation for creating FAISS indexes from protein embeddings.
    
    This class provides methods to:
    1. Load embeddings from TSV/CSV files
    2. Group proteins into families (optional)
    3. Create FAISS indexes for individual families
    4. Create family centroids index for hierarchical search
    5. Save indexes and metadata
    """
    
    def __init__(self, embeddings_file: str, output_dir: str):
        """
        Initialize the FAISS indexing guide.
        
        Args:
            embeddings_file: Path to embeddings file (TSV format: uniprot_id\tembedding)
            output_dir: Directory to save indexes and metadata
        """
        self.embeddings_file = Path(embeddings_file)
        self.output_dir = Path(output_dir)
        self.families_dir = self.output_dir / "families"
        self.metadata_dir = self.output_dir / "metadata"
        self.centroids_dir = self.output_dir / "centroids"
        
        # Create output directories
        self.families_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
        self.centroids_dir.mkdir(parents=True, exist_ok=True)
        
        # Storage for embeddings and metadata
        self.embeddings: Dict[str, np.ndarray] = {}
        self.family_mapping: Dict[str, str] = {}  # protein_id -> family_id
        self.family_proteins: Dict[str, List[str]] = {}  # family_id -> [protein_ids]
        
        if not HAS_FAISS:
            raise ImportError("faiss-cpu is required. Install with: pip install faiss-cpu")
    
    def load_embeddings(self) -> Tuple[Dict[str, np.ndarray], int]:
        """
        Load embeddings from TSV/CSV file.
        
        Expected format:
            uniprot_id\tembedding
            P12345\t0.1,0.2,0.3,...,0.9
        
        Returns:
            Tuple of (embeddings_dict, embedding_dimension)
        """
        logger.info(f"Loading embeddings from {self.embeddings_file}")
        
        if not self.embeddings_file.exists():
            raise FileNotFoundError(f"Embeddings file not found: {self.embeddings_file}")
        
        embeddings = {}
        dim = None
        
        with open(self.embeddings_file, 'r') as f:
            # Auto-detect delimiter
            sample = f.read(1024)
            f.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
            reader = csv.reader(f, dialect)
            
            # Check for header
            header_peek = sample.splitlines()[0].lower()
            has_header = "uniprot_id" in header_peek or "embedding" in header_peek
            
            if has_header:
                next(reader, None)
            
            for row_num, row in enumerate(reader, start=1):
                if not row or len(row) < 2:
                    continue
                
                uniprot_id = row[0].strip()
                emb_str = row[1].strip()
                
                if not uniprot_id or not emb_str:
                    continue
                
                try:
                    # Parse embedding (comma-separated floats)
                    vec = np.fromstring(emb_str, sep=",").astype(np.float32)
                    
                    # Validate dimension consistency
                    if dim is None:
                        dim = vec.shape[0]
                    elif vec.shape[0] != dim:
                        logger.warning(f"Row {row_num}: Inconsistent dimension {vec.shape[0]} (expected {dim})")
                        continue
                    
                    # Check for NaNs
                    if np.isnan(vec).any():
                        logger.warning(f"Row {row_num}: NaN values in embedding for {uniprot_id}")
                        continue
                    
                    embeddings[uniprot_id] = vec
                    
                except Exception as e:
                    logger.warning(f"Row {row_num}: Failed to parse embedding for {uniprot_id}: {e}")
                    continue
        
        logger.info(f"Loaded {len(embeddings)} embeddings with dimension {dim}")
        self.embeddings = embeddings
        return embeddings, dim
    
    def load_family_mapping(self, family_file: Optional[str] = None) -> Dict[str, str]:
        """
        Load protein to family mapping from file or create from metadata.
        
        Args:
            family_file: Optional JSON file with {protein_id: family_id} mapping
        
        Returns:
            Dictionary mapping protein_id to family_id
        """
        if family_file and Path(family_file).exists():
            logger.info(f"Loading family mapping from {family_file}")
            with open(family_file, 'r') as f:
                self.family_mapping = json.load(f)
        else:
            # If no family mapping provided, create default "all" family
            logger.info("No family mapping provided, using single 'all' family")
            self.family_mapping = {pid: "all" for pid in self.embeddings.keys()}
        
        # Build family -> proteins mapping
        self.family_proteins = {}
        for protein_id, family_id in self.family_mapping.items():
            if family_id not in self.family_proteins:
                self.family_proteins[family_id] = []
            self.family_proteins[family_id].append(protein_id)
        
        logger.info(f"Found {len(self.family_proteins)} families")
        for family_id, proteins in self.family_proteins.items():
            logger.info(f"  Family '{family_id}': {len(proteins)} proteins")
        
        return self.family_mapping
    
    def create_family_index(
        self,
        family_id: str,
        protein_ids: List[str],
        index_type: str = "auto"
    ) -> Tuple[str, Dict]:
        """
        Create FAISS index for a specific protein family.
        
        Args:
            family_id: Family identifier
            protein_ids: List of protein IDs in this family
            index_type: Index type ("flat", "ivf", or "auto")
        
        Returns:
            Tuple of (index_file_path, metadata_dict)
        """
        logger.info(f"Creating index for family '{family_id}' with {len(protein_ids)} proteins")
        
        # Get embeddings for this family
        family_embeddings = []
        valid_protein_ids = []
        
        for protein_id in protein_ids:
            if protein_id in self.embeddings:
                family_embeddings.append(self.embeddings[protein_id])
                valid_protein_ids.append(protein_id)
            else:
                logger.warning(f"Protein {protein_id} not found in embeddings, skipping")
        
        if not family_embeddings:
            raise ValueError(f"No valid embeddings found for family '{family_id}'")
        
        # Convert to numpy array
        embeddings_array = np.vstack(family_embeddings).astype(np.float32)
        n_proteins, dim = embeddings_array.shape
        
        logger.info(f"Family '{family_id}': {n_proteins} proteins, dimension {dim}")
        
        # Normalize embeddings for cosine similarity (L2 normalization)
        faiss.normalize_L2(embeddings_array)
        
        # Choose index type based on dataset size
        if index_type == "auto":
            if n_proteins < 10000:
                index_type = "flat"
            else:
                index_type = "ivf"
        
        # Create FAISS index
        if index_type == "flat":
            # Flat index: brute-force search, exact results
            index = faiss.IndexFlatIP(dim)  # Inner product for cosine similarity
            index.add(embeddings_array)
            nlist = None
            logger.info(f"Created Flat index for family '{family_id}'")
            
        elif index_type == "ivf":
            # IVF index: approximate search, faster for large datasets
            # nlist = number of clusters (cells)
            # Rule of thumb: nlist = sqrt(n_proteins), but at least 1 and at most 4096
            nlist = min(max(1, int(np.sqrt(n_proteins))), 4096)
            
            # FAISS requires at least 39 * nlist training points for IVF
            min_training_points = 39 * nlist
            if n_proteins < min_training_points:
                # Adjust nlist to meet training requirements
                nlist = max(1, n_proteins // 39)
                logger.info(f"Adjusted nlist to {nlist} for family '{family_id}' to meet training requirements")
            
            quantizer = faiss.IndexFlatIP(dim)
            index = faiss.IndexIVFFlat(quantizer, dim, nlist, faiss.METRIC_INNER_PRODUCT)
            
            # Train the index
            index.train(embeddings_array)
            index.add(embeddings_array)
            
            # Set number of probes for search (more probes = more accurate but slower)
            index.nprobe = min(10, nlist)  # Search in top 10 cells by default
            
            logger.info(f"Created IVF index for family '{family_id}' (nlist={nlist}, nprobe={index.nprobe})")
        else:
            raise ValueError(f"Unknown index type: {index_type}")
        
        # Save index
        file_safe_family_id = family_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
        index_file = self.families_dir / f"{file_safe_family_id}.faiss"
        faiss.write_index(index, str(index_file))
        logger.info(f"Saved index to {index_file}")
        
        # Create metadata
        metadata = {
            'family_id': family_id,
            'protein_ids': valid_protein_ids,
            'num_proteins': n_proteins,
            'embedding_dim': dim,
            'index_type': index_type,
            'nlist': nlist,
            'metric': 'cosine',
            'index_file': str(index_file.relative_to(self.output_dir))
        }
        
        metadata_file = self.metadata_dir / f"{file_safe_family_id}_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Saved metadata to {metadata_file}")
        
        return str(index_file), metadata
    
    def create_family_centroids_index(self) -> Tuple[str, Dict]:
        """
        Create FAISS index for family centroids.
        
        This index is used for hierarchical search:
        1. Query family centroids to find most similar family
        2. Then search within that family using family-specific index
        
        Returns:
            Tuple of (centroids_index_file, centroids_metadata)
        """
        logger.info("Creating family centroids index")
        
        if not self.family_proteins:
            raise ValueError("No families found. Load family mapping first.")
        
        # Compute centroids for each family
        centroids = []
        family_ids = []
        
        for family_id, protein_ids in self.family_proteins.items():
            # Get embeddings for this family
            family_embeddings = [
                self.embeddings[pid] for pid in protein_ids
                if pid in self.embeddings
            ]
            
            if not family_embeddings:
                logger.warning(f"No embeddings found for family '{family_id}', skipping")
                continue
            
            # Compute centroid (mean of all embeddings in family)
            centroid = np.mean(family_embeddings, axis=0).astype(np.float32)
            centroids.append(centroid)
            family_ids.append(family_id)
        
        if not centroids:
            raise ValueError("No valid family centroids computed")
        
        # Convert to numpy array
        centroids_array = np.vstack(centroids).astype(np.float32)
        n_families, dim = centroids_array.shape
        
        logger.info(f"Computed {n_families} family centroids with dimension {dim}")
        
        # Normalize for cosine similarity
        faiss.normalize_L2(centroids_array)
        
        # Create flat index for centroids (small number of families, exact search is fine)
        index = faiss.IndexFlatIP(dim)
        index.add(centroids_array)
        
        # Save centroids index
        centroids_index_file = self.centroids_dir / "family_centroids.faiss"
        faiss.write_index(index, str(centroids_index_file))
        logger.info(f"Saved centroids index to {centroids_index_file}")
        
        # Save centroids metadata
        centroids_metadata = {
            'family_ids': family_ids,
            'num_families': n_families,
            'embedding_dim': dim,
            'index_type': 'flat',
            'metric': 'cosine',
            'index_file': str(centroids_index_file.relative_to(self.output_dir))
        }
        
        centroids_metadata_file = self.centroids_dir / "family_centroids_metadata.json"
        with open(centroids_metadata_file, 'w') as f:
            json.dump(centroids_metadata, f, indent=2)
        
        logger.info(f"Saved centroids metadata to {centroids_metadata_file}")
        
        return str(centroids_index_file), centroids_metadata
    
    def create_all_indexes(
        self,
        family_mapping_file: Optional[str] = None,
        index_type: str = "auto"
    ) -> Dict[str, Dict]:
        """
        Create all FAISS indexes (family indexes + centroids index).
        
        Args:
            family_mapping_file: Optional file with protein to family mapping
            index_type: Index type for family indexes ("flat", "ivf", or "auto")
        
        Returns:
            Dictionary mapping family_id to metadata
        """
        # Load embeddings
        self.load_embeddings()
        
        # Load family mapping
        self.load_family_mapping(family_mapping_file)
        
        # Create indexes for each family
        family_indexes = {}
        for family_id, protein_ids in self.family_proteins.items():
            try:
                index_file, metadata = self.create_family_index(
                    family_id, protein_ids, index_type=index_type
                )
                family_indexes[family_id] = metadata
            except Exception as e:
                logger.error(f"Failed to create index for family '{family_id}': {e}")
                continue
        
        # Create family centroids index
        try:
            centroids_file, centroids_metadata = self.create_family_centroids_index()
            family_indexes['_centroids'] = centroids_metadata
        except Exception as e:
            logger.error(f"Failed to create centroids index: {e}")
        
        # Save family mapping file
        mapping_file = self.output_dir / "family_mapping.json"
        with open(mapping_file, 'w') as f:
            json.dump(self.family_mapping, f, indent=2)
        logger.info(f"Saved family mapping to {mapping_file}")
        
        # Save summary
        summary = {
            'total_proteins': len(self.embeddings),
            'total_families': len(self.family_proteins),
            'families': list(self.family_proteins.keys()),
            'indexes_created': len([k for k in family_indexes.keys() if k != '_centroids']),
            'embedding_dim': list(self.embeddings.values())[0].shape[0] if self.embeddings else 0
        }
        
        summary_file = self.output_dir / "indexing_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Saved indexing summary to {summary_file}")
        
        logger.info("=" * 60)
        logger.info("Indexing Complete!")
        logger.info(f"  Total proteins: {summary['total_proteins']}")
        logger.info(f"  Total families: {summary['total_families']}")
        logger.info(f"  Indexes created: {summary['indexes_created']}")
        logger.info(f"  Embedding dimension: {summary['embedding_dim']}")
        logger.info("=" * 60)
        
        return family_indexes


def main():
    """Main function to run the FAISS indexing guide."""
    parser = argparse.ArgumentParser(
        description="Create FAISS indexes for protein embeddings with family-based subsetting",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create indexes for all proteins (single family)
  python test/scripts/create_faiss_indexing_guide.py \\
      --embeddings data/embeddings/embeddings.tsv \\
      --output data/indexes
  
  # Create indexes with family grouping
  python test/scripts/create_faiss_indexing_guide.py \\
      --embeddings data/embeddings/embeddings.tsv \\
      --output data/indexes \\
      --families data/families/family_mapping.json
  
  # Create index for specific family only
  python test/scripts/create_faiss_indexing_guide.py \\
      --embeddings data/embeddings/embeddings.tsv \\
      --output data/indexes \\
      --family "kinase_enzyme"
  
  # Force IVF index type for all families
  python test/scripts/create_faiss_indexing_guide.py \\
      --embeddings data/embeddings/embeddings.tsv \\
      --output data/indexes \\
      --index-type ivf
        """
    )
    
    parser.add_argument(
        '--embeddings',
        type=str,
        required=True,
        help='Path to embeddings file (TSV format: uniprot_id\\tembedding)'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output directory for indexes and metadata'
    )
    parser.add_argument(
        '--families',
        type=str,
        default=None,
        help='Path to family mapping file (JSON: {protein_id: family_id})'
    )
    parser.add_argument(
        '--family',
        type=str,
        default=None,
        help='Create index for specific family only'
    )
    parser.add_argument(
        '--index-type',
        type=str,
        choices=['flat', 'ivf', 'auto'],
        default='auto',
        help='Index type: flat (exact), ivf (approximate), or auto (choose based on size)'
    )
    
    args = parser.parse_args()
    
    # Check if FAISS is available
    if not HAS_FAISS:
        logger.error("faiss-cpu is required but not installed.")
        logger.error("Install with: pip install faiss-cpu")
        return 1
    
    try:
        # Create guide instance
        guide = FAISSIndexingGuide(args.embeddings, args.output)
        
        if args.family:
            # Create index for specific family only
            guide.load_embeddings()
            guide.load_family_mapping(args.families)
            
            if args.family not in guide.family_proteins:
                logger.error(f"Family '{args.family}' not found in family mapping")
                logger.info(f"Available families: {list(guide.family_proteins.keys())}")
                return 1
            
            protein_ids = guide.family_proteins[args.family]
            index_file, metadata = guide.create_family_index(
                args.family, protein_ids, index_type=args.index_type
            )
            logger.info(f"Successfully created index for family '{args.family}'")
        else:
            # Create all indexes
            guide.create_all_indexes(
                family_mapping_file=args.families,
                index_type=args.index_type
            )
        
        return 0
        
    except Exception as e:
        logger.error(f"Failed to create indexes: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

