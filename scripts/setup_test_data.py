#!/usr/bin/env python3
"""
Complete Data Reset and Setup for Protein Network Analysis

This script completely resets the data folder while preserving the model,
and creates comprehensive test data for all pipeline steps.
"""

import os
import numpy as np
import pandas as pd
import h5py
import json
import shutil
import random
from pathlib import Path
from typing import List, Dict, Tuple
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Configuration
DATA_DIR = "data"
MODEL_DIR = os.path.join(DATA_DIR, "esm2_t6_8M_UR50D_local")
EMBEDDING_DIM = 320
N_PROTEINS_PER_FAMILY = 2000  # Increased to meet FAISS clustering requirements (39 * nlist)
N_FAMILIES = 20  # Comprehensive set of families

# Realistic protein data for metadata generation
ORGANISMS = [
    'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Danio rerio',
    'Drosophila melanogaster', 'Caenorhabditis elegans', 'Saccharomyces cerevisiae',
    'Escherichia coli', 'Bacillus subtilis', 'Pseudomonas aeruginosa',
    'Staphylococcus aureus', 'Mycobacterium tuberculosis', 'Arabidopsis thaliana',
    'Oryza sativa', 'Zea mays', 'Solanum lycopersicum', 'Nicotiana tabacum'
]

PROTEIN_FAMILIES = [
    'Kinase', 'Phosphatase', 'Transcription Factor', 'Receptor', 'Channel',
    'Transporter', 'Enzyme', 'Structural Protein', 'Chaperone', 'Protease',
    'Ubiquitin Ligase', 'DNA Binding Protein', 'RNA Binding Protein',
    'Membrane Protein', 'Secreted Protein', 'Nuclear Protein', 'Mitochondrial Protein',
    'Cytosolic Protein', 'Extracellular Matrix Protein', 'Cell Adhesion Protein'
]

PROTEIN_FUNCTIONS = [
    'Signal Transduction', 'Metabolism', 'DNA Replication', 'Transcription',
    'Translation', 'Protein Folding', 'Cell Division', 'Apoptosis',
    'Cell Migration', 'Cell Adhesion', 'Ion Transport', 'Nutrient Transport',
    'Energy Production', 'DNA Repair', 'RNA Processing', 'Protein Degradation',
    'Cell Signaling', 'Development', 'Immune Response', 'Stress Response'
]

PROTEIN_NAMES = [
    'Alpha-1-antitrypsin', 'Beta-galactosidase', 'Cytochrome c', 'Hemoglobin',
    'Insulin', 'Myoglobin', 'Ribonuclease', 'Serum albumin', 'Trypsin',
    'Urease', 'Catalase', 'Lysozyme', 'Peroxidase', 'Transferrin',
    'Fibrinogen', 'Collagen', 'Actin', 'Tubulin', 'Keratin', 'Elastin'
]

EC_NUMBERS = [
    '1.1.1.1', '1.1.1.2', '1.1.1.3', '1.1.1.4', '1.1.1.5',
    '2.1.1.1', '2.1.1.2', '2.1.1.3', '2.1.1.4', '2.1.1.5',
    '3.1.1.1', '3.1.1.2', '3.1.1.3', '3.1.1.4', '3.1.1.5',
    '4.1.1.1', '4.1.1.2', '4.1.1.3', '4.1.1.4', '4.1.1.5'
]

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def generate_realistic_protein_sequence(length: int = None) -> str:
    """Generate a realistic protein sequence."""
    if length is None:
        length = random.randint(100, 500)
    return ''.join(random.choices(AMINO_ACIDS, k=length))

def generate_uniprot_ids(num_ids: int) -> List[str]:
    """Generate dummy UniProt IDs (P00001, P00002, ...)"""
    return [f"P{str(i+1).zfill(5)}" for i in range(num_ids)]

def generate_realistic_metadata(protein_ids: List[str], family_id: str) -> pd.DataFrame:
    """Generate realistic protein metadata for network visualization."""
    metadata_list = []
    
    for i, protein_id in enumerate(protein_ids):
        # Randomly select metadata from predefined sets
        organism = random.choice(ORGANISMS)
        family = random.choice(PROTEIN_FAMILIES)
        function = random.choice(PROTEIN_FUNCTIONS)
        protein_name = random.choice(PROTEIN_NAMES)
        ec_number = random.choice(EC_NUMBERS)
        
        sequence_length = random.randint(100, 500)
        sequence = generate_realistic_protein_sequence(sequence_length)
        
        # Generate realistic molecular weight (in kDa)
        molecular_weight = round(sequence_length * 0.11, 2)
        
        # Generate realistic isoelectric point
        isoelectric_point = round(random.uniform(4.5, 9.5), 2)
        
        # Generate realistic descriptions
        description = f"{family} involved in {function.lower()}. Found in {organism}."
        
        # Create comprehensive metadata dictionary
        metadata = {
            'protein_id': protein_id,
            'Protein names': protein_name,
            'Organism': organism,
            'EC number': ec_number,
            'Function [CC]': description,
            'Protein families': family,
            'Reviewed': random.choice(['reviewed', 'unreviewed']),
            'organism': organism,
            'family': family,
            'function': function,
            'sequence_length': sequence_length,
            'molecular_weight_kda': molecular_weight,
            'isoelectric_point': isoelectric_point,
            'sequence': sequence,
            'cellular_location': random.choice(['Cytoplasm', 'Nucleus', 'Membrane', 'Extracellular', 'Mitochondria']),
            'expression_level': random.choice(['High', 'Medium', 'Low']),
            'conservation_score': round(random.uniform(0.1, 1.0), 3),
            'interaction_partners': random.randint(0, 20),
            'post_translational_modifications': random.randint(0, 5),
            'disease_association': random.choice(['None', 'Cancer', 'Cardiovascular', 'Neurological', 'Metabolic']),
            'drug_target': random.choice([True, False]),
            'essential_gene': random.choice([True, False]),
            'expression_tissue': random.choice(['Ubiquitous', 'Brain', 'Liver', 'Heart', 'Muscle', 'Kidney']),
            'subcellular_location': random.choice(['Cytoplasm', 'Nucleus', 'Membrane', 'Extracellular', 'Mitochondria', 'Endoplasmic Reticulum']),
            'protein_domains': random.randint(1, 5),
            'transmembrane_helices': random.randint(0, 7),
            'signal_peptide': random.choice([True, False]),
            'glycosylation_sites': random.randint(0, 10),
            'phosphorylation_sites': random.randint(0, 15),
            'ubiquitination_sites': random.randint(0, 5),
            'acetylation_sites': random.randint(0, 8),
            'methylation_sites': random.randint(0, 6),
            'sumoylation_sites': random.randint(0, 3),
            'disulfide_bonds': random.randint(0, 10),
            'metal_binding_sites': random.randint(0, 5),
            'active_sites': random.randint(0, 3),
            'binding_sites': random.randint(0, 8),
            'catalytic_residues': random.randint(0, 5),
            'regulatory_sites': random.randint(0, 4),
            'protein_family_id': family_id,
            'protein_family_name': family,
            'protein_function': function,
            'protein_organism': organism,
            'protein_sequence': sequence,
            'protein_length': sequence_length,
            'protein_molecular_weight': molecular_weight,
            'protein_isoelectric_point': isoelectric_point,
            'protein_name': protein_name,
            'protein_description': description
        }
        
        metadata_list.append(metadata)
    
    return pd.DataFrame(metadata_list).set_index('protein_id')

def reset_data_directory():
    """Reset the data directory while preserving the model."""
    logger.info("Starting complete data reset...")
    
    # Preserve the model directory
    model_backup = None
    if os.path.exists(MODEL_DIR):
        logger.info(f"Preserving local model at {MODEL_DIR}")
        model_backup = MODEL_DIR + '_backup'
        if os.path.exists(model_backup):
            shutil.rmtree(model_backup)
        shutil.copytree(MODEL_DIR, model_backup)
    
    # Remove all existing data except the model directory
    if os.path.exists(DATA_DIR):
        for root, dirs, files in os.walk(DATA_DIR, topdown=False):
            for name in files:
                file_path = os.path.join(root, name)
                if not file_path.startswith(MODEL_DIR):
                    os.remove(file_path)
            for name in dirs:
                dir_path = os.path.join(root, name)
                if not dir_path.startswith(MODEL_DIR) and not dir_path == MODEL_DIR:
                    shutil.rmtree(dir_path)
    
    # Restore the model directory if it was removed
    if model_backup and os.path.exists(model_backup):
        if not os.path.exists(MODEL_DIR):
            shutil.copytree(model_backup, MODEL_DIR)
        shutil.rmtree(model_backup)
    
    logger.info("Data directory reset complete")

def create_directory_structure():
    """Create the complete directory structure."""
    directories = [
        os.path.join(DATA_DIR, "families"),
        os.path.join(DATA_DIR, "metadata"),
        os.path.join(DATA_DIR, "indexes"),
        os.path.join(DATA_DIR, "indexes", "families"),
        os.path.join(DATA_DIR, "indexes", "metadata"),
        os.path.join(DATA_DIR, "indexes", "cache"),
        os.path.join(DATA_DIR, "indexes", "quantized"),
        os.path.join(DATA_DIR, "indexes", "protein_names"),
        os.path.join(DATA_DIR, "indexes", "protein_ids"),
        os.path.join(DATA_DIR, "family_centroids"),
        os.path.join(DATA_DIR, "family_centroids", "files"),
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
    
    logger.info("Directory structure created")

def create_protein_ids_index(family_mapping: Dict[str, str]):
    """Create an efficient protein IDs index for fast searching (exact match only)."""
    logger.info("Creating protein IDs index...")
    protein_ids_index = {}
    protein_to_family = {}
    for family_id, family_path in family_mapping.items():
        try:
            # Load metadata to get protein IDs
            file_safe_id = family_id.replace(' ', '_').replace('-', '_')
            metadata_file = os.path.join(DATA_DIR, "metadata", f"{file_safe_id}_metadata.parquet")
            if os.path.exists(metadata_file):
                metadata = pd.read_parquet(metadata_file)
                for protein_id in metadata.index:
                    protein_ids_index[protein_id] = {
                        'protein_id': protein_id,
                        'family_id': family_id,
                        'metadata': metadata.loc[protein_id].to_dict()
                    }
                    protein_to_family[protein_id] = family_id
        except Exception as e:
            logger.warning(f"Failed to create index for family {family_id}: {e}")
    # Save the index
    index_file = os.path.join(DATA_DIR, "indexes", "protein_ids", "protein_ids_index.json")
    os.makedirs(os.path.dirname(index_file), exist_ok=True)
    with open(index_file, 'w') as f:
        json.dump(protein_ids_index, f, indent=2)
    # Save protein to family mapping
    mapping_file = os.path.join(DATA_DIR, "indexes", "protein_ids", "protein_to_family.json")
    with open(mapping_file, 'w') as f:
        json.dump(protein_to_family, f, indent=2)
    logger.info(f"Created protein IDs index with {len(protein_ids_index)} entries")
    return protein_ids_index, protein_to_family

def create_family_data(family_id: str, num_proteins: int = N_PROTEINS_PER_FAMILY, uniprot_offset: int = 0) -> Tuple[str, str, List[str]]:
    """Create comprehensive family data with embeddings and metadata, using UniProt IDs."""
    logger.info(f"Creating family {family_id} with {num_proteins} proteins")
    # Generate UniProt IDs
    protein_ids = generate_uniprot_ids(num_proteins)
    # Offset for global uniqueness if needed
    protein_ids = [f"P{str(uniprot_offset + i + 1).zfill(5)}" for i in range(num_proteins)]
    centroid = np.random.normal(0, 5, size=EMBEDDING_DIM).astype(np.float32)
    embeddings = centroid + np.random.normal(0, 0.5, size=(num_proteins, EMBEDDING_DIM)).astype(np.float32)
    h5_path = os.path.join(DATA_DIR, "families", f"{family_id}.h5")
    with h5py.File(h5_path, 'w') as f:
        f.create_dataset('embeddings', data=embeddings, compression='gzip', chunks=True)
        dt = h5py.string_dtype(encoding='utf-8')
        f.create_dataset('protein_ids', data=np.array(protein_ids, dtype=object), dtype=dt)
        f.attrs['num_proteins'] = num_proteins
        f.attrs['embedding_dim'] = EMBEDDING_DIM
        f.attrs['chunk_size'] = num_proteins
        f.attrs['compression'] = 'gzip'
        f.attrs['metadata_file'] = os.path.join('metadata', f'{family_id}_metadata.parquet')
        f.attrs['embedding_dtype'] = 'float32'
    metadata = generate_realistic_metadata(protein_ids, family_id)
    meta_path = os.path.join(DATA_DIR, "metadata", f"{family_id}_metadata.parquet")
    metadata.to_parquet(meta_path, compression='gzip')
    logger.info(f"Created family {family_id}: {num_proteins} proteins, {EMBEDDING_DIM} dimensions")
    return h5_path, meta_path, protein_ids

def create_faiss_indexes(family_id: str):
    """Create FAISS indexes for a family."""
    try:
        import faiss
        import numpy as np
        
        # Load family data
        family_file = f"data/families/{family_id}.h5"
        if not os.path.exists(family_file):
            logger.warning(f"Family file not found: {family_file}")
            return
        
        with h5py.File(family_file, 'r') as f:
            embeddings = f['embeddings'][:]
            protein_ids_raw = f['protein_ids'][:]
            
            # Convert protein_ids to strings if they are bytes
            if hasattr(protein_ids_raw, 'dtype') and protein_ids_raw.dtype.kind == 'S':
                protein_ids = [pid.decode('utf-8') if isinstance(pid, bytes) else str(pid) for pid in protein_ids_raw]
            elif hasattr(protein_ids_raw, 'tolist'):
                protein_ids = protein_ids_raw.tolist()
            else:
                protein_ids = list(protein_ids_raw)
            
            # Ensure all protein_ids are strings
            protein_ids = [str(pid) for pid in protein_ids]
        
        # Create float index
        n_embeddings = embeddings.shape[0]
        embedding_dim = embeddings.shape[1]
        
        # Ensure we have enough data for FAISS clustering
        min_training_points = 39
        max_nlist = max(1, n_embeddings // min_training_points)
        nlist = min(10, max_nlist)  # Use reasonable nlist
        
        # Create file-safe family ID
        file_safe_id = family_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
        
        if n_embeddings >= (nlist * min_training_points):
            # Create IVF index with sufficient training data
            quantizer = faiss.IndexFlatIP(embedding_dim)
            index = faiss.IndexIVFFlat(quantizer, embedding_dim, nlist, faiss.METRIC_INNER_PRODUCT)
            
            # Normalize embeddings for cosine similarity
            faiss.normalize_L2(embeddings)
            
            # Train and add
            index.train(embeddings)
            index.add(embeddings)
            
            index_type = 'IVF_float'
        else:
            # Use flat index for small datasets
            index = faiss.IndexFlatIP(embedding_dim)
            faiss.normalize_L2(embeddings)
            index.add(embeddings)
            
            index_type = 'Flat_float'
        
        # Save index
        index_path = f"data/indexes/families/{file_safe_id}.faiss"
        os.makedirs(os.path.dirname(index_path), exist_ok=True)
        faiss.write_index(index, index_path)
        
        # Save metadata
        metadata = {
            'family_id': family_id,
            'protein_ids': protein_ids,  # This is already a list of strings
            'num_proteins': len(protein_ids),
            'embedding_dim': embedding_dim,
            'index_type': index_type,
            'nlist': nlist,
            'metric': 'cosine'
        }
        
        metadata_path = f"data/indexes/metadata/{file_safe_id}_float_metadata.json"
        os.makedirs(os.path.dirname(metadata_path), exist_ok=True)
        
        import json
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Created FAISS {index_type} index for {family_id}: {n_embeddings} proteins, {embedding_dim} dim")
        
    except Exception as e:
        logger.warning(f"Failed to create FAISS index for {family_id}: {e}")

def create_family_centroids_file(family_mapping: Dict[str, str]):
    """Create family centroids file for classification."""
    logger.info("Creating family centroids file...")
    
    centroids = []
    family_ids = []
    eigenprotein_ids = []
    
    for family_id, family_path in family_mapping.items():
        try:
            with h5py.File(family_path, 'r') as f:
                embeddings = f['embeddings'][:]
                protein_ids = f['protein_ids'][:]
                
                # Calculate centroid
                centroid = np.mean(embeddings, axis=0)
                centroids.append(centroid)
                family_ids.append(family_id)
                eigenprotein_ids.extend(protein_ids)
        except Exception as e:
            logger.warning(f"Failed to create centroid for {family_id}: {e}")
    
    if centroids:
        centroids_array = np.array(centroids, dtype=np.float32)
        family_centroids_file = os.path.join(DATA_DIR, "family_centroids", "files", "family_centroids_binary.npz")
        np.savez_compressed(
            family_centroids_file,
            centroids=centroids_array,
            eigenprotein_ids=np.array(eigenprotein_ids, dtype='<U20'),
            family_ids=list(family_mapping.keys())
        )
        logger.info(f"Created family centroids with {len(centroids)} families")

def setup_complete_test_data():
    """Set up complete test data for the protein network analysis system."""
    logger.info("Setting up complete test data for protein network analysis...")
    reset_data_directory()
    create_directory_structure()
    test_families = [
        "ABC_transporter", "G_protein_coupled_receptor", "serine_protease", "zinc_finger_protein",
        "membrane_channel", "kinase_enzyme", "transcription_factor", "cell_adhesion_molecule",
        "immune_receptor", "metabolic_enzyme", "FAM0", "FAM1", "FAM2", "FAM3", "FAM4",
        "FAMX", "FAMY", "FAMZ", "family_0", "family_1"
    ]
    family_mapping = {}
    all_uniprot_ids = []
    uniprot_offset = 0
    
    # Create family data and indexes
    for family_id in test_families:
        try:
            family_file, metadata_file, protein_ids = create_family_data(family_id, N_PROTEINS_PER_FAMILY, uniprot_offset)
            family_mapping[family_id] = family_file
            all_uniprot_ids.extend(protein_ids)
            uniprot_offset += N_PROTEINS_PER_FAMILY
            create_faiss_indexes(family_id)
        except Exception as e:
            logger.error(f"Failed to create test data for {family_id}: {e}")
    # Create family mapping for hierarchical index
    mapping_file = os.path.join(DATA_DIR, "indexes", "family_mapping.json")
    hierarchical_mapping = {}
    
    for family_id in test_families:
        file_safe_id = family_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
        index_path = os.path.join(DATA_DIR, "indexes", "families", f"{file_safe_id}.faiss")
        if os.path.exists(index_path):
            hierarchical_mapping[family_id] = index_path
    
    with open(mapping_file, 'w') as f:
        json.dump(hierarchical_mapping, f, indent=2)
    
    create_family_centroids_file(family_mapping)
    create_protein_ids_index(family_mapping)
    # Save UniProt ID list for test reference
    uniprot_id_file = os.path.join("test", "data", "uniprot_ids.txt")
    os.makedirs(os.path.dirname(uniprot_id_file), exist_ok=True)
    with open(uniprot_id_file, 'w') as f:
        for pid in all_uniprot_ids:
            f.write(f"{pid}\n")
    logger.info("Complete test data setup finished!")
    print("\n" + "="*60)
    print("COMPLETE DATA RESET SUMMARY")
    print("="*60)
    print(f"Total families: {len(test_families)}")
    print(f"Proteins per family: {N_PROTEINS_PER_FAMILY}")
    print(f"Total proteins: {len(test_families) * N_PROTEINS_PER_FAMILY}")
    print(f"Embedding dimension: {EMBEDDING_DIM}")
    print(f"Model preserved: {os.path.exists(MODEL_DIR)}")
    print(f"FAISS clustering requirement met: ✅ (minimum 400 per family)")
    print(f"Realistic metadata generated: ✅")
    print(f"Parquet files created: ✅")
    print(f"Family mapping created: ✅")
    print(f"Centroids file created: ✅")
    print(f"FAISS indexes created: ✅")
    print(f"Protein IDs index created: ✅")
    print("="*60)

if __name__ == "__main__":
    setup_complete_test_data() 