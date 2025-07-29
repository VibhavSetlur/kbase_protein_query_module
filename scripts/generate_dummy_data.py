import os
import numpy as np
import pandas as pd
import h5py
import faiss
import json
import pyarrow as pa
import pyarrow.parquet as pq
import shutil

# Configurable parameters
N_FAMILIES = 50
N_PROTEINS_PER_FAMILY = 500
MODEL_NAME = 'esm2_t6_8M_UR50D'  # Smallest ESM2 model
EMBEDDING_DIM = 320  # Embedding dimension for esm2_t6_8M_UR50D (must be multiple of 8)
DATA_DIR = 'data'

np.random.seed(42)

families_dir = os.path.join(DATA_DIR, 'families')
metadata_dir = os.path.join(DATA_DIR, 'metadata')
indexes_fam_dir = os.path.join(DATA_DIR, 'indexes', 'families')
indexes_meta_dir = os.path.join(DATA_DIR, 'indexes', 'metadata')
indexes_dir = os.path.join(DATA_DIR, 'indexes')
cache_dir = os.path.join(DATA_DIR, 'cache')

# Remove all existing data to avoid duplicates
if os.path.exists(DATA_DIR):
    for root, dirs, files in os.walk(DATA_DIR, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            shutil.rmtree(os.path.join(root, name))

# Create directories
for d in [families_dir, metadata_dir, indexes_fam_dir, indexes_meta_dir, indexes_dir, cache_dir]:
    os.makedirs(d, exist_ok=True)

family_mapping = {}

centroids = []
eigenprotein_ids = []
family_ids = []

# Generate distinct centroids for each family
FAMILY_CENTROIDS = np.random.normal(0, 5, size=(N_FAMILIES, EMBEDDING_DIM)).astype(np.float32)
FAMILY_SPREAD = 0.5  # Standard deviation for within-family variation

for fam in range(N_FAMILIES):
    family_id = f'family_{fam}'
    # Generate protein IDs - use both UniProt format and dummy format for testing
    if fam < 25:  # First half use UniProt-like format
        protein_ids = [f'P{fam:05d}{i:03d}' for i in range(N_PROTEINS_PER_FAMILY)]
    else:  # Second half use dummy format
        protein_ids = [f'{family_id}_prot_{i}' for i in range(N_PROTEINS_PER_FAMILY)]
    # Generate a unique centroid for this family
    centroid = FAMILY_CENTROIDS[fam]
    # Generate embeddings clustered around the centroid
    embeddings = centroid + np.random.normal(0, FAMILY_SPREAD, size=(N_PROTEINS_PER_FAMILY, EMBEDDING_DIM)).astype(np.float32)
    # Save HDF5 file
    h5_path = os.path.join(families_dir, f'{family_id}.h5')
    with h5py.File(h5_path, 'w') as f:
        f.create_dataset('embeddings', data=embeddings, compression='gzip', chunks=True)
        dt = h5py.string_dtype(encoding='utf-8')
        f.create_dataset('protein_ids', data=np.array(protein_ids, dtype=object), dtype=dt)
        f.attrs['num_proteins'] = N_PROTEINS_PER_FAMILY
        f.attrs['embedding_dim'] = EMBEDDING_DIM
        f.attrs['chunk_size'] = N_PROTEINS_PER_FAMILY
        f.attrs['compression'] = 'gzip'
        f.attrs['metadata_file'] = os.path.join('metadata', f'{family_id}_metadata.parquet')
        f.attrs['embedding_dtype'] = 'float32'
    # Create metadata DataFrame
    metadata = pd.DataFrame({
        'protein_id': protein_ids,
        'protein_name': [f'Protein {i}' for i in range(N_PROTEINS_PER_FAMILY)],
        'description': [f'Description for {pid}' for pid in protein_ids],
        'family_id': family_id
    }).set_index('protein_id')
    # Save as Parquet
    meta_path = os.path.join(metadata_dir, f'{family_id}_metadata.parquet')
    table = pa.Table.from_pandas(metadata)
    pq.write_table(table, meta_path)
    # Build FAISS IVF float index
    dimension = EMBEDDING_DIM
    nlist = min(100, N_PROTEINS_PER_FAMILY // 10) or 1
    quantizer = faiss.IndexFlatL2(dimension)
    index = faiss.IndexIVFFlat(quantizer, dimension, nlist, faiss.METRIC_L2)
    index.train(embeddings)
    index.add(embeddings)
    faiss_path = os.path.join(indexes_fam_dir, f'{family_id}.faiss')
    print(f"Writing float index for {family_id}: type={type(index)}, is_trained={getattr(index, 'is_trained', 'N/A')}, ntotal={getattr(index, 'ntotal', 'N/A')}")
    if hasattr(index, 'is_trained') and not index.is_trained:
        raise RuntimeError(f"Index for {family_id} is not trained!")
    if getattr(index, 'ntotal', 1) == 0:
        raise RuntimeError(f"Index for {family_id} is empty!")
    try:
        faiss.write_index(index, str(faiss_path))
    except Exception as e:
        print(f"Failed to write float index for {family_id}: {e}")
        raise
    # Update family mapping
    family_mapping[family_id] = f'families/{family_id}.faiss'
    # Compute centroid (mean of float32)
    centroid_actual = embeddings.mean(axis=0)
    centroid_bin = (centroid_actual > 0).astype(np.uint8)
    centroid_bin_packed = np.packbits(centroid_bin)
    # Find eigenprotein (closest to centroid in L2)
    dists = np.linalg.norm(embeddings - centroid_actual, axis=1)
    medoid_idx = np.argmin(dists)
    eigenprotein_id = protein_ids[medoid_idx]
    centroids.append(centroid_bin_packed)
    eigenprotein_ids.append(eigenprotein_id)
    family_ids.append(family_id)

# Save family mapping JSON
mapping_path = os.path.join(indexes_dir, 'family_mapping.json')
with open(mapping_path, 'w') as f:
    json.dump(family_mapping, f, indent=2)

# Save centroids and eigenprotein IDs
centroids = np.stack(centroids)
np.savez_compressed(os.path.join(DATA_DIR, 'family_centroids_binary.npz'),
         family_ids=np.array(family_ids),
         centroids=centroids,
         eigenprotein_ids=np.array(eigenprotein_ids))
print(f"Saved family centroids and eigenprotein IDs to '{os.path.join(DATA_DIR, 'family_centroids_binary.npz')}'")

print(f"Dummy data generated in '{DATA_DIR}' for {N_FAMILIES} families.") 