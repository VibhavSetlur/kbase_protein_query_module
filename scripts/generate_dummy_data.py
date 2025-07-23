import os
import numpy as np
import pandas as pd
import h5py
import faiss
import json
import pyarrow as pa
import pyarrow.parquet as pq

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

# Create directories
for d in [families_dir, metadata_dir, indexes_fam_dir, indexes_meta_dir, indexes_dir, cache_dir]:
    os.makedirs(d, exist_ok=True)

family_mapping = {}

centroids = []
eigenprotein_ids = []
family_ids = []

for fam in range(N_FAMILIES):
    family_id = f'family_{fam}'
    # Generate protein IDs
    protein_ids = [f'{family_id}_prot_{i}' for i in range(N_PROTEINS_PER_FAMILY)]
    # Generate random binary embeddings (np.uint8, shape [N, D//8])
    embeddings = np.random.randint(0, 256, size=(N_PROTEINS_PER_FAMILY, EMBEDDING_DIM // 8), dtype=np.uint8)
    # Save HDF5 file
    h5_path = os.path.join(families_dir, f'{family_id}.h5')
    with h5py.File(h5_path, 'w') as f:
        f.create_dataset('embeddings', data=embeddings, compression='gzip', chunks=True)
        dt = h5py.string_dtype(encoding='utf-8')
        f.create_dataset('protein_ids', data=np.array(protein_ids, dtype=object), dtype=dt)
        f.attrs['num_proteins'] = N_PROTEINS_PER_FAMILY
        f.attrs['embedding_dim'] = EMBEDDING_DIM // 8
        f.attrs['chunk_size'] = N_PROTEINS_PER_FAMILY
        f.attrs['compression'] = 'gzip'
        f.attrs['metadata_file'] = os.path.join('metadata', f'{family_id}_metadata.parquet')
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
    # Build FAISS IVF binary index (use flat for small families)
    dimension = EMBEDDING_DIM
    nlist = min(100, N_PROTEINS_PER_FAMILY // 10) or 1
    # FAISS IVF requires at least 39*nlist training points
    if N_PROTEINS_PER_FAMILY < 39 * nlist:
        # Too few points for IVF, use flat index
        index = faiss.IndexBinaryFlat(dimension)
        index.add(embeddings)
    else:
        quantizer = faiss.IndexBinaryFlat(dimension)
        index = faiss.IndexBinaryIVF(quantizer, dimension, nlist)
        index.train(embeddings)
        index.add(embeddings)
    faiss_path = os.path.join(indexes_fam_dir, f'{family_id}.faissbin')
    print(f"Writing index for {family_id}: type={type(index)}, is_trained={getattr(index, 'is_trained', 'N/A')}, ntotal={getattr(index, 'ntotal', 'N/A')}")
    if hasattr(index, 'is_trained') and not index.is_trained:
        raise RuntimeError(f"Index for {family_id} is not trained!")
    if getattr(index, 'ntotal', 1) == 0:
        raise RuntimeError(f"Index for {family_id} is empty!")
    try:
        faiss.write_index_binary(index, str(faiss_path))
    except Exception as e:
        print(f"Failed to write index for {family_id}: {e}")
        raise
    # Update family mapping
    family_mapping[family_id] = f'families/{family_id}.faissbin'
    # Compute centroid (bitwise majority)
    centroid = np.unpackbits(embeddings, axis=1).mean(axis=0) > 0.5
    centroid = centroid.astype(np.uint8)
    # Find eigenprotein (closest to centroid in Hamming distance)
    unpacked = np.unpackbits(embeddings, axis=1)
    dists = np.sum(unpacked != centroid, axis=1)
    medoid_idx = np.argmin(dists)
    eigenprotein_id = protein_ids[medoid_idx]
    centroids.append(centroid)
    eigenprotein_ids.append(eigenprotein_id)
    family_ids.append(family_id)

# Save family mapping JSON
mapping_path = os.path.join(indexes_dir, 'family_mapping.json')
with open(mapping_path, 'w') as f:
    json.dump(family_mapping, f, indent=2)

# Save centroids and eigenprotein IDs
centroids = np.stack(centroids)
np.savez(os.path.join(DATA_DIR, 'family_centroids.npz'),
         family_ids=np.array(family_ids),
         centroids=centroids,
         eigenprotein_ids=np.array(eigenprotein_ids))
print(f"Saved family centroids and eigenprotein IDs to '{os.path.join(DATA_DIR, 'family_centroids.npz')}'")

print(f"Dummy data generated in '{DATA_DIR}' for {N_FAMILIES} families.") 