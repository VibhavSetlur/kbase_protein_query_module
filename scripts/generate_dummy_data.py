import os
import numpy as np
import pandas as pd
import h5py
import faiss
import json
import pyarrow as pa
import pyarrow.parquet as pq

# Configurable parameters
N_FAMILIES = 3
N_PROTEINS_PER_FAMILY = 10
MODEL_NAME = 'esm2_t6_8M_UR50D'  # Smallest ESM2 model
EMBEDDING_DIM = 320  # Embedding dimension for esm2_t6_8M_UR50D
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

for fam in range(N_FAMILIES):
    family_id = f'family_{fam}'
    # Generate protein IDs
    protein_ids = [f'{family_id}_prot_{i}' for i in range(N_PROTEINS_PER_FAMILY)]
    # Generate random embeddings
    embeddings = np.random.randn(N_PROTEINS_PER_FAMILY, EMBEDDING_DIM).astype(np.float32)
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
    # Build FAISS index
    index = faiss.IndexFlatL2(EMBEDDING_DIM)
    index.add(embeddings)
    faiss_path = os.path.join(indexes_fam_dir, f'{family_id}.faiss')
    faiss.write_index(index, faiss_path)
    # Update family mapping
    family_mapping[family_id] = f'families/{family_id}.faiss'

# Save family mapping JSON
mapping_path = os.path.join(indexes_dir, 'family_mapping.json')
with open(mapping_path, 'w') as f:
    json.dump(family_mapping, f, indent=2)

print(f"Dummy data generated in '{DATA_DIR}' for {N_FAMILIES} families.") 