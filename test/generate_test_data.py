"""
Generate synthetic data for pipeline (simple and deterministic).

Creates:
- data/embeddings/embeddings.tsv (uniprot_id \t comma-separated floats)
- data/indexes/id_to_row.json (uniprot_id -> row index)
- data/indexes/binary_embeddings.npy (uint8 matrix: value>0 -> 1 else 0)

Idempotent: overwrite existing files.
"""

import os
import argparse
import shutil
import json
from typing import List, Dict

import numpy as np
import sys

# Ensure project lib/ is on path for module imports when run directly
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir))
LIB_SRC_DIR = os.path.join(PROJECT_ROOT, "lib", "kbase_protein_query_module", "src")
if LIB_SRC_DIR not in sys.path:
    sys.path.insert(0, LIB_SRC_DIR)


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def write_embeddings(out_dir: str, n: int, dim: int, seed: int) -> str:
    rng = np.random.RandomState(seed)
    ids: List[str] = []
    vecs: List[np.ndarray] = []
    for i in range(n):
        pid = f"uniprot_{i:04d}"
        vec = rng.normal(0.0, 1.0, size=(dim,)).astype(np.float32)
        ids.append(pid)
        vecs.append(vec)

    emb_dir = os.path.join(out_dir, "embeddings")
    ensure_dir(emb_dir)
    emb_path = os.path.join(emb_dir, "embeddings.tsv")
    with open(emb_path, "w") as f:
        f.write("uniprot_id\tembedding\n")
        for pid, v in zip(ids, vecs):
            f.write(f"{pid}\t{','.join(f'{x:.6f}' for x in v.tolist())}\n")
    return emb_path


def write_indexes(out_dir: str, ids: List[str], vecs: List[np.ndarray]) -> Dict[str, str]:
    idx_dir = os.path.join(out_dir, "indexes")
    ensure_dir(idx_dir)
    # id_to_row
    id_map = {pid: i for i, pid in enumerate(ids)}
    id_map_path = os.path.join(idx_dir, "id_to_row.json")
    with open(id_map_path, "w") as f:
        json.dump(id_map, f)
    # binary embeddings (value>0 -> 1 else 0)
    mat = np.vstack(vecs).astype(np.float32)
    binary = (mat > 0).astype(np.uint8)
    bin_path = os.path.join(idx_dir, "binary_embeddings.npy")
    np.save(bin_path, binary)
    return {"id_to_row": id_map_path, "binary_embeddings": bin_path}


def build_simple_indexes(out_dir: str, embeddings_file: str) -> Dict[str, str]:
    # Parse embeddings file directly to avoid package init
    ids: List[str] = []
    vecs: List[np.ndarray] = []
    with open(embeddings_file, "r") as f:
        _ = f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ids.append(parts[0])
            vecs.append(np.fromstring(parts[1], sep=",").astype(np.float32))
    return write_indexes(out_dir, ids, vecs)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--dim", type=int, default=320)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--out_dir", default="data")
    args = parser.parse_args()

    ensure_dir(args.out_dir)
    ensure_dir(os.path.join(args.out_dir, "indexes"))
    ensure_dir(os.path.join(args.out_dir, "embeddings"))

    emb_path = write_embeddings(args.out_dir, args.n, args.dim, args.seed)
    idx_paths = build_simple_indexes(args.out_dir, emb_path)

    print("embeddings_count:", args.n)
    print("embeddings_file:", emb_path)
    print("id_to_row:", idx_paths.get("id_to_row"))
    print("binary_embeddings:", idx_paths.get("binary_embeddings"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


