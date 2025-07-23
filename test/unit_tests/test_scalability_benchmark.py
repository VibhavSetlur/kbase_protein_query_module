import pytest
import numpy as np
import time
import tempfile
import os
import faiss
import matplotlib.pyplot as plt

@pytest.mark.skipif(not hasattr(faiss, 'IndexBinaryIVF'), reason="FAISS binary not available")
def test_faiss_binary_vs_float_scalability():
    # Settings
    SIZES = [1000, 5000, 10000, 20000, 50000, 100000]
    D = 320  # Embedding dimension (must be multiple of 8)
    N_QUERIES = 100
    np.random.seed(42)
    build_times_bin, build_times_float = [], []
    query_times_bin, query_times_float = [], []
    sizes_bin, sizes_float = [], []
    for N in SIZES:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Generate data
            xb_bin = np.random.randint(0, 256, size=(N, D // 8), dtype=np.uint8)
            xb_float = np.random.randn(N, D).astype(np.float32)
            xq_bin = np.random.randint(0, 256, size=(N_QUERIES, D // 8), dtype=np.uint8)
            xq_float = np.random.randn(N_QUERIES, D).astype(np.float32)
            nlist = min(100, N // 10) or 1
            # --- Binary index ---
            t0 = time.time()
            quantizer_bin = faiss.IndexBinaryFlat(D)
            index_bin = faiss.IndexBinaryIVF(quantizer_bin, D, nlist)
            index_bin.train(xb_bin)
            index_bin.add(xb_bin)
            build_time_bin = time.time() - t0
            bin_path = os.path.join(tmpdir, f"bin_{N}.faissbin")
            faiss.write_index(index_bin, bin_path)
            size_bin = os.path.getsize(bin_path) / 1024**2
            # Query
            t0 = time.time()
            for i in range(N_QUERIES):
                index_bin.search(xq_bin[i:i+1], 10)
            query_time_bin = (time.time() - t0) / N_QUERIES
            # --- Float index ---
            t0 = time.time()
            quantizer_float = faiss.IndexFlatL2(D)
            index_float = faiss.IndexIVFFlat(quantizer_float, D, nlist, faiss.METRIC_L2)
            index_float.train(xb_float)
            index_float.add(xb_float)
            build_time_float = time.time() - t0
            float_path = os.path.join(tmpdir, f"float_{N}.faiss")
            faiss.write_index(index_float, float_path)
            size_float = os.path.getsize(float_path) / 1024**2
            # Query
            t0 = time.time()
            for i in range(N_QUERIES):
                index_float.search(xq_float[i:i+1], 10)
            query_time_float = (time.time() - t0) / N_QUERIES
            # Record
            build_times_bin.append(build_time_bin)
            build_times_float.append(build_time_float)
            query_times_bin.append(query_time_bin)
            query_times_float.append(query_time_float)
            sizes_bin.append(size_bin)
            sizes_float.append(size_float)
            print(f"N={N:6d} | Bin: build {build_time_bin:.2f}s, query {query_time_bin*1e3:.2f}ms, size {size_bin:.2f}MB | Float: build {build_time_float:.2f}s, query {query_time_float*1e3:.2f}ms, size {size_float:.2f}MB")
    # Plot
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    axs[0].plot(SIZES, build_times_bin, 'o-', label='Binary')
    axs[0].plot(SIZES, build_times_float, 's-', label='Float')
    axs[0].set_title('Index Build Time (s)')
    axs[0].set_xlabel('Dataset Size')
    axs[0].set_ylabel('Time (s)')
    axs[0].legend()
    axs[1].plot(SIZES, query_times_bin, 'o-', label='Binary')
    axs[1].plot(SIZES, query_times_float, 's-', label='Float')
    axs[1].set_title('Query Time (ms)')
    axs[1].set_xlabel('Dataset Size')
    axs[1].set_ylabel('Avg Query Time (ms)')
    axs[1].legend()
    axs[2].plot(SIZES, sizes_bin, 'o-', label='Binary')
    axs[2].plot(SIZES, sizes_float, 's-', label='Float')
    axs[2].set_title('Index File Size (MB)')
    axs[2].set_xlabel('Dataset Size')
    axs[2].set_ylabel('File Size (MB)')
    axs[2].legend()
    plt.tight_layout()
    out_path = os.path.abspath('scalability_benchmark.png')
    plt.savefig(out_path)
    print(f"Scalability plot saved to {out_path}")
    # Print summary
    print("\nSummary Table:")
    print(f"{'N':>8} | {'Bin Build(s)':>10} | {'Float Build(s)':>12} | {'Bin Q(ms)':>9} | {'Float Q(ms)':>11} | {'Bin Size(MB)':>12} | {'Float Size(MB)':>14}")
    for i, N in enumerate(SIZES):
        print(f"{N:8d} | {build_times_bin[i]:10.2f} | {build_times_float[i]:12.2f} | {query_times_bin[i]*1e3:9.2f} | {query_times_float[i]*1e3:11.2f} | {sizes_bin[i]:12.2f} | {sizes_float[i]:14.2f}") 