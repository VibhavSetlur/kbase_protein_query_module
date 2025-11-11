FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer

# -----------------------------------------
# Install system dependencies for scientific computing
# -----------------------------------------

# Install required system dependencies for scientific computing
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

# Install Python dependencies (no cache to keep layers small)
# Install torch CPU-only first (smaller image size for containers, ~500MB vs ~2GB)
RUN pip install --upgrade --no-cache-dir pip setuptools wheel \
 && pip install --no-cache-dir --index-url https://download.pytorch.org/whl/cpu torch>=2.0.0 \
 && pip install --no-cache-dir -r requirements.txt

RUN make all

# Generate test data if it doesn't exist (embeddings and indexes)
# This ensures the data files are available in the Docker image
# Create data directory structure first
RUN mkdir -p /kb/module/data/embeddings /kb/module/data/indexes

# Generate test data - this requires internet access to download models and fetch from UniProt
# If data doesn't exist, generate it. If generation fails, the build will fail.
RUN if [ ! -f /kb/module/data/embeddings/embeddings.tsv ]; then \
    echo "Generating test data (this may take several minutes)..." && \
    python test/scripts/generate_test_data.py \
        --pfam_id PF00096 \
        --max_proteins 50 \
        --reviewed_only \
        --out_dir /kb/module/data || \
    (echo "ERROR: Failed to generate test data. Check logs above for details." && exit 1); \
    else \
    echo "Test data already exists, skipping generation"; \
    fi

# Verify data files exist - fail build if they don't exist
RUN if [ ! -f /kb/module/data/embeddings/embeddings.tsv ]; then \
    echo "ERROR: embeddings.tsv not found. Build cannot continue." && exit 1; \
    else \
    echo "SUCCESS: embeddings.tsv found at /kb/module/data/embeddings/embeddings.tsv" && \
    ls -lh /kb/module/data/embeddings/embeddings.tsv && \
    echo "File size: $(du -h /kb/module/data/embeddings/embeddings.tsv | cut -f1)"; \
    fi

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
