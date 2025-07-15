FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    torch \
    transformers \
    tqdm \
    h5py \
    pandas \
    numpy \
    scikit-learn \
    pyyaml \
    jinja2 \
    networkx \
    plotly \
    faiss-cpu \
    annoy \
    biopython

# Copy the module code
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

# Build the module
RUN make all

# Run tests (fail build if tests fail)
RUN ./scripts/run_tests.sh

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
