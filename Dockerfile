FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# -----------------------------------------
# Fix for EOL Debian Stretch (cover all possible URLs)
RUN sed -i 's|http://deb.debian.org/debian|http://archive.debian.org/debian|g' /etc/apt/sources.list && \
    sed -i 's|http://httpredir.debian.org/debian|http://archive.debian.org/debian|g' /etc/apt/sources.list && \
    sed -i 's|http://security.debian.org|http://archive.debian.org|g' /etc/apt/sources.list && \
    sed -i '/stretch-updates/d' /etc/apt/sources.list && \
    sed -i '/stretch\/updates/d' /etc/apt/sources.list && \
    echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99no-check-valid-until
# -----------------------------------------

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --allow-unauthenticated --no-install-recommends \
    build-essential \
    git \
    libhdf5-dev \
    wget \
    unzip \
    libpcre3-dev \
    libssl-dev \
    autoconf \
    automake \
    libtool \
    pkg-config \
    cmake \
    bison \
    && rm -rf /var/lib/apt/lists/*


# Install SWIG
RUN cd /tmp && \
    wget https://github.com/swig/swig/archive/refs/tags/v4.0.2.zip && \
    unzip v4.0.2.zip && \
    cd swig-4.0.2 && \
    ./autogen.sh && \
    ./configure && \
    make -j4 && \
    make install && \
    cd / && rm -rf /tmp/*

# Upgrade pip and install build tools
RUN pip install --upgrade pip==21.3.1 setuptools==59.5.0 wheel setuptools-rust

# Install required dependencies
RUN pip install --no-cache-dir \

    # Deep Learning & Transformers
    torch==1.10.2 \
    transformers==4.31.0 \
    sentence-transformers==2.2.0 \

    # Data Science & Analysis
    numpy==1.21.0 \
    pandas==1.3.0 \
    scipy==1.7.3 \
    scikit-learn==1.0.2 \


    # Graphs & Networks
    networkx==2.6.0 \

    # Embedding Indexes
    faiss-cpu==1.7.0 \
    annoy==1.17.0 \

    # Plotting & Visualization
    matplotlib==3.5.0 \
    seaborn==0.11.0 \
    plotly==5.0.0 \
    kaleido==0.2.1 \

    # Data Storage & Formats
    h5py==3.6.0 \
    zarr==2.10.3 \
    tables==3.7.0 \
    pyarrow==6.0.0 \

    # Bioinformatics
    biopython==1.79 \

    # Other
    joblib==1.1.0 \
    tqdm==4.62.0 \
    pyyaml==6.0

# Copy the module code
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

# Build the module
RUN make all

# Run tests (fail build if tests fail)
RUN echo 'export PYTHONPATH=$(cd $script_dir/../lib && pwd):$PYTHONPATH' >> test/run_tests.sh
RUN echo 'cd $script_dir/../test' >> test/run_tests.sh
RUN ./scripts/run_tests.sh

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
