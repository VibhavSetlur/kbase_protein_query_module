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

# Install Python dependencies from requirements.txt
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Copy the module code
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

# Build the module
RUN make all


ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
