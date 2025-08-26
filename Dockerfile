FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer

# -----------------------------------------
# Install system dependencies for scientific computing
# -----------------------------------------

# Install required system dependencies for scientific computing
RUN apt-get update && apt-get install -y \
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

# Install Python dependencies
RUN pip install --upgrade pip setuptools wheel
RUN pip install -r requirements.txt

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
