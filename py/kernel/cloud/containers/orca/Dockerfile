FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    wget \
    curl \
    gfortran \
    libopenmpi-dev \
    openmpi-bin \
    libscalapack-mpi-dev \
    liblapack-dev \
    libblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /app

# Install Python dependencies
RUN pip3 install \
    ase \
    numpy \
    scipy \
    matplotlib \
    google-cloud-storage \
    psutil

# Copy Orca binary and basis sets
COPY orca /usr/local/orca/
COPY basis_sets/ /usr/local/orca/basis/

# Set environment variables
ENV ORCA_ROOT=/usr/local/orca
ENV PATH="/usr/local/orca:$PATH"
ENV LD_LIBRARY_PATH="/usr/local/orca:$LD_LIBRARY_PATH"
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

# Copy calculation scripts
COPY calculate.py /app/
COPY utils.py /app/

# Make Orca executable
RUN chmod +x /usr/local/orca/orca*

# Default command
CMD ["python3", "calculate.py"]