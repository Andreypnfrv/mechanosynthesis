FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    wget \
    curl \
    gfortran \
    cmake \
    build-essential \
    libopenmpi-dev \
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

# Copy DFTB+ binary and slakos files
COPY dftb+ /usr/local/bin/dftb+
COPY slakos/ /app/slakos/

# Set environment variables
ENV SKDIR=/app/slakos
ENV PATH="/usr/local/bin:$PATH"
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

# Copy calculation scripts
COPY calculate.py /app/
COPY utils.py /app/

# Make executable
RUN chmod +x /usr/local/bin/dftb+

# Default command
CMD ["python3", "calculate.py"]