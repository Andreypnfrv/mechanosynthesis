#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
import subprocess
from pathlib import Path
from google.cloud import storage
from ase.io import read, write
import numpy as np

def download_from_bucket(bucket_name, source_blob_name, destination_file):
    """Download file from GCS bucket"""
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file)
    print(f"Downloaded {source_blob_name} to {destination_file}")

def upload_to_bucket(bucket_name, source_file, destination_blob_name):
    """Upload file to GCS bucket"""
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_filename(source_file)
    print(f"Uploaded {source_file} to {destination_blob_name}")

def create_orca_input(atoms, method="B3LYP", basis="def2-SVP", calc_type="opt", nprocs=4):
    """Create Orca input file"""
    
    # Job specifications
    if calc_type == "opt":
        keywords = f"! {method} {basis} Opt"
    elif calc_type == "hessian":
        keywords = f"! {method} {basis} Freq"
    elif calc_type == "single":
        keywords = f"! {method} {basis}"
    
    # Add parallelization
    if nprocs > 1:
        keywords += f" PAL{nprocs}"
    
    input_text = f"""{keywords}

%maxcore 2000

* xyz 0 1
"""
    
    # Add coordinates
    for atom in atoms:
        pos = atom.position
        input_text += f"{atom.symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n"
    
    input_text += "*\n"
    
    return input_text

def run_orca_calculation(input_file, output_dir, method="B3LYP", basis="def2-SVP", calc_type="opt", nprocs=4):
    """Run Orca calculation"""
    
    # Read input structure
    atoms = read(input_file)
    
    # Create Orca input
    orca_input = create_orca_input(atoms, method, basis, calc_type, nprocs)
    
    # Write input file
    input_name = "calculation"
    with open(f"{input_name}.inp", "w") as f:
        f.write(orca_input)
    
    # Run Orca
    print(f"Starting Orca calculation: {calc_type}")
    try:
        result = subprocess.run(
            ["orca", f"{input_name}.inp"],
            capture_output=True,
            text=True,
            timeout=3600 * 6  # 6 hour timeout
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"Orca failed: {result.stderr}")
            
        print("Orca calculation completed successfully")
        
    except subprocess.TimeoutExpired:
        raise RuntimeError("Orca calculation timed out after 6 hours")
    
    # Process outputs
    output_files = [
        f"{input_name}.out",
        f"{input_name}.xyz",
        f"{input_name}.gbw",
        f"{input_name}.engrad"
    ]
    
    # Add frequency-specific files
    if calc_type == "hessian":
        output_files.extend([
            f"{input_name}.hess",
            f"{input_name}_property.txt"
        ])
    
    # Copy outputs to output directory
    for filename in output_files:
        if os.path.exists(filename):
            subprocess.run(["cp", filename, str(output_dir / filename)])
    
    # Extract final structure if optimization
    if calc_type == "opt" and os.path.exists(f"{input_name}.xyz"):
        final_atoms = read(f"{input_name}.xyz")
        write(output_dir / "optimized.xyz", final_atoms)
        write(output_dir / "optimized.extxyz", final_atoms)
    
    # Extract energy
    if os.path.exists(f"{input_name}.out"):
        with open(f"{input_name}.out", "r") as f:
            content = f.read()
            
        # Find final energy
        energy_lines = [line for line in content.split('\n') if 'FINAL SINGLE POINT ENERGY' in line]
        if energy_lines:
            energy = float(energy_lines[-1].split()[-1])
            with open(output_dir / "energy.txt", "w") as f:
                f.write(f"Final Energy: {energy:.8f} Hartree\n")

def main():
    parser = argparse.ArgumentParser(description="Orca Cloud Calculator")
    parser.add_argument("--bucket", required=True, help="GCS bucket name")
    parser.add_argument("--input-path", required=True, help="Input file path in bucket")
    parser.add_argument("--output-path", required=True, help="Output directory path in bucket")
    parser.add_argument("--method", default="B3LYP", help="DFT method")
    parser.add_argument("--basis", default="def2-SVP", help="Basis set")
    parser.add_argument("--calc-type", default="opt", choices=["opt", "hessian", "single"])
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs to use")
    
    args = parser.parse_args()
    
    # Force single CPU on certain systems to avoid MPI issues
    if os.uname().sysname == "Darwin":
        args.cpus = 1
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            os.chdir(temp_dir)
            temp_path = Path(temp_dir)
            
            # Download input file
            input_file = temp_path / "input.xyz"
            download_from_bucket(args.bucket, args.input_path, str(input_file))
            
            # Create output directory
            output_dir = temp_path / "outputs"
            output_dir.mkdir()
            
            # Run calculation
            run_orca_calculation(
                input_file, output_dir, 
                args.method, args.basis, args.calc_type, args.cpus
            )
            
            # Upload all output files
            for output_file in output_dir.iterdir():
                if output_file.is_file():
                    upload_to_bucket(
                        args.bucket, 
                        str(output_file), 
                        f"{args.output_path}/{output_file.name}"
                    )
            
            print("Calculation completed successfully!")
            
    except Exception as e:
        print(f"Error during calculation: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()