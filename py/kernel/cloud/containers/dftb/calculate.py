#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
from pathlib import Path
from google.cloud import storage
from ase.io import read, write
from ase.calculators.dftb import Dftb
from ase.optimize import BFGS
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

def run_dftb_calculation(input_file, output_dir, method="relax"):
    """Run DFTB+ calculation"""
    
    # Read input structure
    atoms = read(input_file)
    
    # Setup DFTB+ calculator
    calc = Dftb(
        label='dftb',
        Hamiltonian_SCC='Yes',
        Hamiltonian_SCCTolerance=1e-8,
        Hamiltonian_MaxAngularMomentum_='',
        Hamiltonian_MaxAngularMomentum_C='"p"',
        Hamiltonian_MaxAngularMomentum_H='"s"',
        Hamiltonian_MaxAngularMomentum_O='"p"',
        Hamiltonian_MaxAngularMomentum_N='"p"',
    )
    
    atoms.set_calculator(calc)
    
    if method == "relax":
        # Geometry optimization
        optimizer = BFGS(atoms)
        
        # Trajectory storage
        traj_file = output_dir / "trajectory.traj"
        optimizer.attach(lambda: write(traj_file, atoms, append=True))
        
        # Run optimization
        optimizer.run(fmax=0.05)
        
        # Save final structure
        write(output_dir / "optimized.xyz", atoms)
        write(output_dir / "optimized.extxyz", atoms)
        
    elif method == "hessian":
        # Calculate Hessian matrix
        from ase.vibrations import Vibrations
        
        vib = Vibrations(atoms)
        vib.run()
        vib.summary()
        
        # Save frequencies
        with open(output_dir / "frequencies.txt", "w") as f:
            f.write("Vibrational frequencies (cm-1):\n")
            for freq in vib.get_frequencies():
                f.write(f"{freq:.2f}\n")
    
    elif method == "single":
        # Single point calculation
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        
        with open(output_dir / "energy.txt", "w") as f:
            f.write(f"Energy: {energy:.6f} eV\n")
            f.write(f"Forces:\n{forces}\n")
    
    print(f"DFTB+ calculation completed: {method}")

def main():
    parser = argparse.ArgumentParser(description="DFTB+ Cloud Calculator")
    parser.add_argument("--bucket", required=True, help="GCS bucket name")
    parser.add_argument("--input-path", required=True, help="Input file path in bucket")
    parser.add_argument("--output-path", required=True, help="Output directory path in bucket")
    parser.add_argument("--method", default="relax", choices=["relax", "hessian", "single"])
    parser.add_argument("--cpus", type=int, default=os.cpu_count(), help="Number of CPUs to use")
    
    args = parser.parse_args()
    
    # Set number of threads
    os.environ["OMP_NUM_THREADS"] = str(args.cpus)
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Download input file
            input_file = temp_path / "input.xyz"
            download_from_bucket(args.bucket, args.input_path, str(input_file))
            
            # Create output directory
            output_dir = temp_path / "outputs"
            output_dir.mkdir()
            
            # Run calculation
            run_dftb_calculation(input_file, output_dir, args.method)
            
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