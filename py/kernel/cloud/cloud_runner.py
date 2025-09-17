#!/usr/bin/env python3
"""
Cloud DFT calculation runner with automatic download and cleanup
"""
import os
import sys
import json
import time
import argparse
import requests
from pathlib import Path
from google.cloud import storage
from google.cloud import batch_v1

class CloudDFTRunner:
    def __init__(self, project_id, region='us-central1'):
        self.project_id = project_id
        self.region = region
        self.storage_client = storage.Client()
        self.batch_client = batch_v1.BatchServiceClient()
        
    def upload_input_file(self, bucket_name, local_file, remote_path):
        """Upload input file to GCS bucket"""
        bucket = self.storage_client.bucket(bucket_name)
        blob = bucket.blob(remote_path)
        blob.upload_from_filename(local_file)
        print(f"Uploaded {local_file} to gs://{bucket_name}/{remote_path}")
        
    def trigger_calculation(self, function_url, **kwargs):
        """Trigger cloud calculation via HTTP function"""
        response = requests.post(function_url, json=kwargs)
        response.raise_for_status()
        return response.json()
        
    def wait_for_job_completion(self, job_name, poll_interval=60):
        """Wait for batch job to complete"""
        print(f"Waiting for job {job_name} to complete...")
        
        job_path = f"projects/{self.project_id}/locations/{self.region}/jobs/{job_name}"
        
        while True:
            try:
                job = self.batch_client.get_job(name=job_path)
                state = job.status.state
                
                if state == batch_v1.JobStatus.State.SUCCEEDED:
                    print("✅ Job completed successfully!")
                    return True
                elif state == batch_v1.JobStatus.State.FAILED:
                    print("❌ Job failed!")
                    # Print failure details
                    for task_group in job.status.task_groups.values():
                        for task in task_group.task_status:
                            if task.state == batch_v1.TaskStatus.State.FAILED:
                                print(f"Task failed: {task}")
                    return False
                elif state == batch_v1.JobStatus.State.RUNNING:
                    print("🔄 Job is running...")
                else:
                    print(f"Job status: {state}")
                    
                time.sleep(poll_interval)
                
            except Exception as e:
                print(f"Error checking job status: {e}")
                time.sleep(poll_interval)
                
    def download_results(self, bucket_name, remote_path, local_dir):
        """Download all results from GCS bucket"""
        local_path = Path(local_dir)
        local_path.mkdir(parents=True, exist_ok=True)
        
        bucket = self.storage_client.bucket(bucket_name)
        blobs = bucket.list_blobs(prefix=remote_path)
        
        downloaded_files = []
        for blob in blobs:
            if blob.name.endswith('/'):  # Skip directories
                continue
                
            # Create local file path
            relative_path = blob.name[len(remote_path):].lstrip('/')
            local_file = local_path / relative_path
            local_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Download file
            blob.download_to_filename(str(local_file))
            downloaded_files.append(str(local_file))
            print(f"Downloaded: {blob.name} -> {local_file}")
            
        return downloaded_files
        
    def cleanup_cloud_resources(self, bucket_name, remote_path, job_name):
        """Clean up cloud resources to minimize costs"""
        try:
            # Delete input/output files from bucket
            bucket = self.storage_client.bucket(bucket_name)
            blobs = list(bucket.list_blobs(prefix=remote_path))
            
            for blob in blobs:
                blob.delete()
            print(f"🗑️  Cleaned up {len(blobs)} files from bucket")
            
            # Delete batch job
            job_path = f"projects/{self.project_id}/locations/{self.region}/jobs/{job_name}"
            self.batch_client.delete_job(name=job_path)
            print(f"🗑️  Deleted batch job: {job_name}")
            
        except Exception as e:
            print(f"Warning: Cleanup failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Run DFT calculations in Google Cloud")
    parser.add_argument("--project", required=True, help="GCP project ID")
    parser.add_argument("--bucket", required=True, help="GCS bucket name")
    parser.add_argument("--function-url", required=True, help="Cloud Function trigger URL")
    parser.add_argument("--input", required=True, help="Local input file path")
    parser.add_argument("--output-dir", default="./cloud_results", help="Local output directory")
    parser.add_argument("--calculator", default="dftb", choices=["dftb", "orca"])
    parser.add_argument("--method", default="relax", help="Calculation method")
    parser.add_argument("--cpu-count", type=int, default=64, help="Number of CPUs")
    parser.add_argument("--memory", type=int, default=128, help="Memory in GB")
    parser.add_argument("--cleanup", action="store_true", help="Clean up cloud resources after download")
    parser.add_argument("--poll-interval", type=int, default=60, help="Job status polling interval (seconds)")
    
    args = parser.parse_args()
    
    # Initialize runner
    runner = CloudDFTRunner(args.project)
    
    # Generate unique paths
    input_name = Path(args.input).name
    timestamp = int(time.time())
    remote_input_path = f"inputs/{timestamp}_{input_name}"
    remote_output_path = f"outputs/{timestamp}"
    
    try:
        # Step 1: Upload input file
        print("📤 Uploading input file...")
        runner.upload_input_file(args.bucket, args.input, remote_input_path)
        
        # Step 2: Trigger calculation
        print("🚀 Starting cloud calculation...")
        response = runner.trigger_calculation(
            args.function_url,
            bucket=args.bucket,
            input_path=remote_input_path,
            output_path=remote_output_path,
            calculator=args.calculator,
            method=args.method,
            cpu_count=args.cpu_count,
            memory=args.memory
        )
        
        job_id = response['job_id']
        print(f"✅ Job started: {job_id}")
        
        # Step 3: Wait for completion
        success = runner.wait_for_job_completion(job_id, args.poll_interval)
        if not success:
            print("❌ Calculation failed!")
            sys.exit(1)
            
        # Step 4: Download results
        print("📥 Downloading results...")
        downloaded_files = runner.download_results(
            args.bucket, remote_output_path, args.output_dir
        )
        
        print(f"✅ Downloaded {len(downloaded_files)} files to {args.output_dir}")
        
        # Step 5: Cleanup (optional)
        if args.cleanup:
            print("🧹 Cleaning up cloud resources...")
            runner.cleanup_cloud_resources(args.bucket, remote_input_path.split('/')[0], job_id)
            runner.cleanup_cloud_resources(args.bucket, remote_output_path.split('/')[0], job_id)
            
        print("🎉 Cloud calculation completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()