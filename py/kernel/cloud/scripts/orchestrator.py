import json
import os
from google.cloud import batch_v1
from google.cloud import storage
from google.cloud import functions_v1
import functions_framework

@functions_framework.http
def trigger_calculation(request):
    """Cloud Function to orchestrate DFT calculations"""
    
    # Parse request
    request_json = request.get_json(silent=True)
    if not request_json:
        return {"error": "No JSON payload provided"}, 400
    
    # Required parameters
    bucket_name = request_json.get('bucket')
    input_path = request_json.get('input_path')
    output_path = request_json.get('output_path', 'outputs')
    calculator = request_json.get('calculator', 'dftb')  # dftb or orca
    method = request_json.get('method', 'relax')
    cpu_count = request_json.get('cpu_count', 64)
    memory = request_json.get('memory', 128)
    
    if not bucket_name or not input_path:
        return {"error": "bucket and input_path are required"}, 400
    
    # Project settings
    project_id = os.environ.get('GCP_PROJECT')
    region = os.environ.get('FUNCTION_REGION', 'us-central1')
    
    try:
        # Create batch job
        job_id = create_batch_job(
            project_id, region, bucket_name, input_path, output_path,
            calculator, method, cpu_count, memory
        )
        
        return {
            "status": "success",
            "job_id": job_id,
            "message": f"Started {calculator} calculation job"
        }
        
    except Exception as e:
        return {"error": str(e)}, 500

def create_batch_job(project_id, region, bucket_name, input_path, output_path, 
                    calculator, method, cpu_count, memory):
    """Create Google Batch job for DFT calculation"""
    
    client = batch_v1.BatchServiceClient()
    
    # Container image
    if calculator == 'dftb':
        image = f"gcr.io/{project_id}/dft-runners/dftb:latest"
        args = [
            "python3", "calculate.py",
            "--bucket", bucket_name,
            "--input-path", input_path,
            "--output-path", output_path,
            "--method", method,
            "--cpus", str(cpu_count)
        ]
    elif calculator == 'orca':
        image = f"gcr.io/{project_id}/dft-runners/orca:latest"
        method_parts = method.split('-')
        calc_type = method_parts[0] if len(method_parts) > 0 else 'opt'
        dft_method = method_parts[1] if len(method_parts) > 1 else 'B3LYP'
        basis = method_parts[2] if len(method_parts) > 2 else 'def2-SVP'
        
        args = [
            "python3", "calculate.py",
            "--bucket", bucket_name,
            "--input-path", input_path,
            "--output-path", output_path,
            "--calc-type", calc_type,
            "--method", dft_method,
            "--basis", basis,
            "--cpus", str(cpu_count)
        ]
    else:
        raise ValueError(f"Unknown calculator: {calculator}")
    
    # Job configuration
    runnable = batch_v1.Runnable()
    runnable.container = batch_v1.Runnable.Container()
    runnable.container.image_uri = image
    runnable.container.commands = args
    
    # Task specification
    task_spec = batch_v1.TaskSpec()
    task_spec.runnables = [runnable]
    task_spec.compute_resource = batch_v1.ComputeResource()
    task_spec.compute_resource.cpu_milli = cpu_count * 1000
    task_spec.compute_resource.memory_mib = memory * 1024
    task_spec.max_retry_count = 1
    task_spec.max_run_duration = "21600s"  # 6 hours
    
    # Allocation policy (preemptible)
    allocation_policy = batch_v1.AllocationPolicy()
    instance_policy = batch_v1.AllocationPolicy.InstancePolicy()
    instance_policy.machine_type = f"c2-standard-{min(cpu_count, 60)}"
    instance_policy.provisioning_model = batch_v1.AllocationPolicy.ProvisioningModel.PREEMPTIBLE
    
    policy_template = batch_v1.AllocationPolicy.InstancePolicyOrTemplate()
    policy_template.policy = instance_policy
    allocation_policy.instances = [policy_template]
    
    # Job specification
    job_spec = batch_v1.JobSpec()
    job_spec.task_groups = [batch_v1.TaskGroup()]
    job_spec.task_groups[0].task_spec = task_spec
    job_spec.task_groups[0].task_count = 1
    job_spec.allocation_policy = allocation_policy
    
    # Create job
    job = batch_v1.Job()
    job.spec = job_spec
    
    # Job name
    job_name = f"dft-{calculator}-{input_path.replace('/', '-').replace('.', '-')}"[:63]
    
    # Submit job
    request = batch_v1.CreateJobRequest()
    request.parent = f"projects/{project_id}/locations/{region}"
    request.job_id = job_name
    request.job = job
    
    operation = client.create_job(request=request)
    
    print(f"Created batch job: {job_name}")
    return job_name

@functions_framework.http
def cleanup_resources(request):
    """Clean up old batch jobs and temporary files"""
    
    project_id = os.environ.get('GCP_PROJECT')
    region = os.environ.get('FUNCTION_REGION', 'us-central1')
    
    try:
        # Clean up old batch jobs
        client = batch_v1.BatchServiceClient()
        parent = f"projects/{project_id}/locations/{region}"
        
        for job in client.list_jobs(parent=parent):
            # Delete completed jobs older than 1 day
            if job.status.state in [batch_v1.JobStatus.State.SUCCEEDED, 
                                   batch_v1.JobStatus.State.FAILED]:
                client.delete_job(name=job.name)
                print(f"Deleted old job: {job.name}")
        
        return {"status": "cleanup completed"}
        
    except Exception as e:
        return {"error": str(e)}, 500