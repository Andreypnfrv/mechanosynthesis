output "bucket_name" {
  description = "Name of the storage bucket"
  value       = google_storage_bucket.dft_storage.name
}

output "function_uri" {
  description = "URI of the orchestrator function"
  value       = google_cloudfunctions2_function.dft_orchestrator.service_config[0].uri
}

output "artifact_registry" {
  description = "Artifact Registry repository"
  value       = google_artifact_registry_repository.dft_repo.name
}

output "service_account_email" {
  description = "Service account email for batch jobs"
  value       = google_service_account.dft_compute.email
}

output "network_name" {
  description = "Network name for batch jobs"
  value       = google_compute_network.dft_network.name
}

output "subnet_name" {
  description = "Subnet name for batch jobs"
  value       = google_compute_subnetwork.dft_subnet.name
}