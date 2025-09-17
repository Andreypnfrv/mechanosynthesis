terraform {
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 4.0"
    }
  }
}

provider "google" {
  project = var.project_id
  region  = var.region
  zone    = var.zone
}

# Storage bucket for input/output files
resource "google_storage_bucket" "dft_storage" {
  name     = var.bucket_name
  location = var.region
  
  lifecycle_rule {
    action {
      type = "Delete"
    }
    condition {
      age = 30
    }
  }
}

# Container Registry for Docker images
resource "google_artifact_registry_repository" "dft_repo" {
  location      = var.region
  repository_id = "dft-runners"
  description   = "Docker repository for DFT calculation containers"
  format        = "DOCKER"
}

# Service Account for Batch Jobs
resource "google_service_account" "dft_compute" {
  account_id   = "dft-compute"
  display_name = "DFT Computation Service Account"
}

# IAM bindings
resource "google_project_iam_member" "storage_admin" {
  project = var.project_id
  role    = "roles/storage.admin"
  member  = "serviceAccount:${google_service_account.dft_compute.email}"
}

resource "google_project_iam_member" "batch_agent" {
  project = var.project_id
  role    = "roles/batch.agentReporter"
  member  = "serviceAccount:${google_service_account.dft_compute.email}"
}

# Network for batch jobs
resource "google_compute_network" "dft_network" {
  name                    = "dft-network"
  auto_create_subnetworks = false
}

resource "google_compute_subnetwork" "dft_subnet" {
  name          = "dft-subnet"
  ip_cidr_range = "10.0.0.0/16"
  region        = var.region
  network       = google_compute_network.dft_network.id
}

# Firewall rules
resource "google_compute_firewall" "allow_internal" {
  name    = "dft-allow-internal"
  network = google_compute_network.dft_network.name

  allow {
    protocol = "icmp"
  }

  allow {
    protocol = "tcp"
    ports    = ["0-65535"]
  }

  allow {
    protocol = "udp"
    ports    = ["0-65535"]
  }

  source_ranges = ["10.0.0.0/16"]
}

# Cloud Function for job orchestration
resource "google_storage_bucket" "function_source" {
  name     = "${var.bucket_name}-functions"
  location = var.region
}

resource "google_storage_bucket_object" "function_zip" {
  name   = "orchestrator.zip"
  bucket = google_storage_bucket.function_source.name
  source = "../scripts/orchestrator.zip"
}

resource "google_cloudfunctions2_function" "dft_orchestrator" {
  name        = "dft-orchestrator"
  location    = var.region
  description = "Orchestrates DFT calculations"

  build_config {
    runtime     = "python39"
    entry_point = "trigger_calculation"
    source {
      storage_source {
        bucket = google_storage_bucket.function_source.name
        object = google_storage_bucket_object.function_zip.name
      }
    }
  }

  service_config {
    max_instance_count = 10
    available_memory   = "256M"
    timeout_seconds    = 300
    service_account_email = google_service_account.dft_compute.email
  }
}

# Cloud Scheduler for automatic cleanup
resource "google_cloud_scheduler_job" "cleanup_job" {
  name             = "dft-cleanup"
  description      = "Clean up old batch jobs and files"
  schedule         = "0 2 * * *"
  time_zone        = "UTC"
  region           = var.region

  http_target {
    http_method = "POST"
    uri         = google_cloudfunctions2_function.dft_orchestrator.service_config[0].uri
    
    body = base64encode(jsonencode({
      action = "cleanup"
    }))
  }
}