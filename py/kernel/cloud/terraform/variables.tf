variable "project_id" {
  description = "GCP Project ID"
  type        = string
}

variable "region" {
  description = "GCP Region"
  type        = string
  default     = "us-central1"
}

variable "zone" {
  description = "GCP Zone"
  type        = string
  default     = "us-central1-a"
}

variable "cpu_count" {
  description = "Number of CPUs for calculation"
  type        = number
  default     = 64
}

variable "memory" {
  description = "Memory in GB"
  type        = number
  default     = 128
}

variable "preemptible" {
  description = "Use preemptible instances"
  type        = bool
  default     = true
}

variable "bucket_name" {
  description = "Storage bucket for input/output files"
  type        = string
}