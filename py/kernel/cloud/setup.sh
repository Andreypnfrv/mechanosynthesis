#!/bin/bash
# Setup script for Google Cloud DFT infrastructure

set -e

# Configuration
PROJECT_ID=${1:-"your-project-id"}
REGION=${2:-"us-central1"}
BUCKET_NAME="${PROJECT_ID}-dft-calculations"

echo "🚀 Setting up Google Cloud DFT infrastructure..."
echo "Project: $PROJECT_ID"
echo "Region: $REGION"
echo "Bucket: $BUCKET_NAME"

# Check if gcloud is installed
if ! command -v gcloud &> /dev/null; then
    echo "❌ gcloud CLI not found. Please install it first."
    exit 1
fi

# Check if terraform is installed
if ! command -v terraform &> /dev/null; then
    echo "❌ terraform not found. Please install it first."
    exit 1
fi

# Set project
echo "📋 Setting GCP project..."
gcloud config set project $PROJECT_ID

# Enable required APIs
echo "🔧 Enabling required APIs..."
gcloud services enable \
    cloudbuild.googleapis.com \
    batch.googleapis.com \
    storage.googleapis.com \
    cloudfunctions.googleapis.com \
    artifactregistry.googleapis.com \
    compute.googleapis.com

# Initialize Terraform
echo "🏗️  Initializing Terraform..."
cd terraform
terraform init

# Create terraform.tfvars
cat > terraform.tfvars << EOF
project_id  = "$PROJECT_ID"
region      = "$REGION"
bucket_name = "$BUCKET_NAME"
EOF

# Plan and apply infrastructure
echo "📋 Planning Terraform deployment..."
terraform plan

echo "🚀 Deploying infrastructure..."
terraform apply -auto-approve

# Get outputs
FUNCTION_URL=$(terraform output -raw function_uri)
BUCKET_NAME=$(terraform output -raw bucket_name)

cd ..

# Build containers
echo "🐳 Building Docker containers..."
gcloud builds submit --config=cloudbuild.yaml .

# Create orchestrator function zip
echo "📦 Packaging orchestrator function..."
cd scripts
zip orchestrator.zip orchestrator.py requirements.txt
cd ..

echo "✅ Setup completed successfully!"
echo ""
echo "📝 Next steps:"
echo "1. Copy your DFTB+ binary to containers/dftb/dftb+"
echo "2. Copy your slakos files to containers/dftb/slakos/"
echo "3. Copy your Orca binary to containers/orca/orca/"
echo "4. Copy basis sets to containers/orca/basis_sets/"
echo "5. Rebuild containers: gcloud builds submit --config=cloudbuild.yaml ."
echo ""
echo "🚀 Usage example:"
echo "python3 cloud_runner.py \\"
echo "  --project $PROJECT_ID \\"
echo "  --bucket $BUCKET_NAME \\"
echo "  --function-url $FUNCTION_URL \\"
echo "  --input molecule.xyz \\"
echo "  --calculator dftb \\"
echo "  --method relax \\"
echo "  --cleanup"