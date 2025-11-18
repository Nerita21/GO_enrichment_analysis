#!/usr/bin/env bash
# Simple build helper to create the per-language images used by Nextflow
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

echo "Building Python image: go-enrichment-python:latest"
docker build -f docker/python.Dockerfile -t go-enrichment-python:latest .

echo "Building R image: go-enrichment-r:latest"
docker build -f docker/r.Dockerfile -t go-enrichment-r:latest .

echo "Built images: go-enrichment-python:latest, go-enrichment-r:latest"
