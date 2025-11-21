# Python-only image for GO enrichment pipeline (changed to bullseye for stability)
FROM python:3.11-slim-bullseye

WORKDIR /app

# System deps needed by Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    libxml2-dev \
    libffi-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy and install Python requirements
COPY packages/requirements/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r /app/requirements.txt

# Create results folder
RUN mkdir -p /app/results/Python_based /app/results/comparison /app/data/processed

WORKDIR /app

ENV PYTHONUNBUFFERED=1
ENV PATH="/app/scripts:${PATH}"

ENTRYPOINT ["/bin/bash", "-c"]

