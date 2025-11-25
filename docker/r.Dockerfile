# R-only image for GO enrichment pipeline
FROM rocker/tidyverse:latest

WORKDIR /app

# System deps for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy R requirements installer and run it (expects packages/requirements/install_requirements.R to exist)
COPY packages/requirements/install_requirements.R /app/install_requirements.R
RUN Rscript /app/install_requirements.R || true

# Create expected directories
RUN mkdir -p /app/results/R_based /app/results/comparison /app/data/processed

WORKDIR /app

ENV R_LIBS_USER=/usr/local/lib/R/site-library

ENTRYPOINT ["/bin/bash"]
