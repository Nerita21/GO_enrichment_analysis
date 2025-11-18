# GO Enrichment Analysis: Clustering Method Comparison

[![GitHub](https://img.shields.io/badge/GitHub-Nerita21/GO_enrichment_analysis-blue?logo=github)](https://github.com/Nerita21/GO_enrichment_analysis)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Version: 2.0.0](https://img.shields.io/badge/Version-2.0.0-blue.svg)
![Status: Active Development](https://img.shields.io/badge/Status-Active-brightgreen.svg)

## Overview

This project provides a **comprehensive framework for comparing different approaches to Gene Ontology (GO) term clustering** from enrichment analysis results. The pipeline evaluates multiple semantic similarity methods and clustering strategies across Python and R implementations.

### Key Features

- 🔄 **Parallel execution** of multiple clustering methods
- 📊 **Systematic comparison** with efficiency metrics
- 🐍 **Python pathway:** g:Profiler → GOATOOLS (Lin similarity)
- 📚 **R pathway:** clusterProfiler → GOSemSim (Wang & Lin similarity)
- ⚙️ **Custom clustering** alternatives to built-in methods
- 🐳 **Docker support** for reproducible environments
- 📈 **Professional reporting** with HTML/JSON outputs
- 📋 **Nextflow orchestration** for scalable execution

## Biological Context

miRNAs are critical post-transcriptional regulators of gene expression. When analyzing predicted or validated miRNA target genes, performing GO enrichment analysis reveals the biological processes, molecular functions, and cellular components they regulate.

However, GO enrichment results often contain many redundant terms describing highly related concepts. This project systematically compares methods that address this redundancy through semantic clustering, helping researchers identify the **major functional themes** regulated by their miRNA of interest.

## Workflows Compared

### R-based Approaches (clusterProfiler + GOSemSim)

| Method | Similarity | Approach | Use Case |
|--------|-----------|----------|----------|
| **Wang** | Graph-based topology | Hierarchical clustering + simplify() | Emphasizes structural relationships in GO DAG |
| **Lin** | Information Content | Hierarchical clustering + simplify() | Emphasizes term specificity and informativeness |
| **Custom** | Hybrid (Rel) + p-value | Custom hierarchical + representative selection | Transparent alternative to simplify() |

### Python-based Approaches (GOATOOLS)

| Method | Similarity | Approach | Use Case |
|--------|-----------|----------|----------|
| **Lin** | Information Content | Hierarchical clustering with IC threshold | MICA-based specificity ranking |
| **Custom** | Lin + p-value weight | Custom hierarchical with p-value weighting | Balances similarity and statistical significance |

### Key Differences Summary

| Aspect | R (clusterProfiler) | Python (GOATOOLS) |
|--------|-------------------|-------------------|
| **Enrichment Tool** | clusterProfiler | g:Profiler |
| **Similarity Methods** | Wang, Lin, Rel (Hybrid) | Lin (MICA-based) |
| **Clustering Algorithm** | Hierarchical + simplify() | Hierarchical + custom selection |
| **Labeling Strategy** | Common ancestor or top term | MICA-based representative |
| **Strength** | Built-in, well-tested | Fine-grained control, transparent |
| **Computational Speed** | Fast (pre-computed) | Moderate (IC-based calculation) |
| **Customization** | Limited without code changes | Highly flexible |

