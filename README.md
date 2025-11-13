# Comparison of GO Term Clustering Approaches for miRNA Target Functional Analysis

## Overview
This project compares two complementary approaches for summarizing Gene Ontology (GO) enrichment results derived from predicted or validated miRNA target genes.
The goal is to reduce redundancy among enriched GO terms and to highlight biologically meaningful clusters that describe the main biological processes, molecular functions, and cellular components regulated by the studied miRNA (e.g. hsa-miR-107).

### Two analytical frameworks are compared:
- An R-based approach using the clusterProfiler ecosystem.
- A Python-based approach built upon GOATOOLS and custom information content (IC)–based clustering.

## Aim
To systematically evaluate and compare semantic clustering strategies for interpreting GO enrichment results from miRNA target genes, focusing on:
- How redundant GO terms can be grouped into functional clusters.
- How different methods (R vs Python) prioritize and label key biological processes.
- Which clustering representation best supports biological interpretation of miRNA regulatory roles.

## Methods Summary
| Step                         | R (clusterProfiler)                                                        | Python (GOATOOLS + custom clustering)                               |
| ---------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| **Input**                    | Enrichment results from `enrichGO()` or `gseGO()` using miRNA target genes | GO term list + ontology DAG (`obodag`) + p-values                   |
| **Semantic Similarity**      | `pairwise_termsim()` (Wang or Lin)                                         | Lin similarity via MICA (most informative common ancestor)          |
| **Clustering**               | Hierarchical clustering (`hclust`, average linkage)                        | Hierarchical clustering per namespace (`linkage`, `fcluster`)       |
| **Cluster Labeling**         | Ancestor GO term (via `simplify()` or common node)                         | MICA term with lowest p-value                                       |
| **Representative Selection** | Based on p.adjust or common ancestor                                       | Based on IC and lowest p-value                                      |
| **Cluster Prioritization**   | By significance (`p.adjust`)                                               | By cluster size (>5 terms) or significance                          |
| **Visualization**            | Barplot, dotplot, or `emapplot`                                            | Heatmap and dendrogram                                              |
| **Output**                   | Simplified enrichResult table                                              | Clustered DataFrame with namespace, size, representative GO ID/name |

## Biological Context
miRNAs play a key role in post-transcriptional regulation of gene expression.
Identifying functional clusters among enriched GO terms helps summarize major biological themes affected by miRNA targets.

## Comparison Focus
- Which approach yields clearer, more interpretable clusters?
- Are biologically similar terms grouped consistently between R and Python?
- Do clusters correspond to known biological functions of the studied miRNA?
- How does significance (p.adjust) relate to cluster size and semantic cohesion?

