# This script is used to test clusterprofiler package on one miRNA (hsa-mir-107) predicted targets with significant contextscore (<(-2), non-targetscan scale!)


# Load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install(c("GO.db", "AnnotationDbi", "org.Hs.eg.db"))

library("clusterProfiler", "AnnotationDbi")
library("org.Hs.eg.db")

# Load example dataset
gene_set_hsa_mir_107 <- read.csv("miRNA_target_genes_ex.csv", header = TRUE, stringsAsFactors = FALSE
)

# Exact genes
gene_list <- gene_set_hsa_mir_107$GeneSymbol

GO_result <- enrichGO(gene= gene_list, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

GO_result_df <- as.data.frame(GO_result)

barplot(GO_result, showCategory = 15)
cnetplot(GO_result)
