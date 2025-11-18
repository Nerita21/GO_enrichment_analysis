# Function to analyse miRNA target genes for GO and KEGG enrichment

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(dplyr)


# Get Arguments

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]

out_GO <- args[2]        # e.g., "GO_result.tsv"
out_GO_lin <- args[3]    # e.g., "GO_result_simplified_Lin.tsv"
out_GO_wang <- args[4]   # e.g., "GO_result_simplified_Wang.tsv"
out_KEGG <- args[5]      # e.g., "Kegg_result.tsv"

semData <- godata('org.Hs.eg.db', ont="BP")

analyse_miRNA_targets <- function(input_file) {

    ############################################################################
    # LOADING DATA

    # Load datasets
    mirna_target_gene_set <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)

    #############################################################
    # GENELIST

    # obtaining the gene list for enrichment
    gene_list <- unique(given_mirna_target_gene_set$GeneSymbol)
    gene_list <- gene_list[!is.na(gene_list)]

    # assigning entrez id
    entrez_df <- bitr(given_mirna_target_gene_set$GeneSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    # adding entrez to weights
    entrez_list <- unique(entrez_df$ENTREZID)
    entrez_list <- entrez_list[!is.na(entrez_list)]

    ####################################################################
    # RUN ORA

    set.seed(1234)

    # run GO enrichment
    GO_result <- clusterProfiler::enrichGO(gene= gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.05)
    

    # run simplify to remove redundant terms
    GO_result_simplyfied_Lin <- simplify(GO_result, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Lin", semData = semData)
    GO_result_simplyfied_Wang <- simplify(GO_result, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = semData)

    # run kegg on shortlisted
    Kegg_result <- clusterProfiler::enrichKEGG(gene= entrez_list, organism = 'hsa', minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod= "BH")
    

    return(list(GO_result = GO_result,
                GO_result_simplyfied_Lin = GO_result_simplyfied_Lin,
                GO_result_simplyfied_Wang = GO_result_simplyfied_Wang,
                Kegg_result = Kegg_result))
}

# Run the analysis
enrich_files <- analyse_miRNA_targets(input_file)

# Save the files to results
write.table(as.data.frame(enrich_files$GO_result), file = out_GO, sep = "\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(enrich_files$GO_result_simplyfied_Lin), file = out_GO_lin, sep = "\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(enrich_files$GO_result_simplyfied_Wang), file = out_GO_wang, sep = "\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(enrich_files$Kegg_result), file = out_KEGG, sep = "\t", row.names=FALSE, quote=FALSE)