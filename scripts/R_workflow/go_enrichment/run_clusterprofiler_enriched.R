# Function to analyse miRNA target genes for GO and KEGG enrichment

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(dplyr)
library(dotenv)

# Load environment variables from paths.env
load_dot_env(file = "paths.env")
output_dir <- Sys.getenv("R_RESULT")
input_file <- Sys.getenv("DATA_RAW")
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

# Save the files to results
enrich_files <- analyse_miRNA_targets(input_file)
for (name in names(enrich_files)) {
    output_path <- file.path(output_dir, paste0(name, "_enrichment_results.tsv"))
    write.table(as.data.frame(enrich_files[[name]]), file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
}
