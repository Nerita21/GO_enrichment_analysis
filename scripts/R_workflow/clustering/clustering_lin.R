#!/usr/bin/env Rscript

##########################################################
# R-based GO Term Clustering (Lin Similarity Method)
# 
# This script performs GO term clustering using Lin's method
# which is based on information content of terms.
# Lin similarity emphasizes specificity and information gain.
#
# Input: Enrichment results from clusterProfiler
# Output: Clustered GO terms with Lin-based semantic similarity
##########################################################

# Load libraries
suppressMessages({
    library("clusterProfiler")
    library("AnnotationDbi")
    library("GO.db")
    library("org.Hs.eg.db")
    library("dplyr")
    library("ggplot2")
    library("enrichplot")
    library("GOSemSim")
    library("dotenv")
})

# Load environment variables
load_dot_env(file = "paths.env")
output_dir <- Sys.getenv("R_RESULT")
if (output_dir == "") {
    output_dir <- "./results/R_based"
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("[Lin] Output directory: ", output_dir)

# --- 1. Calculate semantic similarity (Lin method) ---
calculate_sim_clustering_lin <- function(GO_result) {
    message("[Lin] Computing pairwise semantic similarity (Lin method)...")
    
    # Lin similarity uses information content (IC) to measure term specificity
    # It's based on the ratio of IC of the MICA to the average IC of the two terms
    GO_result_sim <- pairwise_termsim(GO_result, method = "Lin")
    
    message("[Lin] Calculated pairwise semantic similarity between GO terms...")
    GO_result_matrix <- GO_result_sim@termsim
    GO_dist_matrix <- 1 - GO_result_matrix
    
    # Hierarchical clustering (average linkage)
    GO_hclustered <- hclust(as.dist(GO_dist_matrix), method = "average")
    
    message("[Lin] Performed hierarchical clustering of GO terms")
    return(list(
        sim_result = GO_result_sim,
        hclust = GO_hclustered,
        dist_matrix = GO_dist_matrix
    ))
}

# --- 2. Create simplified cluster representation ---
simplify_clusters_lin <- function(GO_result) {
    message("[Lin] Simplifying GO results...")
    
    # Use built-in simplify with Lin method
    GO_result_simplified <- simplify(GO_result, cutoff = 0.7, by = "p.adjust", select_fun = min)
    
    # Convert to data frame for output
    result_df <- as.data.frame(GO_result_simplified@result)
    result_df <- result_df %>%
        select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
        mutate(method = "Lin", cluster_type = "simplified")
    
    message("[Lin] Created simplified cluster representation with ", nrow(result_df), " terms")
    return(result_df)
}

# --- 3. Main execution ---
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) < 1) {
        stop("Usage: Rscript clustering_lin.R <enrichment_results.rds> [output_prefix]")
    }
    
    enrichment_file <- args[1]
    output_prefix <- if (length(args) > 1) args[2] else "go_enrichment"
    
    # Load enrichment results (saved from run_clusterprofiler_enriched.R)
    if (!file.exists(enrichment_file)) {
        stop("Enrichment file not found: ", enrichment_file)
    }
    
    message("[Lin] Loading enrichment results from: ", enrichment_file)
    GO_result <- readRDS(enrichment_file)
    
    # Verify object is enrichResult
    if (!("enrichResult" %in% class(GO_result))) {
        stop("Input must be an enrichResult object")
    }
    
    # Run clustering with Lin method
    clustering_results <- calculate_sim_clustering_lin(GO_result)
    
    # Generate simplified representation
    simplified_df <- simplify_clusters_lin(GO_result)
    
    # Save results
    output_file <- file.path(output_dir, paste0(output_prefix, "_lin_clusters.tsv"))
    write.table(
        simplified_df,
        file = output_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    message("[Lin] Saved Lin-method clusters to: ", output_file)
    
    # Create visualization
    plot_file <- file.path(output_dir, paste0(output_prefix, "_lin_plot.png"))
    png(plot_file, width = 1200, height = 800, res = 150)
    
    tryCatch({
        # Barplot of top terms by p-adjust
        p <- barplot(GO_result, showCategory = 20, title = "GO Enrichment (Lin Method)")
        print(p)
    }, error = function(e) {
        message("[Lin] Warning: Could not create barplot - ", e$message)
    })
    
    dev.off()
    message("[Lin] Saved plot to: ", plot_file)
    
    # Return success
    message("[Lin] Clustering complete!")
    invisible(list(clusters = simplified_df, clustering = clustering_results))
}

# Execute
tryCatch({
    main()
}, error = function(e) {
    message("[Lin] ERROR: ", e$message)
    quit(status = 1)
})
