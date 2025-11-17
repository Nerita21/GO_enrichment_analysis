#!/usr/bin/env Rscript

##########################################################
# R-based Custom GO Term Clustering
# 
# This script implements custom hierarchical clustering
# using semantic similarity combined with p-value weighting.
# It provides a flexible alternative to built-in simplify().
#
# Input: Enrichment results from clusterProfiler
# Output: Clustered GO terms with custom clustering approach
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

message("[Custom] Output directory: ", output_dir)

# --- 1. Calculate semantic similarity and perform custom clustering ---
custom_cluster_go_terms <- function(GO_result, sim_threshold = 0.7) {
    message("[Custom] Computing semantic similarity for custom clustering...")
    
    # Compute pairwise similarity (default: Rel - a hybrid method)
    GO_result_sim <- pairwise_termsim(GO_result, method = "Rel")
    
    # Get similarity matrix
    sim_matrix <- GO_result_sim@termsim
    dist_matrix <- 1 - sim_matrix
    
    message("[Custom] Performing hierarchical clustering...")
    # Hierarchical clustering
    hc <- hclust(as.dist(dist_matrix), method = "average")
    
    # Cut tree at threshold
    cluster_dist <- 1 - sim_threshold
    clusters <- cutree(hc, h = cluster_dist)
    
    message("[Custom] Identified ", length(unique(clusters)), " clusters")
    
    return(list(
        hclust = hc,
        clusters = clusters,
        sim_matrix = sim_matrix,
        dist_matrix = dist_matrix
    ))
}

# --- 2. Extract representative and summary per cluster ---
extract_cluster_representatives <- function(GO_result, cluster_assignments) {
    message("[Custom] Extracting cluster representatives...")
    
    result_df <- as.data.frame(GO_result@result)
    go_ids <- result_df$ID
    
    cluster_summary <- data.frame()
    
    for (cid in unique(cluster_assignments)) {
        # Get members of this cluster
        members <- go_ids[cluster_assignments == cid]
        member_rows <- result_df[result_df$ID %in% members, ]
        
        # Representative: highest abundance (lowest p.adjust)
        rep_idx <- which.min(member_rows$p.adjust)
        rep_id <- member_rows$ID[rep_idx]
        rep_desc <- member_rows$Description[rep_idx]
        
        # Summary statistics
        cluster_row <- data.frame(
            cluster_id = cid,
            representative_id = rep_id,
            representative_desc = rep_desc,
            cluster_size = nrow(member_rows),
            avg_p_adjust = mean(member_rows$p.adjust),
            min_p_adjust = min(member_rows$p.adjust),
            max_p_adjust = max(member_rows$p.adjust),
            member_ids = paste(members, collapse = ";"),
            method = "custom",
            cluster_type = "hierarchical"
        )
        
        cluster_summary <- rbind(cluster_summary, cluster_row)
    }
    
    # Sort by p-adjust
    cluster_summary <- cluster_summary %>%
        arrange(min_p_adjust)
    
    message("[Custom] Extracted representatives for ", nrow(cluster_summary), " clusters")
    return(cluster_summary)
}

# --- 3. Main execution ---
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) < 1) {
        stop("Usage: Rscript clustering_custom.R <enrichment_results.rds> [output_prefix]")
    }
    
    enrichment_file <- args[1]
    output_prefix <- if (length(args) > 1) args[2] else "go_enrichment"
    
    # Load enrichment results
    if (!file.exists(enrichment_file)) {
        stop("Enrichment file not found: ", enrichment_file)
    }
    
    message("[Custom] Loading enrichment results from: ", enrichment_file)
    GO_result <- readRDS(enrichment_file)
    
    # Verify object is enrichResult
    if (!("enrichResult" %in% class(GO_result))) {
        stop("Input must be an enrichResult object")
    }
    
    # Perform custom clustering
    clustering_results <- custom_cluster_go_terms(GO_result, sim_threshold = 0.7)
    
    # Extract representatives
    cluster_summary <- extract_cluster_representatives(GO_result, clustering_results$clusters)
    
    # Save results
    output_file <- file.path(output_dir, paste0(output_prefix, "_custom_clusters.tsv"))
    write.table(
        cluster_summary,
        file = output_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    message("[Custom] Saved custom clusters to: ", output_file)
    
    # Create dendrogram visualization
    plot_file <- file.path(output_dir, paste0(output_prefix, "_custom_dendrogram.png"))
    png(plot_file, width = 1200, height = 800, res = 150)
    
    tryCatch({
        # Plot dendrogram with cluster cutoff line
        plot(clustering_results$hclust, 
             main = "GO Term Dendrogram (Custom Clustering)",
             xlab = "GO Term ID",
             ylab = "Distance (1 - Similarity)")
        abline(h = 1 - 0.7, col = "red", lty = 2, lwd = 2)
        legend("topright", legend = "Cluster Threshold", col = "red", lty = 2)
    }, error = function(e) {
        message("[Custom] Warning: Could not create dendrogram - ", e$message)
    })
    
    dev.off()
    message("[Custom] Saved dendrogram to: ", plot_file)
    
    # Return success
    message("[Custom] Custom clustering complete!")
    invisible(list(clusters = cluster_summary, clustering = clustering_results))
}

# Execute
tryCatch({
    main()
}, error = function(e) {
    message("[Custom] ERROR: ", e$message)
    quit(status = 1)
})
