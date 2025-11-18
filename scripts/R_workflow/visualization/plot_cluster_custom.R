# load libraries

# get arguments
args <- commandArgs(trailingOnly=TRUE)
clustering_table_file <- args[1]  # e.g., "GO_wang_clustering.tsv"
out_barplot <- args[2] # e.g., "GO_clusters_barplot.png"
out_dotplot <- args[3] # e.g., "GO_clusters_dotplot.png"


barplot_for_clusters <- function(go_hcluster_table) {
    # Compute cluster size for each node
    go_hcluster_table$cluster_size <- sapply(go_hcluster_table$id, function(cid) length(get_leaves(cid, go_hcluster_table)))

    # Select only internal nodes (clusters) with labels
    clusters <- go_hcluster_table[!is.na(go_hcluster_table$label) & go_hcluster_table$name == "", ]

    # Filter significant clusters (optional, e.g., padjust < 0.05)
    clusters <- clusters[clusters$padjust < 0.05, ]

    # Order by significance (lowest padjust first)
    clusters <- clusters[order(clusters$padjust), ]

    cluster_top15 <- clusters %>%
    slice_head(n = 15)

    library(ggplot2)
    library(ggsci)
    cluserplot <- ggplot(hcluster_top10, aes(x = reorder(label,-log10(padjust)), y = -log10(padjust), fill = label)) +
    geom_col(show.legend = FALSE) + 
    labs(x = "GO terms", y = "Padjust", title = "Enriched GO terms by clusters") +
    coord_flip() +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(cluserplot)
}

dotplot_for_clusters <- function(go_hcluster_table) {
    # Compute cluster size for each node
    go_hcluster_table$cluster_size <- sapply(go_hcluster_table$id, function(cid) length(get_leaves(cid, go_hcluster_table)))

    # Select only internal nodes (clusters) with labels
    clusters <- go_hcluster_table[!is.na(go_hcluster_table$label) & go_hcluster_table$name == "", ]

    # Filter significant clusters (optional, e.g., padjust < 0.05)
    clusters <- clusters[clusters$padjust < 0.05, ]

    # Order by significance (lowest padjust first)
    clusters <- clusters[order(clusters$padjust), ]

    cluster_top15 <- clusters %>%
    slice_head(n = 15)

    library(ggplot2)
    library(ggsci)
    clusterdotplot <- ggplot(cluster_top15, aes(x = reorder(label,-log10(padjust)), y = -log10(padjust), size = cluster_size, color = padjust)) +
    geom_point() + 
    scale_color_viridis_c() +
    labs(x = "GO terms", y = "Padjust", title = "Enriched GO terms by clusters") +
    coord_flip() +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(clusterdotplot)
}

 # Save the plots
go_hcluster_table <- read.delim(clustering_table_file, header = TRUE, stringsAsFactors = FALSE)
barplot_fig <- barplot_for_clusters(go_hcluster_table)
dotplot_fig <- dotplot_for_clusters(go_hcluster_table)
ggsave(out_barplot, plot = barplot_fig, width = 8, height = 6)
ggsave(out_dotplot, plot = dotplot_fig, width = 8, height = 6)