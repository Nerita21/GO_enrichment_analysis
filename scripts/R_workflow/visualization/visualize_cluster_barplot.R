# load libraries

# load path
load_dot_env(file = "paths.env")
output_dir <- Sys.getenv("R_RESULT")

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
    # Save the plot
    output_path <- file.path(output_dir, "GO_clusters_barplot.png")
    ggsave(output_path, plot = cluserplot, width = 8, height = 6)
}

