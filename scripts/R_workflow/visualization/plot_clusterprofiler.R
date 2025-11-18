# Load required libraries
library(clusterProfiler)
library(ggplot2)
library(enrichplot)

# Get Arguments
args <- commandArgs(trailingOnly = TRUE)
GO_result <- args[1]  # This should be an object of class enrichResult from clusterProfiler
out_dotplot <- args[2]  # e.g., "GO_dotplot.png"
out_barplot <- args[3]  # e.g., "GO_barplot.png"
out_cnetplot <- args[4]  # e.g., "GO_cnetplot.png"
out_heatplot <- args[5]  # e.g., "GO_heatplot.png

# Generate visualizations
dotplot_fig <- dotplot(GO_result, showCategory=15) + ggtitle("Dotplot of Enriched GO Terms")
barplot_fig <- barplot(GO_result, showCategory=15) + ggtitle("Barplot of Enriched GO Terms")
cnetplot_fig <- cnetplot(GO_result, showCategory=15, foldChange=NULL) + ggtitle("Cnetplot of Enriched GO Terms")
heatplot_fig <- heatplot(GO_result, showCategory=15, foldChange=NULL) + ggtitle("Heatplot of Enriched GO Terms")

# Save plots
ggsave(out_dotplot, plot = dotplot_fig, width = 8, height = 6)
ggsave(out_barplot, plot = barplot_fig, width = 8, height = 6)
ggsave(out_cnetplot, plot = cnetplot_fig, width = 8, height = 6)
ggsave(out_heatplot, plot = heatplot_fig, width = 8, height = 6)