// R workflow step 3 (substeps lin and wang paths)
// Nextflow module for plotting GO terms with custom clusters

// Step 3a: Plotting custom Wang clusters
process plot_cluster_custom_wang {
    tag "R_Plotting_Custom_wang"
    label "R"
    publishDir "output_dir_r/plots", mode: 'copy'

    input:
        tuple path(wang_clusters), val("r_wang") from labeled

    output:
        path("*_wang_dotplot.png"), emit: plt_dotplot_wang
        path("*_wang_barplot.png"), emit: plt_barplot_wang
        path("*_wang_dendrogram.png"), emit: plt_dendrogram_wang
        path("*_wang_heatmap.png"), emit: plt_heatmap_wang

    script:
    """
    Rscript ${projectDir}/scripts/R_workflow/visualization/plot_clusterprofiler.R \
        ${wang_clusters}
    """
}

// Step 3b: Plotting custom Lin clusters
process plot_cluster_custom_lin {
    tag "R_Plotting_Custom_lin"
    label "R"
    publishDir "output_dir_r/plots", mode: 'copy'

    input:
        tuple path(lin_clusters), val("r_lin") from labeled

    output:
        path("*_lin_dotplot.png"), emit: plt_dotplot_lin
        path("*_lin_barplot.png"), emit: plt_barplot_lin
        path("*_lin_dendrogram.png"), emit: plt_dendrogram_lin
        path("*_lin_heatmap.png"), emit: plt_heatmap_lin

    script:
    """
    Rscript ${projectDir}/scripts/R_workflow/visualization/plot_clusterprofiler.R \
        ${lin_clusters}
    """
}
