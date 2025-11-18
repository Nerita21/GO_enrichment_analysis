// R workflow step 3 (substeps lin and wang paths)
// Nextflow module for plotting GO terms with clusterprofiler default

// Step 3a: Plotting ClusterProfiler results (Wang)
process plot_clusterprofiler_wang {
    tag "R_Plotting_ClusterProfiler_wang"
    label "R"
    publishDir "output_dir_r/plots", mode: 'copy'

    input:
        tuple val(baseName), path("GO_result.tsv"), path("GO_result_simplified_Wang.tsv") from enrichedR

    output:
        path("${baseName}_wang_dotplot.png"), emit: plt_dotplot_wang
        path("${baseName}_wang_barplot.png"), emit: plt_barplot_wang
        path("${baseName}_wang_dendrogram.png"), emit: plt_dendrogram_wang
        path("${baseName}_wang_heatmap.png"), emit: plt_heatmap_wang

    script:
    """
    Rscript ${projectDir}/scripts/R_workflow/visualization/plot_clusterprofiler.R \
        GO_result_simplified_Wang.tsv ${baseName}
    """
}

// Step 3b: Plotting ClusterProfiler results (Lin)
process plot_clusterprofiler_lin {
    tag "R_Plotting_ClusterProfiler_lin"
    label "R"
    publishDir "output_dir_r/plots", mode: 'copy'

    input:
        tuple val(baseName), path("GO_result_simplified_Lin.tsv") from enrichedR

    output:
        path("${baseName}_lin_dotplot.png"), emit: plt_dotplot_lin
        path("${baseName}_lin_barplot.png"), emit: plt_barplot_lin
        path("${baseName}_lin_dendrogram.png"), emit: plt_dendrogram_lin
        path("${baseName}_lin_heatmap.png"), emit: plt_heatmap_lin

    script:
    """
    Rscript ${projectDir}/scripts/R_workflow/visualization/plot_clusterprofiler.R \
        GO_result_simplified_Lin.tsv ${baseName}
    """
}