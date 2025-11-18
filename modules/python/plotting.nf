// Python workflow step 3
// Nextflow module for plotting GO terms
process plot_clusters {
    tag "Python_Plotting"
    label "python"
    publishDir "${projectDir}/results/Python_based/plots/${input_name}", mode: 'copy'
    
    input:
        tuple path(_lin_clusters), path(_lin_plot_json), val(input_name) from clustered
        path(enrichment_results) from enriched
    output:
        path("*_lin_barplot.png"), emit: plt_barblot
        path("*_lin_dotplot.png"), emit: plt_dotplot
        path("*_lin_heatmap.png"), emit: plt_heatmap
        path("*_lin_dendrogram.png"), emit: plt_dendrogram
    
    script:
    """
    python ${projectDir}/scripts/Python_workflow/visualization/plotting.py \
        --clusters ${_lin_clusters} \
        --enrichment ${enrichment_results} \
        --plot-data-json ${_lin_plot_json} 
    """
}