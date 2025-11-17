// Python workflow step 3
// Nextflow module for plotting GO terms
process plot_clusters {
    tag "Python_Plotting"
    label "python"
    publishDir "${projectDir}/results/Python_based/plots/${enriched_file.baseName}", mode: 'copy'
    
    input:
        path plot_json  // receives whatever matches *_lin_plot_data.json (by emit: plot)
    
    output:
        path("*_lin_dendrogram.png"), emit: plt_dendrogram
        path("*_lin_heatmap.png"), emit: plt_heatmap
    
    script:
    """
    python ${projectDir}/scripts/Python_workflow/visualization/plotting.py \
        --input-file ${enrichment_results} 
    """
}