// R workflow step 2 (substep lin path)
// Nextflow module for clustering GO terms using Lin similarity
process cluster_r_lin {
    tag "R_Lin_Clustering"
    label "R"
    publishDir "output_dir_r", mode: 'copy'

    input:
        tuple val(baseName), path("GO_result.tsv"), path("GO_result_simplified_Lin.tsv"), path("GO_result_simplified_Wang.tsv"), path("Kegg_result.tsv") from enrichedR
        path(go_obo) optional true
        path(gaf_file) optional true

    output:
        path("${baseName}_lin_clusters.tsv"), emit: clusters
        tuple path("${baseName}_lin_clusters.tsv"), val("r_lin") emit: labeled

    script:
    """
    Rscript ${projectDir}/scripts/R_workflow/clustering/clustering_goterms_lin.R \
        GO_result.tsv ${baseName}
    """
}