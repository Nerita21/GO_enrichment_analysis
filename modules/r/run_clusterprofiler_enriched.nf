// R workflow step 1
// Nextflow module for enriching GO terms
process run_clusterprofiler {
    tag "clusterProfiler enrichment"
    label "R"
    publishDir "output_dir_r", mode: 'copy'

    input:
        path input_tsv

    output:
        tuple val(input_tsv.baseName), path("GO_result.tsv"), path("GO_result_simplified_Lin.tsv"), path("GO_result_simplified_Wang.tsv"), path("Kegg_result.tsv") emit: enrichedR

    script:
    """
    Rscript ${projectDir}/scripts/R_workflow/go_enrichment/run_clusterprofiler_enriched.R \
        ${input_tsv}
    """
}
