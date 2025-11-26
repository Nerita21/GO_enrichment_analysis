#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ========== PARAMETERS ==========

params.input_tsv = "${projectDir}/data/raw/target_genes_hsa_mir_107.tsv"
params.go_obo = "${projectDir}/data/background/go-basic.obo"
params.gaf_file = null  // optional: GAF annotation file for IC calculation
params.output_dir_py = "${projectDir}/results/Python_based"
params.output_dir_r = "${projectDir}/results/R_based"
params.comparison_dir = "${projectDir}/results/comparison"
params.input_enrich_r = "${launchDir}/modules/r/run_clusterprofiler_enriched.nf"

// Global parameters for harmonized execution across all methods
params.p_cutoff = 0.05
params.min_genes_per_term = 1
params.similarity_threshold = 0.3
params.image_python = 'go-enrichment-python:latest'
params.image_r = 'go-enrichment-r:latest'

// Channel definition
ch_input_genes = file(params.input_tsv)
ch_go_obo = file(params.go_obo)
ch_gaf_file = params.gaf_file ? file(params.gaf_file) : null

log.info "Paths: ${params.input_tsv} and ${params.output_dir_r} and ${params.input_enrich_r} and Project: ${projectDir} and launch: ${launchDir}"

include { run_clusterprofiler } from "${params.input_enrich_r}"

params.input_tsv_abs = "${launchDir}/data/raw/target_genes_hsa_mir_107.tsv"

log.info "Using input TSV: ${params.input_tsv_abs}"

// ========== WORKFLOW DEFINITION ==========
workflow {
    main:
        log.info """
        ================================================
        GO Enrichment Analysis Pipeline (v2.0)
        ================================================
        Input genes       : ${params.input_tsv}
        GO OBO file       : ${params.go_obo}
        P-value cutoff    : ${params.p_cutoff}
        Sim. threshold    : ${params.similarity_threshold}
        Python image      : ${params.image_python}
        R image           : ${params.image_r}
        ================================================
        """

        // Step 1: Run clusterProfiler enrichment
        r_results = run_clusterprofiler(ch_input_genes)
}
