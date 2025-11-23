#!/usr/bin/env nextflow

// Python workflow step 1
// Nextflow process to run g:Profiler enrichment analysis 
process run_gprofiler {
    tag "Python_gprofiler"
    label "python"
    container 'go-enrichment-python:latest'

    publishDir "output_dir_py", mode: 'copy'

    input:
        path input_tsv

    output:
        path("*_gprofiler_enriched.tsv"), emit: enrichedPython

    script:
    """
    python3 ${launchDir}/scripts/Python_workflow/go_enrichment/run_gprofiler_enriched.py \
        ${input_tsv} \
         ${input_tsv.baseName}_gprofiler_enriched.tsv
    """
}

// Input: tsv with miRNA targets (header is "GeneSymbol" in example)
// Output: tsv with enriched GO terms (flat list)
