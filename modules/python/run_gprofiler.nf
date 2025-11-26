#!/usr/bin/env nextflow

// Python workflow step 1
// Nextflow process to run g:Profiler enrichment analysis 
process run_gprofiler {
    tag "pythonGoProfiler"
    label "python"

    publishDir params.output_dir_py, mode: 'copy'

    input:
        path(input_tsv_abs)

    output:
         tuple val(input_tsv_abs.baseName), path("gprofiler_enriched.tsv"), emit: pyenriched

    script:
    """
    python3 ${launchDir}/scripts/Python_workflow/go_enrichment/run_gprofiler_enriched.py \
        ${input_tsv_abs} \
        "gprofiler_enriched.tsv"

    """
}

// Input: tsv with miRNA targets (header is "GeneSymbol" in example)
// Output: tsv with enriched GO terms (flat list)

