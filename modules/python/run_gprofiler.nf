// Python workflow step 1
// Nextflow process to run g:Profiler enrichment analysis 
process run_gprofiler {
    tag "Python_gprofiler"
    label "python"

    publishDir "output_dir_py", mode: 'copy'

    input:
        path input_tsv

    output:
        path("output/*.tsv"), emit: enrichedPython

    script:
    """
    mkdir -p output
    python ${projectDir}/scripts/Python_workflow/go_enrichment/run_gprofiler_enriched.py \
        ${input_tsv} output
    """
}

// Input: tsv with miRNA targets (header is "GeneSymbol" in example)
// Output: tsv with enriched GO terms (flat list)
