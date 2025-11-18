// Python workflow step 2
// Nextflow module for clustering GO terms
process cluster_python_lin {
    tag "Python_Lin_Clustering"
    label "python"
    publishDir "${projectDir}/results/Python_based", mode: 'copy'
    
    input:
        path(enrichment_results)
    
    output:
        tuple path("*_lin_clusters.tsv"), path("*_lin_plot_data.json"), val(enrichment_results.baseName), emit: clustered
    
    script:
    """
    python ${projectDir}/scripts/Python_workflow/clustering/clustering_goterms.py \
        --input-file ${enrichment_results} \
        --out-prefix lin 
    """
}

// Input: tsv with enriched GO terms (from g:Profiler) 

// Optional parameters  (othewise default values used):
//   --sim-threshold : Similarity threshold (Lin) for cutting dendrogram; clusters are formed where sim >= threshold (default: 0.3)
//   --out-prefix   : additional prefix for output files (default "")
//   --gaf         : GAF file for IC calculation (default: human GAF from GO)
//   --taxon-id    : NCBI taxon ID for IC calculation (default 9606 - human)
//   --exclude-evidence : evidence codes to exclude from IC calculation (default: IEA)
//   --min-term-size    : minimum term size for clustering (default 20)
//   --min-intersection : minimum intersection size for clustering (default 1)
//   --pvalue   : Max p_value to keep (default: 0.05)

// Output: tsv with clustered GO terms (Lin similarity)