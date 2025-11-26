#!/usr/bin/env python3

import sys
import pandas as pd
from gprofiler import GProfiler
from pathlib import Path

# Run g:Profiler enrichment analysis and map miRNAs to enriched GO terms
# Input: TSV file with columns "GeneSymbol" and "miRNA"
# Output: TSV file with enriched GO terms, associated genes, miRNAs, and miRNA counts

def main():
    if len(sys.argv) < 2:
        print("Usage: python run_gprofiler_enriched.py <input_file> <output_file> [organism]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    organism = sys.argv[3] if len(sys.argv) > 3 else "hsapiens"

    gene_column = "GeneSymbol"
    # mirna_column = "miRNA"

    # Load input data
    df = pd.read_csv(input_file, sep="\t")
    gene_map = df[gene_column].dropna()
    gene_list = gene_map.str.strip().unique().tolist()

    """ if gene_column not in df.columns or mirna_column not in df.columns:
        raise ValueError(f"Input file must contain columns '{gene_column}' and '{mirna_column}'.")

    gene_mirna_map = df[[gene_column, mirna_column]].dropna()
    gene_list = gene_mirna_map[gene_column].unique().tolist()

    print("gene_mirna_map preview:")
    print(gene_mirna_map.head(10))
 """
    # Run g:Profiler
    gp = GProfiler(return_dataframe=True, )
    results = gp.profile(
        organism=organism,
        query=gene_list,
        sources=["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP"],
        no_evidences=False
    )

    if results.empty:
        print("No enrichment results found.")
        # Write an empty results file with headers so the pipeline can proceed
        empty_cols = [
            "GO_term", "Genes", 
            "native", "description", "p_value", "term_size", "intersection_size"
        ]
        pd.DataFrame(columns=empty_cols).to_csv(output_file, sep="\t", index=False)
        return

    print("g:Profiler results preview:")
    print(results.head(10))


    """ # Build gene-to-miRNA mapping
    gene_to_mirnas = gene_mirna_map.groupby(gene_column)[mirna_column].apply(set).to_dict()

    def get_mirna_set(genes):
        all_mirnas = set()
        for g in genes:
            all_mirnas.update(gene_to_mirnas.get(g, []))
        return sorted(all_mirnas)

    def join_list(items):
        return ",".join(items) if items else ""

    # Use correct column name: "intersection"
    if "intersections" not in results.columns:
        print("Expected 'intersections' column not found in g:Profiler output.")
        print("Available columns:", results.columns.tolist())
        sys.exit(1)

    results["Genes"] = results["intersections"].apply(lambda genes: join_list(sorted(genes)))
    results["miRNAs"] = results["intersections"].apply(lambda genes: join_list(get_mirna_set(genes)))
    results["miRNA_count"] = results["intersections"].apply(lambda genes: len(get_mirna_set(genes)))
    results["GO_term"] = results.get("name", results.get("term_name", "UNKNOWN"))

    # Select columns if they exist
    outcols = ["GO_term", "Genes", "miRNAs", "miRNA_count", "native", "description", "p_value", "term_size", "intersection_size"]
    outcols = [col for col in outcols if col in results.columns]
 """
    # Select columns if they exist
    outcols = ["GO_term", "Genes", "native", "description", "p_value", "term_size", "intersection_size"]
    outcols = [col for col in outcols if col in results.columns]

    # save output
    results[outcols].to_csv(output_file, sep="\t", index=False)
    print(f"GO enrichment saved to: {output_file}")


if __name__ == "__main__":
    main()

