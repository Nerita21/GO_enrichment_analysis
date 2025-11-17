import sys
import pandas as pd
from gprofiler import GProfiler
from pathlib import Path

# Run g:Profiler enrichment analysis and map miRNAs to enriched GO terms
# Input: TSV file with columns "GeneSymbol" and "miRNA"
# Output: TSV file with enriched GO terms, associated genes, miRNAs, and miRNA counts

def main():
    if len(sys.argv) < 3:
        print("Usage: python run_gprofiler_enriched.py <input_file> <output_dir> [organism]")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    organism = sys.argv[3] if len(sys.argv) > 3 else "hsapiens"

    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{input_file.stem}_gprofiler_enriched.tsv"

    gene_column = "GeneSymbol"
    mirna_column = "miRNA"

    # Load input data
    df = pd.read_csv(input_file, sep="\t")
    if gene_column not in df.columns or mirna_column not in df.columns:
        raise ValueError(f"Input file must contain columns '{gene_column}' and '{mirna_column}'.")

    gene_mirna_map = df[[gene_column, mirna_column]].dropna()
    gene_list = gene_mirna_map[gene_column].unique().tolist()

    print("gene_mirna_map preview:")
    print(gene_mirna_map.head(10))

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
        sys.exit(0)

    print("g:Profiler results preview:")
    print(results.head(10))


    # Build gene-to-miRNA mapping
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

    results[outcols].to_csv(output_file, index=False)
    print(f"GO enrichment with miRNAs saved to: {output_file}")


if __name__ == "__main__":
    main()

