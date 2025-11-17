#!/usr/bin/env python3

# Cluster GO terms by semantic similarity using GOATOOLS IC/resnik (Lin similarity).
# Inputs: enrichment TSV, GO OBO, and optional GOA GAF for background IC.
# Usage: python cluster_go_terms.py -i result/some/path/enrichment.tsv --gaf goa_human.gaf.gz


import argparse
import gzip
import json
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


# Try to import goatools; give clear error if not available.
try:
    from goatools.obo_parser import GODag
    # TermCounts is used to compute IC in typical goatools workflows
    from goatools.semantic import TermCounts, get_info_content, resnik_sim
except Exception as e:
    raise ImportError(
        "This script requires goatools. Install with `pip install goatools`.\n"
        f"Original import error: {e}"
    )

# -----------------------
# ARGPARSE
# -----------------------
def get_args():
    p = argparse.ArgumentParser(
        description=(
            "Cluster GO terms by semantic similarity using GOATOOLS IC/resnik (Lin). "
            "Inputs: enrichment TSV, GO OBO, and optional GOA GAF for background IC."
        )
    )
    p.add_argument("--input-file", "-i", required=True, help="Path to enrichment TSV input")
    p.add_argument(
        "--gaf",
        default="goa_human.gaf.gz",
        help="Path to GOA GAF file for background IC (default: goa_human.gaf.gz)",
    )
    p.add_argument("--taxon", default="9606", help="NCBI taxon id filter for GAF (default: 9606)")
    p.add_argument(
        "--exclude-evidence",
        default="IEA,ND",
        help="Comma-separated evidence codes to exclude from background IC (default: IEA,ND)",
    )
    p.add_argument("--min-term-size", type=int, default=20, help="Min term_size to keep (default: 20)")
    p.add_argument(
        "--min-intersection", type=int, default=5, help="Min intersection_size to keep (default: 5)"
    )
    p.add_argument("--pvalue", type=float, default=0.05, help="Max p_value to keep (default: 0.05)")
    p.add_argument(
        "--sim-threshold",
        type=float,
        default=0.3,
        help=(
            "Similarity threshold (Lin) for cutting dendrogram; clusters are formed where sim >= threshold (default: 0.3)"
        ),
    )
    p.add_argument("--species", default="human", help="Species for default GAF mapping (human|mouse).")
    p.add_argument("--out-prefix", default="", help="Optional prefix for output filenames.")
    return p.parse_args()

# -----------------------
# DATA LOADING HELPERS
# -----------------------

def read_enrichment_table(path):
    """Read enrichment table as TSV and enforce required columns."""
    df = pd.read_csv(path)
    required = {"native", "miRNAs", "Genes", "term_size", "intersection_size", "p_value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in input file: {sorted(missing)}")
    # Clean
    df = df.dropna(subset=["native", "miRNAs", "Genes"]).copy()
    df = df.drop_duplicates(subset=["miRNAs", "native"])  # reduce redundancy
    return df


def filter_enrichment(df, min_term_size, min_intersection, pvalue):
    return df[
        (df["term_size"] >= min_term_size)
        & (df["intersection_size"] >= min_intersection)
        & (df["p_value"] < pvalue)
    ].copy()


def build_gene2gos_from_dataset(df):
    """Build gene->set(GO) mapping from enrichment table's Genes/native columns."""
    mapping = defaultdict(set)
    for _, row in df.iterrows():
        go_id = str(row["native"]).strip()
        genes = [g.strip() for g in str(row["Genes"]).split(",") if g and g.strip().lower() != "nan"]
        for g in genes:
            mapping[g].add(go_id)
    return dict(mapping)


def load_gene2gos_from_gaf(path, obodag, taxon_filter="9606", exclude_evidence=("IEA", "ND")):
    """
    Minimal, robust GAF 2.2 reader for background IC:
    - Skips header lines (!)
    - Skips NOT qualifiers
    - Filters by taxon id substring
    - Excludes given evidence codes
    - Keeps only GO IDs present in the provided obodag

    Returns: dict[str, set[str]] mapping DB_Object_ID -> set(GO_ID)
    """
    mapping = defaultdict(set)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line[0] == "!":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            db_object_id = parts[1]
            qualifier = parts[3]
            go_id = parts[4]
            evidence = parts[6]
            taxon_col = parts[12]
            # Filters
            if "NOT" in qualifier:
                continue
            if taxon_filter and (taxon_filter not in taxon_col):
                continue
            if exclude_evidence and evidence in exclude_evidence:
                continue
            if go_id not in obodag:
                continue
            mapping[db_object_id].add(go_id)
    return dict(mapping)


# -----------------------
# SEMANTIC UTILS (GOATOOLS-BASED)
# -----------------------


def build_caches(obodag):
    ancestors_cache = {go_id: set(obodag[go_id].get_all_parents()) | {go_id} for go_id in obodag}
    depth_cache = {go_id: obodag[go_id].depth for go_id in obodag}
    return ancestors_cache, depth_cache


def mica_for_terms(terms, obodag, tc, ancestors_cache, pvalue_map):
    """Return MICA term id and its IC for the list of terms (same namespace), else (None, 0)."""
    valid = [t for t in terms if t in obodag]
    if not valid:
        return None, 0.0
    common = ancestors_cache[valid[0]].copy()
    for t in valid[1:]:
        common.intersection_update(ancestors_cache[t])
        if not common:
            return None, 0.0
    # Choose ancestor with lowest p-value, breaking ties with highest IC
    mica_id = None
    mica_ic = -1.0
    mica_pvalue = float('inf')  # Start with the highest possible p-value
    for a in common:
        ic_a = get_info_content(a, tc)
        pvalue_a = pvalue_map.get(a, float('inf'))  # Default to inf if p-value is not available
        if pvalue_a < mica_pvalue or (pvalue_a == mica_pvalue and ic_a > mica_ic):
            mica_pvalue = pvalue_a
            mica_ic = ic_a
            mica_id = a
    if mica_id is None:
        return None, 0.0
    return mica_id, mica_ic


def lin_similarity(t1, t2, obodag, tc):
    """Compute Lin similarity in [0,1] using GOATOOLS resnik + IC.

    Lin(t1,t2) = 2 * IC(MICA) / (IC(t1) + IC(t2))
    If denominator is 0, returns 0.
    """
    if t1 not in obodag or t2 not in obodag:
        return 0.0
    ic1 = get_info_content(t1, tc)
    ic2 = get_info_content(t2, tc)
    denom = ic1 + ic2
    if denom <= 0:
        return 0.0
    # resnik_sim returns IC(MICA)
    ic_mica = resnik_sim(t1, t2, obodag, tc)
    if ic_mica <= 0:
        return 0.0
    return max(0.0, min(1.0, (2.0 * ic_mica) / denom))


# -----------------------
# CLUSTERING PER NAMESPACE
# -----------------------

def cluster_terms_by_namespace(all_terms, obodag, tc, sim_threshold, pvalue_map):
    """Cluster and plot per GO namespace. Returns a DataFrame of clusters across namespaces."""
    # Organize terms by namespace
    ns_to_terms = defaultdict(list)
    for t in all_terms:
        if t in obodag:
            ns = obodag[t].namespace
            ns_to_terms[ns].append(t)

    all_clusters = []
    plot_data = {}

    for ns, terms in ns_to_terms.items():
        terms = sorted(set(terms))
        n = len(terms)
        if n < 2:
            # Singletons: still record cluster of size 1
            for t in terms:
                all_clusters.append(
                    {
                        "namespace": ns,
                        "cluster_id": 1,
                        "cluster_label": f"{obodag[t].name} ({t}, n=1)",
                        "rep_go_id": t,
                        "rep_go_name": obodag[t].name,
                        "cluster_size": 1,
                    }
                )
            continue

        # Precompute IC for diagonal and speed
        ic_map = {t: get_info_content(t, tc) for t in terms}

        # Build similarity (Lin) matrix
        sim = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    sim[i, j] = 1.0 if ic_map[terms[i]] > 0 else 0.0
                else:
                    s = lin_similarity(terms[i], terms[j], obodag, tc)
                    sim[i, j] = sim[j, i] = s

        # Convert to distance in [0,1]
        dist = 1.0 - sim
        np.fill_diagonal(dist, 0.0)

        # Hierarchical clustering (average linkage)
        Z = linkage(squareform(dist), method="average")
        dist_threshold = 1.0 - sim_threshold
        labels = fcluster(Z, t=dist_threshold, criterion="distance")

        # Map cluster id -> terms
        cid_to_terms = defaultdict(list)
        for idx, cid in enumerate(labels):
            cid_to_terms[cid].append(terms[idx])

        # Label and representative per cluster
        for cid, cterms in cid_to_terms.items():
            # Label using MICA of the cluster
            mica_id, _ = mica_for_terms(cterms, obodag, tc, ancestors_cache, pvalue_map)
            if mica_id is not None:
                label = f"{obodag[mica_id].name} ({mica_id}, n={len(cterms)})"
            else:
                label = f"Unlabeled (n={len(cterms)})"

            # Representative term (medoid)
            idxs = [terms.index(t) for t in cterms]
            subD = dist[np.ix_(idxs, idxs)]
            medoid_local = int(np.argmin(subD.sum(axis=0)))
            medoid_idx = idxs[medoid_local]
            rep = terms[medoid_idx]

            all_clusters.append(
                {
                    "namespace": ns,
                    "cluster_id": int(cid),
                    "cluster_label": label,
                    "rep_go_id": rep,
                    "rep_go_name": obodag[rep].name,
                    "cluster_size": len(cterms),
                }
            )

        plot_data = {
            ns: {
                "terms": terms,
                "sim_matrix": sim,
                "dist_matrix": dist,
                "linkage": Z,
                "dist_threshold": dist_threshold,
            }
        }

        # Aggregate clusters
        df_clusters = pd.DataFrame(all_clusters)
        df_all_clusters = pd.DataFrame(all_clusters)
        if not df_clusters.empty:
            df_clusters = df_clusters.sort_values(["namespace", "cluster_size"], ascending=[True, False]).reset_index(drop=True)
            df_clusters = df_clusters[df_clusters["cluster_size"] > 5].copy()
            print("Clusters processed.")

        return df_clusters, df_all_clusters, plot_data


# -----------------------
# MAIN
# -----------------------

def main():
    args = get_args()

    # Load enrichment and filter
    df = read_enrichment_table(args.input_file)
    df_filt = filter_enrichment(df, args.min_term_size, args.min_intersection, args.pvalue)

    # Load GO ontology
    obodag = GODag(args.obo)

    # Build dataset gene->GO mapping and collect unique GO terms present
    gene2gos = build_gene2gos_from_dataset(df_filt)
    all_terms = sorted({t for gos in gene2gos.values() for t in gos if t in obodag})

    # Background TermCounts using GAF if available
    exclude_evi = tuple(e.strip() for e in args.exclude_evidence.split(",") if e.strip())
    try:
        if os.path.exists(args.gaf):
            bg_gene2gos = load_gene2gos_from_gaf(args.gaf, obodag, taxon_filter=args.taxon, exclude_evidence=exclude_evi)
            if not bg_gene2gos:
                raise ValueError("No annotations loaded from GAF for background IC")
            tc_background = TermCounts(obodag, bg_gene2gos)
            print(f"Loaded {len(bg_gene2gos)} background genes from GAF for IC.")
        else:
            raise FileNotFoundError(f"{args.gaf} not found")
    except Exception as e:
        print(f"WARNING: Could not load GOA background ({e}). Falling back to dataset-derived counts.")
        tc_background = TermCounts(obodag, gene2gos)

    # Build caches shared across steps
    global ancestors_cache
    ancestors_cache, depth_cache = build_caches(obodag)

    if not all_terms:
        print("No valid GO terms after filtering; exiting.")
        return
    
    # Create pvalue_map from the filtered enrichment DataFrame
    pvalue_map = dict(zip(df_filt['native'], df_filt['p_value']))

    # Cluster per namespace using GOATOOLS IC/similarities
    df_clusters, df_all_clusters, plot_data  = cluster_terms_by_namespace(
        all_terms,
        obodag,
        tc_background,
        sim_threshold=args.sim_threshold,
        out_prefix=args.out_prefix,
        show_plots=not args.no_show,
        pvalue_map=pvalue_map
    )

    # Save clusters table
    out_tsv = f"{args.out_prefix}_clustered.tsv"
    df_clusters.to_csv(out_tsv, sep="\t", index=False)
    print(f"Clusters written: {out_tsv}")

    # Save plot data to JSON
    out_json = f"{args.out_prefix}_clustering_plot_data.json"
    serializable_plot_data = {}
    for ns, data in plot_data.items():
        serializable_plot_data[ns] = {
            "terms": data["terms"],
            "sim_matrix": data["sim_matrix"].tolist(),
            "dist_matrix": data["dist_matrix"].tolist(),
            "linkage": data["linkage"].tolist(),
            "dist_threshold": data["dist_threshold"],
        }
    with open(out_json, "w") as f:
        json.dump(serializable_plot_data, f, indent=2)
    print(f"Clustering plot data written: {out_json}")


    # Optionally show (if not disabled)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()