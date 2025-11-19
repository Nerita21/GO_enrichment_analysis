#!/usr/bin/env python3
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
            "Cluster GO terms by semantic similarity using either Lin (IC-based) or Wang (graph-based). "
            "Inputs: enrichment TSV, GO OBO, and optional GOA GAF for background IC."
        )
    )
    p.add_argument("--input-file", "-i", required=True, help="Path to enrichment TSV input")
    p.add_argument("--obo", default="go-basic.obo", help="Path to GO OBO file (default: go-basic.obo)")
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
    p.add_argument("--min-intersection", type=int, default=5, help="Min intersection_size to keep (default: 5)")
    p.add_argument("--pvalue", type=float, default=0.05, help="Max p_value to keep (default: 0.05)")
    p.add_argument(
        "--sim-threshold",
        type=float,
        default=0.3,
        help="Similarity threshold for cutting dendrogram; clusters are formed where sim >= threshold (default: 0.3)",
    )
    p.add_argument("--sim-method", choices=["lin", "wang"], default="wang",
                   help="Similarity method to use: 'lin' (IC-based) or 'wang' (graph-based). Default: lin")
    p.add_argument("--species", default="human", help="Species for default GAF mapping (human|mouse).")
    p.add_argument("--out-prefix", default="goterm", help="Optional prefix for output filenames.")
    p.add_argument("--no-show", action="store_true", help="Do not call plt.show(); only save figures to disk")
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
# SEMANTIC SIMILARITY (GOATOOLS-BASED)
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
# WANG SEMANTIC SIMILARITY
# -----------------------

def get_wang_scores(go_id, obodag):
    """
    Compute Wang semantic contribution scores for a GO term.
    Returns a dict: ancestor_term -> score
    """
    if go_id not in obodag:
        return {}

    term_scores = {}
    queue = [(go_id, 1.0)]  # (current_term, score)

    while queue:
        term, score = queue.pop()
        # Keep maximum score if term already visited
        if term in term_scores:
            if score <= term_scores[term]:
                continue
        term_scores[term] = score

        # Each parent gets score * 0.8
        for parent in obodag[term].parents:
            queue.append((parent.id, score * 0.8))

    return term_scores


def wang_similarity(t1, t2, obodag):
    """
    Compute Wang semantic similarity between two GO terms.
    """
    if t1 not in obodag or t2 not in obodag:
        return 0.0

    S1 = get_wang_scores(t1, obodag)
    S2 = get_wang_scores(t2, obodag)

    # Intersection and union of ancestor sets
    common = set(S1).intersection(S2)
    if not common:
        return 0.0

    numer = sum(min(S1[a], S2[a]) for a in common)
    denom = sum(S1.values()) + sum(S2.values())

    if denom == 0:
        return 0.0

    return 2 * numer / denom


# -----------------------
# CLUSTERING PER NAMESPACE
# -----------------------

def cluster_terms_by_namespace(
    all_terms,
    obodag,
    tc,
    sim_threshold,
    pvalue_map,
    out_prefix="goterm",
    sim_method="lin",
):
    """
    Cluster and plot per GO namespace. Returns (df_clusters, df_all_clusters, plot_data).
    sim_method: "lin" or "wang"
    """
    if sim_method not in ("lin", "wang"):
        raise ValueError("sim_method must be 'lin' or 'wang'")

    # Organize terms by namespace
    ns_to_terms = defaultdict(list)
    for t in all_terms:
        if t in obodag:
            ns = obodag[t].namespace
            ns_to_terms[ns].append(t)

    all_clusters = []
    plot_data = {}

    # Precompute things for Wang if needed (cache scores per term)
    wang_score_cache = {}
    if sim_method == "wang":
        for t in all_terms:
            wang_score_cache[t] = get_wang_scores(t, obodag)

    # Iterate namespaces
    for ns, terms in ns_to_terms.items():
        terms = sorted(set(terms))
        n = len(terms)

        if n < 2:
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
            # still record an empty plot entry
            plot_data[ns] = {"terms": terms, "sim_matrix": np.zeros((len(terms), len(terms))).tolist(), "dist_matrix": [], "linkage": [], "dist_threshold": 1.0}
            continue

        # For Lin only: precompute IC map for speed (and to set diagonal conditionally)
        ic_map = None
        if sim_method == "lin":
            ic_map = {t: get_info_content(t, tc) for t in terms}

        # Build similarity matrix
        sim = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    sim[i, j] = sim[j, i] = 1.0
                    continue
                t1, t2 = terms[i], terms[j]
                if sim_method == "lin":
                    s = lin_similarity(t1, t2, obodag, tc)
                else:
                    # use cached wang scores if available
                    S1 = wang_score_cache.get(t1) or get_wang_scores(t1, obodag)
                    S2 = wang_score_cache.get(t2) or get_wang_scores(t2, obodag)
                    # compute Wang similarity using the same logic as wang_similarity(...)
                    common = set(S1).intersection(S2)
                    if not common:
                        s = 0.0
                    else:
                        numer = sum(min(S1[a], S2[a]) for a in common)
                        denom = sum(S1.values()) + sum(S2.values())
                        s = 0.0 if denom == 0 else 2.0 * numer / denom
                sim[i, j] = sim[j, i] = max(0.0, min(1.0, float(s)))

        # Distance matrix
        dist = 1.0 - sim
        np.fill_diagonal(dist, 0.0)

        # Hierarchical clustering
        Z = linkage(squareform(dist), method="average")
        dist_threshold = 1.0 - sim_threshold
        labels = fcluster(Z, t=dist_threshold, criterion="distance")

        # Map cluster id -> terms
        cid_to_terms = defaultdict(list)
        for idx, cid in enumerate(labels):
            cid_to_terms[cid].append(terms[idx])

        # Label and representative per cluster
        for cid, cterms in cid_to_terms.items():
            # Representative term (medoid)
            idxs = [terms.index(t) for t in cterms]
            subD = dist[np.ix_(idxs, idxs)]
            medoid_local = int(np.argmin(subD.sum(axis=0)))
            rep = cterms[medoid_local]

            # Label
            if sim_method == "lin":
                mica_id, _ = mica_for_terms(cterms, obodag, tc, ancestors_cache, pvalue_map)
                if mica_id is not None:
                    label = f"{obodag[mica_id].name} ({mica_id}, n={len(cterms)})"
                else:
                    label = f"Unlabeled (n={len(cterms)})"
            else:
                # choose the term with highest total wang ancestor score within the cluster
                score_map = {}
                for t in cterms:
                    s = wang_score_cache.get(t)
                    if s is None:
                        s = get_wang_scores(t, obodag)
                    score_map[t] = sum(s.values())
                best = max(score_map, key=score_map.get)
                label = f"{obodag[best].name} ({best}, n={len(cterms)})"

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

        # Save plotting data for this namespace
        plot_data[ns] = {
            "terms": terms,
            "sim_matrix": sim,
            "dist_matrix": dist,
            "linkage": Z,
            "dist_threshold": dist_threshold,
        }

    # after loop: aggregate and return
    df_all_clusters = pd.DataFrame(all_clusters)
    df_clusters = df_all_clusters.sort_values(["namespace", "cluster_size"], ascending=[True, False]).reset_index(drop=True)
    # optional filter for cluster size > 5 (preserve older behavior)
    df_clusters = df_clusters[df_clusters["cluster_size"] > 5].copy()

    print("Clusters processed for all namespaces.")
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
    df_clusters, df_all_clusters, plot_data = cluster_terms_by_namespace(
        all_terms,
        obodag,
        tc_background,
        sim_threshold=args.sim_threshold,
        pvalue_map=pvalue_map,
        out_prefix=args.out_prefix,
        sim_method=args.sim_method,
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