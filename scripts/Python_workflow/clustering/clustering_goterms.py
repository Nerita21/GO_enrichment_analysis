

# Cluster GO terms by semantic similarity using GOATOOLS IC/resnik (Lin similarity).
# Inputs: enrichment TSV, GO OBO, and optional GOA GAF for background IC.
# Usage: python cluster_go_terms.py -i result/some/path/enrichment.tsv --gaf goa_human.gaf.gz


from pathlib import Path
import argparse
import os
import gzip
from collections import defaultdict
import warnings

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

# Try to import goatools; give clear error if not available.
try:
    from goatools.obo_parser import GODag
    # TermCounts is used to compute IC in typical goatools workflows
    from goatools.semantic import TermCounts
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
    p.add_argument("--no-show", action="store_true", help="Do not show plots interactively.")
    return p.parse_args()

# -----------------------
# PATH RESOLUTION
# -----------------------
def resolve_paths(input_file: str, species: str):
    """Given the input file, infer the OBO and output directory locations.

    Expected structure (relative to input file's parent):
        data/background/go-basic.obo
        result/Python_based/
    """
    GAF_MAP = {
        "human": "goa_human.gaf.gz"
    }

    # species → GAF file
    gaf_file = GAF_MAP.get(species.lower())
    if gaf_file is None:
        raise ValueError(f"Unsupported species: {species}")

    input_path = Path(input_file).resolve()
    base = input_path.parent.parent.parent   # parent folder of input file (as your original layout assumed)

    paths = {
        "input": input_path,
        "obo": base / "data" / "background" / "go-basic.obo",
        "output_dir": base / "result" / "Python_based",
        "default_gaf": base / "data" / "background" / gaf_file,
    }

    # create output directory if missing
    paths["output_dir"].mkdir(parents=True, exist_ok=True)

    return paths

def make_output_filename(input_file: Path, output_dir: Path, suffix="clustered", prefix=""):
    stem = input_file.stem
    ext = "".join(input_file.suffixes) or ".tsv"
    new_name = f"{prefix + '.' if prefix else ''}{stem}.{suffix}{ext}"
    return output_dir / new_name

# -----------------------
# DATA LOADING HELPERS
# -----------------------
def read_enrichment_table(path):
    """Read enrichment table as TSV and enforce required columns."""
    # Allow reading TSV or CSV detected by extension
    path = Path(path)
    if path.suffix.lower() in (".csv",):
        df = pd.read_csv(path)
    else:
        # assume TSV by default
        df = pd.read_csv(path, sep="\t")
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
        (df["term_size"].astype(float) >= min_term_size)
        & (df["intersection_size"].astype(float) >= min_intersection)
        & (df["p_value"].astype(float) < pvalue)
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
    opener = gzip.open if str(path).endswith(".gz") else open
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
            # Only keep GO IDs known in the obodag
            if go_id not in obodag:
                continue
            mapping[db_object_id].add(go_id)
    return dict(mapping)

# -----------------------
# SEMANTIC HELPERS (IC / RESNIK / LIN)
# -----------------------
def build_caches(obodag):
    # ancestors_cache: term -> set(all parents + self)
    ancestors_cache = {go_id: set(obodag[go_id].get_all_parents()) | {go_id} for go_id in obodag}
    # depth_cache maybe used later
    depth_cache = {go_id: getattr(obodag[go_id], "depth", 0) for go_id in obodag}
    return ancestors_cache, depth_cache

def get_info_content(go_id, tc: TermCounts):
    """Return info content using TermCounts.get_info_content if available."""
    try:
        # goatools TermCounts has get_info_content or get_ic depending on version
        if hasattr(tc, "get_info_content"):
            return float(tc.get_info_content(go_id))
        elif hasattr(tc, "get_ic"):
            return float(tc.get_ic(go_id))
        else:
            # As fallback, attempt attribute access
            return float(tc.info_content.get(go_id, 0.0))
    except Exception:
        # If something odd occurs, return 0.0 as safe fallback
        return 0.0

def resnik_sim(t1, t2, obodag, tc: TermCounts):
    """Return IC of MICA using goatools' helpers if available; otherwise compute via common ancestors."""
    # Prefer goatools semantic.resnik_sim if present
    try:
        # Some versions of goatools have semantic.resnik_sim
        from goatools.semantic import resnik_sim as goat_resnik
        return float(goat_resnik(t1, t2, obodag, tc))
    except Exception:
        # Fallback: compute MICA as ancestor with max IC
        # Build sets of ancestors (including self)
        anc1 = set(obodag[t1].get_all_parents()) | {t1}
        anc2 = set(obodag[t2].get_all_parents()) | {t2}
        common = anc1.intersection(anc2)
        if not common:
            return 0.0
        # choose ancestor with maximum IC
        ic_vals = [(get_info_content(a, tc), a) for a in common]
        ic_vals.sort(reverse=True)
        return float(ic_vals[0][0])

def lin_similarity(t1, t2, obodag, tc: TermCounts):
    """Compute Lin similarity in [0,1] using Resnik IC.

    Lin(t1,t2) = 2 * IC(MICA) / (IC(t1) + IC(t2))
    """
    if t1 not in obodag or t2 not in obodag:
        return 0.0
    ic1 = get_info_content(t1, tc)
    ic2 = get_info_content(t2, tc)
    denom = ic1 + ic2
    if denom <= 0:
        return 0.0
    ic_mica = resnik_sim(t1, t2, obodag, tc)
    if ic_mica <= 0:
        return 0.0
    sim = (2.0 * ic_mica) / denom
    # numerical stability
    return float(max(0.0, min(1.0, sim)))

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
    # Choose ancestor with lowest p-value (if present), tie-break with highest IC
    mica_id = None
    mica_ic = -1.0
    mica_pvalue = float('inf')
    for a in common:
        ic_a = get_info_content(a, tc)
        pvalue_a = pvalue_map.get(a, float('inf'))
        if (pvalue_a < mica_pvalue) or (pvalue_a == mica_pvalue and ic_a > mica_ic):
            mica_pvalue = pvalue_a
            mica_ic = ic_a
            mica_id = a
    if mica_id is None:
        return None, 0.0
    return mica_id, mica_ic

# -----------------------
# CLUSTERING PER NAMESPACE
# -----------------------
def cluster_terms_by_namespace(all_terms, obodag, tc, ancestors_cache, sim_threshold, pvalue_map, min_cluster_size=1):
    """Cluster and return clusters across namespaces.

    Returns two DataFrames:
      - df_clusters: clusters filtered by min_cluster_size (e.g. > 5)
      - df_all_clusters: all clusters (including small ones)
    """
    # Organize terms by namespace
    ns_to_terms = defaultdict(list)
    for t in all_terms:
        if t in obodag:
            ns = obodag[t].namespace
            ns_to_terms[ns].append(t)

    all_clusters = []
    cluster_counter = 0
    for ns, terms in ns_to_terms.items():
        terms = sorted(set(terms))
        n = len(terms)
        if n == 0:
            continue
        if n == 1:
            cluster_counter += 1
            t = terms[0]
            all_clusters.append(
                {
                    "namespace": ns,
                    "cluster_id": cluster_counter,
                    "cluster_label": f"{obodag[t].name} ({t}, n=1)",
                    "rep_go_id": t,
                    "rep_go_name": obodag[t].name,
                    "cluster_size": 1,
                }
            )
            continue

        # Precompute IC for speed
        ic_map = {t: get_info_content(t, tc) for t in terms}

        # Build similarity (Lin) matrix
        sim = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    sim[i, j] = 1.0 if ic_map.get(terms[i], 0.0) > 0 else 0.0
                else:
                    s = lin_similarity(terms[i], terms[j], obodag, tc)
                    sim[i, j] = sim[j, i] = s

        # Convert to distance in [0,1]
        dist = 1.0 - sim
        # Hierarchical clustering (average linkage)
        # squareform expects a condensed distance matrix (n*(n-1)/2)
        with np.errstate(invalid="ignore"):
            condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        dist_threshold = 1.0 - sim_threshold
        labels = fcluster(Z, t=dist_threshold, criterion="distance")

        # Map cluster id -> terms
        cid_to_terms = defaultdict(list)
        for idx, cid in enumerate(labels):
            cid_to_terms[cid].append(terms[idx])

        # Label and representative per cluster
        for cid, cterms in cid_to_terms.items():
            cluster_counter += 1
            # Label using MICA of the cluster
            mica_id, _ = mica_for_terms(cterms, obodag, tc, ancestors_cache, pvalue_map)
            if mica_id is not None:
                label = f"{obodag[mica_id].name} ({mica_id}, n={len(cterms)})"
            else:
                label = f"Unlabeled (n={len(cterms)})"

            # Representative term (medoid) using sum of distances
            idxs = [terms.index(t) for t in cterms]
            subD = dist[np.ix_(idxs, idxs)]
            # sum of distances per row -> choose minimal
            medoid_local = int(np.argmin(subD.sum(axis=0)))
            medoid_idx = idxs[medoid_local]
            rep = terms[medoid_idx]

            all_clusters.append(
                {
                    "namespace": ns,
                    "cluster_id": int(cluster_counter),
                    "cluster_label": label,
                    "rep_go_id": rep,
                    "rep_go_name": obodag[rep].name,
                    "cluster_size": len(cterms),
                }
            )

    df_all_clusters = pd.DataFrame(all_clusters)
    # make a filtered view for larger clusters ( > min_cluster_size )
    df_clusters = df_all_clusters[df_all_clusters["cluster_size"] >= min_cluster_size].copy()
    if not df_clusters.empty:
        df_clusters = df_clusters.sort_values(["namespace", "cluster_size"], ascending=[True, False]).reset_index(drop=True)
    return df_clusters, df_all_clusters

# -----------------------
# MAIN
# -----------------------
def main():
    args = get_args()
    paths = resolve_paths(args.input_file, args.species)

    # Decide which GAF to use: user-provided or inferred default
    gaf_path = args.gaf if os.path.exists(args.gaf) else (paths["default_gaf"] if os.path.exists(paths["default_gaf"]) else None)

    # Load enrichment and filter
    df = read_enrichment_table(args.input_file)
    df_filt = filter_enrichment(df, args.min_term_size, args.min_intersection, args.pvalue)

    if df_filt.empty:
        print("No terms remaining after filtering. Exiting.")
        return

    # Load GO ontology (GODag)
    obo_path = paths["obo"]
    if not obo_path.exists():
        raise FileNotFoundError(f"OBO file not found at {obo_path}")
    obodag = GODag(str(obo_path))

    # Build dataset gene->GO mapping and collect unique GO terms present
    gene2gos = build_gene2gos_from_dataset(df_filt)
    all_terms = sorted({t for gos in gene2gos.values() for t in gos if t in obodag})

    # Background TermCounts using GAF if available
    exclude_evi = tuple(e.strip() for e in args.exclude_evidence.split(",") if e.strip())
    try:
        tc_background = None
        if gaf_path:
            bg_gene2gos = load_gene2gos_from_gaf(str(gaf_path), obodag, taxon_filter=args.taxon, exclude_evidence=exclude_evi)
            if bg_gene2gos:
                tc_background = TermCounts(obodag, bg_gene2gos)
                print(f"Loaded {len(bg_gene2gos)} background genes from GAF for IC ({gaf_path}).")
        if tc_background is None:
            # fall back to dataset-derived counts
            tc_background = TermCounts(obodag, gene2gos)
            print("Falling back to dataset-derived TermCounts for IC (GAF not available or produced no annotations).")
    except Exception as e:
        warnings.warn(f"Could not build TermCounts using GAF: {e}. Falling back to dataset counts.")
        tc_background = TermCounts(obodag, gene2gos)

    # Build caches shared across steps
    ancestors_cache, depth_cache = build_caches(obodag)

    if not all_terms:
        print("No valid GO terms after filtering; exiting.")
        return

    # Create pvalue_map from the filtered enrichment DataFrame
    pvalue_map = dict(zip(df_filt['native'].astype(str), df_filt['p_value'].astype(float)))

    # Cluster per namespace using GOATOOLS IC/similarities
    df_clusters, df_all_clusters = cluster_terms_by_namespace(
        all_terms,
        obodag,
        tc_background,
        ancestors_cache,
        sim_threshold=args.sim_threshold,
        pvalue_map=pvalue_map,
        min_cluster_size=6  # original code filtered cluster_size > 5
    )

    # Save clusters table
    out_tsv = make_output_filename(paths["input"], paths["output_dir"], suffix="clusters", prefix=args.out_prefix)
    df_clusters.to_csv(out_tsv, sep="\t", index=False)
    print(f"Clusters written: {out_tsv}")

    # Optionally show (if not disabled)
    if not args.no_show and not df_all_clusters.empty:
        # Simple barplot of cluster sizes per namespace (one plot)
        try:
            ax = df_all_clusters.groupby("namespace")["cluster_size"].sum().plot(kind="bar")
            ax.set_title("Aggregate cluster sizes by GO namespace")
            ax.set_ylabel("Sum of cluster sizes")
            plt.tight_layout()
            plt.show()
        except Exception as e:
            warnings.warn(f"Could not render plot: {e}")

if __name__ == "__main__":
    main()
