#!/usr/bin/env python3
"""
Custom GO Term Clustering using GOATOOLS

This module provides custom hierarchical clustering with semantic similarity
combined with p-value weighting, as an alternative to built-in GOATOOLS simplify.

Usage:
    python clustering_custom.py --input-file enrichment.tsv [options]
"""

import argparse
import gzip
from collections import defaultdict
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

try:
    from goatools.obo_parser import GODag
    from goatools.semantic import TermCounts, get_info_content, resnik_sim
except ImportError:
    raise ImportError(
        "This script requires goatools. Install with: pip install goatools"
    )


# --- ARGUMENT PARSING ---
def get_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Custom GO term clustering with semantic similarity and p-value weighting"
    )
    parser.add_argument("--input-file", "-i", required=True,
                        help="Path to enrichment TSV file")
    parser.add_argument("--go-obo", default="data/background/go-basic.obo",
                        help="Path to GO OBO file (default: data/background/go-basic.obo)")
    parser.add_argument("--gaf", default=None,
                        help="Path to GOA GAF file for background IC computation")
    parser.add_argument("--sim-threshold", type=float, default=0.6,
                        help="Similarity threshold for clustering (default: 0.6)")
    parser.add_argument("--min-cluster-size", type=int, default=2,
                        help="Minimum cluster size to report (default: 2)")
    parser.add_argument("--pvalue-weight", type=float, default=0.3,
                        help="Weight for p-value in representative selection (default: 0.3)")
    parser.add_argument("--out-prefix", default="custom",
                        help="Output filename prefix (default: custom)")
    parser.add_argument("--no-show", action="store_true",
                        help="Do not display plots interactively")
    
    return parser.parse_args()


# --- DATA LOADING ---
def load_go_obo(obo_path):
    """Load GO DAG from OBO file."""
    print(f"Loading GO OBO from {obo_path}...")
    if not Path(obo_path).exists():
        raise FileNotFoundError(f"OBO file not found: {obo_path}")
    return GODag(str(obo_path))


def load_enrichment_results(tsv_path):
    """Load enrichment results from TSV file."""
    print(f"Loading enrichment results from {tsv_path}...")
    df = pd.read_csv(tsv_path, sep="\t")
    required = {"native", "miRNAs", "Genes", "p_value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df


def load_background_gaf(gaf_path, obodag, taxon="9606"):
    """Load background gene annotations from GAF file for IC computation."""
    if not gaf_path or not Path(gaf_path).exists():
        print("Warning: GAF file not provided or not found. Using dataset-derived IC.")
        return None
    
    print(f"Loading background GAF from {gaf_path}...")
    gene2gos = defaultdict(set)
    opener = gzip.open if str(gaf_path).endswith(".gz") else open
    
    with opener(gaf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line[0] == "!":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            
            db_object_id = parts[1]
            go_id = parts[4]
            taxon_col = parts[12]
            
            if "NOT" in parts[3]:
                continue
            if taxon and taxon not in taxon_col:
                continue
            if go_id not in obodag:
                continue
            
            gene2gos[db_object_id].add(go_id)
    
    if gene2gos:
        return TermCounts(obodag, dict(gene2gos))
    return None


# --- SIMILARITY COMPUTATION ---
def lin_similarity(t1, t2, obodag, tc):
    """Compute Lin similarity between two GO terms."""
    if t1 not in obodag or t2 not in obodag:
        return 0.0
    
    ic1 = float(tc.get_ic(t1)) if hasattr(tc, 'get_ic') else float(tc.get_info_content(t1))
    ic2 = float(tc.get_ic(t2)) if hasattr(tc, 'get_ic') else float(tc.get_info_content(t2))
    
    denom = ic1 + ic2
    if denom <= 0:
        return 0.0
    
    # Compute IC of MICA
    anc1 = set(obodag[t1].get_all_parents()) | {t1}
    anc2 = set(obodag[t2].get_all_parents()) | {t2}
    common = anc1.intersection(anc2)
    
    if not common:
        return 0.0
    
    ic_mica = max([
        float(tc.get_ic(a)) if hasattr(tc, 'get_ic') else float(tc.get_info_content(a))
        for a in common
    ])
    
    if ic_mica <= 0:
        return 0.0
    
    sim = (2.0 * ic_mica) / denom
    return float(max(0.0, min(1.0, sim)))


def compute_similarity_matrix(go_ids, obodag, tc):
    """Compute similarity matrix for GO IDs."""
    n = len(go_ids)
    sim_matrix = np.zeros((n, n), dtype=float)
    
    print(f"Computing similarity matrix for {n} GO terms...")
    for i in range(n):
        for j in range(i, n):
            if i == j:
                sim_matrix[i, j] = 1.0
            else:
                s = lin_similarity(go_ids[i], go_ids[j], obodag, tc)
                sim_matrix[i, j] = sim_matrix[j, i] = s
    
    return sim_matrix


# --- CLUSTERING ---
def cluster_terms(go_ids, sim_matrix, sim_threshold=0.6):
    """Perform hierarchical clustering of GO terms."""
    dist_matrix = 1.0 - sim_matrix
    
    print(f"Performing hierarchical clustering (threshold={sim_threshold})...")
    condensed = squareform(dist_matrix, checks=False)
    Z = linkage(condensed, method="average")
    
    # Cut dendrogram at distance threshold
    dist_threshold = 1.0 - sim_threshold
    cluster_labels = fcluster(Z, t=dist_threshold, criterion="distance")
    
    return Z, cluster_labels


def extract_cluster_info(go_ids, cluster_labels, pvalue_map, obodag, 
                         pvalue_weight=0.3, min_cluster_size=2):
    """Extract cluster representatives and summary statistics."""
    clusters_info = []
    
    for cluster_id in np.unique(cluster_labels):
        members = [go_ids[i] for i in range(len(go_ids)) if cluster_labels[i] == cluster_id]
        
        if len(members) < min_cluster_size:
            continue
        
        # Select representative based on p-value and centrality
        pvalues = np.array([pvalue_map.get(m, 1.0) for m in members])
        
        # Score = p-value (lower is better) + IC-based centrality
        if len(members) > 1:
            ics = np.array([
                float(getattr(obodag[m], '_ic', 0)) for m in members
            ])
            # Normalize
            scores = pvalue_weight * (pvalues / (pvalues.max() + 1e-6)) + \
                     (1 - pvalue_weight) * (1 - ics / (ics.max() + 1e-6))
            rep_idx = np.argmin(scores)
        else:
            rep_idx = 0
        
        rep_id = members[rep_idx]
        rep_name = obodag[rep_id].name if rep_id in obodag else rep_id
        
        clusters_info.append({
            'cluster_id': int(cluster_id),
            'representative_id': rep_id,
            'representative_name': rep_name,
            'cluster_size': len(members),
            'member_count': len(members),
            'min_pvalue': pvalues.min(),
            'mean_pvalue': pvalues.mean(),
            'member_ids': ';'.join(members),
            'method': 'custom',
            'cluster_type': 'hierarchical_with_pvalue'
        })
    
    return pd.DataFrame(clusters_info).sort_values('min_pvalue')


# --- VISUALIZATION ---
def plot_dendrogram(Z, go_ids, output_path, sim_threshold=0.6):
    """Plot hierarchical clustering dendrogram."""
    fig, ax = plt.subplots(figsize=(14, 8))
    
    dendrogram(
        Z,
        labels=go_ids[:len(set(Z[:, 0].astype(int)) | set(Z[:, 1].astype(int)))],
        ax=ax,
        orientation='right'
    )
    
    ax.set_title("GO Term Dendrogram (Custom Clustering)", fontsize=14, fontweight='bold')
    ax.set_xlabel("Distance (1 - Lin Similarity)", fontsize=12)
    ax.set_ylabel("GO Terms", fontsize=12)
    
    # Draw threshold line
    ax.axvline(x=1-sim_threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold ({sim_threshold})')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Dendrogram saved to {output_path}")
    plt.close()


# --- MAIN ---
def main():
    args = get_args()
    
    # Load data
    obodag = load_go_obo(args.go_obo)
    df = load_enrichment_results(args.input_file)
    
    # Extract unique GO IDs and p-values
    go_ids = sorted(df['native'].unique())
    pvalue_map = dict(zip(df['native'], df['p_value']))
    
    # Load background IC (optional)
    tc = load_background_gaf(args.gaf, obodag)
    if tc is None:
        # Use dataset-derived IC
        gene2gos = defaultdict(set)
        for _, row in df.iterrows():
            genes = [g.strip() for g in str(row['Genes']).split(',') if g]
            for g in genes:
                gene2gos[g].add(str(row['native']))
        tc = TermCounts(obodag, dict(gene2gos))
    
    # Compute similarity matrix
    sim_matrix = compute_similarity_matrix(go_ids, obodag, tc)
    
    # Cluster terms
    Z, cluster_labels = cluster_terms(go_ids, sim_matrix, args.sim_threshold)
    
    # Extract cluster information
    clusters_df = extract_cluster_info(
        go_ids, cluster_labels, pvalue_map, obodag,
        pvalue_weight=args.pvalue_weight,
        min_cluster_size=args.min_cluster_size
    )
    
    # Save results
    output_file = f"{args.out_prefix}_custom_clusters.tsv"
    clusters_df.to_csv(output_file, sep="\t", index=False)
    print(f"Clusters saved to {output_file}")
    
    # Visualize
    plot_file = f"{args.out_prefix}_custom_dendrogram.png"
    plot_dendrogram(Z, go_ids, plot_file, args.sim_threshold)
    
    print(f"\nClustering complete!")
    print(f"  GO terms: {len(go_ids)}")
    print(f"  Clusters: {len(clusters_df)}")
    print(f"  Avg cluster size: {clusters_df['cluster_size'].mean():.2f}")


if __name__ == "__main__":
    main()
