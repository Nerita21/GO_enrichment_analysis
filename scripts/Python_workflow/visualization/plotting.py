#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram
from __future__ import annotations
import argparse
import sys
import textwrap

# ------------------------------------------------------------
# ARGUMENT PARSING
# ------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot GO term clustering results: dotplot, heatmap, dendrogram."
    )
    parser.add_argument(
        "--clusters",
        required=True,
        help="Path to cluster table TSV with columns: rep_go_id, cluster_label, cluster_size.",
    )
    parser.add_argument(
        "--enrichment",
        required=True,
        help="Path to enrichment table TSV with columns: native, p_value.",
    )
    parser.add_argument(
        "--plot_data_json",
        required=True,
        help="Json file with precomputed data for plotting (similarity matrix, terms, linkage matrix, obodag, etc.).",
    )
    parser.add_argument("--out-prefix", default="", help="Optional prefix for output filenames.")
    return parser.parse_args()

# ------------------------------------------------------------
# LOADING AND MANIPULATING DATA
# ------------------------------------------------------------

def load_json(path):
    """Load JSON file containing precomputed plot data."""
    try:
        import json
        with open(path, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        sys.exit(f"Failed to read {path}: {e}")

def load_table(path: str) -> pd.DataFrame:
    """Load TSV with at least columns: rep_go_name, cluster_size."""
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.exit(f"Failed to read {path}: {e}")

    required_cols = {"rep_go_name", "cluster_size"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(
            "Missing required columns: "
            + ", ".join(sorted(missing))
            + f". Columns found: {list(df.columns)}"
        )

    # Ensure numeric cluster_size and drop invalid rows
    df["cluster_size"] = pd.to_numeric(df["cluster_size"], errors="coerce")
    df = df.dropna(subset=["rep_go_name", "cluster_size"]).copy()
    return df

def wrap_label(s: str, width: int = 25) -> str:
    """Wrap a label string to multiple lines for readability."""
    if pd.isna(s):
        return ""
    return textwrap.fill(str(s), width=width)

# ------------------------------------------------------------
# DOTPLOT
# ------------------------------------------------------------
def dotplot_clusters(df_all_clusters, df_enrich):
    """
    Dotplot of representative GO terms colored by cluster label.

    Parameters
    ----------
    df_all_clusters : DataFrame
        Cluster table containing rep_go_id, cluster_label, cluster_size.
    df_enrich : DataFrame
        Enrichment table containing native, p_value.
    """

    df_plot = df_all_clusters.merge(
        df_enrich[['native', 'p_value']],
        left_on="rep_go_id",
        right_on="native",
        how="left"
    )

    # p-value cleaning
    df_plot["p_value"] = pd.to_numeric(df_plot["p_value"], errors="coerce")
    df_plot = df_plot.dropna(subset=["p_value"])
    df_plot.loc[df_plot["p_value"] <= 0, "p_value"] = 1e-300

    # Case: No data → return empty fig with a message
    if df_plot.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No p-values available", ha="center", va="center")
        ax.axis("off")
        return fig
       

    df_plot["neglog10_p"] = -np.log10(df_plot["p_value"])

    # Plot
    fig_h = max(6, min(14, len(df_plot["cluster_label"].unique()) * 0.35))
    fig, ax = plt.subplots(figsize=(10, fig_h))

    sns.scatterplot(
        data=df_plot,
        x="neglog10_p",
        y="cluster_label",
        size="cluster_size",
        hue="cluster_label",
        palette="tab20",
        sizes=(40, 400),
        alpha=0.85,
        edgecolor="white",
        linewidth=0.5,
        legend="brief",
        ax=ax
    )

    ax.set_xlabel("-log10(p-value)")
    ax.set_ylabel("Cluster Label")
    ax.set_title("GO Term Clusters Dotplot")
    fig.tight_layout()

    return fig
   

# ------------------------------------------------------------
# HEATMAP
# ------------------------------------------------------------
def heatmap_similarity(sim_matrix, terms, namespace):
    """
    Heatmap of GO term Lin similarity.

    Parameters
    ----------
    sim_matrix : 2D numpy array or DataFrame
        Symmetric similarity matrix.
    terms : list
        Term IDs in order of the matrix.
    namespace : str
        BP / MF / CC indicator.
    """
    # Case: No data → return empty fig with a message
    if sim_matrix.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No p-values available", ha="center", va="center")
        ax.axis("off")
        return fig
    
    n = len(terms)

    fig, ax = plt.subplots(
        figsize=(max(6, min(16, n * 0.3)),
                 max(5, min(14, n * 0.3)))
    )

    sns.heatmap(
        sim_matrix,
        xticklabels=terms,
        yticklabels=terms,
        cmap="viridis",
        cbar=True,
        ax=ax
    )

    ax.set_title(f"GO Term Lin Similarity Heatmap [{namespace}] (n={n})")
    fig.tight_layout()

    return fig



# ------------------------------------------------------------
# DENDROGRAM
# ------------------------------------------------------------
def dendrogram_terms(Z, terms, obodag, namespace, dist_threshold, sim_threshold):
    """
    GO term dendrogram with term names.

    Parameters
    ----------
    Z : linkage matrix
        Hierarchical clustering linkage matrix.
    terms : list
        GO term IDs corresponding to leaves.
    obodag : dict-like GO DAG object
        Must support obodag[term].name to retrieve GO term names.
    namespace : str
        BP / MF / CC indicator.
    dist_threshold : float
        Threshold for dendrogram line coloring.
    sim_threshold : float
        Lin similarity threshold (displayed on title).
    """
    # Case: No data → return empty fig with a message
    if terms.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No p-values available", ha="center", va="center")
        ax.axis("off")
        return fig

    labels = [obodag[t].name for t in terms]

    fig, ax = plt.subplots(figsize=(12, 6))

    dendrogram(
        Z,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=dist_threshold,
        ax=ax
    )

    ax.set_title(f"GO Terms Dendrogram [{namespace}] (Lin sim > {sim_threshold})")
    fig.tight_layout()

    return fig

# ------------------------------------------------------------
# BARPLOT
# ------------------------------------------------------------
essential_palette = "viridis"

def bar_plot(
    df: pd.DataFrame,
    title: str | None = None,
    top_n: int | None = None,
    sort_desc: bool = True,
    label_wrap: int = 25,
    rotate: int = 45,
):
    """Create a bar plot for cluster sizes. Returns the figure."""

    # Sort and optionally truncate to top N
    df_sorted = df.sort_values("cluster_size", ascending=not sort_desc).copy()
    if top_n is not None:
        df_sorted = df_sorted.head(top_n)

    # Wrap labels for readability
    df_sorted["rep_go_name_wrapped"] = df_sorted["rep_go_name"].apply(
        lambda s: wrap_label(s, width=label_wrap)
    )
    order = df_sorted["rep_go_name_wrapped"].tolist()

    # Dynamic figure width
    n = len(df_sorted)
    width = max(8, min(40, 0.55 * n))
    height = 6

    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=False)

    sns.barplot(
        data=df_sorted,
        x="rep_go_name_wrapped",
        y="cluster_size",
        order=order,
        palette=essential_palette,
        ax=ax,
    )

    ax.set_xlabel("Representative GO term (wrapped)")
    ax.set_ylabel("Cluster size")
    if title:
        ax.set_title(title)

    # Rotate tick labels
    for tick in ax.get_xticklabels():
        tick.set_rotation(rotate)
        tick.set_horizontalalignment("right")

    # Annotate bar values
    try:
        if ax.containers:
            ax.bar_label(ax.containers[0], fmt="%d", padding=2, fontsize=9)
    except Exception:
        pass

    fig.tight_layout()
    return fig

# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():
    args = parse_args()

    # Load data
    df_clusters = load_table(args.clusters)
    df_enrich = load_table(args.enrichment)
    df_plot_data = load_json(args.plot_data_json) 

    # Create plots
    barplot_fig = bar_plot(
        df_clusters,
        title="GO Term Cluster Sizes",
        top_n=None,
        sort_desc=True,
        label_wrap=25,
        rotate=45,
    )
    dotplot_fig = dotplot_clusters(df_clusters, df_enrich)
    heatmap_fig = heatmap_similarity(
        df_plot_data["similarity_matrix"],
        df_plot_data["terms"],
        df_plot_data["namespace"]
    )
    dendrogram_fig = dendrogram_terms(
        df_plot_data["linkage_matrix"],
        df_plot_data["terms"],
        df_plot_data["obodag"],
        df_plot_data["namespace"],
        df_plot_data["dist_threshold"],
        df_plot_data["sim_threshold"]
    )
    # Save plots
    barplot_fig.savefig(args.barplot, dpi=300)
    plt.close(barplot_fig)

    dotplot_fig.savefig(args.dotplot, dpi=300)
    plt.close(dotplot_fig)

    heatmap_fig.savefig(args.heatmap, dpi=300)
    plt.close(heatmap_fig)

    dendrogram_fig.savefig(args.dendrogram, dpi=300)
    plt.close(dendrogram_fig)

if __name__ == "__main__":
    main()