# plotting.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram
from dotenv import load_dotenv
import os

load_dotenv(file = "paths.env")
output_dir = os.getenv("PY_RESULT")

# ------------------------------------------------------------
# DOTPLOT
# ------------------------------------------------------------
def dotplot_clusters(df_all_clusters, df_enrich, out_prefix, image_suffix):
    """
    Dotplot of representative GO terms colored by cluster label.

    Parameters
    ----------
    df_all_clusters : DataFrame
        Cluster table containing rep_go_id, cluster_label, cluster_size.
    df_enrich : DataFrame
        Enrichment table containing native, p_value.
    out_prefix : str
        Path prefix for saved plot.
    image_suffix : str
        Identifier added to output file name.
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

    if df_plot.empty:
        plt.figure(figsize=(8, 3))
        plt.text(0.5, 0.5, "No p-values available", ha="center", va="center")
        plt.axis("off")
        fname = f"{out_prefix}_dotplot_{image_suffix}.png"
        plt.savefig(fname, dpi=300, bbox_inches="tight")
        plt.close()
        return fname

    df_plot["neglog10_p"] = -np.log10(df_plot["p_value"])

    # Plot
    fig_h = max(6, min(14, len(df_plot["cluster_label"].unique()) * 0.35))
    plt.figure(figsize=(10, fig_h))

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
    )

    plt.xlabel("-log10(p-value)")
    plt.ylabel("Cluster Label")
    plt.title("GO Term Clusters Dotplot")
    plt.tight_layout()

    fname = f"{out_prefix}_dotplot_{image_suffix}.png"
    plt.savefig(fname, dpi=300)
    plt.close()
    return fname


# ------------------------------------------------------------
# HEATMAP
# ------------------------------------------------------------
def heatmap_similarity(sim_matrix, terms, out_prefix, namespace, image_suffix):
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
    out_prefix : str
        Output file prefix.
    """

    n = len(terms)

    plt.figure(figsize=(max(6, min(16, n * 0.3)), max(5, min(14, n * 0.3))))
    sns.heatmap(sim_matrix, xticklabels=terms, yticklabels=terms, cmap="viridis", cbar=True)
    plt.title(f"GO Term Lin Similarity Heatmap [{namespace}] (n={n})")
    plt.tight_layout()

    fname = f"{out_prefix}_heatmap_{image_suffix}_{namespace}.png"
    plt.savefig(fname, dpi=300)
    plt.close()
    return fname


# ------------------------------------------------------------
# DENDROGRAM
# ------------------------------------------------------------
def dendrogram_terms(Z, terms, obodag, out_prefix, namespace, dist_threshold, sim_threshold,
                     image_suffix):
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
    out_prefix : str
        Output file prefix.
    dist_threshold : float
        Threshold for dendrogram line coloring.
    sim_threshold : float
        Lin similarity threshold (displayed on title).
    """

    labels = [obodag[t].name for t in terms]

    plt.figure(figsize=(12, 6))
    dendrogram(
        Z,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=dist_threshold,
    )

    plt.title(f"GO Terms Dendrogram [{namespace}] (Lin sim > {sim_threshold})")
    plt.tight_layout()

    fname = f"{out_prefix}_dendrogram_{image_suffix}_{namespace}.png"
    plt.savefig(fname, dpi=300)
    plt.close()
    return fname

# ------------------------------------------------------------
# Save figures
