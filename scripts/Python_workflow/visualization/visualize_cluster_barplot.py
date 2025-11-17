#!/usr/bin/env python3
"""
Visualize GO clustering results as a bar plot.
- x-axis: rep_go_name (wrapped for readability)
- y-axis: cluster_size

Usage examples:
  python3 Python_GSEA/visualize_cluster_barplot.py \
      --input Python_GSEA/Results/clustering_goterms_4_hsa_miR_107_clusters_hsa_miR_107.tsv \
      --output Python_GSEA/Results/clustering_goterms_4_hsa_miR_107_cluster_barplot.png

Optional flags:
  --top 10           Show only top 10 clusters by size
  --ascending        Sort ascending (default is descending)
  --wrap 25          Wrap GO term labels to given width (chars)
  --rotate 45        Rotate x tick labels by degrees
  --title "..."       Custom plot title

How to run:
python3 Python_GSEA/visualize_cluster_barplot.py --input Python_GSEA/Results/clustering_goterms_4_hsa_miR_107_clusters_hsa_miR_107.tsv --output Python_GSEA/Results/clustering_goterms_4_hsa_miR_107_cluster_barplot.png --title "GO cluster sizes: hsa-miR-107"
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys
import textwrap

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



def wrap_label(s: str, width: int = 25) -> str:
    """Wrap a label string to multiple lines for readability."""
    if pd.isna(s):
        return ""
    return textwrap.fill(str(s), width=width)


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


essential_palette = "viridis"


def make_plot(
    df: pd.DataFrame,
    out_path: str | Path,
    title: str | None = None,
    top_n: int | None = None,
    sort_desc: bool = True,
    label_wrap: int = 25,
    rotate: int = 45,
) -> None:
    """Create and save the bar plot."""
    # Sort and optionally truncate to top N
    df_sorted = df.sort_values("cluster_size", ascending=not sort_desc).copy()
    if top_n is not None:
        df_sorted = df_sorted.head(top_n)

    # Wrap labels for readability
    df_sorted["rep_go_name_wrapped"] = df_sorted["rep_go_name"].apply(lambda s: wrap_label(s, width=label_wrap))
    order = df_sorted["rep_go_name_wrapped"].tolist()

    # Dynamic figure width based on number of bars
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

    # Rotate x tick labels for readability
    for tick in ax.get_xticklabels():
        tick.set_rotation(rotate)
        tick.set_horizontalalignment("right")

    # Annotate bar values
    try:
        if ax.containers:
            ax.bar_label(ax.containers[0], fmt="%d", padding=2, fontsize=9)
    except Exception:
        pass

    plt.tight_layout()

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Barplot of cluster_size by rep_go_name from GO clustering TSV."
    )
    default_in = (
        "Python_GSEA/Results/"
        "clustering_goterms_4_hsa_miR_107_clusters_hsa_miR_107.tsv"
    )
    default_out = (
        "Python_GSEA/Results/"
        "clustering_goterms_4_hsa_miR_107_cluster_barplot.png"
    )

    parser.add_argument(
        "-i",
        "--input",
        default=default_in,
        help="Path to TSV file with columns rep_go_name and cluster_size.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=default_out,
        help="Output image file (PNG/SVG/PDF).",
    )
    parser.add_argument(
        "--title",
        default="GO cluster sizes by representative term",
        help="Plot title",
    )
    parser.add_argument(
        "--top", type=int, default=None, help="Show only top N clusters."
    )
    parser.add_argument(
        "--ascending",
        action="store_true",
        help="Sort ascending instead of descending.",
    )
    parser.add_argument(
        "--wrap", type=int, default=25, help="Wrap width for GO term labels."
    )
    parser.add_argument(
        "--rotate", type=int, default=45, help="Rotation angle for x tick labels."
    )

    args = parser.parse_args()

    df = load_table(args.input)
    make_plot(
        df,
        args.output,
        args.title,
        args.top,
        sort_desc=(not args.ascending),
        label_wrap=args.wrap,
        rotate=args.rotate,
    )


if __name__ == "__main__":
    main()
