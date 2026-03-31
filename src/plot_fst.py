#!/usr/bin/env python3
"""Render pairwise FST outputs as a clustered heatmap."""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.config import get_breed_labels, get_breed_order, get_breed_palette, load_config
from src.utils import read_fst_summary


SUMMARY_BASENAME = "graega_ldpruned_fst_wc"


def build_matrix(summary: pd.DataFrame, value_col: str, breed_order: list[str]) -> pd.DataFrame:
    """Convert PLINK's long pairwise summary into a symmetric distance matrix."""
    populations = set(summary["POP1"]) | set(summary["POP2"])
    order = [breed for breed in breed_order if breed in populations]
    matrix = pd.DataFrame(np.nan, index=order, columns=order, dtype=float)

    for row in summary.itertuples(index=False):
        matrix.loc[row.POP1, row.POP2] = getattr(row, value_col)
        matrix.loc[row.POP2, row.POP1] = getattr(row, value_col)

    for idx in range(len(matrix)):
        matrix.iat[idx, idx] = 0.0
    return matrix


def build_linkage(matrix: pd.DataFrame):
    """Construct a hierarchical clustering from the symmetric FST matrix."""
    condensed = squareform(matrix.to_numpy(), checks=False)
    return linkage(condensed, method="average")


def main(config_path: str | None = None):
    cfg = load_config(config_path)
    results = Path(os.environ.get("ENIPRO_RUN_DIR") or cfg["results_dir"])
    fst_dir = results / "fst"
    fig_dir = results / "figures" / "fst"
    fig_dir.mkdir(parents=True, exist_ok=True)

    summary_path = fst_dir / f"{SUMMARY_BASENAME}.fst.summary"
    summary = read_fst_summary(summary_path)

    value_col = "WC_FST" if "WC_FST" in summary.columns else "HUDSON_FST"
    breed_order = get_breed_order(cfg)
    breed_labels = get_breed_labels(cfg)
    breed_palette = get_breed_palette(cfg)
    matrix = build_matrix(summary, value_col, breed_order)

    if matrix.isna().any().any():
        missing = int(matrix.isna().sum().sum())
        raise ValueError(f"FST matrix contains {missing} missing pairwise values")

    matrix_path = fst_dir / f"{SUMMARY_BASENAME}.matrix.tsv"
    matrix.to_csv(matrix_path, sep="\t", float_format="%.6f")

    label_map = {breed: breed_labels.get(breed, breed) for breed in matrix.index}
    plot_matrix = matrix.rename(index=label_map, columns=label_map)
    plot_palette = {label_map[breed]: breed_palette[breed] for breed in matrix.index}
    label_order = [label_map[breed] for breed in matrix.index]
    color_strip = pd.Series(label_order, index=label_order).map(plot_palette)
    linkage_matrix = build_linkage(matrix)

    sns.set_theme(style="white")
    vmax = float(np.nanmax(matrix.to_numpy()))
    cluster_grid = sns.clustermap(
        plot_matrix,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        row_colors=color_strip,
        col_colors=color_strip,
        cmap="mako",
        vmin=0,
        vmax=vmax if vmax > 0 else 1,
        annot=True,
        fmt=".3f",
        linewidths=0.5,
        linecolor="white",
        dendrogram_ratio=(0.14, 0.14),
        colors_ratio=(0.04, 0.04),
        figsize=(10, 9.5),
        cbar_kws={"label": "Pairwise FST"},
    )
    cluster_grid.ax_heatmap.set_xlabel("")
    cluster_grid.ax_heatmap.set_ylabel("")
    cluster_grid.ax_heatmap.set_xticklabels(
        cluster_grid.ax_heatmap.get_xticklabels(),
        rotation=45,
        ha="right",
    )
    cluster_grid.ax_heatmap.set_yticklabels(
        cluster_grid.ax_heatmap.get_yticklabels(),
        rotation=0,
    )
    cluster_grid.fig.suptitle(
        "Pairwise FST Clustered Heatmap\nWeir-Cockerham, LD-pruned SNP set",
        y=1.02,
    )
    cluster_grid.fig.subplots_adjust(top=0.94)
    cluster_grid.savefig(fig_dir / "fst_distance_matrix.png", dpi=300, bbox_inches="tight")
    cluster_grid.savefig(fig_dir / "fst_distance_matrix.pdf", bbox_inches="tight")
    plt.close(cluster_grid.fig)

    mask = ~np.eye(len(matrix), dtype=bool)
    off_diag = matrix.where(mask)
    min_idx = off_diag.stack().idxmin()
    max_idx = off_diag.stack().idxmax()
    min_val = off_diag.loc[min_idx]
    max_val = off_diag.loc[max_idx]

    print(f"FST summary: {summary_path}")
    print(f"Matrix written to {matrix_path}")
    print(f"Figures written to {fig_dir}/")
    print(f"Closest pair: {min_idx[0]} vs {min_idx[1]} ({min_val:.4f})")
    print(f"Most differentiated pair: {max_idx[0]} vs {max_idx[1]} ({max_val:.4f})")


if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(config_path)
