#!/usr/bin/env python3
"""Render pairwise FST outputs as a clustered heatmap."""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.config import get_breed_labels, get_breed_order, load_config
from src.utils import read_fam, read_fst_summary


SUMMARY_BASENAME = "graega_ldpruned_fst_wc"
PRUNED_FAM_BASENAME = "graega_ldpruned.fam"


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


def get_leaf_order(linkage_matrix, labels: list[str]) -> list[str]:
    """Return labels in dendrogram leaf order."""
    dendro = dendrogram(linkage_matrix, no_plot=True)
    return [labels[idx] for idx in dendro["leaves"]]


def get_breed_counts(results: Path, breeds: list[str]) -> dict[str, int]:
    """Count samples per breed from the LD-pruned FAM file when available."""
    fam_path = results / "pca" / "pruned_snps" / PRUNED_FAM_BASENAME
    if not fam_path.exists():
        return {breed: 0 for breed in breeds}

    fam = read_fam(fam_path)
    counts = fam["FID"].value_counts()
    return {breed: int(counts.get(breed, 0)) for breed in breeds}


def format_top_labels(
    breeds: list[str],
    breed_labels: dict[str, str],
    breed_counts: dict[str, int],
) -> list[str]:
    """Build top-axis labels with sample sizes."""
    labels = []
    for breed in breeds:
        label = breed_labels.get(breed, breed)
        count = breed_counts.get(breed, 0)
        if count:
            label = f"{label}\n(n={count})"
        labels.append(label)
    return labels


def format_side_labels(breeds: list[str], breed_labels: dict[str, str]) -> list[str]:
    """Build right-axis labels for row lookup."""
    return [breed_labels.get(breed, breed) for breed in breeds]


def draw_top_dendrogram(ax, linkage_matrix, n_items: int) -> None:
    """Draw a top dendrogram aligned to heatmap cell centers."""
    dendro = dendrogram(linkage_matrix, no_plot=True)
    for xs, ys in zip(dendro["icoord"], dendro["dcoord"]):
        x_coords = (np.asarray(xs) - 5.0) / 10.0 + 0.5
        ax.plot(x_coords, ys, color="#8f8f8f", linewidth=0.9)

    max_height = max(max(ys) for ys in dendro["dcoord"])
    ax.set_xlim(0, n_items)
    ax.set_ylim(0, max_height * 1.06 if max_height > 0 else 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def annotate_upper_triangle(ax, matrix: pd.DataFrame, vmax: float, fontsize: float) -> None:
    """Add compact value labels only for visible upper-triangle cells."""
    values = matrix.to_numpy()
    threshold = vmax * 0.58 if vmax > 0 else 0

    for row_idx in range(values.shape[0]):
        for col_idx in range(row_idx + 1, values.shape[1]):
            value = values[row_idx, col_idx]
            color = "white" if value >= threshold else "#4d4d4d"
            ax.text(
                col_idx + 0.5,
                row_idx + 0.5,
                f"{value:.3f}",
                ha="center",
                va="center",
                fontsize=fontsize,
                color=color,
            )


def position_colorbar(fig, ax_heatmap, ax_cbar) -> None:
    """Place the colorbar after layout is resolved."""
    heat_pos = ax_heatmap.get_position()
    cbar_width = 0.018
    cbar_height = heat_pos.height * 0.42
    cbar_x = min(0.96 - cbar_width, heat_pos.x1 + 0.11)
    cbar_y = heat_pos.y0 + (heat_pos.height - cbar_height) / 2
    ax_cbar.set_position([cbar_x, cbar_y, cbar_width, cbar_height])


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
    matrix = build_matrix(summary, value_col, breed_order)

    if matrix.isna().any().any():
        missing = int(matrix.isna().sum().sum())
        raise ValueError(f"FST matrix contains {missing} missing pairwise values")

    linkage_matrix = build_linkage(matrix)
    leaf_order = get_leaf_order(linkage_matrix, matrix.index.tolist())
    stats_matrix = matrix.loc[leaf_order, leaf_order]

    matrix_path = fst_dir / f"{SUMMARY_BASENAME}.matrix.tsv"
    stats_matrix.to_csv(matrix_path, sep="\t", float_format="%.6f")

    breed_counts = get_breed_counts(results, leaf_order)
    top_labels = format_top_labels(leaf_order, breed_labels, breed_counts)
    side_labels = format_side_labels(leaf_order, breed_labels)

    n_breeds = len(stats_matrix)
    heatmap_size = max(4.8, 0.50 * n_breeds)
    fig_width = heatmap_size + 2.9
    fig_height = heatmap_size + 1.9
    annot_size = max(4.6, 6.0 - 0.10 * n_breeds)
    top_tick_size = max(6.1, 6.9 - 0.08 * n_breeds)
    side_tick_size = max(6.6, 7.4 - 0.06 * n_breeds)

    sns.set_theme(style="white", context="paper")
    vmax = float(np.nanmax(stats_matrix.to_numpy()))
    mask = np.tril(np.ones(stats_matrix.shape, dtype=bool), k=0)
    cmap = sns.color_palette("YlOrRd", as_cmap=True)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1.0, heatmap_size], hspace=0.06)
    ax_dendrogram = fig.add_subplot(gs[0, 0])
    ax_heatmap = fig.add_subplot(gs[1, 0])
    ax_cbar = fig.add_axes([0.88, 0.28, 0.02, 0.34])
    fig.subplots_adjust(left=0.08, right=0.77, top=0.89, bottom=0.10)

    draw_top_dendrogram(ax_dendrogram, linkage_matrix, n_breeds)
    sns.heatmap(
        stats_matrix,
        mask=mask,
        ax=ax_heatmap,
        cmap=cmap,
        vmin=0,
        vmax=vmax if vmax > 0 else 1,
        square=True,
        linewidths=0.6,
        linecolor="#ffffff",
        cbar_ax=ax_cbar,
        cbar_kws={"label": "Pairwise FST"},
    )
    ax_heatmap.set_facecolor("#ffffff")
    annotate_upper_triangle(ax_heatmap, stats_matrix, vmax, annot_size)

    ax_heatmap.set_xlabel("")
    ax_heatmap.set_ylabel("")
    ax_heatmap.tick_params(
        axis="x",
        top=True,
        bottom=False,
        labeltop=True,
        labelbottom=False,
        length=0,
        pad=2,
    )
    ax_heatmap.tick_params(
        axis="y",
        left=False,
        right=True,
        labelleft=False,
        labelright=True,
        length=0,
        pad=4,
    )
    ax_heatmap.set_xticklabels(top_labels, rotation=90, ha="center", va="bottom", fontsize=top_tick_size)
    ax_heatmap.set_yticklabels(side_labels, rotation=0, fontsize=side_tick_size)
    for spine in ax_heatmap.spines.values():
        spine.set_visible(False)

    position_colorbar(fig, ax_heatmap, ax_cbar)
    ax_cbar.yaxis.set_major_locator(MaxNLocator(4))
    ax_cbar.tick_params(labelsize=top_tick_size - 0.2, length=0)
    ax_cbar.set_ylabel("Pairwise FST", fontsize=top_tick_size)

    fig.suptitle("Pairwise FST Clustered Heatmap", y=0.975, fontsize=13, fontweight="semibold")
    fig.text(
        0.5,
        0.935,
        "Weir-Cockerham estimator, LD-pruned SNP set, upper triangle shown",
        ha="center",
        fontsize=8,
        color="#666666",
    )
    fig.text(
        0.08,
        0.035,
        "Average-linkage clustering of pairwise FST; the dendrogram is descriptive, not phylogenetic.",
        fontsize=7.4,
        color="#666666",
    )
    fig.savefig(fig_dir / "fst_distance_matrix.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    mask = ~np.eye(len(stats_matrix), dtype=bool)
    off_diag = stats_matrix.where(mask)
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
