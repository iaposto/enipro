#!/usr/bin/env python3
"""PCA visualization: scree plot, PC scatter plots, and pairplot."""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.config import (
    get_breed_labels,
    get_breed_order,
    get_breed_palette,
    load_config,
)
from src.utils import read_eigenvec, read_eigenval, read_fam


def main(config_path: str | None = None):
    cfg = load_config(config_path)
    results = Path(os.environ.get("ENIPRO_RUN_DIR") or cfg["results_dir"])
    pca_dir = results / "pca" / "eigenvec_eigenval"
    fig_dir = results / "figures" / "pca"
    fig_dir.mkdir(parents=True, exist_ok=True)

    breed_labels = get_breed_labels(cfg)
    breed_palette = get_breed_palette(cfg)

    # Load PCA results
    eigenvec = read_eigenvec(pca_dir / "graega_pca.eigenvec")
    eigenval = read_eigenval(pca_dir / "graega_pca.eigenval")

    # Proportion of variance explained
    total_var = sum(eigenval)
    pve = [ev / total_var * 100 for ev in eigenval]

    # Add breed info
    eigenvec["breed"] = eigenvec["FID"]

    breed_order = [breed for breed in get_breed_order(cfg) if breed in eigenvec["breed"].unique()]
    palette = {breed: breed_palette[breed] for breed in breed_order}

    sns.set_style("whitegrid")

    # --- 1. Scree plot ---
    fig, ax = plt.subplots(figsize=(10, 4))
    x = range(1, len(pve) + 1)
    ax.bar(x, pve, color="#377eb8", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Variance Explained (%)")
    ax.set_title("PCA Scree Plot")
    ax.set_xticks(list(x))
    # Annotate top PCs
    for i in range(min(5, len(pve))):
        ax.text(i + 1, pve[i] + 0.2, f"{pve[i]:.1f}%", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "scree_plot.png", dpi=300)
    plt.close(fig)

    # --- Helper for PC scatter ---
    def plot_pcs(pc_x: int, pc_y: int, filename: str):
        col_x = f"PC{pc_x}"
        col_y = f"PC{pc_y}"
        fig, ax = plt.subplots(figsize=(8, 7))
        for breed in breed_order:
            mask = eigenvec["breed"] == breed
            ax.scatter(
                eigenvec.loc[mask, col_x],
                eigenvec.loc[mask, col_y],
                c=palette[breed],
                label=f"{breed} ({breed_labels.get(breed, breed)})",
                s=30, alpha=0.8, edgecolors="white", linewidth=0.3,
            )
        ax.set_xlabel(f"PC{pc_x} ({pve[pc_x - 1]:.1f}%)")
        ax.set_ylabel(f"PC{pc_y} ({pve[pc_y - 1]:.1f}%)")
        ax.set_title(f"PCA: PC{pc_x} vs PC{pc_y}")
        ax.legend(
            bbox_to_anchor=(1.02, 1), loc="upper left",
            fontsize=8, frameon=True,
        )
        fig.tight_layout()
        fig.savefig(fig_dir / f"{filename}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # --- 2-4. PC scatter plots ---
    plot_pcs(1, 2, "pca_pc1_pc2")
    plot_pcs(1, 3, "pca_pc1_pc3")
    plot_pcs(2, 3, "pca_pc2_pc3")

    # --- 5. Pairplot of PC1-PC4 ---
    pc_cols = ["PC1", "PC2", "PC3", "PC4"]
    plot_df = eigenvec[["breed"] + pc_cols].copy()
    plot_df["breed"] = pd.Categorical(plot_df["breed"], categories=breed_order, ordered=True)

    g = sns.pairplot(
        plot_df, hue="breed", vars=pc_cols,
        palette=palette, diag_kind="kde",
        plot_kws={"s": 20, "alpha": 0.7, "edgecolor": "white", "linewidth": 0.2},
        height=2.2,
    )
    g.figure.suptitle("PCA Pairplot (PC1–PC4)", y=1.01)
    g.savefig(fig_dir / "pca_pairplot_4pc.png", dpi=300, bbox_inches="tight")
    plt.close(g.figure)

    print(f"PCA figures saved to {fig_dir}/")
    print(f"Variance explained by first 5 PCs: {', '.join(f'{v:.1f}%' for v in pve[:5])}")


if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(config_path)
