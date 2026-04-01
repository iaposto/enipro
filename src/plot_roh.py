#!/usr/bin/env python3
"""ROH visualization: length classes and FROH."""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.config import get_breed_labels, get_breed_order, get_breed_palette, load_config
from src.plotting import plot_breed_boxplot


CLASS_ORDER = ["short (1-5 Mb)", "medium (5-10 Mb)", "long (>10 Mb)"]
CLASS_COLORS = {
    "short (1-5 Mb)": "#4daf4a",
    "medium (5-10 Mb)": "#ff7f00",
    "long (>10 Mb)": "#e41a1c",
}


def main(config_path: str | None = None):
    cfg = load_config(config_path)
    results = Path(os.environ.get("ENIPRO_RUN_DIR") or cfg["results_dir"])
    roh_dir = results / "roh"
    fig_dir = results / "figures" / "roh"
    fig_dir.mkdir(parents=True, exist_ok=True)

    froh = pd.read_csv(roh_dir / "stats" / "froh_per_individual.tsv", sep="\t", dtype={"breed": str, "IID": str})
    class_summary = pd.read_csv(
        roh_dir / "stats" / "roh_length_class_per_individual.tsv",
        sep="\t",
        dtype={"breed": str, "IID": str, "length_class": str},
    )

    breed_labels = get_breed_labels(cfg)
    breed_order = [breed for breed in get_breed_order(cfg) if breed in froh["breed"].unique()]
    breed_palette = get_breed_palette(cfg)

    sns.set_style("whitegrid")

    # Mean ROH count per individual in each length class, by breed.
    breed_class_mean = (
        class_summary.groupby(["breed", "length_class"], sort=False)["n_roh"]
        .mean()
        .reset_index()
    )
    pivot = (
        breed_class_mean.pivot(index="breed", columns="length_class", values="n_roh")
        .reindex(breed_order)
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for length_class in CLASS_ORDER:
        values = pivot[length_class].to_numpy() if length_class in pivot.columns else np.zeros(len(pivot))
        ax.bar(x, values, bottom=bottom, color=CLASS_COLORS[length_class], width=0.75, label=length_class)
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels([breed_labels.get(breed, breed) for breed in pivot.index], rotation=45, ha="right")
    ax.set_ylabel("Mean ROH count per individual")
    ax.set_title("ROH Length Classes by Breed")
    ax.legend(title="ROH length")
    fig.tight_layout()
    fig.savefig(fig_dir / "roh_length_classes_by_breed.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # FROH distribution by breed.
    froh["breed"] = pd.Categorical(froh["breed"], categories=breed_order, ordered=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    plot_breed_boxplot(
        data=froh,
        x="breed",
        y="froh",
        breed_order=breed_order,
        breed_palette=breed_palette,
        breed_labels=breed_labels,
        ax=ax,
    )
    ax.set_xlabel("Breed")
    ax.set_ylabel("FROH")
    ax.set_title("Genomic Inbreeding (FROH) by Breed")
    fig.tight_layout()
    fig.savefig(fig_dir / "froh_by_breed.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"ROH figures saved to {fig_dir}/")


if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(config_path)
