#!/usr/bin/env python3
"""ADMIXTURE visualization: CV error curve and Q-matrix stacked bar plots."""

import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.config import load_config, get_breed_colors, get_breed_labels
from src.utils import read_fam


# K values to show in Q-matrix panel (around the CV minimum)
PLOT_K_VALUES = [5, 6, 7, 8, 10]


def parse_cv_errors(admix_dir: Path) -> dict[int, float]:
    """Read CV error from each log{K}.out file."""
    cv = {}
    for log_file in admix_dir.glob("log*.out"):
        m = re.search(r"CV error \(K=(\d+)\): ([0-9.]+)", log_file.read_text())
        if m:
            cv[int(m.group(1))] = float(m.group(2))
    return dict(sorted(cv.items()))


def plot_cv_error(cv: dict[int, float], fig_dir: Path):
    """Line plot of CV error vs K, highlighting the minimum."""
    ks = list(cv.keys())
    errs = list(cv.values())
    best_k = min(cv, key=cv.get)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(ks, errs, "o-", color="#377eb8", linewidth=1.5, markersize=5)
    ax.axvline(best_k, color="#e41a1c", linestyle="--", linewidth=1, label=f"K={best_k} (min CV={cv[best_k]:.5f})")
    ax.set_xlabel("K (number of ancestral populations)")
    ax.set_ylabel("Cross-validation error (10-fold)")
    ax.set_title("ADMIXTURE Cross-Validation Error")
    ax.set_xticks(ks)
    ax.tick_params(axis="x", labelsize=7)
    ax.legend()
    plt.tight_layout()
    out = fig_dir / "admixture_cv_error.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")
    return best_k


def plot_q_matrix(
    q_file: Path,
    fam: pd.DataFrame,
    breed_colors: dict,
    breed_labels: dict,
    K: int,
    fig_dir: Path,
):
    """Stacked bar chart for one K value, samples sorted by breed."""
    q = pd.read_csv(q_file, sep=r"\s+", header=None)
    q.columns = [f"C{i+1}" for i in range(K)]
    df = pd.concat([fam[["FID", "IID"]].reset_index(drop=True), q], axis=1)
    df = df.sort_values("FID").reset_index(drop=True)

    # Assign consistent colors: components sorted by dominant breed
    breed_order = sorted(breed_colors.keys())
    component_cols = [f"C{i+1}" for i in range(K)]

    # Use a qualitative palette for components (independent of breed colors)
    cmap = plt.get_cmap("tab20")
    comp_colors = [cmap(i / max(K - 1, 1)) for i in range(K)]

    # Breed boundary positions
    breed_groups = df.groupby("FID", sort=False)
    breed_sizes = df["FID"].value_counts()[breed_order]
    breed_sizes = breed_sizes[[b for b in breed_order if b in df["FID"].values]]

    fig, ax = plt.subplots(figsize=(14, 3))

    bottoms = np.zeros(len(df))
    for i, col in enumerate(component_cols):
        ax.bar(range(len(df)), df[col].values, bottom=bottoms,
               color=comp_colors[i], width=1.0, linewidth=0)
        bottoms += df[col].values

    # Breed boundary lines and labels
    pos = 0
    for breed in breed_order:
        if breed not in breed_sizes.index:
            continue
        n = breed_sizes[breed]
        if pos > 0:
            ax.axvline(pos - 0.5, color="white", linewidth=0.8)
        ax.text(pos + n / 2 - 0.5, 1.02, breed_labels.get(breed, breed),
                ha="center", va="bottom", fontsize=7, rotation=45)
        pos += n

    ax.set_xlim(-0.5, len(df) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("Ancestry proportion")
    ax.set_xticks([])
    ax.set_title(f"ADMIXTURE  K={K}")

    plt.tight_layout()
    out = fig_dir / f"admixture_K{K}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_q_panel(
    admix_dir: Path,
    fam: pd.DataFrame,
    breed_colors: dict,
    breed_labels: dict,
    k_values: list[int],
    fig_dir: Path,
):
    """Multi-row panel of Q-matrix plots for several K values."""
    n = len(k_values)
    fig, axes = plt.subplots(n, 1, figsize=(14, 2.5 * n))
    if n == 1:
        axes = [axes]

    breed_order = sorted(breed_colors.keys())
    cmap = plt.get_cmap("tab20")

    for ax, K in zip(axes, k_values):
        q_file = admix_dir / f"graega_ldpruned.{K}.Q"
        if not q_file.exists():
            print(f"WARNING: {q_file} not found, skipping K={K}")
            ax.set_visible(False)
            continue

        q = pd.read_csv(q_file, sep=r"\s+", header=None)
        q.columns = [f"C{i+1}" for i in range(K)]
        df = pd.concat([fam[["FID", "IID"]].reset_index(drop=True), q], axis=1)
        df = df.sort_values("FID").reset_index(drop=True)

        comp_colors = [cmap(i / max(K - 1, 1)) for i in range(K)]
        component_cols = [f"C{i+1}" for i in range(K)]
        breed_sizes = df["FID"].value_counts()
        breed_sizes = pd.Series(
            {b: breed_sizes.get(b, 0) for b in breed_order if b in df["FID"].values}
        )

        bottoms = np.zeros(len(df))
        for i, col in enumerate(component_cols):
            ax.bar(range(len(df)), df[col].values, bottom=bottoms,
                   color=comp_colors[i], width=1.0, linewidth=0)
            bottoms += df[col].values

        pos = 0
        for breed, n_breed in breed_sizes.items():
            if pos > 0:
                ax.axvline(pos - 0.5, color="white", linewidth=0.8)
            ax.text(pos + n_breed / 2 - 0.5, 1.03,
                    breed_labels.get(breed, breed),
                    ha="center", va="bottom", fontsize=7)
            pos += n_breed

        ax.set_xlim(-0.5, len(df) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        ax.set_ylabel(f"K={K}", fontsize=9)
        ax.set_xticks([])

    plt.suptitle("ADMIXTURE ancestry proportions", y=1.01, fontsize=11)
    plt.tight_layout()
    out = fig_dir / "admixture_panel.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def main(config_path: str | None = None):
    cfg = load_config(config_path)
    results = Path(os.environ.get("ENIPRO_RUN_DIR") or cfg["results_dir"])
    admix_dir = results / "admixture"
    fam_path = results / "pca" / "pruned_snps" / "graega_ldpruned.fam"
    fig_dir = results / "figures" / "admixture"
    fig_dir.mkdir(parents=True, exist_ok=True)

    breed_colors = get_breed_colors(cfg)
    breed_labels = get_breed_labels(cfg)
    fam = read_fam(fam_path)

    # 1. CV error curve
    cv = parse_cv_errors(admix_dir)
    if not cv:
        print("ERROR: no CV errors found in log files", file=sys.stderr)
        sys.exit(1)
    best_k = plot_cv_error(cv, fig_dir)
    print(f"Optimal K by CV: {best_k}  (CV error={cv[best_k]:.5f})")

    # 2. Individual Q-matrix plots for K values around the minimum
    available_ks = [k for k in PLOT_K_VALUES if (admix_dir / f"graega_ldpruned.{k}.Q").exists()]
    for K in available_ks:
        plot_q_matrix(admix_dir / f"graega_ldpruned.{K}.Q", fam,
                      breed_colors, breed_labels, K, fig_dir)

    # 3. Multi-K panel
    plot_q_panel(admix_dir, fam, breed_colors, breed_labels, available_ks, fig_dir)

    print(f"\nAll figures saved to {fig_dir}")


if __name__ == "__main__":
    main()
