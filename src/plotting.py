"""Shared plotting helpers for consistent breed-based figures."""

from __future__ import annotations

import seaborn as sns


def apply_breed_xticklabels(ax, breed_order: list[str], breed_labels: dict[str, str]) -> None:
    """Set breed tick labels in a consistent order and style."""
    ticks = ax.get_xticks()
    labels = [breed_labels.get(breed, breed) for breed in breed_order]
    ax.set_xticks(ticks, labels, rotation=45, ha="right")


def plot_breed_boxplot(
    *,
    data,
    x: str,
    y: str,
    breed_order: list[str],
    breed_palette: dict[str, str],
    breed_labels: dict[str, str],
    ax,
    width: float = 0.6,
) -> None:
    """Draw a consistently styled breed-colored boxplot."""
    palette = [breed_palette[breed] for breed in breed_order]
    sns.boxplot(
        data=data,
        x=x,
        y=y,
        hue=x,
        order=breed_order,
        hue_order=breed_order,
        palette=palette,
        dodge=False,
        width=width,
        linewidth=1.0,
        fliersize=2,
        legend=False,
        saturation=1,
        ax=ax,
    )
    apply_breed_xticklabels(ax, breed_order, breed_labels)
