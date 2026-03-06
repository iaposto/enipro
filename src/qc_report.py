#!/usr/bin/env python3
"""Generate QC summary report with tables and diagnostic plots."""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Allow running as script from repo root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.config import load_config, get_breed_colors
from src.utils import read_fam, read_het, read_smiss, read_vmiss, read_freq


def main(config_path: str | None = None):
    cfg = load_config(config_path)
    results = Path(cfg["results_dir"])
    stats_dir = results / "qc" / "stats"
    fig_dir = results / "figures" / "qc"
    fig_dir.mkdir(parents=True, exist_ok=True)

    breed_colors = get_breed_colors(cfg)

    # --- Load data ---
    fam = read_fam(results / "qc" / "final" / "graega_qc.fam")
    fam.rename(columns={"FID": "breed"}, inplace=True)

    smiss = read_smiss(stats_dir / "missing.smiss")
    smiss.rename(columns={"#FID": "breed"}, inplace=True)

    vmiss = read_vmiss(stats_dir / "missing.vmiss")

    het = read_het(stats_dir / "het.het")

    freq = read_freq(stats_dir / "freq.afreq")

    # --- Textual summary ---
    breed_counts = fam["breed"].value_counts().sort_index()

    # Count variants at each QC step by reading BIM files
    steps = {
        "Raw input": Path(cfg["input_prefix"]).with_suffix(".bim"),
        "Step 0 (autosomal)": results / "qc" / "step0_autosomal" / "graega_autosomal.bim",
        "Step 1 (SNP QC)": results / "qc" / "step1_snp_qc" / "graega_snpqc.bim",
        "Step 2 (sample QC)": results / "qc" / "step2_sample_qc" / "graega_sampleqc.bim",
        "Final (post-relatedness)": results / "qc" / "final" / "graega_qc.bim",
    }

    print("\n=== QC Summary Report ===\n")
    print("Variant counts by step:")
    for label, bim_path in steps.items():
        if bim_path.exists():
            n = sum(1 for _ in open(bim_path))
            print(f"  {label:30s} {n:>6d}")

    # Sample counts at each step
    fam_steps = {
        "Raw input": Path(cfg["input_prefix"]).with_suffix(".fam"),
        "Step 2 (sample QC)": results / "qc" / "step2_sample_qc" / "graega_sampleqc.fam",
        "Final": results / "qc" / "final" / "graega_qc.fam",
    }
    print("\nSample counts by step:")
    for label, fam_path in fam_steps.items():
        if fam_path.exists():
            n = sum(1 for _ in open(fam_path))
            print(f"  {label:30s} {n:>6d}")

    print("\nPer-breed sample counts (final):")
    for breed, count in breed_counts.items():
        print(f"  {breed:6s} {count:>4d}")

    # Heterozygosity summary
    het["O_HET"] = (het["N(NM)"] - het["O(HOM)"]) / het["N(NM)"]
    het_by_breed = het.merge(fam[["breed", "IID"]], on="IID")

    print("\nPer-breed heterozygosity (observed):")
    for breed, grp in het_by_breed.groupby("breed"):
        print(f"  {breed:6s}  mean={grp['O_HET'].mean():.4f}  sd={grp['O_HET'].std():.4f}")

    print(f"\nMean inbreeding coefficient F: {het['F'].mean():.4f} (sd={het['F'].std():.4f})")

    # --- Plots ---
    sns.set_style("whitegrid")

    # 1. MAF distribution
    fig, ax = plt.subplots(figsize=(8, 4))
    alt_freq = freq["ALT_FREQS"].values
    maf = np.minimum(alt_freq, 1 - alt_freq)
    ax.hist(maf, bins=50, color="#377eb8", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Minor Allele Frequency")
    ax.set_ylabel("Number of SNPs")
    ax.set_title("MAF Distribution (post-QC)")
    fig.tight_layout()
    fig.savefig(fig_dir / "maf_distribution.png", dpi=300)
    fig.savefig(fig_dir / "maf_distribution.pdf")
    plt.close(fig)

    # 2. Per-sample missingness
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(smiss["F_MISS"].values, bins=50, color="#4daf4a", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Missing Rate")
    ax.set_ylabel("Number of Samples")
    ax.set_title("Per-Sample Missing Rate (post-QC)")
    fig.tight_layout()
    fig.savefig(fig_dir / "sample_missingness.png", dpi=300)
    fig.savefig(fig_dir / "sample_missingness.pdf")
    plt.close(fig)

    # 3. Per-SNP missingness
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(vmiss["F_MISS"].values, bins=50, color="#e41a1c", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Missing Rate")
    ax.set_ylabel("Number of SNPs")
    ax.set_title("Per-SNP Missing Rate (post-QC)")
    fig.tight_layout()
    fig.savefig(fig_dir / "snp_missingness.png", dpi=300)
    fig.savefig(fig_dir / "snp_missingness.pdf")
    plt.close(fig)

    # 4. Heterozygosity by breed
    breed_order = sorted(breed_colors.keys())
    palette = [breed_colors[b] for b in breed_order]

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(
        data=het_by_breed, x="breed", y="O_HET",
        hue="breed", order=breed_order, palette=palette, legend=False, ax=ax,
    )
    ax.set_xlabel("Breed")
    ax.set_ylabel("Observed Heterozygosity")
    ax.set_title("Observed Heterozygosity by Breed")
    fig.tight_layout()
    fig.savefig(fig_dir / "het_by_breed.png", dpi=300)
    fig.savefig(fig_dir / "het_by_breed.pdf")
    plt.close(fig)

    # 5. Inbreeding F by breed
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(
        data=het_by_breed, x="breed", y="F",
        hue="breed", order=breed_order, palette=palette, legend=False, ax=ax,
    )
    ax.set_xlabel("Breed")
    ax.set_ylabel("Inbreeding Coefficient (F)")
    ax.set_title("Inbreeding Coefficient by Breed")
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(fig_dir / "inbreeding_by_breed.png", dpi=300)
    fig.savefig(fig_dir / "inbreeding_by_breed.pdf")
    plt.close(fig)

    print(f"\nQC figures saved to {fig_dir}/")


if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(config_path)
