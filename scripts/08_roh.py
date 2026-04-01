#!/usr/bin/env python3
"""Step 8: prepare ROH input, run detectRUNS, and generate ROH figures."""

from __future__ import annotations

import sys

from _common import (
    copy_plink_prefix,
    count_king_removed,
    ensure_dir,
    line_count,
    load_context,
    run_command,
)


def main() -> int:
    """Run the ROH-specific QC path and downstream summaries."""
    ctx = load_context()
    prep_dir = ctx.run_dir / "roh" / "input_prep"
    step0_dir = ensure_dir(prep_dir / "step0_autosomal")
    step1_dir = ensure_dir(prep_dir / "step1_snp_qc_nomaf")
    step2_dir = ensure_dir(prep_dir / "step2_sample_qc")
    step3_dir = ensure_dir(prep_dir / "step3_relatedness")
    input_dir = ensure_dir(ctx.run_dir / "roh" / "input")

    print("=========================================")
    print(" enipro ROH pipeline")
    print(f" Run: {ctx.run_stamp}")
    print(f" Output: {ctx.run_dir}")
    print("=========================================")

    print("=== ROH prep 0: autosomal biallelic SNPs ===")
    run_command(
        [
            "plink2",
            "--bfile",
            str(ctx.input_prefix),
            "--chr-set",
            *ctx.chr_set_args,
            "--chr",
            ctx.autosomes,
            "--snps-only",
            "just-acgt",
            "--max-alleles",
            "2",
            "--make-bed",
            "--out",
            str(step0_dir / "graega_roh_autosomal"),
        ],
        env=ctx.env,
    )

    print("=== ROH prep 1: SNP QC without MAF filtering ===")
    run_command(
        [
            "plink2",
            "--bfile",
            str(step0_dir / "graega_roh_autosomal"),
            "--chr-set",
            *ctx.chr_set_args,
            "--geno",
            ctx.qc_geno,
            "--hwe",
            ctx.qc_hwe,
            "--make-bed",
            "--out",
            str(step1_dir / "graega_roh_snpqc"),
        ],
        env=ctx.env,
    )

    print("=== ROH prep 2: sample QC ===")
    run_command(
        [
            "plink2",
            "--bfile",
            str(step1_dir / "graega_roh_snpqc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--mind",
            ctx.qc_mind,
            "--make-bed",
            "--out",
            str(step2_dir / "graega_roh_sampleqc"),
        ],
        env=ctx.env,
    )

    print("=== ROH prep 3: relatedness filtering ===")
    run_command(
        [
            "plink2",
            "--bfile",
            str(step2_dir / "graega_roh_sampleqc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--indep-pairwise",
            ctx.ld_window,
            ctx.ld_step,
            ctx.ld_r2,
            "--out",
            str(step3_dir / "ld_prune"),
        ],
        env=ctx.env,
    )

    run_command(
        [
            "plink2",
            "--bfile",
            str(step2_dir / "graega_roh_sampleqc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--extract",
            str(step3_dir / "ld_prune.prune.in"),
            "--king-cutoff",
            ctx.qc_king,
            "--out",
            str(step3_dir / "king"),
        ],
        env=ctx.env,
    )

    removed_path = step3_dir / "king.king.cutoff.out.id"
    n_removed = count_king_removed(removed_path)
    roh_input_prefix = input_dir / "graega_roh_input"
    if n_removed > 0:
        print(f"Removing {n_removed} related sample(s) from the unpruned ROH input")
        run_command(
            [
                "plink2",
                "--bfile",
                str(step2_dir / "graega_roh_sampleqc"),
                "--chr-set",
                *ctx.chr_set_args,
                "--remove",
                str(removed_path),
                "--make-bed",
                "--out",
                str(roh_input_prefix),
            ],
            env=ctx.env,
        )
    else:
        print("No related samples found for ROH input")
        copy_plink_prefix(step2_dir / "graega_roh_sampleqc", roh_input_prefix)

    print("=== ROH prep 4: PED/MAP export for detectRUNS ===")
    run_command(
        [
            "plink",
            "--bfile",
            str(roh_input_prefix),
            "--chr-set",
            *ctx.chr_set_args,
            "--recode",
            "--out",
            str(roh_input_prefix),
        ],
        env=ctx.env,
    )

    n_variants = line_count(roh_input_prefix.with_suffix(".bim"))
    n_samples = line_count(roh_input_prefix.with_suffix(".fam"))
    print(f"ROH input dataset: {n_variants} variants, {n_samples} samples")

    print("=== Step 8a: detectRUNS ROH scan ===")
    run_command(
        [
            "Rscript",
            str(ctx.repo_dir / "src" / "run_roh.R"),
            str(ctx.config_path),
        ],
        env=ctx.env,
    )

    print("=== Step 8b: ROH plots ===")
    run_command(
        [
            sys.executable,
            str(ctx.repo_dir / "src" / "plot_roh.py"),
            str(ctx.config_path),
        ],
        env=ctx.env,
    )

    print("ROH analysis complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
