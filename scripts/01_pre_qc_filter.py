#!/usr/bin/env python3
"""Step 0: filter to autosomal biallelic SNPs only."""

from __future__ import annotations

from _common import ensure_dir, line_count, load_context, run_command


def main() -> int:
    """Run the pre-QC autosomal SNP filter."""
    ctx = load_context()
    outdir = ensure_dir(ctx.run_dir / "qc" / "step0_autosomal")
    output_prefix = outdir / "graega_autosomal"

    print("=== Step 0: Pre-QC filter — autosomal biallelic SNPs ===")
    print(f"Input: {ctx.input_prefix}")

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
            str(output_prefix),
        ],
        env=ctx.env,
    )

    n_variants = line_count(output_prefix.with_suffix(".bim"))
    n_samples = line_count(output_prefix.with_suffix(".fam"))
    print(f"Output: {n_variants} variants, {n_samples} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
