#!/usr/bin/env python3
"""Step 1: SNP-level QC."""

from __future__ import annotations

from _common import ensure_dir, line_count, load_context, run_command


def main() -> int:
    """Apply SNP missingness, MAF, and HWE filters."""
    ctx = load_context()
    indir = ctx.run_dir / "qc" / "step0_autosomal"
    outdir = ensure_dir(ctx.run_dir / "qc" / "step1_snp_qc")
    output_prefix = outdir / "graega_snpqc"

    print("=== Step 1: SNP-level QC ===")
    print(f"Thresholds: --geno {ctx.qc_geno} --maf {ctx.qc_maf} --hwe {ctx.qc_hwe}")

    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_autosomal"),
            "--chr-set",
            *ctx.chr_set_args,
            "--geno",
            ctx.qc_geno,
            "--maf",
            ctx.qc_maf,
            "--hwe",
            ctx.qc_hwe,
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
