#!/usr/bin/env python3
"""Step 2: sample-level QC."""

from __future__ import annotations

from _common import ensure_dir, line_count, load_context, run_command


def main() -> int:
    """Apply per-sample missingness filtering."""
    ctx = load_context()
    indir = ctx.run_dir / "qc" / "step1_snp_qc"
    outdir = ensure_dir(ctx.run_dir / "qc" / "step2_sample_qc")
    output_prefix = outdir / "graega_sampleqc"

    print("=== Step 2: Sample-level QC ===")
    print(f"Threshold: --mind {ctx.qc_mind}")

    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_snpqc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--mind",
            ctx.qc_mind,
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
