#!/usr/bin/env python3
"""Step 3: relatedness check and final QC dataset creation."""

from __future__ import annotations

from _common import (
    copy_plink_prefix,
    count_king_removed,
    ensure_dir,
    line_count,
    load_context,
    run_command,
)


def main() -> int:
    """LD prune, run KING, and remove related individuals if needed."""
    ctx = load_context()
    indir = ctx.run_dir / "qc" / "step2_sample_qc"
    outdir = ensure_dir(ctx.run_dir / "qc" / "step3_relatedness")
    final_dir = ensure_dir(ctx.run_dir / "qc" / "final")
    prune_prefix = outdir / "ld_prune"
    king_prefix = outdir / "king"
    final_prefix = final_dir / "graega_qc"

    print("=== Step 3: Relatedness check ===")

    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_sampleqc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--indep-pairwise",
            ctx.ld_window,
            ctx.ld_step,
            ctx.ld_r2,
            "--out",
            str(prune_prefix),
        ],
        env=ctx.env,
    )

    n_pruned_in = line_count(prune_prefix.with_suffix(".prune.in"))
    print(f"LD-pruned SNP set: {n_pruned_in} variants")

    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_sampleqc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--extract",
            str(prune_prefix.with_suffix(".prune.in")),
            "--king-cutoff",
            ctx.qc_king,
            "--out",
            str(king_prefix),
        ],
        env=ctx.env,
    )

    removed_path = outdir / "king.king.cutoff.out.id"
    n_removed = count_king_removed(removed_path)
    if n_removed > 0:
        print(f"Removing {n_removed} related sample(s)")
        run_command(
            [
                "plink2",
                "--bfile",
                str(indir / "graega_sampleqc"),
                "--chr-set",
                *ctx.chr_set_args,
                "--remove",
                str(removed_path),
                "--make-bed",
                "--out",
                str(final_prefix),
            ],
            env=ctx.env,
        )
    else:
        print("No related samples found — copying to final")
        copy_plink_prefix(indir / "graega_sampleqc", final_prefix)

    n_variants = line_count(final_prefix.with_suffix(".bim"))
    n_samples = line_count(final_prefix.with_suffix(".fam"))
    print(f"Final QC dataset: {n_variants} variants, {n_samples} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
