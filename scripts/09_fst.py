#!/usr/bin/env python3
"""Step 9: pairwise FST on the LD-pruned dataset."""

from __future__ import annotations

import sys

from _common import ensure_dir, load_context, run_command


def main() -> int:
    """Compute pairwise FST and render the summary matrix."""
    ctx = load_context()
    indir = ctx.run_dir / "pca" / "pruned_snps"
    outdir = ensure_dir(ctx.run_dir / "fst")
    input_prefix = indir / "graega_ldpruned"
    output_prefix = outdir / "graega_ldpruned_fst_wc"

    if not input_prefix.with_suffix(".bed").exists():
        print(f"ERROR: {input_prefix}.bed not found.", file=sys.stderr)
        print(
            "Run scripts/06_pca.py first, or export RUN_STAMP for an existing run directory.",
            file=sys.stderr,
        )
        return 1

    print("=== Step 9: Pairwise FST on LD-pruned data ===")
    run_command(
        [
            "plink2",
            "--bfile",
            str(input_prefix),
            "--chr-set",
            *ctx.chr_set_args,
            "--family",
            "BREED",
            "--fst",
            "BREED",
            f"method=wc",
            f"blocksize={ctx.fst_blocksize}",
            "--out",
            str(output_prefix),
        ],
        env=ctx.env,
    )

    run_command(
        [
            "Rscript",
            str(ctx.repo_dir / "src" / "plot_fst.R"),
            str(ctx.config_path),
        ],
        env=ctx.env,
    )

    print(f"FST summary: {output_prefix}.fst.summary")
    print(f"FST matrix:  {output_prefix}.matrix.tsv")
    print(f"Figures:     {ctx.run_dir / 'figures' / 'fst'}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
