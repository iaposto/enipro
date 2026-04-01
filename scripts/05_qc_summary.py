#!/usr/bin/env python3
"""Step 4: QC summary statistics."""

from __future__ import annotations

from _common import ensure_dir, load_context, run_command


def main() -> int:
    """Generate missingness, heterozygosity, and frequency summaries."""
    ctx = load_context()
    indir = ctx.run_dir / "qc" / "final"
    outdir = ensure_dir(ctx.run_dir / "qc" / "stats")

    print("=== Step 4: QC summary statistics ===")

    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_qc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--missing",
            "--out",
            str(outdir / "missing"),
        ],
        env=ctx.env,
    )

    run_command(
        [
            "plink",
            "--bfile",
            str(indir / "graega_qc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--het",
            "--out",
            str(outdir / "het"),
        ],
        env=ctx.env,
    )

    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_qc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--freq",
            "--out",
            str(outdir / "freq"),
        ],
        env=ctx.env,
    )

    print(f"Summary statistics written to {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
