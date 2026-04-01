#!/usr/bin/env python3
"""Steps 5 and 6: LD pruning plus PCA."""

from __future__ import annotations

from _common import ensure_dir, line_count, load_context, run_command


def main() -> int:
    """Generate the LD-pruned PCA input and compute principal components."""
    ctx = load_context()
    indir = ctx.run_dir / "qc" / "final"
    prune_dir = ensure_dir(ctx.run_dir / "pca" / "pruned_snps")
    pca_dir = ensure_dir(ctx.run_dir / "pca" / "eigenvec_eigenval")
    prune_prefix = prune_dir / "ld_prune"

    print("=== Step 5: LD pruning for PCA ===")

    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_qc"),
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
    print(f"LD-pruned SNP set for PCA: {n_pruned_in} variants")

    print("=== Step 5b: Save LD-pruned PLINK binary (for downstream analyses, e.g. ADMIXTURE) ===")
    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_qc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--extract",
            str(prune_prefix.with_suffix(".prune.in")),
            "--make-bed",
            "--out",
            str(prune_dir / "graega_ldpruned"),
        ],
        env=ctx.env,
    )

    print("=== Step 6: PCA computation ===")
    run_command(
        [
            "plink2",
            "--bfile",
            str(indir / "graega_qc"),
            "--chr-set",
            *ctx.chr_set_args,
            "--extract",
            str(prune_prefix.with_suffix(".prune.in")),
            "--pca",
            ctx.n_pcs,
            "--out",
            str(pca_dir / "graega_pca"),
        ],
        env=ctx.env,
    )

    print(f"PCA complete: {ctx.n_pcs} PCs computed")
    print(f"Eigenvectors: {pca_dir / 'graega_pca.eigenvec'}")
    print(f"Eigenvalues:  {pca_dir / 'graega_pca.eigenval'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
