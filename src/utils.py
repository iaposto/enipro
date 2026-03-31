"""Utility functions for parsing PLINK output files."""

from pathlib import Path

import pandas as pd


def read_fam(path: str | Path) -> pd.DataFrame:
    """Parse a PLINK .fam file into a DataFrame."""
    return pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "father", "mother", "sex", "pheno"],
        dtype={"FID": str, "IID": str},
    )


def read_bim(path: str | Path) -> pd.DataFrame:
    """Parse a PLINK .bim file into a DataFrame."""
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "snp_id", "cm", "bp", "a1", "a2"],
        dtype={"chr": str, "snp_id": str},
    )


def read_eigenvec(path: str | Path) -> pd.DataFrame:
    """Parse PLINK2 .eigenvec (with header) into a DataFrame."""
    df = pd.read_csv(path, sep="\t", dtype={"#FID": str, "IID": str})
    df.rename(columns={"#FID": "FID"}, inplace=True)
    return df


def read_eigenval(path: str | Path) -> list[float]:
    """Parse PLINK2 .eigenval file into a list of eigenvalues."""
    with open(path) as f:
        return [float(line.strip()) for line in f if line.strip()]


def read_het(path: str | Path) -> pd.DataFrame:
    """Parse PLINK 1.9 .het output."""
    return pd.read_csv(
        path,
        sep=r"\s+",
        dtype={"FID": str, "IID": str},
    )


def read_smiss(path: str | Path) -> pd.DataFrame:
    """Parse PLINK2 .smiss (per-sample missingness)."""
    return pd.read_csv(path, sep="\t", dtype={"#FID": str, "IID": str})


def read_vmiss(path: str | Path) -> pd.DataFrame:
    """Parse PLINK2 .vmiss (per-variant missingness)."""
    return pd.read_csv(path, sep="\t", dtype={"#CHROM": str})


def read_freq(path: str | Path) -> pd.DataFrame:
    """Parse PLINK2 .afreq file."""
    return pd.read_csv(path, sep="\t", dtype={"#CHROM": str})


def read_fst_summary(path: str | Path) -> pd.DataFrame:
    """Parse PLINK2 .fst.summary output."""
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={"#POP1": str, "POP1": str, "POP2": str},
    )
    df.rename(columns=lambda col: col.lstrip("#"), inplace=True)
    return df
