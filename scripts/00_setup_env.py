#!/usr/bin/env python3
"""Install all pipeline dependencies into the enipro micromamba environment."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


FALLBACK_MAMBA = Path("/home/i/iapostof/y/micromamba")

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(line_buffering=True)


def resolve_mamba() -> str:
    """Locate the micromamba executable."""
    env_value = os.environ.get("MAMBA_EXE")
    if env_value:
        return env_value

    discovered = shutil.which("micromamba")
    if discovered:
        return discovered

    if FALLBACK_MAMBA.is_file() and os.access(FALLBACK_MAMBA, os.X_OK):
        return str(FALLBACK_MAMBA)

    raise SystemExit("micromamba executable not found")


def run(command: list[str], *, capture_output: bool = False) -> subprocess.CompletedProcess[str]:
    """Run a subprocess and optionally capture its stdout."""
    return subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=capture_output,
    )


def print_first_line(command: list[str]) -> None:
    """Run a command and print only the first stdout line."""
    completed = run(command, capture_output=True)
    first_line = completed.stdout.splitlines()[0]
    print(first_line)


def main() -> int:
    """Install packages and print tool versions."""
    mamba_bin = resolve_mamba()

    run(
        [
            mamba_bin,
            "install",
            "-n",
            "enipro",
            "-c",
            "conda-forge",
            "-c",
            "bioconda",
            "--channel-priority",
            "flexible",
            "-y",
            "plink2",
            "plink",
            "python=3.11",
            "pandas",
            "numpy",
            "matplotlib-base",
            "seaborn",
            "scipy",
            "pyyaml",
            "r-base",
            "r-yaml",
            "r-ggplot2",
            "r-data.table",
            "r-rcpp",
            "r-plyr",
            "r-iterators",
            "r-itertools",
            "r-gridextra",
            "r-reshape2",
        ]
    )

    print("Environment setup complete.")

    run(
        [
            mamba_bin,
            "run",
            "-n",
            "enipro",
            "Rscript",
            "-e",
            (
                "if (!requireNamespace('detectRUNS', quietly = TRUE)) "
                "install.packages('detectRUNS', repos = 'https://cran.r-project.org')"
            ),
        ]
    )

    print_first_line([mamba_bin, "run", "-n", "enipro", "plink2", "--version"])
    print_first_line([mamba_bin, "run", "-n", "enipro", "plink", "--version"])
    print_first_line([mamba_bin, "run", "-n", "enipro", "python", "--version"])
    print_first_line([mamba_bin, "run", "-n", "enipro", "R", "--version"])
    print_first_line(
        [
            mamba_bin,
            "run",
            "-n",
            "enipro",
            "Rscript",
            "-e",
            "cat(as.character(packageVersion('detectRUNS')), '\\n')",
        ]
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
