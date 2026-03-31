#!/usr/bin/env bash
# Install all dependencies into the enipro micromamba environment.
set -euo pipefail

MAMBA_BIN="${MAMBA_EXE:-}"
if [[ -z "${MAMBA_BIN}" ]]; then
    MAMBA_BIN="$(command -v micromamba || true)"
fi
if [[ -z "${MAMBA_BIN}" && -x "/home/i/iapostof/y/micromamba" ]]; then
    MAMBA_BIN="/home/i/iapostof/y/micromamba"
fi
if [[ -z "${MAMBA_BIN}" ]]; then
    echo "micromamba executable not found" >&2
    exit 1
fi

"${MAMBA_BIN}" install -n enipro -c conda-forge -c bioconda --channel-priority flexible -y \
    plink2 plink python=3.11 \
    pandas numpy matplotlib-base seaborn scipy pyyaml \
    r-base r-yaml r-ggplot2 r-data.table r-rcpp r-plyr r-iterators r-itertools r-gridextra r-reshape2

echo "Environment setup complete."
"${MAMBA_BIN}" run -n enipro Rscript -e "if (!requireNamespace('detectRUNS', quietly = TRUE)) install.packages('detectRUNS', repos = 'https://cran.r-project.org')"
"${MAMBA_BIN}" run -n enipro plink2 --version
"${MAMBA_BIN}" run -n enipro plink --version | head -1
"${MAMBA_BIN}" run -n enipro python --version
"${MAMBA_BIN}" run -n enipro R --version | head -1
"${MAMBA_BIN}" run -n enipro Rscript -e "cat(as.character(packageVersion('detectRUNS')), '\n')"
