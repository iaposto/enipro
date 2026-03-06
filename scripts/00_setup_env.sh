#!/usr/bin/env bash
# Install all dependencies into the enipro micromamba environment.
set -euo pipefail

micromamba install -n enipro -c conda-forge -c bioconda --channel-priority flexible -y \
    plink2 plink python=3.11 \
    pandas numpy matplotlib-base seaborn scipy pyyaml

echo "Environment setup complete."
micromamba run -n enipro plink2 --version
micromamba run -n enipro plink --version | head -1
micromamba run -n enipro python --version
