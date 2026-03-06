#!/usr/bin/env bash
# Shared variables for all pipeline scripts.
# Sourced by each step script — not run directly.

# Load paths from params.yaml using a simple grep/sed approach
# (avoids Python dependency in shell scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../config/params.yaml"

_yaml_val() {
    grep "^${1}:" "${CONFIG}" | sed "s/^${1}:[[:space:]]*//"
}

INPUT_PREFIX=$(_yaml_val "input_prefix")
RESULTS_DIR=$(_yaml_val "results_dir")
LOG_DIR=$(_yaml_val "log_dir")
CHR_SET=$(_yaml_val "chr_set")
AUTOSOMES=$(_yaml_val "autosomes" | tr -d '"')

# QC thresholds (nested under qc: — parse with indentation)
QC_GENO=$(grep "geno:" "${CONFIG}" | head -1 | awk '{print $2}')
QC_MIND=$(grep "mind:" "${CONFIG}" | head -1 | awk '{print $2}')
QC_MAF=$(grep "maf:" "${CONFIG}" | head -1 | awk '{print $2}')
QC_HWE=$(grep "hwe:" "${CONFIG}" | head -1 | awk '{print $2}')
QC_KING=$(grep "king_cutoff:" "${CONFIG}" | head -1 | awk '{print $2}')

# LD pruning
LD_WINDOW=$(grep "window:" "${CONFIG}" | head -1 | awk '{print $2}')
LD_STEP=$(grep "step:" "${CONFIG}" | head -1 | awk '{print $2}')
LD_R2=$(grep "r2:" "${CONFIG}" | head -1 | awk '{print $2}')

# PCA
N_PCS=$(grep "n_pcs:" "${CONFIG}" | head -1 | awk '{print $2}')

mkdir -p "${LOG_DIR}"
