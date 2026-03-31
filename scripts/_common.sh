#!/usr/bin/env bash
# Shared variables for all pipeline scripts.
# Sourced by each step script — not run directly.

# Load paths from params.yaml using a simple grep/sed approach
# (avoids Python dependency in shell scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../config/params.yaml"

_yaml_val() {
    grep "^${1}:" "${CONFIG}" | sed "s/^${1}:[[:space:]]*//" | sed 's/[[:space:]]*#.*//'
}

INPUT_PREFIX=$(_yaml_val "input_prefix")
RESULTS_DIR=$(_yaml_val "results_dir")
LOG_DIR=$(_yaml_val "log_dir")
CHR_SET=$(_yaml_val "chr_set")
AUTOSOMES=$(_yaml_val "autosomes" | tr -d '"')
NO_XY=$(_yaml_val "no_xy")

# Build --chr-set argument array (supports optional no-xy)
CHR_SET_ARGS=("${CHR_SET}")
[[ "${NO_XY}" == "true" ]] && CHR_SET_ARGS+=("no-xy")

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

# --- Timestamped run directory ---
# If RUN_STAMP was not exported by run_all.sh (i.e., individual script run),
# generate a new stamp so each standalone run also gets its own directory.
if [[ -z "${RUN_STAMP:-}" ]]; then
    export RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
fi

RUN_DIR="${RESULTS_DIR}/${RUN_STAMP}"
mkdir -p "${RUN_DIR}"

# Copy the config used for this run (idempotent — same file on repeated calls)
cp "${CONFIG}" "${RUN_DIR}/params.yaml"

# Export run dir for Python scripts
export ENIPRO_RUN_DIR="${RUN_DIR}"

mkdir -p "${LOG_DIR}"
