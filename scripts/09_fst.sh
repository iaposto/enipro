#!/usr/bin/env bash
# Step 9: Pairwise FST on the LD-pruned dataset, plus matrix plotting.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${SCRIPT_DIR}/_common.sh"

INDIR="${RUN_DIR}/pca/pruned_snps"
OUTDIR="${RUN_DIR}/fst"
INPUT_PREFIX="${INDIR}/graega_ldpruned"
OUTPUT_PREFIX="${OUTDIR}/graega_ldpruned_fst_wc"
mkdir -p "${OUTDIR}"

if [[ ! -f "${INPUT_PREFIX}.bed" ]]; then
    echo "ERROR: ${INPUT_PREFIX}.bed not found." >&2
    echo "Run scripts/06_pca.sh first, or export RUN_STAMP for an existing run directory." >&2
    exit 1
fi

echo "=== Step 9: Pairwise FST on LD-pruned data ==="

plink2 \
    --bfile "${INPUT_PREFIX}" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --family BREED \
    --fst BREED method=wc blocksize="${FST_BLOCKSIZE}" \
    --out "${OUTPUT_PREFIX}"

python "${REPO_DIR}/src/plot_fst.py" "${REPO_DIR}/config/params.yaml"

echo "FST summary: ${OUTPUT_PREFIX}.fst.summary"
echo "FST matrix:  ${OUTPUT_PREFIX}.matrix.tsv"
echo "Figures:     ${RUN_DIR}/figures/fst/"
