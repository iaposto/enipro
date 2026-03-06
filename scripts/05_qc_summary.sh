#!/usr/bin/env bash
# Step 4: Generate QC summary statistics on the final dataset.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/_common.sh"

INDIR="${RESULTS_DIR}/qc/final"
OUTDIR="${RESULTS_DIR}/qc/stats"
mkdir -p "${OUTDIR}"

echo "=== Step 4: QC summary statistics ==="

# Per-sample and per-variant missingness
plink2 \
    --bfile "${INDIR}/graega_qc" \
    --chr-set "${CHR_SET}" \
    --missing \
    --out "${OUTDIR}/missing"

# Heterozygosity and inbreeding coefficient (PLINK 1.9)
plink \
    --bfile "${INDIR}/graega_qc" \
    --chr-set "${CHR_SET}" \
    --het \
    --out "${OUTDIR}/het"

# Allele frequencies
plink2 \
    --bfile "${INDIR}/graega_qc" \
    --chr-set "${CHR_SET}" \
    --freq \
    --out "${OUTDIR}/freq"

echo "Summary statistics written to ${OUTDIR}/"
