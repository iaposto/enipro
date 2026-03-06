#!/usr/bin/env bash
# Step 2: Sample-level QC — per-individual call rate filter.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/_common.sh"

INDIR="${RESULTS_DIR}/qc/step1_snp_qc"
OUTDIR="${RESULTS_DIR}/qc/step2_sample_qc"
mkdir -p "${OUTDIR}"

echo "=== Step 2: Sample-level QC ==="
echo "Threshold: --mind ${QC_MIND}"

plink2 \
    --bfile "${INDIR}/graega_snpqc" \
    --chr-set "${CHR_SET}" \
    --mind "${QC_MIND}" \
    --make-bed \
    --out "${OUTDIR}/graega_sampleqc"

N_VARIANTS=$(wc -l < "${OUTDIR}/graega_sampleqc.bim")
N_SAMPLES=$(wc -l < "${OUTDIR}/graega_sampleqc.fam")
echo "Output: ${N_VARIANTS} variants, ${N_SAMPLES} samples"
