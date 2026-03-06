#!/usr/bin/env bash
# Step 3: Relatedness check — LD prune then KING kinship filter.
# Removes one individual per pair closer than 1st-degree relatives.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/_common.sh"

INDIR="${RESULTS_DIR}/qc/step2_sample_qc"
OUTDIR="${RESULTS_DIR}/qc/step3_relatedness"
FINALDIR="${RESULTS_DIR}/qc/final"
mkdir -p "${OUTDIR}" "${FINALDIR}"

echo "=== Step 3: Relatedness check ==="

# LD pruning for relatedness estimation
plink2 \
    --bfile "${INDIR}/graega_sampleqc" \
    --chr-set "${CHR_SET}" \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${OUTDIR}/ld_prune"

N_PRUNED_IN=$(wc -l < "${OUTDIR}/ld_prune.prune.in")
echo "LD-pruned SNP set: ${N_PRUNED_IN} variants"

# KING relatedness with automatic removal
plink2 \
    --bfile "${INDIR}/graega_sampleqc" \
    --chr-set "${CHR_SET}" \
    --extract "${OUTDIR}/ld_prune.prune.in" \
    --king-cutoff "${QC_KING}" \
    --out "${OUTDIR}/king"

# Check how many samples were flagged for removal
N_REMOVED=0
if [ -f "${OUTDIR}/king.king.cutoff.out.id" ]; then
    # Subtract 1 for header line (#FID IID)
    N_REMOVED=$(tail -n +2 "${OUTDIR}/king.king.cutoff.out.id" | wc -l)
fi

if [ "${N_REMOVED}" -gt 0 ]; then
    echo "Removing ${N_REMOVED} related sample(s)"
    plink2 \
        --bfile "${INDIR}/graega_sampleqc" \
        --chr-set "${CHR_SET}" \
        --remove "${OUTDIR}/king.king.cutoff.out.id" \
        --make-bed \
        --out "${FINALDIR}/graega_qc"
else
    echo "No related samples found — copying to final"
    cp "${INDIR}/graega_sampleqc.bed" "${FINALDIR}/graega_qc.bed"
    cp "${INDIR}/graega_sampleqc.bim" "${FINALDIR}/graega_qc.bim"
    cp "${INDIR}/graega_sampleqc.fam" "${FINALDIR}/graega_qc.fam"
fi

N_VARIANTS=$(wc -l < "${FINALDIR}/graega_qc.bim")
N_SAMPLES=$(wc -l < "${FINALDIR}/graega_qc.fam")
echo "Final QC dataset: ${N_VARIANTS} variants, ${N_SAMPLES} samples"
