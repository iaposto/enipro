#!/usr/bin/env bash
# Step 1: SNP-level QC — call rate, MAF, HWE filters.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/_common.sh"

INDIR="${RUN_DIR}/qc/step0_autosomal"
OUTDIR="${RUN_DIR}/qc/step1_snp_qc"
mkdir -p "${OUTDIR}"

echo "=== Step 1: SNP-level QC ==="
echo "Thresholds: --geno ${QC_GENO} --maf ${QC_MAF} --hwe ${QC_HWE}"

plink2 \
    --bfile "${INDIR}/graega_autosomal" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --geno "${QC_GENO}" \
    --maf "${QC_MAF}" \
    --hwe "${QC_HWE}" \
    --make-bed \
    --out "${OUTDIR}/graega_snpqc"

N_VARIANTS=$(wc -l < "${OUTDIR}/graega_snpqc.bim")
N_SAMPLES=$(wc -l < "${OUTDIR}/graega_snpqc.fam")
echo "Output: ${N_VARIANTS} variants, ${N_SAMPLES} samples"
