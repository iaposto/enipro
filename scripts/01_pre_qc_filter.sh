#!/usr/bin/env bash
# Step 0: Filter to autosomal biallelic SNPs only (chr 1-29, ACGT alleles).
#
# Removes: sex chromosomes, mitochondrial, unplaced, scaffold-mapped markers,
# indels, deletions, and monomorphic sites with missing allele codes.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/_common.sh"

OUTDIR="${RUN_DIR}/qc/step0_autosomal"
mkdir -p "${OUTDIR}"

echo "=== Step 0: Pre-QC filter — autosomal biallelic SNPs ==="
echo "Input: ${INPUT_PREFIX}"

plink2 \
    --bfile "${INPUT_PREFIX}" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --chr "${AUTOSOMES}" \
    --snps-only just-acgt \
    --max-alleles 2 \
    --make-bed \
    --out "${OUTDIR}/graega_autosomal"

N_VARIANTS=$(wc -l < "${OUTDIR}/graega_autosomal.bim")
N_SAMPLES=$(wc -l < "${OUTDIR}/graega_autosomal.fam")
echo "Output: ${N_VARIANTS} variants, ${N_SAMPLES} samples"
