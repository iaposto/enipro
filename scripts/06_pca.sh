#!/usr/bin/env bash
# Steps 5-6: LD prune and compute PCA on the QC'd dataset.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/_common.sh"

INDIR="${RESULTS_DIR}/qc/final"
PRUNEDIR="${RESULTS_DIR}/pca/pruned_snps"
PCADIR="${RESULTS_DIR}/pca/eigenvec_eigenval"
mkdir -p "${PRUNEDIR}" "${PCADIR}"

echo "=== Step 5: LD pruning for PCA ==="

plink2 \
    --bfile "${INDIR}/graega_qc" \
    --chr-set "${CHR_SET}" \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${PRUNEDIR}/ld_prune"

N_PRUNED_IN=$(wc -l < "${PRUNEDIR}/ld_prune.prune.in")
echo "LD-pruned SNP set for PCA: ${N_PRUNED_IN} variants"

echo "=== Step 6: PCA computation ==="

plink2 \
    --bfile "${INDIR}/graega_qc" \
    --chr-set "${CHR_SET}" \
    --extract "${PRUNEDIR}/ld_prune.prune.in" \
    --pca "${N_PCS}" \
    --out "${PCADIR}/graega_pca"

echo "PCA complete: ${N_PCS} PCs computed"
echo "Eigenvectors: ${PCADIR}/graega_pca.eigenvec"
echo "Eigenvalues:  ${PCADIR}/graega_pca.eigenval"
