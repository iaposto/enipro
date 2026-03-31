#!/usr/bin/env bash
# Step 8: Prepare ROH input, run detectRUNS, and generate ROH figures.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${SCRIPT_DIR}/_common.sh"

PREPDIR="${RUN_DIR}/roh/input_prep"
STEP0DIR="${PREPDIR}/step0_autosomal"
STEP1DIR="${PREPDIR}/step1_snp_qc_nomaf"
STEP2DIR="${PREPDIR}/step2_sample_qc"
STEP3DIR="${PREPDIR}/step3_relatedness"
INPUTDIR="${RUN_DIR}/roh/input"

mkdir -p "${STEP0DIR}" "${STEP1DIR}" "${STEP2DIR}" "${STEP3DIR}" "${INPUTDIR}"

echo "========================================="
echo " enipro ROH pipeline"
echo " Run: ${RUN_STAMP}"
echo " Output: ${RUN_DIR}"
echo "========================================="

echo "=== ROH prep 0: autosomal biallelic SNPs ==="
plink2 \
    --bfile "${INPUT_PREFIX}" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --chr "${AUTOSOMES}" \
    --snps-only just-acgt \
    --max-alleles 2 \
    --make-bed \
    --out "${STEP0DIR}/graega_roh_autosomal"

echo "=== ROH prep 1: SNP QC without MAF filtering ==="
plink2 \
    --bfile "${STEP0DIR}/graega_roh_autosomal" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --geno "${QC_GENO}" \
    --hwe "${QC_HWE}" \
    --make-bed \
    --out "${STEP1DIR}/graega_roh_snpqc"

echo "=== ROH prep 2: sample QC ==="
plink2 \
    --bfile "${STEP1DIR}/graega_roh_snpqc" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --mind "${QC_MIND}" \
    --make-bed \
    --out "${STEP2DIR}/graega_roh_sampleqc"

echo "=== ROH prep 3: relatedness filtering ==="
plink2 \
    --bfile "${STEP2DIR}/graega_roh_sampleqc" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${STEP3DIR}/ld_prune"

plink2 \
    --bfile "${STEP2DIR}/graega_roh_sampleqc" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --extract "${STEP3DIR}/ld_prune.prune.in" \
    --king-cutoff "${QC_KING}" \
    --out "${STEP3DIR}/king"

N_REMOVED=0
if [ -f "${STEP3DIR}/king.king.cutoff.out.id" ]; then
    N_REMOVED=$(tail -n +2 "${STEP3DIR}/king.king.cutoff.out.id" | wc -l)
fi

if [ "${N_REMOVED}" -gt 0 ]; then
    echo "Removing ${N_REMOVED} related sample(s) from the unpruned ROH input"
    plink2 \
        --bfile "${STEP2DIR}/graega_roh_sampleqc" \
        --chr-set "${CHR_SET_ARGS[@]}" \
        --remove "${STEP3DIR}/king.king.cutoff.out.id" \
        --make-bed \
        --out "${INPUTDIR}/graega_roh_input"
else
    echo "No related samples found for ROH input"
    cp "${STEP2DIR}/graega_roh_sampleqc.bed" "${INPUTDIR}/graega_roh_input.bed"
    cp "${STEP2DIR}/graega_roh_sampleqc.bim" "${INPUTDIR}/graega_roh_input.bim"
    cp "${STEP2DIR}/graega_roh_sampleqc.fam" "${INPUTDIR}/graega_roh_input.fam"
fi

echo "=== ROH prep 4: PED/MAP export for detectRUNS ==="
plink \
    --bfile "${INPUTDIR}/graega_roh_input" \
    --chr-set "${CHR_SET_ARGS[@]}" \
    --recode \
    --out "${INPUTDIR}/graega_roh_input"

N_VARIANTS=$(wc -l < "${INPUTDIR}/graega_roh_input.bim")
N_SAMPLES=$(wc -l < "${INPUTDIR}/graega_roh_input.fam")
echo "ROH input dataset: ${N_VARIANTS} variants, ${N_SAMPLES} samples"

echo "=== Step 8a: detectRUNS ROH scan ==="
Rscript "${REPO_DIR}/src/run_roh.R" "${REPO_DIR}/config/params.yaml"

echo "=== Step 8b: ROH plots ==="
python "${REPO_DIR}/src/plot_roh.py" "${REPO_DIR}/config/params.yaml"

echo "ROH analysis complete."
