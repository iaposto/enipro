#!/usr/bin/env bash
# Master runner: execute the full QC + PCA pipeline.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "========================================="
echo " enipro Phase 1: QC + PCA Pipeline"
echo "========================================="

bash "${SCRIPT_DIR}/01_pre_qc_filter.sh"
echo ""
bash "${SCRIPT_DIR}/02_snp_qc.sh"
echo ""
bash "${SCRIPT_DIR}/03_sample_qc.sh"
echo ""
bash "${SCRIPT_DIR}/04_relatedness.sh"
echo ""
bash "${SCRIPT_DIR}/05_qc_summary.sh"
echo ""
bash "${SCRIPT_DIR}/06_pca.sh"
echo ""

echo "========================================="
echo " Pipeline complete. Generating reports..."
echo "========================================="

python "${REPO_DIR}/src/qc_report.py"
python "${REPO_DIR}/src/plot_pca.py"

echo "All done."
