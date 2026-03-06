# enipro Phase 1: Goat SNP QC + PCA Workflow
SHELL := /bin/bash
SCRIPTS := scripts
ENV := micromamba run -n enipro

.PHONY: all qc pca figures clean help

help:
	@echo "Targets:"
	@echo "  all      — Run full pipeline (QC + PCA + figures)"
	@echo "  qc       — Run QC steps 0–4"
	@echo "  pca      — Run PCA (requires qc)"
	@echo "  figures  — Generate all plots (requires qc + pca)"
	@echo "  clean    — Print instructions for cleaning outputs"

all: qc pca figures

qc:
	$(ENV) bash $(SCRIPTS)/01_pre_qc_filter.sh
	$(ENV) bash $(SCRIPTS)/02_snp_qc.sh
	$(ENV) bash $(SCRIPTS)/03_sample_qc.sh
	$(ENV) bash $(SCRIPTS)/04_relatedness.sh
	$(ENV) bash $(SCRIPTS)/05_qc_summary.sh

pca: qc
	$(ENV) bash $(SCRIPTS)/06_pca.sh

figures: pca
	$(ENV) python src/qc_report.py
	$(ENV) python src/plot_pca.py

clean:
	@echo "To remove all results, run:"
	@echo "  rm -rf /home/i/iapostof/projects/enipro/results /home/i/iapostof/projects/enipro/logs"
