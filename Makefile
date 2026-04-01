# enipro Phase 1: Goat SNP QC + PCA Workflow
SHELL := /bin/bash
SCRIPTS := scripts
ENV := micromamba run -n enipro

.PHONY: all qc pca figures fst roh clean help

help:
	@echo "Targets:"
	@echo "  all      — Run full pipeline (QC + PCA + figures)"
	@echo "  qc       — Run QC steps 0–4"
	@echo "  pca      — Run PCA (requires qc)"
	@echo "  figures  — Generate all plots (requires qc + pca)"
	@echo "  fst      — Run pairwise FST on an existing run with PCA outputs"
	@echo "  roh      — Run ROH detection + plots"
	@echo "  clean    — Print instructions for cleaning outputs"

all: qc pca figures

qc:
	$(ENV) python $(SCRIPTS)/01_pre_qc_filter.py
	$(ENV) python $(SCRIPTS)/02_snp_qc.py
	$(ENV) python $(SCRIPTS)/03_sample_qc.py
	$(ENV) python $(SCRIPTS)/04_relatedness.py
	$(ENV) python $(SCRIPTS)/05_qc_summary.py

pca: qc
	$(ENV) python $(SCRIPTS)/06_pca.py

figures: pca
	$(ENV) python src/qc_report.py
	$(ENV) python src/plot_pca.py

fst:
	$(ENV) python $(SCRIPTS)/09_fst.py

roh:
	$(ENV) python $(SCRIPTS)/08_roh.py

clean:
	@echo "To remove all results, run:"
	@echo "  rm -rf /home/i/iapostof/projects/enipro/results /home/i/iapostof/projects/enipro/logs"
