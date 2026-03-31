# enipro

Bioinformatics pipeline for sheep/goat SNP genotyping analysis, supporting genetic diversity assessment, population structure characterization, and ROH-based genomic inbreeding analysis.

## Dataset

The test dataset (GRAEGA) consists of goat SNP genotypes from 10 breeds (ANG, ARI, BHA, CRO, DAM, IND, MUR, PAG, SER, SKO), genotyped on the Axiom Caprine Genotyping Array.
Data is in PLINK binary format (`.bed`/`.bim`/`.fam`) and mapped to the goat reference genome ARS1/CHI\_v2 (29 autosomes).

Raw data lives in the project directory and is never modified by the pipeline:

```
/home/i/iapostof/projects/enipro/test_data/graega_top_alleles/
```

## Environment

All tools and dependencies are managed in a micromamba environment named `enipro`.

### Setup

```bash
bash scripts/00_setup_env.sh
```

This installs:

| Tool | Purpose |
|---|---|
| **plink2** (v2.0) | Primary tool for QC filtering, LD pruning, PCA, allele frequencies |
| **plink** (v1.9) | Heterozygosity and inbreeding coefficient (`--het`); used because `plink2 --het` was unreliable here |
| **Python 3.11** | Data processing, visualization, reporting |
| **R** | ROH detection with `detectRUNS` |
| **detectRUNS** | Sliding-window ROH calling and ROH-derived summaries |
| pandas, numpy, matplotlib, seaborn, scipy, pyyaml | Python analysis stack |

**Note:** The install command requires `--channel-priority flexible` to resolve cross-channel dependencies between bioconda (plink) and conda-forge (Python/R). `detectRUNS` is installed from CRAN into the same `enipro` environment after the conda packages are resolved.

### Tool choices

- **PLINK 2.0:** PLINK 2.0 handles all QC operations and has a purpose-built `--pca` command.
- **PLINK2 `--pca`:** PLINK2 `--pca` computes eigenvectors/eigenvalues.
- **R + detectRUNS:** ROH detection uses the CRAN `detectRUNS` package because it provides a well-established sliding-window implementation and direct ROH summaries from PLINK PED/MAP input.
- **Python visualization:** PCA, QC, and ROH figures are still generated in Python so the plotting style remains consistent across the project.

## Pipeline

The pipeline is organized into numbered shell scripts. Each script reads parameters from `config/params.yaml` via a shared `scripts/_common.sh` module.

### Run directories

Each pipeline run writes its output to a timestamped subdirectory:

```
results/<YYYYMMDD_HHMMSS>/
├── params.yaml          # copy of config used for this run
├── qc/
│   ├── step0_autosomal/
│   ├── step1_snp_qc/
│   ├── step2_sample_qc/
│   ├── step3_relatedness/
│   ├── final/
│   └── stats/
├── pca/
│   ├── pruned_snps/     # LD-pruned .bed/.bim/.fam also saved here for ADMIXTURE input
│   └── eigenvec_eigenval/
├── admixture/           # ADMIXTURE outputs (copied from HPC after run)
├── roh/
│   ├── input_prep/      # ROH-specific QC path (no MAF filter; LD only for KING relatedness)
│   ├── input/           # Final ROH PLINK + PED/MAP input
│   ├── runs/            # ROH segment calls
│   └── stats/           # FROH tables and summary files
└── figures/
    ├── qc/
    ├── pca/
    ├── admixture/
    └── roh/
```

When `run_all.sh` is invoked it generates a single `RUN_STAMP` shared by all steps and the Python visualization scripts. When an individual step script is run directly, it generates its own stamp. The `params.yaml` snapshot makes every run directory self-contained and reproducible.

### Running

Full pipeline:

```bash
micromamba run -n enipro bash scripts/run_all.sh
```

Or via Make:

```bash
micromamba run -n enipro make all
```

ROH analysis:

```bash
micromamba run -n enipro make roh
```

Individual steps can be run independently:

```bash
micromamba run -n enipro bash scripts/01_pre_qc_filter.sh
```

The ROH workflow is self-contained and can also be run directly:

```bash
micromamba run -n enipro bash scripts/08_roh.sh
```

### Step 0 — Pre-QC filter (`01_pre_qc_filter.sh`)

Filters the raw dataset to retain only autosomal, biallelic SNPs with standard ACGT alleles.

```
plink2 --bfile <input> --chr-set 29 no-xy --chr 1-29 --snps-only just-acgt --max-alleles 2 --make-bed
```

**What is removed:**

| Category | Reason |
|---|---|
| Sex chromosomes (chr 30/31) | All samples have sex coded as unknown (0). |
| Mitochondrial markers (chr 33) | |
| Unplaced SNPs (chr 0) | No genomic position. |
| Scaffold-mapped SNPs (NW\_\*) | On unlocalized scaffolds not assigned to chromosomes |
| Indels, deletions, non-ACGT alleles | Not standard biallelic SNPs; the `--snps-only just-acgt` flag removes multi-character alleles (e.g., TAC), deletion codes (D), and missing allele codes (0) |

**Output:** `results/qc/step0_autosomal/graega_autosomal.{bed,bim,fam}`

### Step 1 — SNP-level QC (`02_snp_qc.sh`)

Applies per-variant quality filters.

```
plink2 --chr-set 29 no-xy --geno 0.10 --maf 0.02 --hwe 1e-6 --make-bed
```

| Filter | Threshold | Rationale |
|---|---|---|
| `--geno 0.10` | Remove SNPs with >10% missing genotypes | Standard call rate threshold (90%) |
| `--maf 0.02` | Remove SNPs with minor allele frequency <2% | Low-MAF variants add noise to PCA and inflate LD estimates. |
| `--hwe 1e-6` | Remove SNPs with Hardy-Weinberg exact test p < 1e-6 | Lenient threshold. |

**Output:** `results/qc/step1_snp_qc/graega_snpqc.{bed,bim,fam}`

### Step 2 — Sample-level QC (`03_sample_qc.sh`)

Applies per-individual call rate filter.

```
plink2 --chr-set 29 --mind 0.10 --make-bed
```

Removes individuals with >10% missing genotypes.

**Output:** `results/qc/step2_sample_qc/graega_sampleqc.{bed,bim,fam}`

### Step 3 — Relatedness check (`04_relatedness.sh`)

Identifies and removes closely related individuals that would bias PCA and downstream analyses.

```
plink2 --indep-pairwise 50 5 0.2 --out ld_prune        # LD pruning
plink2 --extract ld_prune.prune.in --king-cutoff 0.177  # KING kinship
```

**LD pruning** is performed first (window=50 SNPs, step=5, r²=0.2).

**KING kinship** is used instead of PLINK's IBD-based `--genome` because KING is robust to population structure. "In a multi-breed dataset, IBD estimates are inflated by shared ancestry between breeds, leading to false positives. KING correctly handles this."

| KING kinship | Relationship |
|---|---|
| > 0.354 | Duplicate / MZ twin |
| > 0.177 | 1st-degree (parent-offspring, full sibling) |
| > 0.0884 | 2nd-degree (half-sibling, grandparent) |
| > 0.0442 | 3rd-degree |

The threshold is set at **0.177** (1st-degree relatives). PLINK2's `--king-cutoff` automatically selects which individual to remove from each pair (minimizing total removals) and writes two files:

- `.king.cutoff.in.id` — samples to keep
- `.king.cutoff.out.id` — samples to remove

If any samples are flagged, they are removed with `--remove`. If none are flagged, the dataset is copied to the final directory unchanged.

**Output:** `results/qc/final/graega_qc.{bed,bim,fam}` — the final QC'd dataset used for all downstream analyses.

### Step 4 — QC summary statistics (`05_qc_summary.sh`)

Generates descriptive statistics on the final QC'd dataset.

```
plink2 --missing       # per-sample (.smiss) and per-variant (.vmiss) missingness
plink  --het           # observed/expected homozygosity, inbreeding coefficient F # plink2 --het was not working properly so plink was used
plink2 --freq          # allele frequencies (.afreq)
```

The inbreeding coefficient is calculated as:

```
F = (E(HOM) - O(HOM)) / E(HOM)
```

where E(HOM) is the expected number of homozygous genotypes under Hardy-Weinberg equilibrium and O(HOM) is the observed count. Positive *F* indicates excess homozygosity (inbreeding); negative *F* indicates excess heterozygosity (possible admixture or outbreeding).

**Output:** `results/qc/stats/{missing.smiss, missing.vmiss, het.het, freq.afreq}`

### Step 5–6 — PCA (`06_pca.sh`)

LD pruning followed by principal component analysis.

```
plink2 --indep-pairwise 50 5 0.2 --out ld_prune         # Step 5: LD pruning
plink2 --extract ld_prune.prune.in --pca 20 --out pca    # Step 6: PCA
```

**LD pruning for PCA** chromosomal regions with extended LD can inflate/dominate the principal components, distorting the population structure. The same standard parameters are used (50/5/0.2).

**20 PCs are computed.** 

PLINK2 `--pca` is used.

**Output:**

| File | Content |
|---|---|
| `results/pca/pruned_snps/ld_prune.prune.in` | SNP IDs retained after LD pruning |
| `results/pca/pruned_snps/ld_prune.prune.out` | SNP IDs removed by LD pruning |
| `results/pca/eigenvec_eigenval/graega_pca.eigenvec` | Tab-separated file: FID, IID, PC1–PC20 |
| `results/pca/eigenvec_eigenval/graega_pca.eigenval` | One eigenvalue per line (20 values) |

### Step 7 — ADMIXTURE (`07_admixture.sh`)

Model-based ancestry estimation run on an **external HPC cluster** (not part of `run_all.sh`). The input is the LD-pruned binary dataset saved by Step 5b.

**Input:** `results/pca/pruned_snps/graega_ldpruned.{bed,bim,fam}` — copy these to the HPC working directory alongside `07_admixture.sh`.

```bash
# Submitted via SLURM
sbatch scripts/07_admixture.sh
```

The script runs K=2–30 sequentially, each with 10-fold cross-validation:

```bash
admixture --cv=10 -j16 -s 42 graega_ldpruned.bed $K | tee log${K}.out
```

| Parameter | Value | Rationale |
|---|---|---|
| `--cv=10` | 10-fold cross-validation | More precise CV error estimate than default 5-fold |
| `-j16` | 16 threads | Conservative on shared HPC node (64 CPUs available) |
| `-s 42` | Fixed random seed | Reproducibility |
| K=2–30 | Range covers all plausible structures | CV curve identifies the optimal K |

**Optimal K:** determined by the lowest CV error. For the GRAEGA dataset K=7 (CV=0.628), reflecting that the 10 named breeds cluster into 7 genetically distinct ancestral groups (MUR and ANG are each highly distinct; ARI splits into two sub-clusters; BHA/IND/PAG/SER/CRO seem to share a common Greek mainland component; SKO has its own island-derived component).

**Outputs** (copy back from HPC to `results/<RUN_STAMP>/admixture/`):

| File | Content |
|---|---|
| `graega_ldpruned.{K}.Q` | Ancestry proportions per individual (267 rows × K columns) |
| `graega_ldpruned.{K}.P` | Ancestral allele frequencies per SNP (45,192 rows × K columns) |
| `log{K}.out` | ADMIXTURE log including CV error for each K |
| `provenance.txt` | Records which pipeline run's pruned files were used as input |

### Step 8 — ROH detection (`08_roh.sh`)

Runs of homozygosity are called with the CRAN `detectRUNS` package using the **sliding-window** method.

The ROH input follows the same overall QC logic as the main pipeline, but **skips MAF filtering and skips LD pruning in the final ROH dataset**. The ROH-specific preparation path is:

1. autosomal, biallelic, ACGT SNP filter
2. SNP QC with `--geno 0.10` and `--hwe 1e-6` only
3. sample QC with `--mind 0.10`
4. KING relatedness filtering, where LD pruning is used only to estimate kinship and choose samples to remove
5. export the final unpruned ROH dataset to PLINK PED/MAP for `detectRUNS`

This preserves the requested SNP set for ROH calling while still applying the non-MAF QC filters and relatedness filter.

The sliding-window ROH parameters are:

| Parameter | Value |
|---|---|
| Minimum ROH length | 1 Mb |
| Minimum SNPs per ROH | 50 |
| Minimum density | 1 SNP / 70 kb |
| Window size | 50 SNPs |
| Maximum gap between consecutive SNPs | 100 kb |
| Maximum missing SNPs per window / per ROH | 5 |
| Maximum heterozygous SNPs per window / per ROH | 1 |
| Window threshold | 0.05 |

Equivalent `detectRUNS::slidingRUNS.run()` settings:

```r
windowSize   = 50
threshold    = 0.05
minSNP       = 50
maxGap       = 100000
minLengthBps = 1000000
minDensity   = 1 / 70000
maxMissWindow = 5
maxMissRun    = 5
maxOppWindow  = 1
maxOppRun     = 1
```

**FROH** is calculated as:

```text
FROH = ΣLROH / Laut
```

where `ΣLROH` is the total ROH length per individual and `Laut` is the autosomal array length estimated from the maximum mapped position on each autosome in the autosomal step-0 SNP map. This keeps the denominator tied to the genotyping array rather than to the ROH calls themselves.

**Outputs:**

| File | Content |
|---|---|
| `results/<RUN_STAMP>/roh/input/graega_roh_input.{bed,bim,fam,ped,map}` | Final ROH input dataset after non-MAF QC + relatedness filtering |
| `results/<RUN_STAMP>/roh/runs/roh_runs.tsv` | One row per detected ROH segment |
| `results/<RUN_STAMP>/roh/stats/froh_per_individual.tsv` | Per-individual total ROH length and `FROH` |
| `results/<RUN_STAMP>/roh/stats/autosomal_genome_length_by_chr.tsv` | Per-chromosome contribution to `Laut` |
| `results/<RUN_STAMP>/roh/stats/roh_length_class_per_individual.tsv` | ROH counts and total length in the three length classes |
| `results/<RUN_STAMP>/roh/stats/roh_summary.txt` | Text summary of run counts and `FROH` |

### QC report (`src/qc_report.py`)

Parses all QC statistics and produces a text summary and diagnostic plots:

```bash
micromamba run -n enipro python src/qc_report.py
```

**Figures generated** (in `results/figures/qc/`):

| File | Description |
|---|---|
| `maf_distribution.{png,pdf}` | Histogram of minor allele frequencies after QC |
| `sample_missingness.{png,pdf}` | Distribution of per-sample missing rates |
| `snp_missingness.{png,pdf}` | Distribution of per-SNP missing rates |
| `het_by_breed.{png,pdf}` | Boxplot of observed heterozygosity per breed |
| `inbreeding_by_breed.{png,pdf}` | Boxplot of inbreeding coefficient *F* per breed |

The text report prints() variant/sample counts at each QC step, per-breed sample sizes, and per-breed heterozygosity summaries.

### ADMIXTURE visualization (`src/plot_admixture.py`)

Produces CV error and Q-matrix plots from the ADMIXTURE output:

```bash
ENIPRO_RUN_DIR=/home/i/iapostof/projects/enipro/results/<RUN_STAMP> \
  micromamba run -n enipro python src/plot_admixture.py
```

**Figures generated** (in `results/figures/admixture/`):

| File | Description |
|---|---|
| `admixture_cv_error.png` | CV error vs K (K=2–30), red dashed line at optimal K |
| `admixture_K{K}.png` | Stacked bar chart of ancestry proportions for one K value, samples sorted by breed |
| `admixture_panel.png` | Multi-row panel showing K=5,6,7,8,10 side by side for comparison |

### PCA visualization (`src/plot_pca.py`)

Produces PCA plots from the eigenvector/eigenvalue output:

```bash
micromamba run -n enipro python src/plot_pca.py
```

**Figures generated** (in `results/figures/pca/`):

| File | Description |
|---|---|
| `scree_plot.{png,pdf}` | Bar chart of variance explained (%) per PC, with annotations on the top PCs |
| `pca_pc1_pc2.{png,pdf}` | PC1 vs PC2 scatter plot, colored by breed |
| `pca_pc1_pc3.{png,pdf}` | PC1 vs PC3 scatter plot |
| `pca_pc2_pc3.{png,pdf}` | PC2 vs PC3 scatter plot |
| `pca_pairplot_4pc.{png,pdf}` | Pairwise scatter of PC1–PC4 with KDE on the diagonal |

All plots use a consistent breed color palette defined in `config/params.yaml`. Each scatter plot shows axis labels with the proportion of variance explained by that PC.

### ROH execution (`src/run_roh.R`)

Runs `detectRUNS::slidingRUNS.run()` on the ROH-specific PED/MAP files and writes:

- the full ROH table (`roh_runs.tsv`)
- per-individual `FROH`
- per-chromosome autosomal lengths used for `Laut`
- per-individual ROH counts in the short/medium/long categories

Example:

```bash
ENIPRO_RUN_DIR=/home/i/iapostof/projects/enipro/results/<RUN_STAMP> \
  micromamba run -n enipro Rscript src/run_roh.R config/params.yaml
```

### ROH visualization (`src/plot_roh.py`)

Produces ROH summary figures from the `detectRUNS` outputs:

```bash
ENIPRO_RUN_DIR=/home/i/iapostof/projects/enipro/results/<RUN_STAMP> \
  micromamba run -n enipro python src/plot_roh.py config/params.yaml
```

**Figures generated** (in `results/figures/roh/`):

| File | Description |
|---|---|
| `roh_length_classes_by_breed.{png,pdf}` | Mean ROH count per individual, split into short (1-5 Mb), medium (5-10 Mb), and long (>10 Mb) classes |
| `froh_by_breed.{png,pdf}` | Breed-wise distribution of `FROH` |

## Configuration

All parameters are centralized in `config/params.yaml`:

```yaml
# Paths
input_prefix: /home/i/iapostof/projects/enipro/test_data/graega_top_alleles/graega_top_alleles
results_dir:  /home/i/iapostof/projects/enipro/results

# Species
chr_set: 29        # goat: 29 autosomes
autosomes: "1-29"

# QC thresholds
qc:
  geno: 0.10       # max per-SNP missing rate
  mind: 0.10       # max per-sample missing rate
  maf: 0.02        # min minor allele frequency
  hwe: 1e-6        # HWE p-value threshold
  king_cutoff: 0.177  # KING kinship threshold

# LD pruning
ld_prune:
  window: 50
  step: 5
  r2: 0.2

# PCA
pca:
  n_pcs: 20

# ROH
roh:
  min_length_bp: 1000000
  min_snps: 50
  min_density_bp_per_snp: 70000
  window_size: 50
  max_gap_bp: 100000
  max_missing_window: 5
  max_missing_run: 5
  max_heterozygous_window: 1
  max_heterozygous_run: 1
  window_threshold: 0.05

# Breed colors for visualization
breeds:
  ANG: { color: "#e41a1c" }
  ARI: { color: "#377eb8" }
  # ... (10 breeds total)
```

Shell scripts parse this file via `scripts/_common.sh`. Python scripts load it via `src/config.py`.

## Directory layout

```
repos/enipro/                         # Git-tracked code
├── config/
│   └── params.yaml                   # Central configuration
├── scripts/
│   ├── _common.sh                    # Shared variables (sourced by all scripts)
│   ├── 00_setup_env.sh               # Environment setup
│   ├── 01_pre_qc_filter.sh           # Step 0: autosomal biallelic filter
│   ├── 02_snp_qc.sh                  # Step 1: SNP-level QC
│   ├── 03_sample_qc.sh               # Step 2: sample-level QC
│   ├── 04_relatedness.sh             # Step 3: KING relatedness check
│   ├── 05_qc_summary.sh              # Step 4: summary statistics
│   ├── 06_pca.sh                     # Steps 5–6: LD prune + PCA (also saves pruned .bed for ADMIXTURE)
│   ├── 07_admixture.sh               # ADMIXTURE K=2–30 (SLURM script, run on HPC)
│   ├── 08_roh.sh                     # ROH-specific QC path + detectRUNS execution
│   └── run_all.sh                    # Master runner (Steps 0–6 + visualization)
├── src/
│   ├── config.py                     # YAML config loader
│   ├── utils.py                      # PLINK output file parsers
│   ├── qc_report.py                  # QC summary report + diagnostic plots
│   ├── plot_pca.py                   # PCA visualization
│   ├── plot_admixture.py             # ADMIXTURE CV error + Q-matrix plots
│   ├── run_roh.R                     # detectRUNS ROH calling + FROH tables
│   └── plot_roh.py                   # ROH length-class and FROH figures
├── Makefile                          # Workflow orchestration
└── .gitignore

projects/enipro/                      # Data and results (NOT in git)
├── test_data/
│   └── graega_top_alleles/           # Raw input (never modified)
│       ├── graega_top_alleles.bed
│       ├── graega_top_alleles.bim
│       └── graega_top_alleles.fam
├── results/
│   └── <RUN_STAMP>/
│       ├── qc/
│       │   ├── step0_autosomal/      # After pre-QC filter
│       │   ├── step1_snp_qc/         # After SNP QC
│       │   ├── step2_sample_qc/      # After sample QC
│       │   ├── step3_relatedness/    # KING outputs + LD prune lists
│       │   ├── final/                # Final QC'd dataset
│       │   └── stats/                # Missingness, het, freq files
│       ├── pca/
│       │   ├── pruned_snps/          # LD-pruned SNP lists + graega_ldpruned.{bed,bim,fam}
│       │   └── eigenvec_eigenval/    # PCA results
│       ├── admixture/                # ADMIXTURE outputs (copied from HPC)
│       │   ├── graega_ldpruned.{K}.Q # Ancestry fractions for each K
│       │   ├── graega_ldpruned.{K}.P # Ancestral allele frequencies for each K
│       │   ├── log{K}.out            # ADMIXTURE logs with CV errors
│       │   └── provenance.txt        # Records input run stamp
│       ├── roh/
│       │   ├── input_prep/           # ROH-only QC path without MAF filtering
│       │   ├── input/                # Final ROH input PLINK + PED/MAP files
│       │   ├── runs/                 # Detected ROH segments
│       │   └── stats/                # FROH and ROH length-class tables
│       └── figures/
│           ├── qc/                   # QC diagnostic plots
│           ├── pca/                  # PCA scatter and scree plots
│           ├── admixture/            # CV error curve and Q-matrix bar charts
│           └── roh/                  # ROH class and FROH figures
└── logs/
```
