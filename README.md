# enipro

Bioinformatics pipeline for sheep/goat SNP genotyping analysis, supporting genetic diversity assessment, population structure characterization.

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
| **plink** (v1.9) | Heterozygosity and inbreeding coefficient (`--het`) | # Issues with plink2 when using --het.
| **Python 3.11** | Data processing, visualization, reporting |
| pandas, numpy, matplotlib, seaborn, scipy, pyyaml | Python analysis stack |

**Note:** The install command requires `--channel-priority flexible` to resolve cross-channel dependencies between bioconda (plink) and conda-forge (Python).

### Tool choices

- **PLINK 2.0:** PLINK 2.0 handles all QC operations and has a purpose-built `--pca` command.
- **PLINK2 `--pca`:** PLINK2 `--pca` computes eigenvectors/eigenvalues.
- **R not used in Phase 1.** For now, all visualization is handled by Python (matplotlib + seaborn). R may be introduced in later phases for specialized packages (e.g., `detectRUNS` for ROH analysis).

## Pipeline

The pipeline is organized into numbered shell scripts. Each script reads parameters from `config/params.yaml` via a shared `scripts/_common.sh` module.

### Running

Full pipeline:

```bash
micromamba run -n enipro bash scripts/run_all.sh
```

Or via Make:

```bash
micromamba run -n enipro make all
```

Individual steps can be run independently:

```bash
micromamba run -n enipro bash scripts/01_pre_qc_filter.sh
```

### Step 0 — Pre-QC filter (`01_pre_qc_filter.sh`)

Filters the raw dataset to retain only autosomal, biallelic SNPs with standard ACGT alleles.

```
plink2 --bfile <input> --chr-set 29 --chr 1-29 --snps-only just-acgt --max-alleles 2 --make-bed
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
plink2 --chr-set 29 --geno 0.10 --maf 0.01 --hwe 1e-6 --make-bed
```

| Filter | Threshold | Rationale |
|---|---|---|
| `--geno 0.10` | Remove SNPs with >10% missing genotypes | Standard call rate threshold (90%) |
| `--maf 0.01` | Remove SNPs with minor allele frequency <1% | Low-MAF variants add noise to PCA and inflate LD estimates. |
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
  maf: 0.01        # min minor allele frequency
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
│   ├── 06_pca.sh                     # Steps 5–6: LD prune + PCA
│   └── run_all.sh                    # Master runner
├── src/
│   ├── config.py                     # YAML config loader
│   ├── utils.py                      # PLINK output file parsers
│   ├── qc_report.py                  # QC summary report + diagnostic plots
│   └── plot_pca.py                   # PCA visualization
├── Makefile                          # Workflow orchestration
└── .gitignore

projects/enipro/                      # Data and results (NOT in git)
├── test_data/
│   └── graega_top_alleles/           # Raw input (never modified)
│       ├── graega_top_alleles.bed
│       ├── graega_top_alleles.bim
│       └── graega_top_alleles.fam
├── results/
│   ├── qc/
│   │   ├── step0_autosomal/          # After pre-QC filter
│   │   ├── step1_snp_qc/            # After SNP QC
│   │   ├── step2_sample_qc/         # After sample QC
│   │   ├── step3_relatedness/        # KING outputs + LD prune lists
│   │   ├── final/                    # Final QC'd dataset
│   │   └── stats/                    # Missingness, het, freq files
│   ├── pca/
│   │   ├── pruned_snps/              # LD-pruned SNP lists
│   │   └── eigenvec_eigenval/        # PCA results
│   └── figures/
│       ├── qc/                       # QC diagnostic plots
│       └── pca/                      # PCA scatter and scree plots
└── logs/
```


