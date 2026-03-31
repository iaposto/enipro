#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(detectRUNS)
  library(yaml)
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  sub("^--file=", "", file_arg[1])
}

config_path <- local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    normalizePath(args[1], mustWork = TRUE)
  } else {
    script_path <- normalizePath(get_script_path(), mustWork = TRUE)
    normalizePath(file.path(dirname(script_path), "..", "config", "params.yaml"), mustWork = TRUE)
  }
})

cfg <- yaml.load_file(config_path)
run_dir <- Sys.getenv("ENIPRO_RUN_DIR")
if (identical(run_dir, "")) {
  stop("ENIPRO_RUN_DIR is not set. Run this script from scripts/08_roh.sh.")
}
run_dir <- normalizePath(run_dir, mustWork = TRUE)

ped_path <- file.path(run_dir, "roh", "input", "graega_roh_input.ped")
map_path <- file.path(run_dir, "roh", "input", "graega_roh_input.map")
fam_path <- file.path(run_dir, "roh", "input", "graega_roh_input.fam")
array_bim_path <- file.path(run_dir, "roh", "input_prep", "step0_autosomal", "graega_roh_autosomal.bim")
runs_dir <- file.path(run_dir, "roh", "runs")
stats_dir <- file.path(run_dir, "roh", "stats")
dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

roh_cfg <- cfg$roh
runs <- slidingRUNS.run(
  genotypeFile = ped_path,
  mapFile = map_path,
  windowSize = roh_cfg$window_size,
  threshold = roh_cfg$window_threshold,
  minSNP = roh_cfg$min_snps,
  ROHet = FALSE,
  maxOppWindow = roh_cfg$max_heterozygous_window,
  maxMissWindow = roh_cfg$max_missing_window,
  maxGap = roh_cfg$max_gap_bp,
  minLengthBps = roh_cfg$min_length_bp,
  minDensity = 1 / roh_cfg$min_density_bp_per_snp,
  maxOppRun = roh_cfg$max_heterozygous_run,
  maxMissRun = roh_cfg$max_missing_run
)

runs <- as.data.frame(runs, stringsAsFactors = FALSE)
colnames(runs) <- c("breed", "IID", "chromosome", "n_snp", "start_bp", "end_bp", "length_bp")
runs$length_mb <- runs$length_bp / 1e6
runs$length_class <- ifelse(
  runs$length_bp < 5e6,
  "short (1-5 Mb)",
  ifelse(runs$length_bp < 1e7, "medium (5-10 Mb)", "long (>10 Mb)")
)
runs <- runs[order(runs$breed, runs$IID, runs$chromosome, runs$start_bp), ]

fam <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("breed", "IID", "father", "mother", "sex", "pheno")
fam <- fam[, c("breed", "IID")]

array_bim <- read.table(array_bim_path, header = FALSE, stringsAsFactors = FALSE)
colnames(array_bim) <- c("chromosome", "snp_id", "cm", "bp", "a1", "a2")
autosomal_lengths <- aggregate(bp ~ chromosome, data = array_bim, FUN = max)
colnames(autosomal_lengths) <- c("chromosome", "autosomal_length_bp")
autosomal_lengths <- autosomal_lengths[order(as.integer(autosomal_lengths$chromosome)), ]
laut_bp <- sum(autosomal_lengths$autosomal_length_bp)

if (nrow(runs) > 0) {
  roh_totals <- aggregate(length_bp ~ breed + IID, data = runs, FUN = sum)
  colnames(roh_totals)[3] <- "total_roh_bp"
} else {
  roh_totals <- data.frame(
    breed = character(),
    IID = character(),
    total_roh_bp = numeric(),
    stringsAsFactors = FALSE
  )
}
froh <- merge(fam, roh_totals, by = c("breed", "IID"), all.x = TRUE)
froh$total_roh_bp[is.na(froh$total_roh_bp)] <- 0
froh$autosomal_genome_bp <- laut_bp
froh$froh <- froh$total_roh_bp / laut_bp
froh <- froh[order(froh$breed, froh$IID), ]

length_class_levels <- c("short (1-5 Mb)", "medium (5-10 Mb)", "long (>10 Mb)")
class_grid <- merge(
  fam,
  data.frame(length_class = length_class_levels, stringsAsFactors = FALSE),
  by = NULL
)
if (nrow(runs) > 0) {
  class_counts <- aggregate(length_bp ~ breed + IID + length_class, data = runs, FUN = length)
  colnames(class_counts)[4] <- "n_roh"
  class_lengths <- aggregate(length_bp ~ breed + IID + length_class, data = runs, FUN = sum)
  colnames(class_lengths)[4] <- "total_roh_bp"
  class_summary <- merge(class_grid, class_counts, by = c("breed", "IID", "length_class"), all.x = TRUE)
  class_summary <- merge(class_summary, class_lengths, by = c("breed", "IID", "length_class"), all.x = TRUE)
  class_summary$n_roh[is.na(class_summary$n_roh)] <- 0
  class_summary$total_roh_bp[is.na(class_summary$total_roh_bp)] <- 0
} else {
  class_summary <- class_grid
  class_summary$n_roh <- 0
  class_summary$total_roh_bp <- 0
}
class_summary <- class_summary[order(class_summary$breed, class_summary$IID, class_summary$length_class), ]

write.table(runs, file.path(runs_dir, "roh_runs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(froh, file.path(stats_dir, "froh_per_individual.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(autosomal_lengths, file.path(stats_dir, "autosomal_genome_length_by_chr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(class_summary, file.path(stats_dir, "roh_length_class_per_individual.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_lines <- c(
  "detectRUNS ROH summary",
  sprintf("Runs detected: %d", nrow(runs)),
  sprintf("Samples in ROH input: %d", nrow(fam)),
  sprintf("Autosomal genome length used for FROH (bp): %.0f", laut_bp),
  sprintf("Mean FROH: %.6f", mean(froh$froh)),
  sprintf("Median FROH: %.6f", median(froh$froh)),
  sprintf("Max FROH: %.6f", max(froh$froh))
)
writeLines(summary_lines, file.path(stats_dir, "roh_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
