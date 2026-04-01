#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

SUMMARY_BASENAME <- "graega_ldpruned_fst_wc"
PRUNED_FAM_BASENAME <- "graega_ldpruned.fam"

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  sub("^--file=", "", file_arg[1])
}

get_config_path <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && nzchar(args[1])) {
    return(normalizePath(args[1], mustWork = TRUE))
  }

  script_path <- normalizePath(get_script_path(), mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), "..", "config", "params.yaml"), mustWork = TRUE)
}

resolve_run_dir <- function() {
  run_dir <- Sys.getenv("ENIPRO_RUN_DIR")
  if (identical(run_dir, "")) {
    stop("ENIPRO_RUN_DIR is not set. Run this script from scripts/09_fst.py or export it explicitly.")
  }
  normalizePath(run_dir, mustWork = TRUE)
}

get_breed_order <- function(cfg, populations) {
  cfg_order <- names(cfg$breeds)
  extras <- setdiff(sort(populations), cfg_order)
  c(cfg_order[cfg_order %in% populations], extras)
}

get_breed_labels <- function(cfg, breeds) {
  labels <- setNames(breeds, breeds)
  for (breed in breeds) {
    breed_cfg <- cfg$breeds[[breed]]
    if (!is.null(breed_cfg) && !is.null(breed_cfg$label) && nzchar(breed_cfg$label)) {
      labels[[breed]] <- breed_cfg$label
    }
  }
  labels
}

build_matrix <- function(summary, value_col, breed_order) {
  matrix_values <- matrix(
    NA_real_,
    nrow = length(breed_order),
    ncol = length(breed_order),
    dimnames = list(breed_order, breed_order)
  )

  for (idx in seq_len(nrow(summary))) {
    pop1 <- summary$POP1[[idx]]
    pop2 <- summary$POP2[[idx]]
    value <- as.numeric(summary[[value_col]][[idx]])
    matrix_values[pop1, pop2] <- value
    matrix_values[pop2, pop1] <- value
  }

  diag(matrix_values) <- 0

  if (anyNA(matrix_values)) {
    stop(sprintf("FST matrix contains %d missing pairwise values", sum(is.na(matrix_values))))
  }

  matrix_values
}

write_matrix_tsv <- function(matrix_values, output_path) {
  matrix_df <- as.data.frame(matrix_values, check.names = FALSE)
  for (column in names(matrix_df)) {
    matrix_df[[column]] <- sprintf("%.6f", matrix_df[[column]])
  }
  write.table(
    matrix_df,
    file = output_path,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )
}

get_breed_counts <- function(fam_path, breeds) {
  counts <- setNames(rep(0L, length(breeds)), breeds)
  if (!file.exists(fam_path)) {
    return(counts)
  }

  fam <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
  fam_counts <- table(fam[[1]])
  shared <- intersect(names(fam_counts), breeds)
  counts[shared] <- as.integer(fam_counts[shared])
  counts
}

make_column_labels <- function(breeds, breed_labels) {
  vapply(breeds, function(breed) breed_labels[[breed]], character(1))
}

make_row_labels <- function(breeds, breed_labels, breed_counts) {
  vapply(
    breeds,
    function(breed) {
      label <- breed_labels[[breed]]
      count <- breed_counts[[breed]]
      if (!is.null(count) && !is.na(count) && count > 0L) {
        sprintf("%s (n=%d)", label, count)
      } else {
        label
      }
    },
    character(1)
  )
}

palette_values <- function(n = 256L) {
  grDevices::colorRampPalette(
    c("#fff3d6", "#fdd17d", "#fca55d", "#ef6548", "#d7301f", "#7f0000")
  )(n)
}

scale_offdiag_value <- function(value, vmin, vmax) {
  if (isTRUE(all.equal(vmax, vmin))) {
    return(0.5)
  }
  sqrt((value - vmin) / (vmax - vmin))
}

cell_fill <- function(value, palette, vmin, vmax, diagonal_fill) {
  if (isTRUE(all.equal(value, 0))) {
    return(diagonal_fill)
  }
  scaled <- scale_offdiag_value(value, vmin, vmax)
  scaled <- min(max(scaled, 0), 1)
  index <- floor(scaled * (length(palette) - 1L)) + 1L
  palette[[index]]
}

cell_text_color <- function(value, vmin, vmax) {
  if (isTRUE(all.equal(value, 0))) {
    return("#4f4f4f")
  }
  scaled <- scale_offdiag_value(value, vmin, vmax)
  if (scaled >= 0.62) "#ffffff" else "#2b2b2b"
}

draw_top_dendrogram <- function(hc) {
  if (is.null(hc)) {
    plot.new()
    return(invisible(NULL))
  }

  par(mar = c(0.2, 2.4, 0.5, 8.0))
  plot(as.dendrogram(hc), leaflab = "none", axes = FALSE, xlab = "", ylab = "")
  axis(2, las = 1, cex.axis = 0.8, tck = -0.02)
}

draw_heatmap <- function(matrix_values, column_labels, row_labels, palette, vmin, vmax) {
  diagonal_fill <- "#e7e7e7"
  n_items <- nrow(matrix_values)
  row_positions <- rev(seq_len(n_items))

  par(mar = c(6.6, 2.4, 0.2, 8.0), xpd = NA)
  plot.new()
  plot.window(xlim = c(0.5, n_items + 0.5), ylim = c(0.5, n_items + 0.5), asp = 1, xaxs = "i", yaxs = "i")

  for (row_idx in seq_len(n_items)) {
    y <- row_positions[[row_idx]]
    for (col_idx in seq_len(n_items)) {
      value <- matrix_values[row_idx, col_idx]
      rect(
        xleft = col_idx - 0.5,
        ybottom = y - 0.5,
        xright = col_idx + 0.5,
        ytop = y + 0.5,
        col = cell_fill(value, palette, vmin, vmax, diagonal_fill),
        border = "#ffffff",
        lwd = 0.9
      )
      label <- if (row_idx == col_idx) "0" else sprintf("%.3f", value)
      text(
        x = col_idx,
        y = y,
        labels = label,
        cex = 0.75,
        col = cell_text_color(value, vmin, vmax)
      )
    }
  }

  axis(1, at = seq_len(n_items), labels = FALSE, tick = FALSE)
  axis(2, at = row_positions, labels = FALSE, tick = FALSE)

  text(
    x = seq_len(n_items),
    y = rep(-0.22, n_items),
    labels = column_labels,
    srt = 45,
    adj = 1,
    cex = 0.82
  )
  text(
    x = rep(n_items + 0.85, n_items),
    y = row_positions,
    labels = row_labels,
    adj = 0,
    cex = 0.82
  )
}

draw_colorbar <- function(palette, vmin, vmax) {
  par(mar = c(6.6, 0.6, 0.2, 2.8), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  palette_positions <- seq(0, 1, length.out = length(palette) + 1L)
  for (idx in seq_along(palette)) {
    rect(0.25, palette_positions[[idx]], 0.60, palette_positions[[idx + 1L]], col = palette[[idx]], border = NA)
  }
  rect(0.25, 0, 0.60, 1, border = "#666666", lwd = 0.8)

  if (isTRUE(all.equal(vmax, vmin))) {
    axis(4, at = 0.5, labels = sprintf("%.3f", vmin), las = 1, tick = FALSE, line = -0.2, cex.axis = 0.82)
  } else {
    ticks <- unique(c(vmin, pretty(c(vmin, vmax), n = 5), vmax))
    ticks <- ticks[ticks >= vmin & ticks <= vmax]
    tick_positions <- sqrt((ticks - vmin) / (vmax - vmin))
    axis(
      4,
      at = tick_positions,
      labels = sprintf("%.3f", ticks),
      las = 1,
      tick = FALSE,
      line = -0.2,
      cex.axis = 0.82
    )
  }

  text(0.45, 1.05, labels = "Pairwise FST", cex = 0.88)
}

summarize_extremes <- function(matrix_values) {
  off_diag_mask <- row(matrix_values) != col(matrix_values)
  off_diag_values <- matrix_values
  off_diag_values[!off_diag_mask] <- NA_real_

  min_index <- which(off_diag_values == min(off_diag_values, na.rm = TRUE), arr.ind = TRUE)[1, ]
  max_index <- which(off_diag_values == max(off_diag_values, na.rm = TRUE), arr.ind = TRUE)[1, ]

  list(
    min_pair = rownames(matrix_values)[min_index[[1]]],
    min_pair_other = colnames(matrix_values)[min_index[[2]]],
    min_value = matrix_values[min_index[[1]], min_index[[2]]],
    max_pair = rownames(matrix_values)[max_index[[1]]],
    max_pair_other = colnames(matrix_values)[max_index[[2]]],
    max_value = matrix_values[max_index[[1]], max_index[[2]]]
  )
}

config_path <- get_config_path()
cfg <- yaml.load_file(config_path)
run_dir <- resolve_run_dir()
fst_dir <- file.path(run_dir, "fst")
fig_dir <- file.path(run_dir, "figures", "fst")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}

summary_path <- file.path(fst_dir, sprintf("%s.fst.summary", SUMMARY_BASENAME))
summary <- read.table(
  summary_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)
names(summary) <- sub("^#", "", names(summary))
value_col <- if ("WC_FST" %in% names(summary)) {
  "WC_FST"
} else if ("HUDSON_FST" %in% names(summary)) {
  "HUDSON_FST"
} else {
  stop("FST summary must contain either WC_FST or HUDSON_FST.")
}

populations <- unique(c(summary$POP1, summary$POP2))
breed_order <- get_breed_order(cfg, populations)
breed_labels <- get_breed_labels(cfg, breed_order)
matrix_values <- build_matrix(summary, value_col, breed_order)

hc <- if (nrow(matrix_values) >= 2L) hclust(as.dist(matrix_values), method = "average") else NULL
ordered_breeds <- if (is.null(hc)) rownames(matrix_values) else hc$labels[hc$order]
ordered_matrix <- matrix_values[ordered_breeds, ordered_breeds, drop = FALSE]

matrix_path <- file.path(fst_dir, sprintf("%s.matrix.tsv", SUMMARY_BASENAME))
write_matrix_tsv(ordered_matrix, matrix_path)

breed_counts <- get_breed_counts(
  file.path(run_dir, "pca", "pruned_snps", PRUNED_FAM_BASENAME),
  ordered_breeds
)
column_labels <- make_column_labels(ordered_breeds, breed_labels)
row_labels <- make_row_labels(ordered_breeds, breed_labels, breed_counts)

off_diag_values <- ordered_matrix[row(ordered_matrix) != col(ordered_matrix)]
vmin <- min(off_diag_values)
vmax <- max(off_diag_values)
palette <- palette_values()
figure_path <- file.path(fig_dir, "fst_distance_matrix.png")
estimator_label <- if (identical(value_col, "WC_FST")) "Weir-Cockerham" else "Hudson"

png(filename = figure_path, width = 2400, height = 2200, res = 300)
layout(matrix(c(1, 0, 2, 3), nrow = 2, byrow = TRUE), widths = c(6.8, 1.4), heights = c(1.55, 6.8))
par(oma = c(0.6, 0.4, 2.6, 0.4), family = "sans")

draw_top_dendrogram(hc)
draw_heatmap(ordered_matrix, column_labels, row_labels, palette, vmin, vmax)
draw_colorbar(palette, vmin, vmax)

mtext("Pairwise FST Matrix", side = 3, line = 0.4, outer = TRUE, cex = 1.5, font = 2)
dev.off()

extremes <- summarize_extremes(ordered_matrix)
cat(sprintf("FST summary: %s\n", summary_path))
cat(sprintf("Matrix written to %s\n", matrix_path))
cat(sprintf("Figures written to %s/\n", fig_dir))
cat(sprintf("Closest pair: %s vs %s (%.4f)\n", extremes$min_pair, extremes$min_pair_other, extremes$min_value))
cat(sprintf("Most differentiated pair: %s vs %s (%.4f)\n", extremes$max_pair, extremes$max_pair_other, extremes$max_value))
