#!/usr/bin/env Rscript

# Within-month position coefficients (mean_lambda) — Italy screenshot models.
# Input:  results/position_decomposition_it.csv
# Output: results/plots/position_decomposition_it/

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(
    dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
    "v2_paths.R"
  )
} else {
  file.path(getwd(), "v2_paths.R")
}
source(v2.boot.path)

script_root <- script_dir
pos_file <- file.path(script_root, "results", "position_decomposition_it.csv")
out_dir  <- file.path(script_root, "results", "plots", "position_decomposition_it")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pos_file)) stop("Missing: ", pos_file)

pos <- read.csv(pos_file, stringsAsFactors = FALSE)

model_defs <- list(
  list(
    spec_id = "spec_136",
    label = "R-MIDAS + gas + ISTAT BCI",
    lf_vars = c("istat_bci_it")
  ),
  list(
    spec_id = "spec_141",
    label = "R-MIDAS + gas + IPI C.Goods + IPI Energy + ISTAT BCI",
    lf_vars = c("istat_bci_it", "Consumer_Goods", "Energy")
  )
)

pretty_lf <- function(x) {
  switch(x,
         istat_bci_it = "ISTAT BCI",
         Consumer_Goods = "IPI Consumer Goods",
         Energy = "IPI Energy",
         gsub("_", " ", x, fixed = TRUE))
}

palette_cols <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e")

plot_position_profile <- function(wide, coef_cols, title, png_out) {
  png(png_out, width = 1400, height = 800, res = 150)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(5, 5, 4, 2))

  ymat <- as.matrix(wide[, coef_cols, drop = FALSE])
  mode(ymat) <- "numeric"
  leg_labels <- vapply(
    sub("^mean_lambda\\.", "", coef_cols),
    pretty_lf,
    character(1)
  )

  matplot(
    x = wide$position,
    y = ymat,
    type = "l",
    lty = 1,
    lwd = 2,
    col = palette_cols[seq_len(ncol(ymat))],
    xlab = "Within-month position (day 1-28)",
    ylab = expression("Mean " * lambda[i]),
    main = title,
    xaxt = "n",
    xlim = c(1, 28)
  )
  axis(1, at = seq(1, 28, by = 1))
  abline(h = 0, col = "gray50", lty = 2, lwd = 1)
  legend(
    "topright",
    legend = leg_labels,
    col = palette_cols[seq_len(ncol(ymat))],
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.85
  )
}

for (md in model_defs) {
  block <- pos[pos$spec_id == md$spec_id & pos$lf_var %in% md$lf_vars, , drop = FALSE]
  if (nrow(block) == 0L) stop("No rows for ", md$spec_id)

  wide <- reshape(
    block[, c("position", "lf_var", "mean_lambda")],
    direction = "wide",
    idvar = "position",
    timevar = "lf_var",
    v.names = "mean_lambda"
  )
  wide <- wide[order(wide$position), , drop = FALSE]
  coef_cols <- grep("^mean_lambda", names(wide), value = TRUE)

  png_out <- file.path(out_dir, sprintf("position_profile_%s.png", md$spec_id))
  plot_position_profile(
    wide, coef_cols,
    sprintf("Italy — %s (%s)", md$label, md$spec_id),
    png_out
  )
  cat("wrote", png_out, "\n")
}

# Both screenshot models on one figure (all LF series)
all_rows <- list()
for (md in model_defs) {
  block <- pos[pos$spec_id == md$spec_id & pos$lf_var %in% md$lf_vars, , drop = FALSE]
  block$series <- paste0(md$spec_id, ": ", vapply(block$lf_var, pretty_lf, character(1)))
  all_rows[[length(all_rows) + 1L]] <- block[, c("position", "series", "mean_lambda")]
}
long_df <- do.call(rbind, all_rows)
wide_all <- reshape(
  long_df,
  direction = "wide",
  idvar = "position",
  timevar = "series",
  v.names = "mean_lambda"
)
wide_all <- wide_all[order(wide_all$position), , drop = FALSE]
coef_cols_all <- grep("^mean_lambda", names(wide_all), value = TRUE)

png_all <- file.path(out_dir, "position_profile_screenshot_models_it.png")
plot_position_profile(
  wide_all, coef_cols_all,
  "Italy — within-month LF coefficients (screenshot models)",
  png_all
)
cat("wrote", png_all, "\n")
