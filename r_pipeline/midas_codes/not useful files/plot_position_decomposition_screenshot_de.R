#!/usr/bin/env Rscript

# Within-month position profiles (mean λ_i) — Germany "Results I" R-MIDAS specs.
# For each calendar day i in {1..28}, coefficients are averaged across months
# at that position (see position_decomposition_<cc>.csv from legacy short run).
#
# Input:  results/position_decomposition_de.csv  (override: POS_FILE)
# Output: results/plots/position_decomposition_de/

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
pos_file <- Sys.getenv(
  "POS_FILE",
  unset = file.path(script_root, "results", "position_decomposition_de.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(script_root, "results", "plots", "position_decomposition_de")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pos_file)) stop("Missing: ", pos_file)

pos <- read.csv(pos_file, stringsAsFactors = FALSE)

model_defs <- list(
  list(
    spec_id = "spec_239",
    label = "R-MIDAS + gas + brent + PMI MFG + ifo BCI + IPI (C.Goods, Mfg, Energy)",
    lf_vars = c(
      "pmi_mfg_de", "ifo_bci_de", "Consumer_Goods", "Manufacturing", "Energy"
    )
  ),
  list(
    spec_id = "spec_174",
    label = "R-MIDAS + gas + PMI MFG + ifo BCI + IPI (C.Goods, Mfg)",
    lf_vars = c("pmi_mfg_de", "ifo_bci_de", "Consumer_Goods", "Manufacturing")
  ),
  list(
    spec_id = "spec_170",
    label = "R-MIDAS + gas + PMI MFG + ifo BCI + IPI (Mfg)",
    lf_vars = c("pmi_mfg_de", "ifo_bci_de", "Manufacturing")
  ),
  list(
    spec_id = "spec_168",
    label = "R-MIDAS + gas + PMI MFG + ifo BCI",
    lf_vars = c("pmi_mfg_de", "ifo_bci_de")
  )
)

pretty_lf <- function(x) {
  switch(x,
         pmi_mfg_de = "PMI manufacturing DE",
         ifo_bci_de = "ifo BCI DE",
         Consumer_Goods = "IPI consumer goods",
         Manufacturing = "IPI manufacturing",
         Energy = "IPI energy",
         gsub("_", " ", x, fixed = TRUE))
}

palette_cols <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf")

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
    ylab = expression("Mean " * lambda[i] * " (avg. across months)"),
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

map_rows <- list()

for (k in seq_along(model_defs)) {
  md <- model_defs[[k]]
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

  csv_out <- file.path(out_dir, sprintf("position_profile_%s.csv", md$spec_id))
  write.csv(wide, csv_out, row.names = FALSE)

  coef_cols <- grep("^mean_lambda", names(wide), value = TRUE)
  png_out <- file.path(out_dir, sprintf("position_profile_%s.png", md$spec_id))
  plot_position_profile(
    wide,
    coef_cols,
    sprintf("Germany — %s", md$label),
    png_out
  )
  cat("wrote", png_out, "\n")
  cat("wrote", csv_out, "\n")

  map_rows[[k]] <- data.frame(
    spec_id = md$spec_id,
    label = md$label,
    lf_vars = paste(md$lf_vars, collapse = "|"),
    png = png_out,
    csv = csv_out,
    stringsAsFactors = FALSE
  )
}

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

png_all <- file.path(out_dir, "position_profile_screenshot_models_de.png")
plot_position_profile(
  wide_all,
  coef_cols_all,
  "Germany — within-month LF coefficients (Results I R-MIDAS specs)",
  png_all
)
cat("wrote", png_all, "\n")

map_df <- do.call(rbind, map_rows)
write.csv(
  map_df,
  file.path(out_dir, "position_profile_screenshot_map_de.csv"),
  row.names = FALSE
)

cat("Saved plots in: ", out_dir, "\n", sep = "")
