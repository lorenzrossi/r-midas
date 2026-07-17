#!/usr/bin/env Rscript

# Plot Germany forecasts for models in "Results I" screenshot.
# Uses standard R-MIDAS forecasts (r_midas_de_h*.csv) + AR (ar_de_h*.csv).

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(
    dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
    "v2_paths.R"
  )
} else file.path(getwd(), "v2_paths.R")
source(v2.boot.path)

script_root <- script_dir
forecast_dir <- file.path(script_root, "forecasts")
out_dir <- file.path(script_root, "results", "plots", "screenshot_models_de")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

HORIZONS <- c(1, 7, 14, 21, 28)

model_defs <- list(
  list(
    label = "AR(7) + dum",
    family = "AR",
    col = "ar_dum"
  ),
  list(
    label = "AR(7) + dum + gas + brent",
    family = "AR",
    col = "ar_dum_gas_brent"
  ),
  list(
    label = "R-MIDAS + gas + brent + PMI + ifo + IPI (all)",
    family = "R_MIDAS",
    spec_id = "spec_239",
    hf_vars = "gas|brent",
    lf_vars = "pmi_mfg_de|ifo_bci_de|Consumer_Goods|Manufacturing|Energy"
  ),
  list(
    label = "R-MIDAS + gas + PMI + ifo + IPI (C.Goods, Mfg)",
    family = "R_MIDAS",
    spec_id = "spec_174",
    hf_vars = "gas",
    lf_vars = "pmi_mfg_de|ifo_bci_de|Consumer_Goods|Manufacturing"
  ),
  list(
    label = "R-MIDAS + gas + PMI + ifo + IPI (Mfg)",
    family = "R_MIDAS",
    spec_id = "spec_170",
    hf_vars = "gas",
    lf_vars = "pmi_mfg_de|ifo_bci_de|Manufacturing"
  ),
  list(
    label = "R-MIDAS + gas + PMI + ifo",
    family = "R_MIDAS",
    spec_id = "spec_168",
    hf_vars = "gas",
    lf_vars = "pmi_mfg_de|ifo_bci_de"
  )
)

palette_cols <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"
)

plot_from <- as.Date("2019-01-01")
plot_to <- as.Date("2025-12-31")

for (h in HORIZONS) {
  ar_path <- file.path(forecast_dir, sprintf("ar_de_h%02d.csv", h))
  rm_path <- file.path(forecast_dir, sprintf("r_midas_de_h%02d.csv", h))
  if (!file.exists(ar_path) || !file.exists(rm_path)) {
    warning("Missing forecast file for h=", h)
    next
  }

  ar <- read.csv(ar_path, stringsAsFactors = FALSE)
  rm <- read.csv(rm_path, stringsAsFactors = FALSE)
  ar$target_date <- as.Date(ar$target_date)
  rm$target_date <- as.Date(rm$target_date)

  base <- ar[, c("target_date", "y_actual"), drop = FALSE]
  base <- base[base$target_date >= plot_from & base$target_date <= plot_to, , drop = FALSE]
  base <- base[order(base$target_date), , drop = FALSE]

  ymat <- matrix(NA_real_, nrow = nrow(base), ncol = length(model_defs))
  colnames(ymat) <- vapply(model_defs, function(m) m$label, character(1))

  for (j in seq_along(model_defs)) {
    m <- model_defs[[j]]
    if (m$family == "AR") {
      src <- ar
      if (!m$col %in% names(src)) next
      v <- src[[m$col]]
      names(v) <- as.character(src$target_date)
      ymat[, j] <- v[as.character(base$target_date)]
    } else {
      src <- rm
      if (!m$spec_id %in% names(src)) next
      v <- src[[m$spec_id]]
      names(v) <- as.character(src$target_date)
      ymat[, j] <- v[as.character(base$target_date)]
    }
  }

  png(
    file.path(out_dir, sprintf("forecasts_de_h%02d.png", h)),
    width = 1600, height = 900, res = 150
  )
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(5, 5, 4, 1))

  y_all <- c(base$y_actual, ymat)
  y_all <- y_all[is.finite(y_all)]
  y_rng <- range(y_all, na.rm = TRUE)

  plot(
    base$target_date, base$y_actual,
    type = "l", lwd = 2.2, col = "black",
    xlab = "Target date", ylab = "EUR/MWh",
    main = sprintf("Germany — h = %d (2019–2025)", h),
    ylim = y_rng
  )
  for (j in seq_len(ncol(ymat))) {
    lines(
      base$target_date, ymat[, j],
      lwd = 1.6, col = palette_cols[j]
    )
  }
  legend(
    "topright",
    legend = c("Actual", colnames(ymat)),
    col = c("black", palette_cols[seq_len(ncol(ymat))]),
    lty = 1,
    lwd = c(2.2, rep(1.6, ncol(ymat))),
    bty = "n",
    cex = 0.75
  )
  dev.off()
}

# Combined panel: h = 1 only, full sample (lighter)
h <- 1L
ar <- read.csv(file.path(forecast_dir, "ar_de_h01.csv"), stringsAsFactors = FALSE)
rm <- read.csv(file.path(forecast_dir, "r_midas_de_h01.csv"), stringsAsFactors = FALSE)
ar$target_date <- as.Date(ar$target_date)
rm$target_date <- as.Date(rm$target_date)
base <- ar[, c("target_date", "y_actual")]
base <- base[order(base$target_date), ]

ymat <- matrix(NA_real_, nrow = nrow(base), ncol = length(model_defs))
colnames(ymat) <- vapply(model_defs, function(m) m$label, character(1))
for (j in seq_along(model_defs)) {
  m <- model_defs[[j]]
  if (m$family == "AR") {
    ymat[, j] <- ar[[m$col]]
  } else {
    ymat[, j] <- rm[[m$spec_id]]
  }
}

png(file.path(out_dir, "forecasts_de_h01_full_sample.png"), width = 1600, height = 900, res = 150)
par(mar = c(5, 5, 4, 1))
y_all <- c(base$y_actual, ymat)
y_rng <- range(y_all, na.rm = TRUE)
plot(base$target_date, base$y_actual, type = "l", lwd = 2.2, col = "black",
     xlab = "Target date", ylab = "EUR/MWh",
     main = "Germany — h = 1 (full evaluation sample)", ylim = y_rng)
for (j in seq_len(ncol(ymat))) lines(base$target_date, ymat[, j], lwd = 1.4, col = palette_cols[j])
legend("topright", legend = c("Actual", colnames(ymat)),
       col = c("black", palette_cols), lty = 1, lwd = c(2.2, rep(1.4, ncol(ymat))),
       bty = "n", cex = 0.7)
dev.off()

map_df <- do.call(rbind, lapply(model_defs, function(m) {
  data.frame(
    label = m$label,
    family = m$family,
    spec_or_col = if (m$family == "AR") m$col else m$spec_id,
    hf_vars = if (m$family == "AR") NA_character_ else m$hf_vars,
    lf_vars = if (m$family == "AR") NA_character_ else m$lf_vars,
    stringsAsFactors = FALSE
  )
}))
write.csv(map_df, file.path(out_dir, "screenshot_models_map_de.csv"), row.names = FALSE)

cat("Saved plots in: ", out_dir, "\n", sep = "")
