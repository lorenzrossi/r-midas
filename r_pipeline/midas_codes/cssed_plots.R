#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — CUMULATIVE SUM OF SQUARED-ERROR DIFFERENCES (Welch-Goyal 2008)
#
# One panel per country x horizon, x-axis = target date.  A RISING line means
# the model is beating its reference.  Curves:
#   1. ARX daily-only (ar_dum_gas_brent) vs AR benchmark
#   2. real-time macro combination comb_invmse of each family (ARX, R_MIDAS,
#      R_MIDAS_L7) vs AR benchmark
#   3. each family's comb_invmse vs its OWN daily-only gas+brent twin
#      (dashed: macro information alone, lag structure held fixed)
#   4. R_MIDAS_L7 spec_000 (no macro) vs AR benchmark (lag-structure effect)
# Subperiods shaded: pre <= 2020-12-31, crisis 2021-2022, post >= 2023.
# Base graphics only.  Outputs: results/cssed/cssed_<cc>_h<hh>.png and
# results/cssed/cssed_<cc>.csv (long: h, target_date, curve, value).
# The Giacomini-Rossi fluctuation test (results/gr_fluctuation_<cc>.csv)
# provides the formal critical values alongside these plots.
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
            "v2_paths.R")
} else file.path(getwd(), "v2_paths.R")
source(v2.boot.path)
source(file.path(script_dir, "common_utils.R"))

FORECAST_DIR <- file.path(script_dir, "forecasts")
OUT_DIR      <- file.path(script_dir, "results", "cssed")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FAMS <- c(ARX = "arx", R_MIDAS = "r_midas", R_MIDAS_L7 = "r_midas_l7")
COLS <- c(ARX = "#1b9e77", R_MIDAS = "#d95f02", R_MIDAS_L7 = "#7570b3")

read_aligned <- function(fp, key) {
  if (!file.exists(fp)) return(NULL)
  d <- read.csv(fp, stringsAsFactors = FALSE)
  d[match(key, paste(d$origin_date, d$target_date)), , drop = FALSE]
}

run_country <- function(cc) {
  cat(sprintf("\n=== CSSED plots — %s ===\n", toupper(cc)))
  long <- list()
  for (h in HORIZONS) {
    ar_fp <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    if (!file.exists(ar_fp)) next
    ar <- read.csv(ar_fp, stringsAsFactors = FALSE)
    ar$target_date <- as.Date(ar$target_date)
    ar <- ar[order(ar$target_date), ]
    key <- paste(ar$origin_date, ar$target_date)
    y <- ar$y_actual; td <- ar$target_date
    e_b <- ar$ar_dum - y
    curves <- list()
    if ("ar_dum_gas_brent" %in% names(ar))
      curves[["AR + gas + brent vs AR"]] <- list(v = cssed(e_b, ar$ar_dum_gas_brent - y),
                                                 col = "black", lty = 1, lwd = 2)
    for (famlab in names(FAMS)) {
      fam <- FAMS[[famlab]]
      cb <- read_aligned(file.path(FORECAST_DIR, sprintf("comb_%s_%s_h%02d.csv", fam, cc, h)), key)
      fc <- read_aligned(file.path(FORECAST_DIR, sprintf("%s_%s_h%02d.csv", fam, cc, h)), key)
      dp <- file.path(FORECAST_DIR, sprintf("%s_spec_dictionary_%s.csv", fam, cc))
      if (is.null(cb) || is.null(fc) || !file.exists(dp)) next
      dict <- read.csv(dp, stringsAsFactors = FALSE)
      dict$hf_vars[is.na(dict$hf_vars)] <- ""; dict$lf_vars[is.na(dict$lf_vars)] <- ""
      twin_gb <- dict$spec_id[dict$hf_vars == "gas|brent" & dict$lf_vars == ""][1]
      if ("comb_invmse" %in% names(cb)) {
        curves[[sprintf("%s comb_invmse vs AR", famlab)]] <-
          list(v = cssed(e_b, cb$comb_invmse - y), col = COLS[[famlab]], lty = 1, lwd = 2)
        if (!is.na(twin_gb) && twin_gb %in% names(fc))
          curves[[sprintf("%s comb_invmse vs own gas+brent twin", famlab)]] <-
            list(v = cssed(fc[[twin_gb]] - y, cb$comb_invmse - y), col = COLS[[famlab]], lty = 2, lwd = 1.5)
      }
      if (famlab == "R_MIDAS_L7" && "spec_000" %in% names(fc))
        curves[["R_MIDAS_L7 spec_000 (no macro) vs AR"]] <-
          list(v = cssed(e_b, fc$spec_000 - y), col = COLS[[famlab]], lty = 3, lwd = 2)
    }
    if (!length(curves)) next
    for (nm in names(curves))
      long[[length(long) + 1L]] <- data.frame(country = cc, h = h, target_date = td,
                                              curve = nm, value = curves[[nm]]$v,
                                              stringsAsFactors = FALSE)
    png(file.path(OUT_DIR, sprintf("cssed_%s_h%02d.png", cc, h)), width = 1400, height = 800, res = 130)
    rng <- range(unlist(lapply(curves, `[[`, "v")), na.rm = TRUE)
    par(mar = c(4, 5, 3, 1))
    plot(td, curves[[1]]$v, type = "n", ylim = rng, xlab = "target date",
         ylab = "cumulative (e_ref^2 - e_model^2), (EUR/MWh)^2",
         main = sprintf("%s, h = %d: cumulative squared-error differences (rising = model better)",
                        toupper(cc), h))
    usr <- par("usr")
    rect(as.numeric(as.Date("2021-01-01")), usr[3], as.numeric(as.Date("2022-12-31")), usr[4],
         col = "#fde0dd", border = NA)
    rect(as.numeric(as.Date("2023-01-01")), usr[3], usr[2], usr[4], col = "#edf8e9", border = NA)
    abline(h = 0, col = "grey50")
    for (nm in names(curves)) lines(td, curves[[nm]]$v, col = curves[[nm]]$col,
                                    lty = curves[[nm]]$lty, lwd = curves[[nm]]$lwd)
    legend("topleft", legend = names(curves), col = vapply(curves, `[[`, "", "col"),
           lty = vapply(curves, `[[`, 1, "lty"), lwd = vapply(curves, `[[`, 1, "lwd"),
           bty = "n", cex = 0.8)
    dev.off()
    cat(sprintf("  h=%02d: %d curves\n", h, length(curves)))
  }
  if (length(long)) {
    L <- do.call(rbind, long)
    write.csv(L, file.path(OUT_DIR, sprintf("cssed_%s.csv", cc)), row.names = FALSE)
    cat(sprintf("  wrote %s\n", file.path(OUT_DIR, sprintf("cssed_%s.csv", cc))))
  }
}
for (cc in eval_countries()) run_country(cc)
cat("\nCSSED plots done.\n")
