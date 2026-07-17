#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — PER-SPEC FORECAST PLOTS
#
# For every model that has been computed by the active pipeline, generate a
# plot of three time series:
#   * y_actual                           — realised electricity price
#   * y_hat_<spec>                       — forecast from the current model
#   * y_hat_ar_dum                       — AR(1,2,3,7) + dummies benchmark
# All three series are on the €/MWh (back-transformed) scale.
#
# Output is organised as one multi-page PDF per (family, country, horizon)
# combination — one page per spec.  This gives the user a flippable book
# of forecasts without exploding the file count.
#
# Family / file mapping (must match summary_table.R FAMILY_PATTERNS):
#   AR                : ar_<cc>_h<h>.csv               (4 specs incl. benchmark)
#   RU_MIDAS          : ru_midas_<cc>_h<h>.csv          (256 specs)
#   R_MIDAS           : r_midas_<cc>_h<h>.csv           (256 specs)
#   R_MIDAS_extended  : r_midas_extended_<cc>_h<h>.csv  (256 specs)
#
# Parallelism: parallel::mclapply across (family, country, h) jobs.  Each
# worker writes its own PDF — no contention.  N_CORES env var controls
# concurrency.  Expected runtime: ~5-15 minutes for the full set on M4.
#
# Knobs:
#   TOP_N_PER_PDF env var (default unlimited) — restrict each PDF to the
#                                                top-N specs by RMSE.  Useful
#                                                for a preview without
#                                                generating 256-page PDFs.
#   COUNTRY_NAME, HORIZONS — inherited from common_utils.R defaults.
# =============================================================================

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grep("^--file=", args)]
  if (length(farg) > 0) dirname(normalizePath(sub("^--file=", "", farg))) else getwd()
}, error = function(e) getwd())
if (dir.exists(script_dir)) setwd(script_dir)

source(file.path(script_dir, "common_utils.R"))

FORECAST_DIR <- file.path(script_dir, "forecasts")
PLOT_DIR     <- file.path(script_dir, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N_PER_PDF <- suppressWarnings(as.integer(Sys.getenv("TOP_N_PER_PDF", "")))
if (is.na(TOP_N_PER_PDF) || TOP_N_PER_PDF <= 0L) TOP_N_PER_PDF <- .Machine$integer.max

# -----------------------------------------------------------------------------
# Families and the filename patterns they emit
# -----------------------------------------------------------------------------
FAMILIES <- list(
  AR               = list(file = function(cc, h) sprintf("ar_%s_h%02d.csv",                cc, h),
                          dict = NULL),
  RU_MIDAS         = list(file = function(cc, h) sprintf("ru_midas_%s_h%02d.csv",          cc, h),
                          dict = function(cc) sprintf("ru_midas_spec_dictionary_%s.csv",   cc)),
  R_MIDAS          = list(file = function(cc, h) sprintf("r_midas_%s_h%02d.csv",           cc, h),
                          dict = function(cc) sprintf("r_midas_spec_dictionary_%s.csv",    cc)),
  R_MIDAS_extended = list(file = function(cc, h) sprintf("r_midas_extended_%s_h%02d.csv",  cc, h),
                          dict = function(cc) sprintf("r_midas_extended_spec_dictionary_%s.csv", cc))
)

# -----------------------------------------------------------------------------
# Helper: per-spec page plot
# -----------------------------------------------------------------------------
plot_one_spec <- function(target_date, y_actual, y_bench, y_spec,
                          country, h, family, spec_id, hf_vars, lf_vars,
                          rmse_lev, rmse_ratio_lev, rmse_asinh, rmse_ratio_asinh,
                          dm_stat_lev = NA, dm_pval_lev = NA) {
  # Robust upper y-cap for visual readability (allows the 2022-23 spikes to
  # show but caps the y-axis at p99.5 + a margin so non-crisis dynamics are
  # still legible).  Use the union of actual + benchmark + spec.
  yrange <- range(c(y_actual, y_bench, y_spec), na.rm = TRUE, finite = TRUE)
  cap_hi <- quantile(c(y_actual, y_bench, y_spec), 0.995, na.rm = TRUE) * 1.15
  cap_lo <- min(yrange[1], 0)
  ylim   <- c(cap_lo, max(yrange[2], cap_hi, na.rm = TRUE))

  hf_lab <- if (nzchar(hf_vars)) hf_vars else "—"
  lf_lab <- if (nzchar(lf_vars)) lf_vars else "—"

  par(mar = c(4.0, 4.0, 4.6, 1.2))
  plot(target_date, y_actual, type = "l", col = "black", lwd = 1.0,
       xlab = "target date", ylab = "elec price  (€/MWh)",
       ylim = ylim,
       main = sprintf("%s  |  %s  h=%d  |  %s",
                      family, toupper(country), h, spec_id))
  # Subtitle / metrics in mtext line 2 of the main margin
  mtext(sprintf("HF: %s   LF: %s", hf_lab, lf_lab),
        side = 3, line = 1.6, cex = 0.85)
  mtext(sprintf("RMSE €/MWh = %.2f  (ratio %.3f)   |   RMSE asinh = %.4f  (ratio %.3f)   |   DM = %.2f  (p = %.3g)",
                rmse_lev, rmse_ratio_lev, rmse_asinh, rmse_ratio_asinh,
                dm_stat_lev, dm_pval_lev),
        side = 3, line = 0.55, cex = 0.75)
  lines(target_date, y_bench, col = "firebrick", lwd = 1.0, lty = 2)
  lines(target_date, y_spec,  col = "steelblue", lwd = 1.0)
  legend("topleft", inset = c(0.005, 0.005), bty = "n", cex = 0.80,
         legend = c("Actual", "Benchmark  (ar_dum)", sprintf("Model  (%s)", spec_id)),
         col    = c("black",  "firebrick",            "steelblue"),
         lty    = c(1,        2,                      1),
         lwd    = c(1.4,      1.4,                    1.4))
  grid(lty = 3, col = "grey85")
}

# -----------------------------------------------------------------------------
# Build the job list (family, country, h)
# -----------------------------------------------------------------------------
build_jobs <- function() {
  jobs <- list()
  for (fam in names(FAMILIES)) {
    for (cc in c("de", "it")) {
      for (h in HORIZONS) {
        f <- FAMILIES[[fam]]$file(cc, h)
        if (file.exists(file.path(FORECAST_DIR, f))) {
          jobs[[length(jobs) + 1L]] <- list(family = fam, cc = cc, h = h)
        }
      }
    }
  }
  jobs
}

# -----------------------------------------------------------------------------
# One worker = one PDF (family × country × h)
# -----------------------------------------------------------------------------
run_one_job <- function(job) {
  fam <- job$family; cc <- job$cc; h <- job$h
  fmeta <- FAMILIES[[fam]]
  fpath <- file.path(FORECAST_DIR, fmeta$file(cc, h))
  ar_path <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))

  if (!file.exists(fpath) || !file.exists(ar_path))
    return(list(family = fam, cc = cc, h = h, status = "missing"))

  d  <- read.csv(fpath, stringsAsFactors = FALSE)
  ar <- read.csv(ar_path, stringsAsFactors = FALSE)
  d$target_date  <- as.Date(d$target_date)
  ar$target_date <- as.Date(ar$target_date)

  # Align on target_date
  common <- intersect(as.character(ar$target_date), as.character(d$target_date))
  if (length(common) == 0L)
    return(list(family = fam, cc = cc, h = h, status = "no_overlap"))
  d  <- d[match(common, as.character(d$target_date)), , drop = FALSE]
  ar <- ar[match(common, as.character(ar$target_date)), , drop = FALSE]
  ar_dum <- ar$ar_dum
  td     <- ar$target_date
  ya     <- ar$y_actual

  # Identify spec columns: AR family has named cols, MIDAS families have spec_NNN
  spec_cols <- if (fam == "AR") {
    setdiff(colnames(d), c("origin_date", "target_date", "y_actual",
                            "target_pos"))   # target_pos appears in some files
  } else {
    grep("^spec_", colnames(d), value = TRUE)
  }
  if (length(spec_cols) == 0L)
    return(list(family = fam, cc = cc, h = h, status = "no_specs"))

  # Optional dictionary join
  dict <- NULL
  if (!is.null(fmeta$dict)) {
    dp <- file.path(FORECAST_DIR, fmeta$dict(cc))
    if (file.exists(dp)) dict <- read.csv(dp, stringsAsFactors = FALSE)
  }

  # Compute per-spec RMSE for ordering / annotation
  rmse_levs <- sapply(spec_cols, function(sp) {
    e <- ya - d[[sp]]; ok <- is.finite(e); if (!any(ok)) NA_real_ else sqrt(mean(e[ok]^2))
  })
  rmse_asinhs <- sapply(spec_cols, function(sp) {
    e <- asinh(ya) - asinh(d[[sp]]); ok <- is.finite(e)
    if (!any(ok)) NA_real_ else sqrt(mean(e[ok]^2))
  })
  rmse_bench_lev <- sqrt(mean((ya - ar_dum)^2, na.rm = TRUE))
  rmse_bench_asinh <- sqrt(mean((asinh(ya) - asinh(ar_dum))^2, na.rm = TRUE))
  rmse_ratio_levs   <- rmse_levs   / rmse_bench_lev
  rmse_ratio_asinhs <- rmse_asinhs / rmse_bench_asinh

  # DM-HLN per spec (on €/MWh, sqerr loss) for the annotation
  dm <- sapply(spec_cols, function(sp) {
    e_t <- ya - d[[sp]]; e_b <- ya - ar_dum
    res <- dm_hln(e_t, e_b, h, loss = "sqerr")
    c(stat = res$stat, pval = res$pvalue)
  })
  dm_stats <- dm["stat", ]
  dm_pvals <- dm["pval", ]

  # Order specs by RMSE ascending; truncate to TOP_N_PER_PDF if set
  ord <- order(rmse_levs, na.last = TRUE)
  if (TOP_N_PER_PDF < length(ord)) ord <- ord[1:TOP_N_PER_PDF]

  # Write PDF: one page per spec
  out_pdf <- file.path(PLOT_DIR, sprintf("%s_%s_h%02d.pdf",
                                          tolower(fam), cc, h))
  pdf(out_pdf, width = 11, height = 6.2, paper = "special")
  for (idx in ord) {
    sp <- spec_cols[idx]
    hfv <- if (!is.null(dict)) {
      i <- match(sp, dict$spec_id); if (!is.na(i)) dict$hf_vars[i] else ""
    } else ""
    lfv <- if (!is.null(dict)) {
      i <- match(sp, dict$spec_id); if (!is.na(i)) dict$lf_vars[i] else ""
    } else ""
    plot_one_spec(td, ya, ar_dum, d[[sp]],
                  country = cc, h = h, family = fam, spec_id = sp,
                  hf_vars = hfv, lf_vars = lfv,
                  rmse_lev = rmse_levs[idx], rmse_ratio_lev = rmse_ratio_levs[idx],
                  rmse_asinh = rmse_asinhs[idx], rmse_ratio_asinh = rmse_ratio_asinhs[idx],
                  dm_stat_lev = dm_stats[idx], dm_pval_lev = dm_pvals[idx])
  }
  dev.off()
  list(family = fam, cc = cc, h = h, status = "ok",
       n_pages = length(ord), file = out_pdf)
}

# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------
n_cores <- get_n_cores()
jobs <- build_jobs()
cat(sprintf("=== plot_forecasts.R — %d (family, country, h) jobs  (n_cores=%d) ===\n",
            length(jobs), n_cores))

t_start <- Sys.time()
job_results <- parallel_lapply(jobs, run_one_job, n_cores = n_cores)
el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

for (r in job_results) {
  if (is.null(r)) next
  if (r$status == "ok")
    cat(sprintf("  %s %s h=%d → %s  (%d pages)\n",
                r$family, r$cc, r$h, r$file, r$n_pages))
  else
    cat(sprintf("  %s %s h=%d → skipped (%s)\n", r$family, r$cc, r$h, r$status))
}
cat(sprintf("\nplot_forecasts done in %.1fs (%.1fmin) — wrote %d PDFs to %s\n",
            el, el / 60,
            sum(vapply(job_results, function(r) isTRUE(!is.null(r) && r$status == "ok"),
                       logical(1))),
            PLOT_DIR))
