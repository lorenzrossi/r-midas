#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — SUMMARY TABLE FOR R-MIDAS LEGACY ONLY
#
# Same metric block as summary_table.R, but restricted to the R_MIDAS_legacy
# family.  Reads:
#   forecasts/ar_<cc>_h<h>.csv                            (benchmark = ar_dum)
#   forecasts/r_midas_legacy_<cc>_h<h>.csv                (legacy specs)
#   forecasts/r_midas_legacy_spec_dictionary_<cc>.csv     (HF / LF membership)
#
# Writes:
#   results/summary_legacy.csv                            (every legacy spec)
#   results/summary_legacy_top5.csv                       (top 5 per country, h)
#
# Path resolution — IDENTICAL to summary_table.R so this script works in any
# setup where summary_table.R already works.  Priority:
#   1) V2_DATA_DIR env var
#   2) V2_SCRIPT_DIR env var
#   3) directory of this script (derived from commandArgs)
#   4) getwd()
# =============================================================================

code_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grep("^--file=", args)]
  if (length(farg) > 0) {
    raw <- sub("^--file=", "", farg[1])
    raw <- gsub("~+~", " ", raw, fixed = TRUE)
    cand <- suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))
    if (is.na(cand) || !nzchar(cand)) getwd() else dirname(cand)
  } else getwd()
}, error = function(e) getwd())

data_root <- Sys.getenv("V2_DATA_DIR", unset = "")
if (!nzchar(data_root)) data_root <- Sys.getenv("V2_SCRIPT_DIR", unset = "")
if (nzchar(data_root)) {
  data_root <- suppressWarnings(normalizePath(data_root, winslash = "/", mustWork = FALSE))
}
if (!nzchar(data_root) || is.na(data_root) || !dir.exists(data_root)) {
  data_root <- code_dir
}
if (dir.exists(data_root)) setwd(data_root)

# Loud diagnostics so any misconfiguration is visible
cat("code_dir  :", code_dir,  "\n")
cat("data_root :", data_root, "\n")
cat("getwd()   :", getwd(),   "\n\n")

cu <- file.path(code_dir,  "common_utils.R")
if (!file.exists(cu)) cu <- file.path(data_root, "common_utils.R")
if (!file.exists(cu)) {
  stop("Cannot find common_utils.R.  Tried:\n  ",
       file.path(code_dir,  "common_utils.R"), "\n  ",
       file.path(data_root, "common_utils.R"), "\n",
       "Set V2_DATA_DIR to the v2_new_dataset folder that contains forecasts/ and common_utils.R.")
}
cat("sourcing common_utils.R from:", cu, "\n")
source(cu)

FORECAST_DIR <- file.path(data_root, "forecasts")
RESULT_DIR   <- file.path(data_root, "results")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("FORECAST_DIR:", FORECAST_DIR,
    "  (exists =", dir.exists(FORECAST_DIR), ")\n")
cat("RESULT_DIR  :", RESULT_DIR,
    "  (exists =", dir.exists(RESULT_DIR), ")\n\n")

# -----------------------------------------------------------------------------
# Metric blocks (same as summary_table.R)
# -----------------------------------------------------------------------------
.rmse_block <- function(e_t, e_b, h) {
  n <- length(e_t)
  if (n < 5L) return(list(rmse = NA_real_, rmse_se = NA_real_,
                          rmse_ratio = NA_real_,
                          dm_stat = NA_real_, dm_pval = NA_real_))
  rmse   <- sqrt(mean(e_t^2))
  rmse_b <- sqrt(mean(e_b^2))
  v_e2   <- var(e_t^2)
  rmse_se    <- if (rmse > 0) sqrt(v_e2 / (n * 4 * rmse^2)) else NA_real_
  rmse_ratio <- if (rmse_b > 0) rmse / rmse_b else NA_real_
  dm <- dm_hln(e_t, e_b, h, loss = "sqerr")
  list(rmse = rmse, rmse_se = rmse_se, rmse_ratio = rmse_ratio,
       dm_stat = dm$stat, dm_pval = dm$pvalue)
}
.mae_block <- function(e_t, e_b, h) {
  n <- length(e_t)
  if (n < 5L) return(list(mae = NA_real_, mae_se = NA_real_,
                          mae_ratio = NA_real_,
                          dm_stat = NA_real_, dm_pval = NA_real_))
  mae   <- mean(abs(e_t))
  mae_b <- mean(abs(e_b))
  mae_se <- sd(abs(e_t)) / sqrt(n)
  mae_ratio <- if (mae_b > 0) mae / mae_b else NA_real_
  dm <- dm_hln(e_t, e_b, h, loss = "abserr")
  list(mae = mae, mae_se = mae_se, mae_ratio = mae_ratio,
       dm_stat = dm$stat, dm_pval = dm$pvalue)
}
compute_metrics <- function(y_actual, y_hat_test, y_hat_bench, h) {
  ok <- is.finite(y_actual) & is.finite(y_hat_test) & is.finite(y_hat_bench)
  y_a <- y_actual[ok]; y_t <- y_hat_test[ok]; y_b <- y_hat_bench[ok]

  e_t_p <- y_a - y_t; e_b_p <- y_a - y_b
  rmse_lev <- .rmse_block(e_t_p, e_b_p, h)
  mae_lev  <- .mae_block (e_t_p, e_b_p, h)

  ya_a <- asinh(y_a); ya_t <- asinh(y_t); ya_b <- asinh(y_b)
  e_t_a <- ya_a - ya_t; e_b_a <- ya_a - ya_b
  rmse_a <- .rmse_block(e_t_a, e_b_a, h)
  mae_a  <- .mae_block (e_t_a, e_b_a, h)

  sig2_t <- mean(e_t_a^2); sig2_b <- mean(e_b_a^2)
  y_t_bc <- y_t * (1 + sig2_t / 2); y_b_bc <- y_b * (1 + sig2_b / 2)
  e_t_bc <- y_a - y_t_bc; e_b_bc <- y_a - y_b_bc
  rmse_bc <- .rmse_block(e_t_bc, e_b_bc, h)

  list(
    n = length(y_a),
    rmse = rmse_lev$rmse, rmse_se = rmse_lev$rmse_se,
    rmse_ratio = rmse_lev$rmse_ratio,
    dm_stat = rmse_lev$dm_stat, dm_pval = rmse_lev$dm_pval,
    rmse_asinh = rmse_a$rmse, rmse_se_asinh = rmse_a$rmse_se,
    rmse_ratio_asinh = rmse_a$rmse_ratio,
    dm_stat_asinh = rmse_a$dm_stat, dm_pval_asinh = rmse_a$dm_pval,
    mae = mae_lev$mae, mae_se = mae_lev$mae_se,
    mae_ratio = mae_lev$mae_ratio,
    dm_mae_stat = mae_lev$dm_stat, dm_mae_pval = mae_lev$dm_pval,
    mae_asinh = mae_a$mae, mae_se_asinh = mae_a$mae_se,
    mae_ratio_asinh = mae_a$mae_ratio,
    dm_mae_stat_asinh = mae_a$dm_stat, dm_mae_pval_asinh = mae_a$dm_pval,
    rmse_bc = rmse_bc$rmse, rmse_bc_se = rmse_bc$rmse_se,
    rmse_bc_ratio = rmse_bc$rmse_ratio,
    dm_bc_stat = rmse_bc$dm_stat, dm_bc_pval = rmse_bc$dm_pval,
    sigma2_asinh_test = sig2_t, sigma2_asinh_bench = sig2_b
  )
}

# -----------------------------------------------------------------------------
# Main loop — country × horizon, R_MIDAS_legacy specs only
# -----------------------------------------------------------------------------
all_rows <- list()
for (cc in c("de", "it")) {
  for (h in HORIZONS) {
    ar_path  <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    leg_path <- file.path(FORECAST_DIR, sprintf("r_midas_legacy_%s_h%02d.csv", cc, h))
    if (!file.exists(ar_path)) {
      cat(sprintf("[skip] %s h=%02d : missing AR file  (%s)\n",
                  cc, h, basename(ar_path)))
      next
    }
    if (!file.exists(leg_path)) {
      cat(sprintf("[skip] %s h=%02d : missing legacy file (%s)\n",
                  cc, h, basename(leg_path)))
      next
    }
    ar  <- read.csv(ar_path,  stringsAsFactors = FALSE)
    leg <- read.csv(leg_path, stringsAsFactors = FALSE)
    if (!("ar_dum" %in% names(ar))) {
      cat(sprintf("[skip] %s h=%02d : AR file has no 'ar_dum' column\n", cc, h))
      next
    }
    ar$target_date  <- as.Date(ar$target_date)
    leg$target_date <- as.Date(leg$target_date)
    common <- intersect(as.character(ar$target_date),
                        as.character(leg$target_date))
    if (length(common) == 0L) {
      cat(sprintf("[skip] %s h=%02d : no overlap on target_date\n", cc, h))
      next
    }
    ar_a  <- ar [match(common, as.character(ar$target_date)),  , drop = FALSE]
    leg_a <- leg[match(common, as.character(leg$target_date)), , drop = FALSE]
    bench <- ar_a$ar_dum
    spec_cols <- grep("^spec_", colnames(leg_a), value = TRUE)
    cat(sprintf("[ok]  %s h=%02d : %d specs, n=%d common dates\n",
                cc, h, length(spec_cols), length(common)))
    for (sp in spec_cols) {
      m <- compute_metrics(leg_a$y_actual, leg_a[[sp]], bench, h)
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        country = cc, h = h, family = "R_MIDAS_legacy", model = sp,
        n = m$n,
        rmse = m$rmse, rmse_se = m$rmse_se, rmse_ratio = m$rmse_ratio,
        dm_stat = m$dm_stat, dm_pval = m$dm_pval,
        rmse_asinh = m$rmse_asinh, rmse_se_asinh = m$rmse_se_asinh,
        rmse_ratio_asinh = m$rmse_ratio_asinh,
        dm_stat_asinh = m$dm_stat_asinh, dm_pval_asinh = m$dm_pval_asinh,
        mae = m$mae, mae_se = m$mae_se, mae_ratio = m$mae_ratio,
        dm_mae_stat = m$dm_mae_stat, dm_mae_pval = m$dm_mae_pval,
        mae_asinh = m$mae_asinh, mae_se_asinh = m$mae_se_asinh,
        mae_ratio_asinh = m$mae_ratio_asinh,
        dm_mae_stat_asinh = m$dm_mae_stat_asinh,
        dm_mae_pval_asinh = m$dm_mae_pval_asinh,
        rmse_bc = m$rmse_bc, rmse_bc_se = m$rmse_bc_se,
        rmse_bc_ratio = m$rmse_bc_ratio,
        dm_bc_stat = m$dm_bc_stat, dm_bc_pval = m$dm_bc_pval,
        sigma2_asinh_test  = m$sigma2_asinh_test,
        sigma2_asinh_bench = m$sigma2_asinh_bench,
        stringsAsFactors = FALSE
      )
    }
  }
}

cat(sprintf("\ncollected %d rows total\n", length(all_rows)))

if (length(all_rows) == 0L) {
  stop("No legacy rows produced.\n",
       "Check that these files actually exist under FORECAST_DIR:\n",
       "  ar_<cc>_h<h>.csv\n",
       "  r_midas_legacy_<cc>_h<h>.csv\n",
       "If FORECAST_DIR above is wrong, set V2_DATA_DIR to the v2_new_dataset folder.")
}

summary_df <- do.call(rbind, all_rows)

# Attach hf_vars / lf_vars from the legacy dictionary
summary_df$hf_vars <- NA_character_
summary_df$lf_vars <- NA_character_
for (cc in c("de", "it")) {
  dp <- file.path(FORECAST_DIR,
                  sprintf("r_midas_legacy_spec_dictionary_%s.csv", cc))
  if (!file.exists(dp)) next
  d <- read.csv(dp, stringsAsFactors = FALSE)
  sel <- summary_df$country == cc
  idx <- match(summary_df$model[sel], d$spec_id)
  summary_df$hf_vars[sel] <- d$hf_vars[idx]
  summary_df$lf_vars[sel] <- d$lf_vars[idx]
}

# Sort
summary_df <- summary_df[order(summary_df$country, summary_df$h,
                               summary_df$rmse_ratio, na.last = TRUE), ]

# Write the two CSVs
out_path <- file.path(RESULT_DIR, "summary_legacy.csv")
write.csv(summary_df, out_path, row.names = FALSE)
cat(sprintf("\n[write] %s  (rows = %d)\n", out_path, nrow(summary_df)))

top5 <- do.call(rbind, by(summary_df, list(summary_df$country, summary_df$h),
                          function(g) {
                            g <- g[is.finite(g$rmse_ratio), , drop = FALSE]
                            if (nrow(g) == 0L) return(NULL)
                            g[order(g$rmse_ratio)[1:min(5, nrow(g))], , drop = FALSE]
                          }))
top_path <- file.path(RESULT_DIR, "summary_legacy_top5.csv")
write.csv(top5, top_path, row.names = FALSE)
cat(sprintf("[write] %s  (rows = %d)\n", top_path, nrow(top5)))
