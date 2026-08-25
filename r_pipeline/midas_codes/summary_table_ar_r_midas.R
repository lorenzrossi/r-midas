#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — SUMMARY TABLE (AR + R-MIDAS STANDARD ONLY)
#
# Same metrics as summary_table.R, restricted to:
#   - AR       (ar_<cc>_h<h>.csv)
#   - R_MIDAS  (r_midas_<cc>_h<h>.csv)
#
# Outputs:
#   results/summary_table_ar_r_midas.csv
#   results/summary_top5_ar_r_midas_per_country_h.csv
#
# Paths: V2_DATA_DIR or V2_SCRIPT_DIR (see v2_paths.R).
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
            "v2_paths.R")
} else file.path(getwd(), "v2_paths.R")
source(v2.boot.path)

code_dir  <- script_dir
data_root <- Sys.getenv("V2_DATA_DIR", unset = "")
if (!nzchar(data_root)) data_root <- Sys.getenv("V2_SCRIPT_DIR", unset = script_dir)
if (nzchar(data_root)) {
  data_root <- suppressWarnings(normalizePath(data_root, winslash = "/", mustWork = FALSE))
}
if (!nzchar(data_root) || is.na(data_root) || !dir.exists(data_root)) {
  data_root <- code_dir
}
if (dir.exists(data_root)) setwd(data_root)

cu <- file.path(code_dir, "common_utils.R")
if (!file.exists(cu)) cu <- file.path(data_root, "common_utils.R")
if (!file.exists(cu)) {
  stop("Cannot find common_utils.R. Set V2_DATA_DIR to the v2_new_dataset folder.")
}
source(cu)

FORECAST_DIR <- file.path(data_root, "forecasts")
RESULT_DIR   <- file.path(data_root, "results")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

.rmse_block <- function(e_t, e_b, h) {
  n <- length(e_t)
  if (n < 5L) return(list(rmse = NA_real_, rmse_se = NA_real_,
                          rmse_ratio = NA_real_,
                          dm_stat = NA_real_, dm_pval = NA_real_))
  rmse <- sqrt(mean(e_t^2))
  rmse_b <- sqrt(mean(e_b^2))
  v_e2 <- var(e_t^2)
  rmse_se <- if (rmse > 0) sqrt(v_e2 / (n * 4 * rmse^2)) else NA_real_
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

compute_metrics <- function(y_actual, y_hat_test, y_hat_bench, h,
                            y_hat_test_bc = NULL, y_hat_bench_bc = NULL) {
  ok <- is.finite(y_actual) & is.finite(y_hat_test) & is.finite(y_hat_bench)
  y_a <- y_actual[ok]; y_t <- y_hat_test[ok]; y_b <- y_hat_bench[ok]

  e_t_p <- y_a - y_t
  e_b_p <- y_a - y_b
  rmse_lev <- .rmse_block(e_t_p, e_b_p, h)
  mae_lev  <- .mae_block (e_t_p, e_b_p, h)

  ya_a <- asinh(y_a); ya_t <- asinh(y_t); ya_b <- asinh(y_b)
  e_t_a <- ya_a - ya_t
  e_b_a <- ya_a - ya_b
  rmse_a <- .rmse_block(e_t_a, e_b_a, h)
  mae_a  <- .mae_block (e_t_a, e_b_a, h)

  # Feasible bias-corrected forecasts are produced at each origin using
  # that rolling window's residual variance. Never estimate the correction
  # from the full evaluation sample.
  if (is.null(y_hat_test_bc) || is.null(y_hat_bench_bc)) {
    rmse_bc <- list(rmse = NA_real_, rmse_se = NA_real_,
                    rmse_ratio = NA_real_, dm_stat = NA_real_,
                    dm_pval = NA_real_)
  } else {
    ytbc0 <- y_hat_test_bc[ok]; ybbc0 <- y_hat_bench_bc[ok]
    okbc <- is.finite(y_a) & is.finite(ytbc0) & is.finite(ybbc0)
    e_t_bc <- y_a[okbc] - ytbc0[okbc]
    e_b_bc <- y_a[okbc] - ybbc0[okbc]
    rmse_bc <- .rmse_block(e_t_bc, e_b_bc, h)
  }
  sig2_t <- NA_real_
  sig2_b <- NA_real_

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
    sigma2_asinh_test  = sig2_t,
    sigma2_asinh_bench = sig2_b
  )
}

.metrics_row <- function(cc, h, family, sp, m) {
  data.frame(
    country = cc, h = h, family = family, model = sp,
    n = m$n,
    rmse = m$rmse, rmse_se = m$rmse_se,
    rmse_ratio = m$rmse_ratio,
    dm_stat = m$dm_stat, dm_pval = m$dm_pval,
    rmse_asinh = m$rmse_asinh, rmse_se_asinh = m$rmse_se_asinh,
    rmse_ratio_asinh = m$rmse_ratio_asinh,
    dm_stat_asinh = m$dm_stat_asinh, dm_pval_asinh = m$dm_pval_asinh,
    mae = m$mae, mae_se = m$mae_se,
    mae_ratio = m$mae_ratio,
    dm_mae_stat = m$dm_mae_stat, dm_mae_pval = m$dm_mae_pval,
    mae_asinh = m$mae_asinh, mae_se_asinh = m$mae_se_asinh,
    mae_ratio_asinh = m$mae_ratio_asinh,
    dm_mae_stat_asinh = m$dm_mae_stat_asinh, dm_mae_pval_asinh = m$dm_mae_pval_asinh,
    rmse_bc = m$rmse_bc, rmse_bc_se = m$rmse_bc_se,
    rmse_bc_ratio = m$rmse_bc_ratio,
    dm_bc_stat = m$dm_bc_stat, dm_bc_pval = m$dm_bc_pval,
    sigma2_asinh_test  = m$sigma2_asinh_test,
    sigma2_asinh_bench = m$sigma2_asinh_bench,
    stringsAsFactors = FALSE
  )
}

FAMILY_PATTERNS <- list(
  R_MIDAS = function(cc, h) sprintf("r_midas_%s_h%02d.csv", cc, h),
  # Real-time combinations written by combination_subperiod_gr.R (columns
  # comb_eq / comb_trim / comb_invmse_fam / comb_invmse).  Run the
  # combination stage BEFORE this summary so the files exist; the family is
  # skipped silently when they do not.
  COMB    = function(cc, h) sprintf("comb_r_midas_%s_h%02d.csv", cc, h)
)

DICT_FILES <- c(
  R_MIDAS = "r_midas_spec_dictionary_%s.csv"
)

process_country_h <- function(cc, h) {
  ar_path <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
  if (!file.exists(ar_path))
    stop("Missing AR forecasts: ", ar_path)

  ar <- read.csv(ar_path, stringsAsFactors = FALSE)
  ar$target_date <- as.Date(ar$target_date)
  ar$origin_date <- as.Date(ar$origin_date)
  bench <- ar$ar_dum

  rows <- list()

  ar_specs <- setdiff(colnames(ar), c("origin_date", "target_date", "y_actual",
                                      grep("__bc$", colnames(ar), value = TRUE)))
  for (sp in ar_specs) {
    m <- compute_metrics(ar$y_actual, ar[[sp]], bench, h,
                         ar[[paste0(sp, "__bc")]], ar[["ar_dum__bc"]])
    rows[[length(rows) + 1L]] <- .metrics_row(cc, h, "AR", sp, m)
  }

  for (fam in names(FAMILY_PATTERNS)) {
    rel <- FAMILY_PATTERNS[[fam]](cc, h)
    fp <- file.path(FORECAST_DIR, rel)
    if (!file.exists(fp)) next
    mi <- read.csv(fp, stringsAsFactors = FALSE)
    mi$target_date <- as.Date(mi$target_date)
    common <- intersect(as.character(ar$target_date),
                        as.character(mi$target_date))
    if (length(common) == 0L) {
      warning("No overlap for ", fam, " ", cc, " h=", h)
      next
    }
    mi_a <- mi[match(common, as.character(mi$target_date)), , drop = FALSE]
    ar_a <- ar[match(common, as.character(ar$target_date)), , drop = FALSE]
    bench_a <- ar_a$ar_dum
    col_prefix <- if (fam == "COMB") "^comb_" else "^spec_"
    spec_cols <- grep(col_prefix, colnames(mi_a), value = TRUE)
    spec_cols <- spec_cols[!grepl("__bc$", spec_cols)]
    for (sp in spec_cols) {
      m <- compute_metrics(mi_a$y_actual, mi_a[[sp]], bench_a, h,
                           mi_a[[paste0(sp, "__bc")]], ar_a[["ar_dum__bc"]])
      rows[[length(rows) + 1L]] <- .metrics_row(cc, h, fam, sp, m)
    }
  }

  do.call(rbind, rows)
}

all_rows <- list()
for (cc in c("de", "it")) {
  for (h in HORIZONS) {
    cat(sprintf("processing %s  h=%d ...\n", cc, h))
    block <- tryCatch(process_country_h(cc, h),
                      error = function(e) {
                        warning("Skipping ", cc, " h=", h, ": ", conditionMessage(e))
                        NULL
                      })
    if (!is.null(block)) all_rows[[length(all_rows) + 1L]] <- block
  }
}
if (length(all_rows) == 0L) stop("No forecast blocks processed — check forecasts/.")
summary_df <- do.call(rbind, all_rows)

summary_df$hf_vars <- NA_character_
summary_df$lf_vars <- NA_character_
for (cc in c("de", "it")) {
  for (fam in names(DICT_FILES)) {
    dp <- file.path(FORECAST_DIR, sprintf(DICT_FILES[[fam]], cc))
    if (!file.exists(dp)) next
    d <- read.csv(dp, stringsAsFactors = FALSE)
    sel <- summary_df$country == cc & summary_df$family == fam
    if (!any(sel)) next
    idx <- match(summary_df$model[sel], d$spec_id)
    summary_df$hf_vars[sel] <- d$hf_vars[idx]
    summary_df$lf_vars[sel] <- d$lf_vars[idx]
  }
}

summary_df <- summary_df[order(summary_df$country, summary_df$h,
                               summary_df$family, summary_df$rmse_ratio,
                               na.last = TRUE), ]

out_path <- file.path(RESULT_DIR, "summary_table_ar_r_midas.csv")
write.csv(summary_df, out_path, row.names = FALSE)
cat(sprintf("\nwrote %s  (rows = %d)\n", out_path, nrow(summary_df)))

best <- do.call(rbind, by(summary_df, list(summary_df$country, summary_df$h),
                          function(g) {
                            g <- g[is.finite(g$rmse_ratio), , drop = FALSE]
                            if (nrow(g) == 0L) return(NULL)
                            g[order(g$rmse_ratio)[1:min(5, nrow(g))], , drop = FALSE]
                          }))
best_path <- file.path(RESULT_DIR, "summary_top5_ar_r_midas_per_country_h.csv")
write.csv(best, best_path, row.names = FALSE)
cat(sprintf("wrote %s  (rows = %d)\n", best_path, nrow(best)))
