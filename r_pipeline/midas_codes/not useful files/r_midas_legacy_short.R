#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — R-MIDAS LEGACY (UNIFIED, ASYMMETRIC-WINDOW, RAW VARIABLES)
#                                                              (PARALLEL)
#
# ASYMMETRIC-WINDOW VARIANT.  This unified driver uses DIFFERENT rolling
# window sizes per country, chosen so that both countries' first feasible
# forecast origin lands around mid-January 2015 (same evaluation period):
#
#   Germany : 10-year (3650-day) window, full 2005-2025 daily sample
#   Italy   : 2-year (730-day) window, daily sample filtered to >= 2013-01-01
#
# Methodologically asymmetric — kept as a deliberate compromise and must be
# disclosed in any writeup.  Otherwise identical to the symmetric short-
# window variant: raw EUR/MWh y, raw HF/LF, 256-spec grid, 5-start NLS,
# in-window position decomposition of LF coefficients.
#
# Per-country logic is identical to r_midas_legacy_short_de.R and
# r_midas_legacy_short_it.R; this script is just a thin wrapper that loops
# over the COUNTRY_NAME setting.  Outputs use the same per-country
# filenames as the standalone scripts, so it is safe to run this driver
# AFTER (or instead of) running either of the two standalone scripts —
# files will be overwritten cleanly with the same content.
#
# Driver knobs (env vars):
#   COUNTRY_NAME = "Germany" | "Italy" | "ALL"   (default: "ALL")
#   SPEC_LIMIT   = integer                       (default: 256)
#   V2_DATA_DIR  = path to v2_new_dataset folder (if not running from inside)
#
# OUTPUTS (per country cc ∈ {de, it}):
#   forecasts/r_midas_legacy_short_<cc>_h<h>.csv
#   forecasts/r_midas_legacy_short_<cc>_spec_dictionary.csv
#   results/summary_legacy_short_<cc>.csv
#   results/summary_legacy_short_<cc>_top5.csv
#   results/coefficients/r_midas_legacy_short_<cc>.csv
#   results/position_decomposition_<cc>.csv
# =============================================================================

# -----------------------------------------------------------------------------
# Bootstrap
# -----------------------------------------------------------------------------
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

cu <- file.path(code_dir, "common_utils.R")
if (!file.exists(cu)) cu <- file.path(data_root, "common_utils.R")
if (!file.exists(cu)) stop("Cannot find common_utils.R.  Set V2_DATA_DIR.")
source(cu)

DATA_DIR     <- file.path(data_root, "data")
FORECAST_DIR <- file.path(data_root, "forecasts")
RESULT_DIR   <- file.path(data_root, "results")
COEF_DIR     <- file.path(RESULT_DIR, "coefficients")
dir.create(FORECAST_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(COEF_DIR,     showWarnings = FALSE, recursive = TRUE)

cat("code_dir    :", code_dir,  "\n")
cat("data_root   :", data_root, "\n")
cat("DATA_DIR    :", DATA_DIR,
    "  (exists =", dir.exists(DATA_DIR), ")\n\n")

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
# Asymmetric per-country settings: window_days + optional data_start cutoff
# applied to the DAILY series so that both countries' first forecast origin
# lands at roughly the same calendar date (~mid-Jan 2015).
COUNTRY_SETTINGS <- list(
  Germany = list(window_days = 3650L, data_start = NA_character_),         # 10y, full sample
  Italy   = list(window_days = 730L,  data_start = "2013-01-01")           #  2y, 2013-2025
)
Q_AR              <- 14L
K_POS             <- 28L
SPEC_LIMIT        <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L
theta_init   <- c(0, 0)
theta_lower  <- c(-8, -8)
theta_upper  <- c( 8,  0)
optim_method <- "L-BFGS-B"
optim_maxit  <- 3000L
optim_factr  <- 1e7
optim_starts <- 5L

# -----------------------------------------------------------------------------
# Legacy NLS helpers (shared across countries)
# -----------------------------------------------------------------------------
build_weighted_ar <- function(theta, ar_lags) {
  as.numeric(ar_lags %*% exp_almon_weights(ncol(ar_lags), theta))
}
concentrated_sse <- function(theta, y_vec, x_lin, ar_lags) {
  z_ar <- build_weighted_ar(theta, ar_lags)
  ok <- stats::complete.cases(cbind(x_lin, z_ar, y_vec))
  if (sum(ok) < ncol(x_lin) + 2L) return(1e12)
  X_full <- cbind(x_lin[ok, , drop = FALSE], z_ar = z_ar[ok])
  fit <- stats::lm.fit(X_full, y_vec[ok])
  sum(fit$residuals^2)
}
fit_window_model_legacy <- function(y_vec, x_lin, ar_lags, theta_start) {
  starts <- matrix(NA_real_, nrow = optim_starts, ncol = 2L)
  starts[1, ] <- if (length(theta_start) == 2L && all(is.finite(theta_start)))
                   theta_start else theta_init
  if (optim_starts > 1L) {
    for (s in 2:optim_starts)
      starts[s, ] <- runif(2L, min = theta_lower, max = theta_upper)
  }
  best <- NULL; best_val <- Inf
  for (s in 1:optim_starts) {
    cand <- tryCatch(
      stats::optim(par = as.numeric(starts[s, ]),
                   fn = function(th) concentrated_sse(th, y_vec, x_lin, ar_lags),
                   method = optim_method,
                   lower = theta_lower, upper = theta_upper,
                   control = list(maxit = optim_maxit, factr = optim_factr)),
      error = function(e) NULL
    )
    if (is.null(cand)) next
    if (!is.finite(cand$value) || cand$value >= 1e11) next
    if (!(cand$convergence %in% c(0L, 1L))) next
    if (cand$value < best_val) { best_val <- cand$value; best <- cand }
  }
  if (is.null(best)) return(NULL)
  theta_hat <- best$par
  z_ar <- build_weighted_ar(theta_hat, ar_lags)
  ok   <- stats::complete.cases(cbind(x_lin, z_ar, y_vec))
  if (sum(ok) < ncol(x_lin) + 2L) return(NULL)
  X_full <- cbind(x_lin[ok, , drop = FALSE], z_ar = z_ar[ok])
  fit <- stats::lm.fit(X_full, y_vec[ok])
  beta <- fit$coefficients
  if (anyNA(beta)) return(NULL)
  list(theta = theta_hat, beta = beta, beta_names = colnames(X_full),
       weights = exp_almon_weights(ncol(ar_lags), theta_hat),
       n_obs = sum(ok))
}
legacy_iterated_predict <- function(fit, x_lin_row, history_y_aug,
                                    target_dates, holiday_dates,
                                    h_max = H_MAX) {
  beta <- fit$beta; bnms <- fit$beta_names; w <- fit$weights
  Q <- length(w)
  z_ar_idx <- match("z_ar", bnms)
  delta_y  <- beta[z_ar_idx]
  beta_lin <- beta[-z_ar_idx]
  lin_names <- bnms[-z_ar_idx]
  names(beta_lin) <- lin_names
  y_hat_h <- numeric(h_max)
  for (s in 1:h_max) {
    row_lin <- refresh_calendar_row(x_lin_row, lin_names,
                                    target_dates[s], holiday_dates)
    lin_part <- sum(beta_lin * row_lin[lin_names])
    n_hist <- length(history_y_aug)
    ar_idx <- n_hist - (1:Q) + 1L
    if (any(ar_idx < 1L)) {
      y_hat_h[s] <- NA_real_
      history_y_aug <- c(history_y_aug, NA_real_)
      next
    }
    z_y <- sum(w * history_y_aug[ar_idx])
    yh  <- lin_part + delta_y * z_y
    y_hat_h[s] <- yh
    history_y_aug <- c(history_y_aug, yh)
  }
  y_hat_h
}

# -----------------------------------------------------------------------------
# Summary metric blocks (shared)
# -----------------------------------------------------------------------------
.rmse_block <- function(e_t, e_b, h) {
  n_ <- length(e_t)
  if (n_ < 5L) return(list(rmse=NA, rmse_se=NA, rmse_ratio=NA, dm_stat=NA, dm_pval=NA))
  rmse <- sqrt(mean(e_t^2)); rmse_b <- sqrt(mean(e_b^2))
  rmse_se <- if (rmse > 0) sqrt(var(e_t^2) / (n_ * 4 * rmse^2)) else NA
  rmse_ratio <- if (rmse_b > 0) rmse / rmse_b else NA
  dm <- dm_hln(e_t, e_b, h, loss = "sqerr")
  list(rmse=rmse, rmse_se=rmse_se, rmse_ratio=rmse_ratio,
       dm_stat=dm$stat, dm_pval=dm$pvalue)
}
.mae_block <- function(e_t, e_b, h) {
  n_ <- length(e_t)
  if (n_ < 5L) return(list(mae=NA, mae_se=NA, mae_ratio=NA, dm_stat=NA, dm_pval=NA))
  mae <- mean(abs(e_t)); mae_b <- mean(abs(e_b))
  mae_se <- sd(abs(e_t)) / sqrt(n_)
  mae_ratio <- if (mae_b > 0) mae / mae_b else NA
  dm <- dm_hln(e_t, e_b, h, loss = "abserr")
  list(mae=mae, mae_se=mae_se, mae_ratio=mae_ratio,
       dm_stat=dm$stat, dm_pval=dm$pvalue)
}
compute_metrics_local <- function(y_actual, y_hat_test, y_hat_bench, h) {
  ok <- is.finite(y_actual) & is.finite(y_hat_test) & is.finite(y_hat_bench)
  y_a <- y_actual[ok]; y_t <- y_hat_test[ok]; y_b <- y_hat_bench[ok]
  e_t_p <- y_a - y_t; e_b_p <- y_a - y_b
  rmse_lev <- .rmse_block(e_t_p, e_b_p, h)
  mae_lev  <- .mae_block (e_t_p, e_b_p, h)
  ya_a <- asinh(y_a); ya_t <- asinh(y_t); ya_b <- asinh(y_b)
  e_t_a <- ya_a - ya_t; e_b_a <- ya_a - ya_b
  rmse_a <- .rmse_block(e_t_a, e_b_a, h)
  mae_a  <- .mae_block (e_t_a, e_b_a, h)
  list(n = length(y_a),
       rmse = rmse_lev$rmse, rmse_se = rmse_lev$rmse_se,
       rmse_ratio = rmse_lev$rmse_ratio,
       dm_stat = rmse_lev$dm_stat, dm_pval = rmse_lev$dm_pval,
       rmse_asinh = rmse_a$rmse, rmse_ratio_asinh = rmse_a$rmse_ratio,
       dm_stat_asinh = rmse_a$dm_stat, dm_pval_asinh = rmse_a$dm_pval,
       mae = mae_lev$mae, mae_ratio = mae_lev$mae_ratio,
       dm_mae_stat = mae_lev$dm_stat, dm_mae_pval = mae_lev$dm_pval,
       mae_asinh = mae_a$mae, mae_ratio_asinh = mae_a$mae_ratio)
}

# -----------------------------------------------------------------------------
# Per-country runner
# -----------------------------------------------------------------------------
run_country <- function(country) {
  stopifnot(country %in% c("Germany", "Italy"))
  cc <- if (country == "Germany") "de" else "it"

  # Per-country asymmetric settings
  settings    <- COUNTRY_SETTINGS[[country]]
  window_days <- settings$window_days
  data_start  <- if (!is.na(settings$data_start)) as.Date(settings$data_start) else NA

  lf_names <- if (country == "Germany")
    c("pmi_mfg_de", "pmi_serv_de", "ifo_bci_de",
      "Consumer_Goods", "Manufacturing", "Energy")
  else
    c("pmi_mfg_it", "pmi_serv_it", "istat_bci_it",
      "Consumer_Goods", "Manufacturing", "Energy")
  daily_path   <- file.path(DATA_DIR, sprintf("prices_%s_daily.csv", cc))
  monthly_path <- file.path(DATA_DIR, sprintf("dataset_%s_m.csv",   cc))
  if (!file.exists(daily_path))   stop("Missing daily file: ",   daily_path)
  if (!file.exists(monthly_path)) stop("Missing monthly file: ", monthly_path)
  daily   <- read.csv(daily_path,   stringsAsFactors = FALSE)
  monthly <- read.csv(monthly_path, stringsAsFactors = FALSE)
  daily$date   <- as.Date(daily$date)
  monthly$date <- as.Date(monthly$date)
  daily   <- daily[order(daily$date), , drop = FALSE]
  monthly <- monthly[order(monthly$date), , drop = FALSE]
  names(monthly)[names(monthly) == "Electricity_Gas_Steam"] <- "Energy"
  monthly$ref_month <- as.Date(format(monthly$date, "%Y-%m-01"))

  # Asymmetric variant: optional daily-data filter (Italy: >= 2013-01-01)
  if (!is.na(data_start)) {
    n_before <- nrow(daily)
    daily <- daily[daily$date >= data_start, , drop = FALSE]
    cat(sprintf("%s daily filtered to >= %s : %d rows (was %d)\n",
                country, format(data_start), nrow(daily), n_before))
  }
  cat(sprintf("%s window_days = %d\n", country, window_days))
  dates <- daily$date
  n     <- nrow(daily)
  y     <- as.numeric(daily$elec_price)
  gas_lag1   <- c(NA_real_, as.numeric(daily$gas_price)[1:(n - 1L)])
  brent_lag1 <- c(NA_real_, as.numeric(daily$brent_price)[1:(n - 1L)])

  cat(sprintf("\n=== %s : y range [%.2f, %.2f]   mean=%.2f   n=%d ===\n",
              country,
              min(y, na.rm=TRUE), max(y, na.rm=TRUE),
              mean(y, na.rm=TRUE), n))

  ar_lag_mat <- matrix(NA_real_, nrow = n, ncol = Q_AR)
  colnames(ar_lag_mat) <- paste0("y_lag", 1:Q_AR)
  for (j in 1:Q_AR) ar_lag_mat[(j + 1L):n, j] <- y[1:(n - j)]

  lf_pub <- c(pmi_mfg_de = PUB_LAG_SURVEY, pmi_serv_de = PUB_LAG_SURVEY,
              ifo_bci_de = PUB_LAG_SURVEY,
              pmi_mfg_it = PUB_LAG_SURVEY, pmi_serv_it = PUB_LAG_SURVEY,
              istat_bci_it = PUB_LAG_SURVEY,
              Consumer_Goods = PUB_LAG_IPI, Manufacturing = PUB_LAG_IPI,
              Energy = PUB_LAG_IPI)
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(dates, monthly, v, lf_pub[[v]])

  cal_mat   <- build_calendar_dummies(dates, country)
  X_base    <- cbind(intercept = 1, cal_mat)
  base_cols <- colnames(X_base)
  hol_set   <- country_holidays(unique(as.integer(format(dates, "%Y"))), country)
  within_month_pos <- as.integer(format(dates, "%d"))
  within_month_pos[within_month_pos > K_POS] <- NA_integer_

  hf_names   <- c("gas", "brent")
  flag_names <- c(hf_names, lf_names)
  masks <- expand.grid(replicate(length(flag_names), 0:1, simplify = FALSE),
                       KEEP.OUT.ATTRS = FALSE)
  masks <- as.matrix(masks); colnames(masks) <- flag_names
  masks_int <- as.integer(masks %*% (2^(rev(seq_len(ncol(masks))) - 1L)))
  ord <- order(masks_int)
  masks <- masks[ord, , drop = FALSE]
  spec_ids <- sprintf("spec_%03d", 0:(nrow(masks) - 1L))
  rownames(masks) <- spec_ids
  if (SPEC_LIMIT < nrow(masks)) {
    masks <- masks[1:SPEC_LIMIT, , drop = FALSE]
    spec_ids <- spec_ids[1:SPEC_LIMIT]
  }

  dict <- data.frame(
    spec_id = spec_ids,
    n_hf    = rowSums(masks[, hf_names, drop = FALSE]),
    n_lf    = rowSums(masks[, lf_names, drop = FALSE]),
    hf_vars = vapply(seq_len(nrow(masks)), function(i) {
      v <- hf_names[masks[i, hf_names] == 1L]
      if (length(v) == 0L) "" else paste(v, collapse = "|")
    }, character(1)),
    lf_vars = vapply(seq_len(nrow(masks)), function(i) {
      v <- lf_names[masks[i, lf_names] == 1L]
      if (length(v) == 0L) "" else paste(v, collapse = "|")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  dict_path <- file.path(FORECAST_DIR,
                         sprintf("r_midas_legacy_short_%s_spec_dictionary.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)

  origins <- make_origins(n, window_days,
                          pre_lag_max = max(c(AR_LAGS, Q_HF, Q_AR)),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  coef_cols <- c("delta_y", lf_names)
  fit_one_spec <- function(si) {
    set.seed(20240501L + si)
    sid <- spec_ids[si]; mk <- masks[si, ]
    hf_used <- character(0)
    if (mk["gas"]   == 1L) hf_used <- c(hf_used, "gas")
    if (mk["brent"] == 1L) hf_used <- c(hf_used, "brent")
    lf_used <- lf_names[mk[lf_names] == 1L]
    theta_warm <- NULL
    yhat_mat <- matrix(NA_real_, n_orig, H_MAX)
    coef_mat <- matrix(NA_real_, n_orig, length(coef_cols),
                       dimnames = list(NULL, coef_cols))
    pos_sum   <- if (length(lf_used) > 0L)
      matrix(0, length(lf_used), K_POS,
             dimnames = list(lf_used, paste0("D", sprintf("%02d", 1:K_POS))))
    else NULL
    pos_count <- if (length(lf_used) > 0L)
      matrix(0L, length(lf_used), K_POS,
             dimnames = list(lf_used, paste0("D", sprintf("%02d", 1:K_POS))))
    else NULL

    for (oi in seq_along(origins)) {
      t0 <- origins[oi]
      win_lo <- t0 - window_days + 1L
      win_rows <- win_lo:t0
      ywin <- y[win_rows]
      ar_lags_win <- ar_lag_mat[win_rows, , drop = FALSE]
      lin_cols  <- list(X_base[win_rows, , drop = FALSE])
      lin_names <- base_cols
      if (length(lf_used) > 0L) {
        lin_cols[[length(lin_cols) + 1L]] <- lf_mat[win_rows, lf_used, drop = FALSE]
        lin_names <- c(lin_names, lf_used)
      }
      if ("gas" %in% hf_used) {
        lin_cols[[length(lin_cols) + 1L]] <- matrix(gas_lag1[win_rows], ncol = 1,
                                                    dimnames = list(NULL, "gas_lag1"))
        lin_names <- c(lin_names, "gas_lag1")
      }
      if ("brent" %in% hf_used) {
        lin_cols[[length(lin_cols) + 1L]] <- matrix(brent_lag1[win_rows], ncol = 1,
                                                    dimnames = list(NULL, "brent_lag1"))
        lin_names <- c(lin_names, "brent_lag1")
      }
      x_lin_win <- do.call(cbind, lin_cols); colnames(x_lin_win) <- lin_names
      keep_train <- stats::complete.cases(x_lin_win) &
                    apply(ar_lags_win, 1, function(r) all(is.finite(r))) &
                    is.finite(ywin)
      if (sum(keep_train) < ncol(x_lin_win) + 3L) next
      x_lin_fit <- x_lin_win[keep_train, , drop = FALSE]
      ar_lags_fit <- ar_lags_win[keep_train, , drop = FALSE]
      yfit <- ywin[keep_train]
      if (!is.null(theta_warm) && length(theta_warm) != 2L) theta_warm <- NULL
      fit <- fit_window_model_legacy(yfit, x_lin_fit, ar_lags_fit, theta_warm)
      if (is.null(fit)) next
      theta_warm <- fit$theta
      coef_mat[oi, "delta_y"] <- fit$beta["z_ar"]
      if (length(lf_used) > 0L) {
        for (v in lf_used) {
          b <- fit$beta[v]
          if (!is.na(b)) coef_mat[oi, v] <- as.numeric(b)
        }
      }
      # Position decomposition
      if (length(lf_used) > 0L) {
        w_y <- fit$weights
        z_y_win <- as.numeric(ar_lags_win %*% w_y)
        pos_win <- within_month_pos[win_rows]
        lf_pos_blk <- matrix(0, length(win_rows), length(lf_used) * K_POS)
        pos_col_names <- character(length(lf_used) * K_POS)
        idx <- 1L
        for (v in lf_used) for (i in 1:K_POS) {
          lf_pos_blk[, idx] <- lf_mat[win_rows, v] *
                               as.integer(!is.na(pos_win) & pos_win == i)
          pos_col_names[idx] <- sprintf("lfpos__%s__D%02d", v, i)
          idx <- idx + 1L
        }
        colnames(lf_pos_blk) <- pos_col_names
        hf_blk <- list()
        if ("gas" %in% hf_used)   hf_blk[["gas_lag1"]]   <- gas_lag1[win_rows]
        if ("brent" %in% hf_used) hf_blk[["brent_lag1"]] <- brent_lag1[win_rows]
        X_pos <- cbind(X_base[win_rows, , drop = FALSE], z_y = z_y_win,
                       if (length(hf_blk) > 0L) do.call(cbind, hf_blk) else NULL,
                       lf_pos_blk)
        keep_pos <- !is.na(pos_win) & stats::complete.cases(X_pos) & is.finite(ywin)
        if (sum(keep_pos) > ncol(X_pos) + 5L) {
          qrfit <- qr(X_pos[keep_pos, , drop = FALSE])
          if (qrfit$rank >= ncol(X_pos)) {
            beta_pos <- qr.coef(qrfit, ywin[keep_pos])
            names(beta_pos) <- colnames(X_pos)
            if (!anyNA(beta_pos)) {
              for (v in lf_used) for (i in 1:K_POS) {
                cn <- sprintf("lfpos__%s__D%02d", v, i)
                b  <- beta_pos[cn]
                if (!is.na(b) && is.finite(b)) {
                  pos_sum[v, i]   <- pos_sum[v, i] + as.numeric(b)
                  pos_count[v, i] <- pos_count[v, i] + 1L
                }
              }
            }
          }
        }
      }
      # Iterated forecast
      lin_row_parts <- list(X_base[t0 + 1L, , drop = TRUE])
      lin_row_names <- base_cols
      if (length(lf_used) > 0L) {
        lin_row_parts[[length(lin_row_parts) + 1L]] <- lf_mat[t0 + 1L, lf_used, drop = TRUE]
        lin_row_names <- c(lin_row_names, lf_used)
      }
      if ("gas" %in% hf_used) {
        lin_row_parts[[length(lin_row_parts) + 1L]] <- gas_lag1[t0 + 1L]
        lin_row_names <- c(lin_row_names, "gas_lag1")
      }
      if ("brent" %in% hf_used) {
        lin_row_parts[[length(lin_row_parts) + 1L]] <- brent_lag1[t0 + 1L]
        lin_row_names <- c(lin_row_names, "brent_lag1")
      }
      x_lin_row <- unlist(lin_row_parts)
      names(x_lin_row) <- lin_row_names
      if (anyNA(x_lin_row)) next
      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      yhat_mat[oi, ] <- legacy_iterated_predict(
        fit = fit, x_lin_row = x_lin_row,
        history_y_aug = y[1:t0],
        target_dates = target_dates,
        holiday_dates = hol_set,
        h_max = H_MAX
      )
    }
    list(yhat_mat = yhat_mat, coef_mat = coef_mat,
         pos_sum = pos_sum, pos_count = pos_count, lf_used = lf_used)
  }

  n_cores <- get_n_cores()
  cat(sprintf("running %s — %d specs × %d origins on %d cores\n",
              country, length(spec_ids), n_orig, n_cores))
  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec,
                                  n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("fitting %s done in %.1fs (%.2fh)\n", country, el, el / 3600))

  # Forecast CSVs
  for (h in HORIZONS) {
    yhat_h_mat <- vapply(seq_along(spec_ids),
                         function(si) spec_results[[si]]$yhat_mat[, h],
                         numeric(n_orig))
    colnames(yhat_h_mat) <- spec_ids
    out <- data.frame(origin_date = dates[origins],
                      target_date = dates[origins + h],
                      y_actual    = y[origins + h],
                      stringsAsFactors = FALSE)
    out <- cbind(out, as.data.frame(yhat_h_mat, stringsAsFactors = FALSE))
    out <- out[!is.na(out$y_actual), , drop = FALSE]
    fp <- file.path(FORECAST_DIR,
                    sprintf("r_midas_legacy_short_%s_h%02d.csv", cc, h))
    write.csv(out, fp, row.names = FALSE)
    cat(sprintf("wrote %s  (rows = %d)\n", fp, nrow(out)))
  }
  # Rolling pooled λ
  coef_blocks <- vector("list", length(spec_ids))
  for (si in seq_along(spec_ids)) {
    cm <- spec_results[[si]]$coef_mat
    blk <- data.frame(origin_date = dates[origins], spec_id = spec_ids[si],
                      stringsAsFactors = FALSE)
    for (col in coef_cols) blk[[col]] <- cm[, col]
    coef_blocks[[si]] <- blk
  }
  coef_df <- do.call(rbind, coef_blocks)
  all_na <- rowSums(is.na(coef_df[, coef_cols, drop = FALSE])) == length(coef_cols)
  coef_df <- coef_df[!all_na, , drop = FALSE]
  coef_path <- file.path(COEF_DIR, sprintf("r_midas_legacy_short_%s.csv", cc))
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat(sprintf("wrote %s  (rows = %d)\n", coef_path, nrow(coef_df)))
  # Position decomposition
  pos_rows <- list()
  for (si in seq_along(spec_ids)) {
    r <- spec_results[[si]]
    if (is.null(r$pos_sum) || is.null(r$lf_used) || length(r$lf_used) == 0L) next
    for (v in r$lf_used) for (i in 1:K_POS) {
      cnt <- r$pos_count[v, i]
      pos_rows[[length(pos_rows) + 1L]] <- data.frame(
        spec_id = spec_ids[si], lf_var = v, position = i,
        mean_lambda = if (cnt > 0L) r$pos_sum[v, i] / cnt else NA_real_,
        n_origins = cnt, stringsAsFactors = FALSE
      )
    }
  }
  if (length(pos_rows) > 0L) {
    pos_df <- do.call(rbind, pos_rows)
    pos_path <- file.path(RESULT_DIR, sprintf("position_decomposition_%s.csv", cc))
    write.csv(pos_df, pos_path, row.names = FALSE)
    cat(sprintf("wrote %s  (rows = %d)\n", pos_path, nrow(pos_df)))
  }
  # Summary table
  sum_rows <- list()
  for (h in HORIZONS) {
    ar_path  <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    leg_path <- file.path(FORECAST_DIR, sprintf("r_midas_legacy_short_%s_h%02d.csv", cc, h))
    if (!file.exists(ar_path) || !file.exists(leg_path)) next
    ar  <- read.csv(ar_path,  stringsAsFactors = FALSE)
    leg <- read.csv(leg_path, stringsAsFactors = FALSE)
    ar$target_date  <- as.Date(ar$target_date)
    leg$target_date <- as.Date(leg$target_date)
    common <- intersect(as.character(ar$target_date), as.character(leg$target_date))
    if (length(common) == 0L) next
    ar_a  <- ar [match(common, as.character(ar$target_date)),  , drop = FALSE]
    leg_a <- leg[match(common, as.character(leg$target_date)), , drop = FALSE]
    bench <- ar_a$ar_dum
    spec_cols_h <- grep("^spec_", colnames(leg_a), value = TRUE)
    for (sp in spec_cols_h) {
      m <- compute_metrics_local(leg_a$y_actual, leg_a[[sp]], bench, h)
      sum_rows[[length(sum_rows) + 1L]] <- data.frame(
        country = cc, h = h, family = "R_MIDAS_legacy_short", model = sp,
        n = m$n,
        rmse = m$rmse, rmse_se = m$rmse_se, rmse_ratio = m$rmse_ratio,
        dm_stat = m$dm_stat, dm_pval = m$dm_pval,
        rmse_asinh = m$rmse_asinh, rmse_ratio_asinh = m$rmse_ratio_asinh,
        dm_stat_asinh = m$dm_stat_asinh, dm_pval_asinh = m$dm_pval_asinh,
        mae = m$mae, mae_ratio = m$mae_ratio,
        dm_mae_stat = m$dm_mae_stat, dm_mae_pval = m$dm_mae_pval,
        mae_asinh = m$mae_asinh, mae_ratio_asinh = m$mae_ratio_asinh,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(sum_rows) > 0L) {
    summary_df <- do.call(rbind, sum_rows)
    d <- read.csv(dict_path, stringsAsFactors = FALSE)
    idx <- match(summary_df$model, d$spec_id)
    summary_df$hf_vars <- d$hf_vars[idx]
    summary_df$lf_vars <- d$lf_vars[idx]
    summary_df <- summary_df[order(summary_df$country, summary_df$h,
                                   summary_df$rmse_ratio, na.last = TRUE), ]
    sum_path <- file.path(RESULT_DIR, sprintf("summary_legacy_short_%s.csv", cc))
    write.csv(summary_df, sum_path, row.names = FALSE)
    cat(sprintf("wrote %s  (rows = %d)\n", sum_path, nrow(summary_df)))
    top5 <- do.call(rbind, by(summary_df, list(summary_df$country, summary_df$h),
                              function(g) {
                                g <- g[is.finite(g$rmse_ratio), , drop = FALSE]
                                if (nrow(g) == 0L) return(NULL)
                                g[order(g$rmse_ratio)[1:min(5, nrow(g))], , drop = FALSE]
                              }))
    top_path <- file.path(RESULT_DIR, sprintf("summary_legacy_short_%s_top5.csv", cc))
    write.csv(top5, top_path, row.names = FALSE)
    cat(sprintf("wrote %s  (rows = %d)\n", top_path, nrow(top5)))
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------
country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
cat("Will run countries:", paste(countries, collapse = ", "), "\n")
for (cc_ in countries) run_country(cc_)
cat("\nR-MIDAS legacy short-window unified driver done.\n")
