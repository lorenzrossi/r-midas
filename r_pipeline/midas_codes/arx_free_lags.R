#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — ARX FREE-LAG CONTROL FAMILY  (PARALLEL via mclapply)
#
# Same 256 specifications (fuel lag-1 x monthly LF indicators), same
# information set, same rolling window and the same iterated multi-step
# scheme as r_midas.R — but the four AR coefficients on lags {1, 2, 3, 7}
# are FREE and the whole model is estimated by OLS.  Family: "ARX".
#
# Purpose (referee / TODO item 6): with four lags the exp-Almon block saves
# only one parameter (delta, theta1, theta2 instead of four phi_j) and forces
# every implied phi_j to share the sign of delta.  Comparing ARX with R_MIDAS
# spec by spec separates the forecasting value of the monthly indicators
# from the effect of the lag restriction.  Consistency checks: spec_000
# coincides with ar_dum, the gas-only spec with ar_dum_gas, etc.
#
# Design per spec (all linear, OLS on the WINDOW_DAYS rolling window):
#   intercept + 18 calendar dummies + y_lag{1,2,3,7}
#   [+ one release-aware daily column per LF variable]
#   [+ gas_lag1] [+ brent_lag1]
# Multi-step: LF frozen at the origin row t0 (latest released value), fuel
# lag-1 frozen at row t0+1 (price observed at t0), calendar refreshed at
# each target date, AR cells fed with realised / previous-step values.
#
# OUTPUTS:
#   forecasts/arx_<cc>_h<h>.csv             point forecasts (+ __bc twins)
#   forecasts/arx_<cc>_h<h>_sd.csv          predictive sd (asinh scale) per spec
#   forecasts/arx_spec_dictionary_<cc>.csv  spec -> variable membership
#   results/coefficients/arx_<cc>.csv       rolling phi_1, phi_2, phi_3, phi_7,
#                                           sum_phi, min_root_mod, <each LF>
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

DATA_DIR     <- file.path(script_dir, "data")
FORECAST_DIR <- file.path(script_dir, "forecasts")
COEF_DIR     <- file.path(script_dir, "results", "coefficients")
dir.create(FORECAST_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(COEF_DIR,     showWarnings = FALSE, recursive = TRUE)

SPEC_LIMIT <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_arx <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== ARX free-lag control (OLS, AR lags {1,2,3,7} free) — %s (n_cores=%d) ===\n",
              country, n_cores))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  P       <- daily$elec_price              # RAW prices; normalised asinh per window

  ar_mat  <- build_ar_lags(P, AR_LAGS)     # raw lagged prices, transformed per origin
  ar_cols <- colnames(ar_mat)

  gas_lag1 <- c(NA_real_, daily$gas_log[1:(n - 1L)])
  brt_lag1 <- c(NA_real_, daily$brent_log[1:(n - 1L)])

  # LF (daily, release-aware; identical construction to r_midas.R)
  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- LF_PUB_LAG
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(
    dates, monthly, v, lf_pub[[v]],
    locf = (lf_pub[[v]] == PUB_LAG_SURVEY))

  cal_mat   <- build_calendar_dummies(dates, country)
  X_base    <- cbind(intercept = 1, cal_mat)
  base_cols <- colnames(X_base)
  hol_set   <- country_holidays(unique(as.integer(format(dates, "%Y"))), country)

  # Full regressor pool; each spec selects its columns.
  X_all <- cbind(X_base, ar_mat, lf_mat,
                 gas_lag1 = gas_lag1, brent_lag1 = brt_lag1)

  hf_names   <- c("gas", "brent")
  flag_names <- c(hf_names, lf_names)

  # Same spec enumeration / ordering as r_midas.R, so spec ids coincide.
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
    cat(sprintf("  SPEC_LIMIT=%d (truncating from 256)\n", SPEC_LIMIT))
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
  dict_path <- file.path(FORECAST_DIR, sprintf("arx_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS, pre_lag_max = max(AR_LAGS), h_max = H_MAX)
  n_orig  <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]), format(dates[origins[n_orig]])))

  # Normalised-asinh location / scale per origin (same as ar_benchmark.R).
  VP <- vst_params_by_origin(P, origins, WINDOW_DAYS)
  write.csv(data.frame(origin_date = dates[origins], a = VP$a, b = VP$b),
            file.path(FORECAST_DIR, sprintf("vst_params_%s.csv", cc)), row.names = FALSE)

  phi_cols  <- paste0("phi_", AR_LAGS)
  coef_cols <- c(phi_cols, "sum_phi", "min_root_mod", lf_names)

  fit_one_spec <- function(si) {
    mk <- masks[si, ]
    hf_used <- c(if (mk["gas"] == 1L) "gas_lag1", if (mk["brent"] == 1L) "brent_lag1")
    lf_used <- lf_names[mk[lf_names] == 1L]
    cols    <- c(base_cols, ar_cols, lf_used, hf_used)
    ar_idx  <- match(ar_cols, cols)
    lf_idx  <- match(lf_used, cols)

    yhat_mat   <- matrix(NA_real_, n_orig, H_MAX)
    sigma2_mat <- matrix(NA_real_, n_orig, H_MAX)
    coef_mat   <- matrix(NA_real_, n_orig, length(coef_cols),
                         dimnames = list(NULL, coef_cols))

    for (oi in seq_along(origins)) {
      t0     <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      p      <- list(a = VP$a[oi], b = VP$b[oi])
      Xtr    <- X_all[win_lo:t0, cols, drop = FALSE]
      Xtr[, ar_idx] <- vst(Xtr[, ar_idx, drop = FALSE], p)   # lags on the window's scale
      ytr    <- vst(P[win_lo:t0], p)
      ok     <- rows_all_finite(Xtr) & is.finite(ytr)
      if (sum(ok) < ncol(Xtr) + 5L) next
      qrfit  <- qr(Xtr[ok, , drop = FALSE])
      if (qrfit$rank < ncol(Xtr)) next
      beta   <- qr.coef(qrfit, ytr[ok])
      names(beta) <- cols
      res    <- ytr[ok] - drop(Xtr[ok, , drop = FALSE] %*% beta)
      phi    <- as.numeric(beta[ar_cols])
      sigma2_mat[oi, ] <- iterated_error_variance(AR_LAGS, phi, mean(res^2), H_MAX)

      coef_mat[oi, phi_cols]       <- phi
      coef_mat[oi, "sum_phi"]      <- sum(phi)
      coef_mat[oi, "min_root_mod"] <- ar_min_root_modulus(AR_LAGS, phi)
      if (length(lf_used) > 0L) coef_mat[oi, lf_used] <- as.numeric(beta[lf_used])

      # Forecast row at the origin's information set: calendar at t0+1
      # (refreshed per step anyway), LF frozen at row t0, fuel lag-1 at row
      # t0+1 (= price observed at t0); AR cells are filled by the iterator.
      frow <- X_all[t0 + 1L, cols, drop = TRUE]
      if (length(lf_idx) > 0L) frow[lf_idx] <- lf_mat[t0, lf_used]
      frow[ar_idx] <- 0
      if (anyNA(frow)) next

      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      yhat_mat[oi, ] <- iterated_forecast(
        coef          = beta,
        colnames_     = cols,
        history_y_aug = vst(P[1:t0], p),
        frozen_row    = frow,
        ar_lag_idx    = ar_idx,
        ar_lags       = AR_LAGS,
        mon_idx       = NULL,
        target_dates  = target_dates,
        h_max         = H_MAX,
        holiday_dates = hol_set
      )
    }
    list(yhat_mat = yhat_mat, coef_mat = coef_mat, sigma2 = sigma2_mat)
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec, n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.1fmin) — wall-clock\n", el, el / 60))

  for (h in HORIZONS) {
    yhat_h_mat <- vapply(seq_along(spec_ids),
                         function(si) vst_inv(spec_results[[si]]$yhat_mat[, h], VP),
                         numeric(n_orig))
    colnames(yhat_h_mat) <- spec_ids
    yhat_bc_mat <- vapply(seq_along(spec_ids), function(si) {
      mu <- spec_results[[si]]$yhat_mat[, h]
      s2 <- spec_results[[si]]$sigma2[, h]
      vst_inv(mu, VP, s2)                       # mean forecast, EUR/MWh
    }, numeric(n_orig))
    colnames(yhat_bc_mat) <- paste0(spec_ids, "__bc")
    sd_mat <- vapply(seq_along(spec_ids),
                     function(si) sqrt(spec_results[[si]]$sigma2[, h]),
                     numeric(n_orig))
    colnames(sd_mat) <- spec_ids

    out <- data.frame(
      origin_date = dates[origins],
      target_date = dates[origins + h],
      y_actual    = daily$elec_price[origins + h],
      stringsAsFactors = FALSE
    )
    keep_rows <- !is.na(out$y_actual)
    sd_out <- cbind(out[, c("origin_date", "target_date")],
                    as.data.frame(sd_mat, stringsAsFactors = FALSE))
    out <- cbind(out, as.data.frame(yhat_h_mat, stringsAsFactors = FALSE),
                 as.data.frame(yhat_bc_mat, stringsAsFactors = FALSE))
    out    <- out[keep_rows, , drop = FALSE]
    sd_out <- sd_out[keep_rows, , drop = FALSE]
    fp <- file.path(FORECAST_DIR, sprintf("arx_%s_h%02d.csv", cc, h))
    write.csv(out, fp, row.names = FALSE)
    write.csv(sd_out, file.path(FORECAST_DIR, sprintf("arx_%s_h%02d_sd.csv", cc, h)),
              row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", fp, nrow(out)))
  }

  coef_blocks <- vector("list", length(spec_ids))
  for (si in seq_along(spec_ids)) {
    cm  <- spec_results[[si]]$coef_mat
    blk <- data.frame(origin_date = dates[origins], spec_id = spec_ids[si],
                      stringsAsFactors = FALSE)
    for (col in coef_cols) blk[[col]] <- cm[, col]
    coef_blocks[[si]] <- blk
  }
  coef_df <- do.call(rbind, coef_blocks)
  all_na  <- rowSums(is.na(coef_df[, coef_cols, drop = FALSE])) == length(coef_cols)
  coef_df <- coef_df[!all_na, , drop = FALSE]
  coef_path <- file.path(COEF_DIR, sprintf("arx_%s.csv", cc))
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d, cols = %d)\n",
              coef_path, nrow(coef_df), ncol(coef_df)))
  cat(sprintf("  free AR coefficients, median over windows and specs: phi(1,2,3,7) = (%s); sum = %.3f; min root modulus < 1 in %.1f%% of windows\n",
              paste(sprintf("%.3f", vapply(phi_cols, function(cn) median(coef_df[[cn]], na.rm = TRUE), numeric(1))),
                    collapse = ", "),
              median(coef_df$sum_phi, na.rm = TRUE),
              100 * mean(coef_df$min_root_mod < 1, na.rm = TRUE)))
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_arx(cc_)
cat("\nARX free-lag control family done.\n")
