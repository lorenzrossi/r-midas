#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — AR(1,2,3,7) BENCHMARK FAMILY  (PARALLEL via mclapply)
#
# Four specifications per country and per horizon h ∈ {1, 7, 14, 21, 28}:
#   ar_dum            : AR(1,2,3,7) + calendar dummies                              [BENCHMARK]
#   ar_dum_gas        : ar_dum + lag 1 of log(gas)
#   ar_dum_brent      : ar_dum + lag 1 of log(brent)
#   ar_dum_gas_brent  : ar_dum + lag 1 of log(gas) and log(brent)
#
# Calendar control set (built by build_calendar_dummies in common_utils.R):
#   * 6 weekday dummies (Sun = reference)
#   * 11 month-of-year dummies (Jan = reference)
#   * 1 holiday dummy (country-specific national holidays)
#
# Estimation: OLS on the WINDOW_DAYS rolling window.
# Multi-step: iterated, option c.1, with calendar refresh at the target row.
# Forecasts back-transformed via sinh(.) to euros/MWh.
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

DATA_DIR      <- file.path(script_dir, "data")
FORECAST_DIR  <- file.path(script_dir, "forecasts")
dir.create(FORECAST_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# RUNNER FOR ONE COUNTRY
# -----------------------------------------------------------------------------
run_ar_country <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== AR benchmark — %s (n_cores=%d) ===\n", country, n_cores))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  n       <- nrow(daily)
  y       <- daily$y_asinh
  dates   <- daily$date

  ar_mat  <- build_ar_lags(y, AR_LAGS)
  cal_mat <- build_calendar_dummies(dates, country)
  gas_mat <- build_hf_exo_lags(daily$gas_log,   "gas",   1L)
  brt_mat <- build_hf_exo_lags(daily$brent_log, "brent", 1L)

  ar_cols  <- colnames(ar_mat)
  cal_cols <- colnames(cal_mat)
  gas_cols <- colnames(gas_mat)
  brt_cols <- colnames(brt_mat)

  hol_set <- country_holidays(unique(as.integer(format(dates, "%Y"))), country)

  specs <- list(
    ar_dum           = list(extra_cols = character(0)),
    ar_dum_gas       = list(extra_cols = gas_cols),
    ar_dum_brent     = list(extra_cols = brt_cols),
    ar_dum_gas_brent = list(extra_cols = c(gas_cols, brt_cols))
  )

  X_full <- cbind(intercept = 1, cal_mat, ar_mat, gas_mat, brt_mat)
  base_cols <- c("intercept", cal_cols, ar_cols)

  origins <- make_origins(n, WINDOW_DAYS,
                          # Max lag actually used: AR lags {1,2,3,7} and the
                          # fuel lag-1 columns, so 7 pre-sample days suffice
                          # (referee point 11; was 14 via Q_HF).
                          pre_lag_max = max(AR_LAGS),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  # One AR-family specification across all rolling origins.  The design is
  # fully linear (no Almon block), so each window is a single OLS:
  #   asinh(price)_t = intercept + calendar dummies + AR(1,2,3,7) lags
  #                    [+ gas_lag1] [+ brent_lag1] + error
  fit_one_spec <- function(spec_name) {
    cols <- c(base_cols, specs[[spec_name]]$extra_cols)  # columns of THIS spec
    yhat_mat <- matrix(NA_real_, n_orig, H_MAX)
    # Horizon-specific forecast-error variances sigma2_h (one row per origin,
    # one column per step 1..H_MAX): used by the inverse-asinh bias
    # correction, which must grow with the horizon (referee point 3).
    sigma2_mat <- matrix(NA_real_, n_orig, H_MAX)
    for (oi in seq_along(origins)) {
      t0     <- origins[oi]                 # forecast origin (window end)
      win_lo <- t0 - WINDOW_DAYS + 1L      # window start
      Xtr    <- X_full[win_lo:t0, cols, drop = FALSE]  # window design
      ytr    <- y[win_lo:t0]               # window response (asinh scale)
      ok     <- stats::complete.cases(Xtr) & is.finite(ytr)  # drop NA rows
      if (sum(ok) < ncol(Xtr) + 5L) next   # too few observations: skip
      qrfit  <- qr(Xtr[ok, , drop = FALSE])# OLS via QR decomposition
      if (qrfit$rank < ncol(Xtr)) next     # collinear design: skip origin
      beta   <- qr.coef(qrfit, ytr[ok])    # OLS coefficient vector
      names(beta) <- cols
      # In-window one-step residual variance, then propagated to horizon h
      # through the fitted AR dynamics (MA psi-weights): the __bc correction
      # at horizon h uses sigma2_h = sigma2 * sum_{s<h} psi_s^2, NOT the
      # one-step variance at every horizon.
      res <- ytr[ok] - drop(Xtr[ok, , drop = FALSE] %*% beta)
      phi <- as.numeric(beta[ar_cols])       # AR coefficients at lags 1,2,3,7
      sigma2_mat[oi, ] <- iterated_error_variance(AR_LAGS, phi,
                                                  mean(res^2), H_MAX)

      # Iterated multi-step: roll the one-step OLS equation forward h steps,
      # replacing the AR lag cells with realised prices where observed and
      # with previous-step forecasts otherwise; calendar cells refreshed at
      # each true target date; fuel lag-1 cells frozen at the origin value.
      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      frow         <- X_full[t0 + 1L, cols, drop = TRUE]  # frozen base row
      yhat_mat[oi, ] <- iterated_forecast(
        coef          = beta,
        colnames_     = cols,
        history_y_aug = y[1:t0],
        frozen_row    = frow,
        ar_lag_idx    = match(ar_cols, cols),
        ar_lags       = AR_LAGS,
        mon_idx       = NULL,           # refresh_calendar_row uses name patterns
        target_dates  = target_dates,
        h_max         = H_MAX,
        holiday_dates = hol_set
      )
    }
    list(spec = spec_name, yhat_mat = yhat_mat, sigma2 = sigma2_mat)
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(names(specs), fit_one_spec, n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.1fmin)\n", el, el / 60))

  for (h in HORIZONS) {
    out <- data.frame(
      origin_date = dates[origins],
      target_date = dates[origins + h],
      y_actual    = daily$elec_price[origins + h],
      stringsAsFactors = FALSE
    )
    for (r in spec_results) {
      mu_h <- r$yhat_mat[, h]
      out[[r$spec]] <- sinh(mu_h)
      # Bias correction with the h-step error variance (column h).
      out[[paste0(r$spec, "__bc")]] <- sinh(mu_h) * exp(r$sigma2[, h] / 2)
    }
    out <- out[!is.na(out$y_actual), , drop = FALSE]
    fp <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    write.csv(out, fp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", fp, nrow(out)))
  }
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_ar_country(cc_)
cat("\nAR benchmark family done.\n")
