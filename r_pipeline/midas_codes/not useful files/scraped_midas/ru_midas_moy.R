#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — RU-MIDAS (OLS) — SECONDARY VERSION WITH MONTH-OF-YEAR DUMMIES
#
# This is the SECONDARY RU-MIDAS file in the v2 pipeline.  The PRINCIPAL
# RU-MIDAS is `ru_midas.R`, which follows Foroni-Rossini-et-al. (ECB WP 2250)
# and Cimadomo et al. (2020) by interacting every lag with 28 within-month
# position dummies.  This file is kept as a robustness comparison: it uses
# the much simpler month-of-year (Jan-Dec) seasonal control with FREE Weron
# AR(1,2,3,7) and 14 unrestricted lags per HF exogenous variable, in the
# spirit of Foroni-Marcellino-Schumacher (2015, JoE) unrestricted MIDAS.
# Family name in summary_table.csv: RU_MIDAS_moy.
#
# Per country and per horizon h ∈ {1, 7, 14, 21, 28}:
#   Base design (always included):
#     - intercept
#     - 11 month-of-year dummies (Jan = reference)
#     - Weron AR lags: y_{t-1}, y_{t-2}, y_{t-3}, y_{t-7}     (free coefficients)
#   Candidate HF exogenous blocks (each = 14 daily unrestricted lags):
#     - gas_log_lag1..14
#     - brent_log_lag1..14
#   Candidate LF blocks (each = ONE daily column = LOCF of pub-aware monthly):
#     DE: pmi_mfg_de, pmi_serv_de, ifo_bci_de, Consumer_Goods, Manufacturing, Energy
#     IT: pmi_mfg_it, pmi_serv_it, istat_bci_it, Consumer_Goods, Manufacturing, Energy
#     Pub lags: IPI = 2 months, surveys (PMI / IFO / ISTAT) = 1 month.
#
# Sweep ALL 2^(2 HF + 6 LF) = 256 subsets.  Spec 000 is "MIDAS with just
# dummies" (AR + month dummies only).
#
# Parallelism: per-spec workers via parallel::mclapply.  Each worker iterates
# all rolling origins for one spec and returns its yhat matrix.  N_CORES env
# var (default = detectCores()-1) controls concurrency.  At 256 specs, the
# script saturates as many cores as you give it; expected speedup on M4 is
# 6-10x vs serial.
# =============================================================================

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grep("^--file=", args)]
  if (length(farg) > 0) dirname(normalizePath(sub("^--file=", "", farg))) else getwd()
}, error = function(e) getwd())
if (dir.exists(script_dir)) setwd(script_dir)

source(file.path(script_dir, "common_utils.R"))

DATA_DIR     <- file.path(script_dir, "data")
FORECAST_DIR <- file.path(script_dir, "forecasts")
dir.create(FORECAST_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# RUNNER FOR ONE COUNTRY
# -----------------------------------------------------------------------------
run_ru_midas_country <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== RU-MIDAS combos — %s (n_cores=%d) ===\n", country, n_cores))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  y       <- daily$y_asinh

  # ---- HF exogenous lag blocks --------------------------------------------
  gas_mat   <- build_hf_exo_lags(daily$gas_log,   "gas",   Q_HF)
  brent_mat <- build_hf_exo_lags(daily$brent_log, "brent", Q_HF)

  # ---- LF (daily, publication-aware LOCF) ---------------------------------
  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- c(
    pmi_mfg_de  = PUB_LAG_SURVEY, pmi_serv_de  = PUB_LAG_SURVEY, ifo_bci_de   = PUB_LAG_SURVEY,
    pmi_mfg_it  = PUB_LAG_SURVEY, pmi_serv_it  = PUB_LAG_SURVEY, istat_bci_it = PUB_LAG_SURVEY,
    Consumer_Goods = PUB_LAG_IPI, Manufacturing = PUB_LAG_IPI, Energy        = PUB_LAG_IPI
  )
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  for (v in lf_names)
    lf_mat[, v] <- make_daily_lf(dates, monthly, v, lf_pub[[v]])

  # ---- Base design --------------------------------------------------------
  ar_mat   <- build_ar_lags(y, AR_LAGS)
  mon_mat  <- build_month_dummies(dates)
  ar_cols  <- colnames(ar_mat)
  mon_cols <- colnames(mon_mat)

  X_full <- cbind(
    intercept = 1,
    mon_mat,
    ar_mat,
    gas_mat,
    brent_mat,
    lf_mat
  )
  base_cols     <- c("intercept", mon_cols, ar_cols)
  hf_col_groups <- list(gas = colnames(gas_mat), brent = colnames(brent_mat))

  # ---- Combo enumeration --------------------------------------------------
  hf_names   <- c("gas", "brent")
  flag_names <- c(hf_names, lf_names)
  masks <- expand.grid(replicate(length(flag_names), 0:1, simplify = FALSE),
                       KEEP.OUT.ATTRS = FALSE)
  masks <- as.matrix(masks); colnames(masks) <- flag_names
  masks_int <- as.integer(masks %*% (2^(rev(seq_len(ncol(masks))) - 1L)))
  ord <- order(masks_int)
  masks <- masks[ord, , drop = FALSE]
  spec_ids <- sprintf("spec_%03d", 0:(nrow(masks) - 1L))

  spec_cols <- lapply(seq_len(nrow(masks)), function(i) {
    mk <- masks[i, ]
    extra <- character(0)
    for (h_v in hf_names)
      if (mk[h_v] == 1L) extra <- c(extra, hf_col_groups[[h_v]])
    for (l_v in lf_names)
      if (mk[l_v] == 1L) extra <- c(extra, l_v)
    c(base_cols, extra)
  })
  names(spec_cols) <- spec_ids

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
  dict_path <- file.path(FORECAST_DIR, sprintf("ru_midas_moy_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS,
                          pre_lag_max = max(c(AR_LAGS, Q_HF)),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  # ---- Worker: one spec across all origins -------------------------------
  fit_one_spec <- function(si) {
    cols <- spec_cols[[si]]
    yhat_mat <- matrix(NA_real_, n_orig, H_MAX)
    for (oi in seq_along(origins)) {
      t0     <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      Xtr    <- X_full[win_lo:t0, cols, drop = FALSE]
      ytr    <- y[win_lo:t0]
      ok     <- stats::complete.cases(Xtr) & is.finite(ytr)
      if (sum(ok) < ncol(Xtr) + 5L) next
      qrfit  <- qr(Xtr[ok, , drop = FALSE])
      if (qrfit$rank < ncol(Xtr)) next
      beta   <- qr.coef(qrfit, ytr[ok])
      names(beta) <- cols

      frow <- X_full[t0 + 1L, cols, drop = TRUE]
      if (anyNA(frow)) next

      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      yhat_mat[oi, ] <- iterated_forecast(
        coef          = beta,
        colnames_     = cols,
        history_y_aug = y[1:t0],
        frozen_row    = frow,
        ar_lag_idx    = match(ar_cols, cols),
        ar_lags       = AR_LAGS,
        mon_idx       = match(mon_cols, cols),
        target_dates  = target_dates,
        h_max         = H_MAX
      )
    }
    yhat_mat
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec,
                                  n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.1fmin) — wall-clock\n", el, el / 60))

  # ---- Write per-horizon CSVs --------------------------------------------
  for (h in HORIZONS) {
    yhat_h_mat <- vapply(seq_along(spec_ids),
                         function(si) sinh(spec_results[[si]][, h]),
                         numeric(n_orig))
    colnames(yhat_h_mat) <- spec_ids
    out <- data.frame(
      origin_date = dates[origins],
      target_date = dates[origins + h],
      y_actual    = daily$elec_price[origins + h],
      stringsAsFactors = FALSE
    )
    out <- cbind(out, as.data.frame(yhat_h_mat, stringsAsFactors = FALSE))
    out <- out[!is.na(out$y_actual), , drop = FALSE]
    fp  <- file.path(FORECAST_DIR, sprintf("ru_midas_moy_%s_h%02d.csv", cc, h))
    write.csv(out, fp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d, specs = %d)\n",
                fp, nrow(out), length(spec_ids)))
  }
}

# -----------------------------------------------------------------------------
# DRIVER
# -----------------------------------------------------------------------------
country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_ru_midas_country(cc_)
cat("\nRU-MIDAS combos done.\n")
