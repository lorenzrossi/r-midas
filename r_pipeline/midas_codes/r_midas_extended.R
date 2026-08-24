#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — R-MIDAS EXTENDED (POOLED, EXP-ALMON NLS) AT h ∈ {1,7,14,21,28}
#                                                              (PARALLEL)
#
# "Extended" R-MIDAS: exp-Almon weighting applies to BOTH the AR block of y
# AND to each included HF exogenous block (gas, brent).  Family in
# summary_table: R_MIDAS_extended.
#
# Design per spec:
#   Linear part X_lin (estimated by OLS, profiled out at each theta evaluation):
#     - intercept                                                            (always)
#     - 6 weekday dummies (Sun = reference)                                   (always)
#     - 11 month-of-year dummies (Jan = reference)                            (always)
#     - 1 holiday dummy (country-specific national holidays)                  (always)
#     - LF columns, one daily column per LF variable                          (conditional)
#
#   Non-linear part Z_list (Almon-weighted, NLS over theta):
#     - y_lags: exp-Almon over the last Q=14 daily lags of y                  (always)
#     - gas:    exp-Almon over the last Q=14 daily lags of gas_log            (if gas in spec)
#     - brent:  exp-Almon over the last Q=14 daily lags of brent_log          (if brent in spec)
#
# NLS protocol (b): n_starts = 5 + 5 * (K - 1), where K = number of Almon
# blocks for the current spec.  Scaled multi-start (Ghysels-Wright 2009 +
# Bates-Granger) so the probe density per dimension is roughly constant as
# the theta vector grows from 2-D to 4-D to 6-D.
#
# Iterated multi-step (option c.1):
#   * y AR-Almon: history_y_aug fed back through fixed θ_y Almon weights
#   * HF Almon (gas, brent): FROZEN at row (origin + 1) z value
#   * LF block: FROZEN at row (origin + 1)
#   * Calendar dummies (weekday, month, holiday): REFRESHED at target row
#
# OUTPUTS:
#   forecasts/r_midas_extended_<cc>_h<h>.csv              forecasts per horizon
#   forecasts/r_midas_extended_spec_dictionary_<cc>.csv   spec → variable membership
#   results/coefficients/r_midas_extended_<cc>.csv        rolling coefficients
#                                                         (origin_date, spec_id,
#                                                         delta_y, delta_gas,
#                                                         delta_brent, <LF>)
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
# FWL-accelerated concentrated NLS: identical estimator, objective evaluated
# through pre-computed cross products (see fit_r_midas_fast.R header).
source(file.path(script_dir, "fit_r_midas_fast.R"))
fit_r_midas <- fit_r_midas_fast

DATA_DIR     <- file.path(script_dir, "data")
FORECAST_DIR <- file.path(script_dir, "forecasts")
COEF_DIR     <- file.path(script_dir, "results", "coefficients")
dir.create(FORECAST_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(COEF_DIR,     showWarnings = FALSE, recursive = TRUE)

Q_AR       <- Q_HF
SPEC_LIMIT <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L

# With deterministic shape restarts (fit_r_midas_fast) the effective number
# of distinct starts is capped at 6 (warm start + 5 canonical shapes); the
# K-scaling formula is kept for the n_starts argument but no longer triggers
# uniform random draws.
n_starts_for_K <- function(K) 5L + 5L * (K - 1L)

# Exponentially discounted (weighted) least squares, applied IDENTICALLY to
# both countries.  RMIDAS_DISCOUNT = 1 (default) disables it.
DISCOUNT_LAMBDA <- suppressWarnings(as.numeric(Sys.getenv("RMIDAS_DISCOUNT", "1")))
if (!is.finite(DISCOUNT_LAMBDA) || DISCOUNT_LAMBDA <= 0 || DISCOUNT_LAMBDA > 1)
  DISCOUNT_LAMBDA <- 1

run_r_midas_extended_h <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== R-MIDAS extended, multi-h — %s (n_cores=%d, n_starts scales with K: 5/10/15) ===\n",
              country, n_cores))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  y       <- daily$y_asinh

  ar_mat    <- build_hf_exo_lags(y,               "y",     Q_AR)
  gas_mat   <- build_hf_exo_lags(daily$gas_log,   "gas",   Q_HF)
  brent_mat <- build_hf_exo_lags(daily$brent_log, "brent", Q_HF)
  Z_full    <- list(y_lags = ar_mat, gas = gas_mat, brent = brent_mat)

  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- c(
    pmi_mfg_de = PUB_LAG_SURVEY, pmi_serv_de = PUB_LAG_SURVEY, ifo_bci_de = PUB_LAG_SURVEY,
    pmi_mfg_it = PUB_LAG_SURVEY, pmi_serv_it = PUB_LAG_SURVEY, istat_bci_it = PUB_LAG_SURVEY,
    Consumer_Goods = PUB_LAG_IPI, Manufacturing = PUB_LAG_IPI, Energy = PUB_LAG_IPI
  )
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  # Surveys: LOCF over missing monthly releases; IPI: no LOCF.
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(
    dates, monthly, v, lf_pub[[v]],
    locf = (lf_pub[[v]] == PUB_LAG_SURVEY))

  cal_mat   <- build_calendar_dummies(dates, country)
  X_base    <- cbind(intercept = 1, cal_mat)
  base_cols <- colnames(X_base)
  hol_set   <- country_holidays(unique(as.integer(format(dates, "%Y"))), country)

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
  dict_path <- file.path(FORECAST_DIR,
                         sprintf("r_midas_extended_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS,
                          pre_lag_max = max(c(AR_LAGS, Q_HF, Q_AR)),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  coef_cols <- c("delta_y", "delta_gas", "delta_brent", lf_names)

  fit_one_spec <- function(si) {
    set.seed(20240501L + si)
    sid <- spec_ids[si]
    mk  <- masks[si, ]
    hf_used <- character(0)
    if (mk["gas"]   == 1L) hf_used <- c(hf_used, "gas")
    if (mk["brent"] == 1L) hf_used <- c(hf_used, "brent")
    lf_used <- lf_names[mk[lf_names] == 1L]

    K_spec  <- 1L + length(hf_used)
    n_st    <- n_starts_for_K(K_spec)

    theta_warm <- NULL
    yhat_mat   <- matrix(NA_real_, n_orig, H_MAX)
    coef_mat   <- matrix(NA_real_, n_orig, length(coef_cols),
                         dimnames = list(NULL, coef_cols))

    for (oi in seq_along(origins)) {
      t0 <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      win_rows <- win_lo:t0
      ywin <- y[win_rows]

      Z_list <- list(y_lags = ar_mat[win_rows, , drop = FALSE])
      for (nm in hf_used) Z_list[[nm]] <- Z_full[[nm]][win_rows, , drop = FALSE]

      X_lin <- if (length(lf_used) > 0L)
        cbind(X_base[win_rows, , drop = FALSE], lf_mat[win_rows, lf_used, drop = FALSE])
      else
        X_base[win_rows, , drop = FALSE]

      keep_train <- stats::complete.cases(X_lin) &
                    rows_all_finite(Z_list$y_lags) &
                    is.finite(ywin)
      for (nm in hf_used)
        keep_train <- keep_train & rows_all_finite(Z_list[[nm]])
      if (sum(keep_train) < ncol(X_lin) + 2L * (1L + length(hf_used)) + 5L) next

      Xlin_fit <- X_lin[keep_train, , drop = FALSE]
      yfit     <- ywin[keep_train]
      Z_fit    <- lapply(Z_list, function(M) M[keep_train, , drop = FALSE])

      # Discounted (weighted) LS via row scaling (see r_midas.R).
      if (DISCOUNT_LAMBDA < 1) {
        sw <- sqrt(DISCOUNT_LAMBDA ^ (t0 - win_rows))[keep_train]
        Xlin_fit <- Xlin_fit * sw
        yfit     <- yfit * sw
        Z_fit    <- lapply(Z_fit, function(M) M * sw)
      }

      if (!is.null(theta_warm) && length(theta_warm) != 2L * length(Z_fit))
        theta_warm <- NULL

      fit <- fit_r_midas(yfit, Z_fit, Xlin_fit,
                         theta_init  = theta_warm,
                         optim_maxit = 2000L,
                         optim_factr = 5e8,
                         n_starts    = n_st)
      if (is.null(fit)) next
      theta_warm <- as.numeric(t(fit$theta))

      coef_mat[oi, "delta_y"] <- fit$delta["y_lags"]
      if ("gas"   %in% hf_used) coef_mat[oi, "delta_gas"]   <- fit$delta["gas"]
      if ("brent" %in% hf_used) coef_mat[oi, "delta_brent"] <- fit$delta["brent"]
      if (length(lf_used) > 0L) {
        for (v in lf_used) {
          b <- fit$beta_lin[v]
          if (!is.na(b)) coef_mat[oi, v] <- as.numeric(b)
        }
      }

      # Low-frequency regressors must belong to the forecast-origin
      # information set. Target-day calendar controls remain taken from
      # t0 + 1 because they are known in advance.
      lf_forecast <- if (length(lf_used) > 0L)
        lf_mat[t0, lf_used, drop = FALSE]
      else
        NULL

      X_lin_row <- if (length(lf_used) > 0L)
        c(X_base[t0 + 1L, , drop = TRUE],
          drop(lf_forecast))
      else
        X_base[t0 + 1L, , drop = TRUE]
      names(X_lin_row) <- c(base_cols, lf_used)
      if (anyNA(X_lin_row)) next

      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      yhat_mat[oi, ] <- iterated_forecast_r_midas(
        fit_obj         = fit,
        Z_list_full     = Z_full,
        X_lin_row       = X_lin_row,
        X_lin_colnames  = names(X_lin_row),
        history_y_aug   = y[1:t0],
        origin_idx      = t0,
        target_dates    = target_dates,
        ar_block_name   = "y_lags",
        hf_block_names  = hf_used,
        h_max           = H_MAX,
        holiday_dates   = hol_set
      )
    }
    list(yhat_mat = yhat_mat, coef_mat = coef_mat)
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec,
                                  n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.2fh) — wall-clock\n", el, el / 3600))

  for (h in HORIZONS) {
    yhat_h_mat <- vapply(seq_along(spec_ids),
                         function(si) sinh(spec_results[[si]]$yhat_mat[, h]),
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
    fp  <- file.path(FORECAST_DIR,
                     sprintf("r_midas_extended_%s_h%02d.csv", cc, h))
    write.csv(out, fp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", fp, nrow(out)))
  }

  coef_blocks <- vector("list", length(spec_ids))
  for (si in seq_along(spec_ids)) {
    cm <- spec_results[[si]]$coef_mat
    blk <- data.frame(
      origin_date = dates[origins],
      spec_id     = spec_ids[si],
      stringsAsFactors = FALSE
    )
    for (col in coef_cols) blk[[col]] <- cm[, col]
    coef_blocks[[si]] <- blk
  }
  coef_df <- do.call(rbind, coef_blocks)
  all_na <- rowSums(is.na(coef_df[, coef_cols, drop = FALSE])) == length(coef_cols)
  coef_df <- coef_df[!all_na, , drop = FALSE]
  coef_path <- file.path(COEF_DIR, sprintf("r_midas_extended_%s.csv", cc))
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d, cols = %d)\n",
              coef_path, nrow(coef_df), ncol(coef_df)))
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_r_midas_extended_h(cc_)
cat("\nR-MIDAS extended multi-h done.\n")
