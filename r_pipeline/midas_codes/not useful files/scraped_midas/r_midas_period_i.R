#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — PERIOD-SPECIFIC R-MIDAS (FMS 2015 Eq. 8)  (PARALLEL)
#
# For k = 28 within-month positions, a separate R-MIDAS equation per i:
#
#   x_{tau,i} = alpha_i + lambda_i' z^{LF}_{tau - pub_lag}
#                       + delta_{y,i} * sum_j w(j; theta_{y,i}) y_{t-j}
#                       + delta_{g,i} * sum_j w(j; theta_{g,i}) gas_{t-j}    [if gas]
#                       + delta_{b,i} * sum_j w(j; theta_{b,i}) brent_{t-j}  [if brent]
#                       + eps_{tau,i}
#
# Implementation 1, "textbook per-i, no pooled controls".  Days whose strict
# position > k are dropped from estimation.
#
# Iterated multi-step at h ∈ {1, 7, 14, 21, 28}: at step s the equation for
# position i_s := position(origin + s) is invoked.  AR feedback via y_aug;
# HF Almon and LF FROZEN at row (origin+1).
#
# Spec grid: 256 subsets, but SPEC_LIMIT defaults to 16 because 256 * k * n_orig
# NLS fits is genuinely heavy.
#
# Parallelism: per-spec workers via parallel::mclapply.  Each worker iterates
# all rolling origins for ONE spec, maintaining a local per-position warm-start
# theta cache (k * 1 vector per spec).  N_CORES env var controls concurrency.
# Expected speedup on M4 with SPEC_LIMIT=16: full 10x (one core per spec).
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

K_POS      <- 28L
Q_AR       <- Q_HF
SPEC_LIMIT <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "16")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 16L

within_month_position <- function(dates, k = K_POS) {
  pos <- as.integer(format(dates, "%d"))
  ifelse(pos > k, NA_integer_, pos)
}

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_period_i <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== R-MIDAS period-i (FMS Eq.8) — %s (n_cores=%d) ===\n",
              country, n_cores))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  y       <- daily$y_asinh
  pos     <- within_month_position(dates, K_POS)

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
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(dates, monthly, v, lf_pub[[v]])

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
  cat(sprintf("  specs = %d  (SPEC_LIMIT=%d of 256)\n",
              length(spec_ids), SPEC_LIMIT))

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
                         sprintf("r_midas_period_i_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s\n", dict_path))

  origins <- make_origins(n, WINDOW_DAYS,
                          pre_lag_max = max(c(AR_LAGS, Q_HF, Q_AR)),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  # ---- Worker: one spec across all origins and all 28 positions ----------
  fit_one_spec <- function(si) {
    sid <- spec_ids[si]
    mk  <- masks[si, ]
    hf_used <- character(0)
    if (mk["gas"]   == 1L) hf_used <- c(hf_used, "gas")
    if (mk["brent"] == 1L) hf_used <- c(hf_used, "brent")
    lf_used <- lf_names[mk[lf_names] == 1L]

    theta_warm_by_i <- vector("list", K_POS)   # worker-local warm-start cache
    yhat_mat <- matrix(NA_real_, n_orig, H_MAX)

    for (oi in seq_along(origins)) {
      t0 <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      win_rows <- win_lo:t0

      per_i_fits <- vector("list", K_POS)
      for (i in 1:K_POS) {
        rows_i <- win_rows[pos[win_rows] == i & !is.na(pos[win_rows])]
        if (length(rows_i) < 5L * max(length(lf_used) + length(hf_used) + 1L, 4L)) next

        Z_list_i <- list(y_lags = ar_mat[rows_i, , drop = FALSE])
        for (nm in hf_used) Z_list_i[[nm]] <- Z_full[[nm]][rows_i, , drop = FALSE]

        X_lin_i <- if (length(lf_used) > 0L)
          cbind(intercept = 1, lf_mat[rows_i, lf_used, drop = FALSE])
        else
          matrix(1, length(rows_i), 1L, dimnames = list(NULL, "intercept"))

        keep <- stats::complete.cases(X_lin_i) & is.finite(y[rows_i])
        for (nm in names(Z_list_i))
          keep <- keep & apply(Z_list_i[[nm]], 1, function(r) all(is.finite(r)))
        if (sum(keep) < ncol(X_lin_i) + 2L * (1L + length(hf_used)) + 3L) next

        Z_fit  <- lapply(Z_list_i, function(M) M[keep, , drop = FALSE])
        Xfit   <- X_lin_i[keep, , drop = FALSE]
        yfit_i <- y[rows_i][keep]

        theta_init <- theta_warm_by_i[[i]]
        if (!is.null(theta_init) &&
            length(theta_init) != 2L * length(Z_fit)) theta_init <- NULL

        fit_i <- fit_r_midas(yfit_i, Z_fit, Xfit,
                             theta_init  = theta_init,
                             optim_maxit = if (is.null(theta_init)) 600L else 150L)
        if (!is.null(fit_i)) {
          per_i_fits[[i]] <- fit_i
          theta_warm_by_i[[i]] <- as.numeric(t(fit_i$theta))
        }
      }

      frozen_lf <- if (length(lf_used) > 0L)
        lf_mat[t0 + 1L, lf_used, drop = TRUE] else NULL
      if (!is.null(frozen_lf) && anyNA(frozen_lf)) next

      history_y_aug <- y[1:t0]
      h_path <- numeric(H_MAX); ok_path <- logical(H_MAX)
      for (s in 1:H_MAX) {
        target_date <- dates[t0 + s]
        i_s <- within_month_position(target_date, K_POS)
        if (is.na(i_s)) i_s <- K_POS
        fit_s <- per_i_fits[[i_s]]
        if (is.null(fit_s)) { ok_path[s] <- FALSE; break }

        X_lin_row <- c(intercept = 1, if (length(lf_used) > 0L) frozen_lf else numeric(0))
        names(X_lin_row) <- c("intercept", lf_used)
        lin_part <- sum(fit_s$beta_lin[names(X_lin_row)] * X_lin_row)

        Q <- ncol(ar_mat)
        W_y <- exp_almon_weights(Q, fit_s$theta["y_lags", ])
        n_hist <- length(history_y_aug)
        ar_idx <- n_hist - (1:Q) + 1L
        if (any(ar_idx < 1L)) { ok_path[s] <- FALSE; break }
        z_y <- sum(W_y * history_y_aug[ar_idx])

        hf_part <- 0
        for (nm in hf_used) {
          W_nm <- exp_almon_weights(Q, fit_s$theta[nm, ])
          z_nm <- sum(W_nm * Z_full[[nm]][t0 + 1L, ])
          hf_part <- hf_part + fit_s$delta[nm] * z_nm
        }
        yh <- lin_part + fit_s$delta["y_lags"] * z_y + hf_part
        h_path[s] <- yh; ok_path[s] <- TRUE
        history_y_aug <- c(history_y_aug, yh)
      }

      for (h in HORIZONS) {
        if (h <= length(ok_path) && ok_path[h]) {
          yhat_mat[oi, h] <- sinh(h_path[h])
        }
      }
    }
    yhat_mat
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec,
                                  n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.2fh) — wall-clock\n", el, el / 3600))

  for (h in HORIZONS) {
    target_pos_vec <- within_month_position(dates[origins + h], K_POS)
    out <- data.frame(
      origin_date = dates[origins],
      target_date = dates[origins + h],
      target_pos  = target_pos_vec,
      y_actual    = daily$elec_price[origins + h],
      stringsAsFactors = FALSE
    )
    for (si in seq_along(spec_ids)) {
      out[[spec_ids[si]]] <- spec_results[[si]][, h]
    }
    out <- out[!is.na(out$y_actual), , drop = FALSE]
    fp  <- file.path(FORECAST_DIR,
                     sprintf("r_midas_period_i_%s_h%02d.csv", cc, h))
    write.csv(out, fp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", fp, nrow(out)))
  }
}

# -----------------------------------------------------------------------------
# DRIVER
# -----------------------------------------------------------------------------
country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_period_i(cc_)
cat("\nR-MIDAS period-i (FMS Eq.8) done.\n")
