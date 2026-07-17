#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — R-MIDAS (POOLED, EXPONENTIAL-ALMON NLS), h = 1 ONLY (PARALLEL)
#
# Same exp-Almon AR + optional HF Almon (gas, brent) + optional LF block
# structure as r_midas_surveys_h.R, but ONLY the one-step-ahead direct
# forecast at h = 1.  No iterated multi-step is needed.
#
# Spec grid: 2^8 = 256 subsets.  Warm-start theta across rolling origins.
# Parallelism: per-spec workers via parallel::mclapply (N_CORES).
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

Q_AR       <- Q_HF
SPEC_LIMIT <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_r_midas_h1 <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== R-MIDAS pooled, h=1 — %s (n_cores=%d) ===\n",
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

  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- c(
    pmi_mfg_de = PUB_LAG_SURVEY, pmi_serv_de = PUB_LAG_SURVEY, ifo_bci_de = PUB_LAG_SURVEY,
    pmi_mfg_it = PUB_LAG_SURVEY, pmi_serv_it = PUB_LAG_SURVEY, istat_bci_it = PUB_LAG_SURVEY,
    Consumer_Goods = PUB_LAG_IPI, Manufacturing = PUB_LAG_IPI, Energy = PUB_LAG_IPI
  )
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(dates, monthly, v, lf_pub[[v]])

  mon_mat <- build_month_dummies(dates)
  X_base  <- cbind(intercept = 1, mon_mat)
  base_cols <- colnames(X_base)

  hf_names   <- c("gas", "brent")
  flag_names <- c(hf_names, lf_names)
  hf_lookup  <- list(gas = gas_mat, brent = brent_mat)

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
    cat(sprintf("  SPEC_LIMIT=%d\n", SPEC_LIMIT))
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
                         sprintf("r_midas_h1_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS,
                          pre_lag_max = max(c(AR_LAGS, Q_HF, Q_AR)),
                          h_max = 1L)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  # ---- Worker: one spec across all origins -------------------------------
  fit_one_spec <- function(si) {
    sid <- spec_ids[si]
    mk  <- masks[si, ]
    hf_used <- character(0)
    if (mk["gas"]   == 1L) hf_used <- c(hf_used, "gas")
    if (mk["brent"] == 1L) hf_used <- c(hf_used, "brent")
    lf_used <- lf_names[mk[lf_names] == 1L]

    theta_warm <- NULL
    yhat_vec   <- rep(NA_real_, n_orig)

    for (oi in seq_along(origins)) {
      t0 <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      win_rows <- win_lo:t0
      ywin <- y[win_rows]

      Z_list <- list(y_lags = ar_mat[win_rows, , drop = FALSE])
      for (nm in hf_used) Z_list[[nm]] <- hf_lookup[[nm]][win_rows, , drop = FALSE]

      X_lin <- if (length(lf_used) > 0L)
        cbind(X_base[win_rows, , drop = FALSE], lf_mat[win_rows, lf_used, drop = FALSE])
      else
        X_base[win_rows, , drop = FALSE]

      keep_train <- stats::complete.cases(X_lin) &
                    apply(Z_list$y_lags, 1, function(r) all(is.finite(r))) &
                    is.finite(ywin)
      for (nm in hf_used)
        keep_train <- keep_train & apply(Z_list[[nm]], 1, function(r) all(is.finite(r)))
      if (sum(keep_train) < ncol(X_lin) + 2L * (1L + length(hf_used)) + 5L) next

      Xlin_fit <- X_lin[keep_train, , drop = FALSE]
      yfit     <- ywin[keep_train]
      Z_fit    <- lapply(Z_list, function(M) M[keep_train, , drop = FALSE])

      if (!is.null(theta_warm) && length(theta_warm) != 2L * length(Z_fit))
        theta_warm <- NULL

      fit <- fit_r_midas(yfit, Z_fit, Xlin_fit,
                         theta_init  = theta_warm,
                         optim_maxit = if (is.null(theta_warm)) 800L else 200L)
      if (is.null(fit)) next
      theta_warm <- as.numeric(t(fit$theta))

      X_lin_row <- if (length(lf_used) > 0L)
        c(X_base[t0 + 1L, , drop = TRUE],
          lf_mat[t0 + 1L, lf_used, drop = TRUE])
      else
        X_base[t0 + 1L, , drop = TRUE]
      names(X_lin_row) <- c(base_cols, lf_used)
      if (anyNA(X_lin_row)) next

      Q <- ncol(Z_list$y_lags)
      W_y <- exp_almon_weights(Q, fit$theta["y_lags", ])
      z_y <- sum(W_y * y[(t0 - Q + 1L):t0])

      hf_part <- 0
      for (nm in hf_used) {
        W_nm <- exp_almon_weights(Q, fit$theta[nm, ])
        z_nm <- sum(W_nm * hf_lookup[[nm]][t0 + 1L, ])
        hf_part <- hf_part + fit$delta[nm] * z_nm
      }

      lin_part <- sum(fit$beta_lin * X_lin_row)
      yhat_vec[oi] <- lin_part + fit$delta["y_lags"] * z_y + hf_part
    }
    yhat_vec
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec,
                                  n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.2fh) — wall-clock\n", el, el / 3600))

  out <- data.frame(
    origin_date = dates[origins],
    target_date = dates[origins + 1L],
    y_actual    = daily$elec_price[origins + 1L],
    stringsAsFactors = FALSE
  )
  for (si in seq_along(spec_ids)) {
    out[[spec_ids[si]]] <- sinh(spec_results[[si]])
  }
  out <- out[!is.na(out$y_actual), , drop = FALSE]
  fp  <- file.path(FORECAST_DIR, sprintf("r_midas_h1_%s.csv", cc))
  write.csv(out, fp, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d)\n", fp, nrow(out)))
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# DRIVER
# -----------------------------------------------------------------------------
country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_r_midas_h1(cc_)
cat("\nR-MIDAS pooled h=1 done.\n")
