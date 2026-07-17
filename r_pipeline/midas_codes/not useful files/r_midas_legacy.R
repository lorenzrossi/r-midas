#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — R-MIDAS LEGACY (POOLED, EXP-ALMON ON y ONLY)
#                                                              (PARALLEL)
#
# Methodological twin of `r_midas.R` (exp-Almon polynomial on the Q = 14 daily
# AR lags of y; gas_lag1 and brent_lag1 as plain linear single-lag controls;
# LF block enters linearly), implemented in the style of the legacy
# standard_midas/r_midas/r_midas_surveys_h.R script.  Family in summary_table:
# R_MIDAS_legacy.
#
# Differences from r_midas.R that are deliberate (legacy-style):
#   * NLS protocol: 5-start L-BFGS-B per rolling window per spec, optim_maxit
#     = 3000L, optim_factr = 1e7 (default optim tolerance).  The v2 r_midas.R
#     uses the same 5-start protocol but with optim_maxit = 2000L and a
#     looser factr = 5e8.  Both warm-start theta from the previous origin.
#   * Inline concentrated_sse and fit_window_model functions, following the
#     legacy code layout, instead of calling fit_r_midas from common_utils.R.
#   * lm.fit-based OLS at theta_hat (rather than qr.coef) — same numerical
#     answer, slightly different code path.
#
# Identical to r_midas.R in every other respect:
#   * Same data via common_utils::load_country_data (asinh y, log gas, log
#     brent, Δlog IPI, scale PMI/IFO/ISTAT).
#   * Same calendar control set (6 weekday + 11 month + 1 holiday).
#   * Same 256-spec combo grid over {gas, brent} × {pmi_mfg, pmi_serv, bci,
#     Consumer_Goods, Manufacturing, Energy}.
#   * Same publication lags (IPI = 2, surveys = 1) and pre_lag_max.
#   * Same iterated multi-step (option c.1) with calendar refresh at target.
#   * Same per-spec parallel via parallel_lapply.
#   * Same rolling coefficient extraction: results/coefficients/r_midas_legacy_<cc>.csv
#     with columns (origin_date, spec_id, delta_y, <each LF name>).
#
# The intent is cross-implementation validation: if r_midas.R and
# r_midas_legacy.R produce numerically very similar forecasts and similar
# rolling-coefficient time series, the methodology is solid; meaningful
# divergence between the two would indicate an implementation bug worth
# investigating.
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

# Legacy-style NLS knobs.
Q_ar         <- Q_HF              # 14 daily AR lags (legacy lowercase Q_ar)
theta_init   <- c(0, 0)
theta_lower  <- c(-8, -8)
theta_upper  <- c( 8,  0)
optim_method <- "L-BFGS-B"
optim_maxit  <- 3000L              # legacy default
optim_factr  <- 1e7                # optim default tolerance
optim_starts <- 5L

SPEC_LIMIT <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L

# -----------------------------------------------------------------------------
# Legacy-style inline NLS helpers
# -----------------------------------------------------------------------------
# exp_almon_weights from common_utils.R already uses log-sum-exp; we reuse it.
# (The legacy script defined the same function inline; we keep the single
# implementation in common_utils.R to avoid drift.)

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

# Per-window concentrated NLS with multi-start.  Returns NULL on failure.
fit_window_model <- function(y_vec, x_lin, ar_lags, theta_start) {
  starts <- matrix(NA_real_, nrow = optim_starts, ncol = 2L)
  starts[1, ] <- if (length(theta_start) == 2L && all(is.finite(theta_start)))
                   theta_start else theta_init
  if (optim_starts > 1L) {
    for (s in 2:optim_starts) {
      starts[s, ] <- runif(2L, min = theta_lower, max = theta_upper)
    }
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

  # Refit OLS at theta_hat with lm.fit (legacy style)
  theta_hat <- best$par
  z_ar <- build_weighted_ar(theta_hat, ar_lags)
  ok   <- stats::complete.cases(cbind(x_lin, z_ar, y_vec))
  if (sum(ok) < ncol(x_lin) + 2L) return(NULL)
  X_full <- cbind(x_lin[ok, , drop = FALSE], z_ar = z_ar[ok])
  fit <- stats::lm.fit(X_full, y_vec[ok])
  beta <- fit$coefficients
  if (anyNA(beta)) return(NULL)

  list(theta       = theta_hat,
       beta        = beta,           # named: linear cols then z_ar
       beta_names  = colnames(X_full),
       weights     = exp_almon_weights(ncol(ar_lags), theta_hat),
       rss         = sum(fit$residuals^2),
       n_obs       = sum(ok),
       k_params    = ncol(X_full))
}

# -----------------------------------------------------------------------------
# Legacy-style iterated multi-step forecast
# -----------------------------------------------------------------------------
# Mirrors r_midas.R's iterated_forecast_r_midas behaviour but with the
# beta / weights returned by fit_window_model.
#
# Inputs:
#   fit            output of fit_window_model
#   x_lin_row      frozen linear row at (origin + 1), names matching the
#                  beta_names that index into beta
#   history_y_aug  vector of y_asinh through origin (length = t0)
#   target_dates   length-H_MAX Date vector
#   holiday_dates  Date set for the country's national holidays
legacy_iterated_predict <- function(fit, x_lin_row, history_y_aug,
                                    target_dates, holiday_dates,
                                    h_max = H_MAX) {
  beta <- fit$beta
  bnms <- fit$beta_names
  w    <- fit$weights
  Q    <- length(w)
  # Index of the z_ar slope in beta
  z_ar_idx  <- match("z_ar", bnms)
  if (is.na(z_ar_idx))
    stop("fit$beta is missing the z_ar slope coefficient.")
  delta_y   <- beta[z_ar_idx]
  beta_lin  <- beta[-z_ar_idx]
  lin_names <- bnms[-z_ar_idx]
  names(beta_lin) <- lin_names

  y_hat_h <- numeric(h_max)
  for (s in 1:h_max) {
    # Refresh calendar columns on x_lin_row (mon, weekday, holiday)
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
# RUNNER
# -----------------------------------------------------------------------------
run_r_midas_legacy <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== R-MIDAS LEGACY (Almon on y, linear HF) — %s (n_cores=%d, n_starts=%d) ===\n",
              country, n_cores, optim_starts))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  y       <- daily$y_asinh

  # AR lag matrix (legacy: build with loop, Q_ar columns)
  ar_lag_mat <- matrix(NA_real_, nrow = n, ncol = Q_ar)
  colnames(ar_lag_mat) <- paste0("y_lag", 1:Q_ar)
  for (j in 1:Q_ar) ar_lag_mat[(j + 1L):n, j] <- y[1:(n - j)]

  # Single-lag HF exo (gas_log_lag1, brent_log_lag1)
  gas_lag1 <- c(NA_real_, daily$gas_log[1:(n - 1L)])
  brt_lag1 <- c(NA_real_, daily$brent_log[1:(n - 1L)])

  # LF (publication-aware LOCF)
  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- c(
    pmi_mfg_de = PUB_LAG_SURVEY, pmi_serv_de = PUB_LAG_SURVEY, ifo_bci_de = PUB_LAG_SURVEY,
    pmi_mfg_it = PUB_LAG_SURVEY, pmi_serv_it = PUB_LAG_SURVEY, istat_bci_it = PUB_LAG_SURVEY,
    Consumer_Goods = PUB_LAG_IPI, Manufacturing = PUB_LAG_IPI, Energy = PUB_LAG_IPI
  )
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(dates, monthly, v, lf_pub[[v]])

  # Calendar controls: 6 weekday + 11 month + 1 holiday
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
                         sprintf("r_midas_legacy_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS,
                          pre_lag_max = max(c(AR_LAGS, Q_HF, Q_ar)),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  coef_cols <- c("delta_y", lf_names)

  fit_one_spec <- function(si) {
    set.seed(20240501L + si)
    sid <- spec_ids[si]
    mk  <- masks[si, ]
    hf_used <- character(0)
    if (mk["gas"]   == 1L) hf_used <- c(hf_used, "gas")
    if (mk["brent"] == 1L) hf_used <- c(hf_used, "brent")
    lf_used <- lf_names[mk[lf_names] == 1L]

    theta_warm <- NULL
    yhat_mat   <- matrix(NA_real_, n_orig, H_MAX)
    coef_mat   <- matrix(NA_real_, n_orig, length(coef_cols),
                         dimnames = list(NULL, coef_cols))

    for (oi in seq_along(origins)) {
      t0 <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      win_rows <- win_lo:t0
      ywin <- y[win_rows]
      ar_lags_win <- ar_lag_mat[win_rows, , drop = FALSE]

      # Build x_lin: intercept + calendar + LF + HF lag1 cols (legacy style:
      # all linear regressors stacked into one matrix)
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
        lin_cols[[length(lin_cols) + 1L]] <- matrix(brt_lag1[win_rows], ncol = 1,
                                                    dimnames = list(NULL, "brent_lag1"))
        lin_names <- c(lin_names, "brent_lag1")
      }
      x_lin_win <- do.call(cbind, lin_cols)
      colnames(x_lin_win) <- lin_names

      keep_train <- stats::complete.cases(x_lin_win) &
                    apply(ar_lags_win, 1, function(r) all(is.finite(r))) &
                    is.finite(ywin)
      if (sum(keep_train) < ncol(x_lin_win) + 3L) next

      x_lin_fit <- x_lin_win[keep_train, , drop = FALSE]
      ar_lags_fit <- ar_lags_win[keep_train, , drop = FALSE]
      yfit <- ywin[keep_train]

      if (!is.null(theta_warm) && length(theta_warm) != 2L) theta_warm <- NULL

      fit <- fit_window_model(yfit, x_lin_fit, ar_lags_fit, theta_warm)
      if (is.null(fit)) next
      theta_warm <- fit$theta

      # ---- RECORD ROLLING COEFFICIENTS ----
      coef_mat[oi, "delta_y"] <- fit$beta["z_ar"]
      if (length(lf_used) > 0L) {
        for (v in lf_used) {
          b <- fit$beta[v]
          if (!is.na(b)) coef_mat[oi, v] <- as.numeric(b)
        }
      }

      # Frozen linear row at (origin + 1).  Must include all column names of
      # x_lin_win so refresh_calendar_row + legacy_iterated_predict can index
      # by name.
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
        lin_row_parts[[length(lin_row_parts) + 1L]] <- brt_lag1[t0 + 1L]
        lin_row_names <- c(lin_row_names, "brent_lag1")
      }
      x_lin_row <- unlist(lin_row_parts)
      names(x_lin_row) <- lin_row_names
      if (anyNA(x_lin_row)) next

      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      yhat_mat[oi, ] <- legacy_iterated_predict(
        fit            = fit,
        x_lin_row      = x_lin_row,
        history_y_aug  = y[1:t0],
        target_dates   = target_dates,
        holiday_dates  = hol_set,
        h_max          = H_MAX
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
                     sprintf("r_midas_legacy_%s_h%02d.csv", cc, h))
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
  coef_path <- file.path(COEF_DIR, sprintf("r_midas_legacy_%s.csv", cc))
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d, cols = %d)\n",
              coef_path, nrow(coef_df), ncol(coef_df)))
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_r_midas_legacy(cc_)
cat("\nR-MIDAS legacy (Almon on y, linear HF) done.\n")
