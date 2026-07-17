#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — R-MIDAS (POOLED, EXP-ALMON ON y ONLY) AT h ∈ {1,7,14,21,28}
#                                                              (PARALLEL)
#
# Almon polynomial on the AR lag set J = {1, 2, 3, 7} of y (same lags as the
# AR benchmark; the actual j values enter the Almon exponent), with all other
# daily regressors entering linearly at lag 1.  Family name in summary_table:
# R_MIDAS.
#
# Design per spec:
#   Linear part X_lin (estimated by OLS, profiled out at each theta evaluation):
#     - intercept                                                            (always)
#     - 6 weekday dummies (Sun = reference)                                   (always)
#     - 11 month-of-year dummies (Jan = reference)                            (always)
#     - 1 holiday dummy (country-specific national holidays)                  (always)
#     - LF columns, one daily column per LF variable                          (conditional)
#     - gas_log_lag1                                                          (if gas in spec)
#     - brent_log_lag1                                                        (if brent in spec)
#
#   Non-linear part Z_list (Almon-weighted, NLS over theta):
#     - y_lags: exp-Almon over the AR lag set J = {1, 2, 3, 7}                (always)
#       (yesterday, two/three days ago, one week ago; the j values of J enter
#       the exponential-Almon exponent — see model.tex eq. (almon).  Same lag
#       set as the AR benchmark, so R-MIDAS restricts exactly the block the
#       benchmark leaves free.)
#
# NLS protocol: 5-start L-BFGS-B per rolling window per spec. Start 1 is the
# previous origin's warm-start theta (or (0,0) initially); the remaining
# starts are deterministic canonical Almon shapes supplied by
# fit_r_midas_fast.R. The lowest-SSE admissible optimum is retained.
#
# Iterated multi-step (option c.1):
#   * y AR-Almon: history_y_aug fed back through fixed θ_y Almon weights
#   * LF and HF lag1 columns: FROZEN at row (origin + 1)
#   * Calendar dummies (weekday, month, holiday): REFRESHED at target row
#
# OUTPUTS:
#   forecasts/r_midas_<cc>_h<h>.csv              forecasts per horizon
#   forecasts/r_midas_spec_dictionary_<cc>.csv   spec → variable membership
#   results/coefficients/r_midas_<cc>.csv        rolling-window coefficients
#                                                (wide format: origin_date,
#                                                spec_id, delta_y, <each LF>)
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

Q_AR             <- Q_HF
SPEC_LIMIT       <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L
N_STARTS_R_MIDAS <- 5L

# Exponentially discounted (weighted) least squares, applied IDENTICALLY to
# both countries so the methodology stays symmetric.  Rows of the estimation
# window are weighted lambda^(t0 - t) (weight 1 on the most recent day).
# RMIDAS_DISCOUNT = 1 (default) disables it; e.g. 0.995 halves the weight
# roughly every 139 days.
DISCOUNT_LAMBDA <- suppressWarnings(as.numeric(Sys.getenv("RMIDAS_DISCOUNT", "1")))
if (!is.finite(DISCOUNT_LAMBDA) || DISCOUNT_LAMBDA <= 0 || DISCOUNT_LAMBDA > 1)
  DISCOUNT_LAMBDA <- 1

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_r_midas <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== R-MIDAS standard (Almon on y, linear HF) — %s (n_cores=%d, n_starts=%d) ===\n",
              country, n_cores, N_STARTS_R_MIDAS))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  y       <- daily$y_asinh

  # Almon AR block of y over the lag set J = {1, 2, 3, 7} (model.tex)
  ar_mat <- build_ar_lags(y, AR_LAGS)
  Z_full <- list(y_lags = ar_mat)
  AR_LAG_SETS <- list(y_lags = AR_LAGS)

  # Single-lag HF exo
  gas_lag1 <- c(NA_real_, daily$gas_log[1:(n - 1L)])
  brt_lag1 <- c(NA_real_, daily$brent_log[1:(n - 1L)])

  # LF (daily, publication-aware LOCF)
  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- c(
    pmi_mfg_de = PUB_LAG_SURVEY, pmi_serv_de = PUB_LAG_SURVEY, ifo_bci_de = PUB_LAG_SURVEY,
    pmi_mfg_it = PUB_LAG_SURVEY, pmi_serv_it = PUB_LAG_SURVEY, istat_bci_it = PUB_LAG_SURVEY,
    Consumer_Goods = PUB_LAG_IPI, Manufacturing = PUB_LAG_IPI, Energy = PUB_LAG_IPI
  )
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  # Surveys: LOCF over missing monthly releases (original standard_midas
  # convention); IPI: no LOCF (a missing release drops the row).
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(
    dates, monthly, v, lf_pub[[v]],
    locf = (lf_pub[[v]] == PUB_LAG_SURVEY))

  # Calendar controls: 6 weekday + 11 month + 1 holiday = 18 columns
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
                         sprintf("r_midas_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS,
                          pre_lag_max = max(c(AR_LAGS, Q_HF, Q_AR)),
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
    sigma2_vec <- rep(NA_real_, n_orig)

    for (oi in seq_along(origins)) {
      t0 <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      win_rows <- win_lo:t0
      ywin <- y[win_rows]

      Z_list <- list(y_lags = ar_mat[win_rows, , drop = FALSE])

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
      X_lin <- do.call(cbind, lin_cols)
      colnames(X_lin) <- lin_names

      keep_train <- stats::complete.cases(X_lin) &
                    rows_all_finite(Z_list$y_lags) &
                    is.finite(ywin)
      if (sum(keep_train) < ncol(X_lin) + 3L) next

      Xlin_fit <- X_lin[keep_train, , drop = FALSE]
      yfit     <- ywin[keep_train]
      Z_fit    <- list(y_lags = Z_list$y_lags[keep_train, , drop = FALSE])

      # Discounted (weighted) LS via row scaling: WLS = OLS on sqrt-weighted
      # rows.  Coefficients are then used on UNscaled forecast rows as usual.
      if (DISCOUNT_LAMBDA < 1) {
        sw <- sqrt(DISCOUNT_LAMBDA ^ (t0 - win_rows))[keep_train]
        Xlin_fit <- Xlin_fit * sw
        yfit     <- yfit * sw
        Z_fit    <- lapply(Z_fit, function(M) M * sw)
      }

      if (!is.null(theta_warm) && length(theta_warm) != 2L) theta_warm <- NULL

      fit <- fit_r_midas(yfit, Z_fit, Xlin_fit,
                         theta_init  = theta_warm,
                         optim_maxit = 2000L,
                         optim_factr = 5e8,
                         n_starts    = N_STARTS_R_MIDAS,
                         lag_sets    = AR_LAG_SETS)
      if (is.null(fit)) next
      theta_warm <- as.numeric(t(fit$theta))
      sigma2_vec[oi] <- fit$sigma2

      coef_mat[oi, "delta_y"] <- fit$delta["y_lags"]
      if (length(lf_used) > 0L) {
        for (v in lf_used) {
          b <- fit$beta_lin[v]
          if (!is.na(b)) coef_mat[oi, v] <- as.numeric(b)
        }
      }

      lin_row_parts <- list(X_base[t0 + 1L, , drop = TRUE])
      lin_row_names <- base_cols
      if (length(lf_used) > 0L) {
        # Low-frequency information must be known at the forecast origin.
        lf_forecast <- lf_mat[t0, lf_used, drop = FALSE]
        lin_row_parts[[length(lin_row_parts) + 1L]] <- drop(lf_forecast)
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
      X_lin_row <- unlist(lin_row_parts)
      names(X_lin_row) <- lin_row_names
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
        hf_block_names  = character(0),
        h_max           = H_MAX,
        holiday_dates   = hol_set,
        lag_sets        = AR_LAG_SETS
      )
    }
    list(yhat_mat = yhat_mat, coef_mat = coef_mat, sigma2 = sigma2_vec)
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
    # Feasible origin-by-origin inverse-asinh correction using the residual
    # variance estimated inside that origin's rolling estimation window.
    yhat_bc_mat <- vapply(seq_along(spec_ids), function(si) {
      mu <- spec_results[[si]]$yhat_mat[, h]
      s2 <- spec_results[[si]]$sigma2
      sinh(mu) * exp(s2 / 2)
    }, numeric(n_orig))
    colnames(yhat_bc_mat) <- paste0(spec_ids, "__bc")
    out <- data.frame(
      origin_date = dates[origins],
      target_date = dates[origins + h],
      y_actual    = daily$elec_price[origins + h],
      stringsAsFactors = FALSE
    )
    out <- cbind(out, as.data.frame(yhat_h_mat, stringsAsFactors = FALSE),
                 as.data.frame(yhat_bc_mat, stringsAsFactors = FALSE))
    out <- out[!is.na(out$y_actual), , drop = FALSE]
    fp  <- file.path(FORECAST_DIR,
                     sprintf("r_midas_%s_h%02d.csv", cc, h))
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
  coef_path <- file.path(COEF_DIR, sprintf("r_midas_%s.csv", cc))
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d, cols = %d)\n",
              coef_path, nrow(coef_df), ncol(coef_df)))
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_r_midas(cc_)
cat("\nR-MIDAS standard (Almon on y, linear HF) done.\n")
