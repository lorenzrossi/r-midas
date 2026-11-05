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
# NLS protocol: 7-start L-BFGS-B per rolling window per spec. Start 1 is the
# previous origin's warm-start theta (or (0,0) initially); the remaining
# starts are deterministic canonical Almon shapes supplied by
# fit_r_midas_fast.R (flat, two U-shapes, fast/slow decay, hump). The
# lowest-SSE admissible optimum is retained.
# Theta box: [-8,8] x [-8, RMIDAS_THETA2_MAX] (default 8 = theta2 free, so
# U-shaped weight profiles with a weekly revival at lag 7 are admissible;
# RMIDAS_THETA2_MAX=0 reproduces the earlier log-concave restriction).
#
# Iterated multi-step (option c.1):
#   * y AR-Almon: history_y_aug fed back through fixed θ_y Almon weights
#   * LF columns: FROZEN at the ORIGIN row t0 (latest released value);
#     HF fuel lag1 columns: FROZEN at row (origin + 1), i.e. the fuel
#     price observed at the origin
#   * Calendar dummies (weekday, month, holiday): REFRESHED at target row
#
# OUTPUTS:
#   forecasts/r_midas_<cc>_h<h>.csv              forecasts per horizon
#   forecasts/r_midas_spec_dictionary_<cc>.csv   spec → variable membership
#   forecasts/r_midas_<cc>_h<h>_sd.csv           predictive sd (asinh scale)
#                                                of the h-step Gaussian
#                                                density, per spec (input of
#                                                density_eval.R)
#   results/coefficients/r_midas_<cc>.csv        rolling-window coefficients
#                                                (wide format: origin_date,
#                                                spec_id, delta_y, theta1,
#                                                theta2, phi_1 phi_2 phi_3
#                                                phi_7 = implied AR coefs
#                                                delta*b(j;theta),
#                                                min_root_mod = stability
#                                                diagnostic (> 1 = stationary),
#                                                at_bound / flag_delta_small /
#                                                flag_degenerate / convergence
#                                                = optimizer + identification
#                                                flags, <each LF>)
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

# Family parametrisation (r_midas_l7.R sets these before sourcing this file):
#   RMIDAS_LAGS   lag set of the exp-Almon block ("1,2,3,7" = standard R-MIDAS
#                 on the benchmark's gappy set; "1:7" = consecutive lags 1..7,
#                 family R_MIDAS_L7)
#   RMIDAS_PREFIX file prefix of every output (r_midas | r_midas_l7)
#   RMIDAS_FAMILY label used in the summaries
.parse_lags <- function(x) if (grepl(":", x)) eval(parse(text = x)) else
  as.integer(strsplit(x, ",")[[1]])
RM_LAGS   <- as.integer(.parse_lags(Sys.getenv("RMIDAS_LAGS", "1,2,3,7")))
RM_PREFIX <- Sys.getenv("RMIDAS_PREFIX", "r_midas")
RM_FAMILY <- Sys.getenv("RMIDAS_FAMILY", "R_MIDAS")

DATA_DIR     <- file.path(script_dir, "data")
FORECAST_DIR <- file.path(script_dir, "forecasts")
COEF_DIR     <- file.path(script_dir, "results", "coefficients")
dir.create(FORECAST_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(COEF_DIR,     showWarnings = FALSE, recursive = TRUE)

Q_AR             <- Q_HF
SPEC_LIMIT       <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L
# 7 starts = warm start + 6 canonical shapes (flat, two U-shapes, fast and
# slow decay, hump) so that both log-convex restarts are always tried.
N_STARTS_R_MIDAS <- 7L

# Exponentially discounted (weighted) least squares, applied IDENTICALLY to
# both countries so the methodology stays symmetric.  Rows of the estimation
# window are weighted lambda^(t0 - t) (weight 1 on the most recent day).
# RMIDAS_DISCOUNT = 1 (default) disables it; e.g. 0.995 halves the weight
# roughly every 139 days.
# One-step Huber-reweighted estimation (robust to crisis-day outliers):
# RMIDAS_HUBER=1 enables it; RMIDAS_HUBER_K sets the tuning constant
# (default 1.345 = 95% Gaussian efficiency).  Applied identically to both
# countries.  See huber_weights() in common_utils.R.
HUBER_ON <- Sys.getenv("RMIDAS_HUBER", "0") %in% c("1", "TRUE", "true", "yes")
HUBER_K  <- suppressWarnings(as.numeric(Sys.getenv("RMIDAS_HUBER_K", "1.345")))
if (!is.finite(HUBER_K) || HUBER_K <= 0) HUBER_K <- 1.345

DISCOUNT_LAMBDA <- suppressWarnings(as.numeric(Sys.getenv("RMIDAS_DISCOUNT", "1")))
if (!is.finite(DISCOUNT_LAMBDA) || DISCOUNT_LAMBDA <= 0 || DISCOUNT_LAMBDA > 1)
  DISCOUNT_LAMBDA <- 1

# Upper bound on theta2 of the exp-Almon exponent.  DEFAULT 8: the box is
# [-8,8]^2 (bounded for numerical stability only), admitting log-CONVEX
# (U-shaped) profiles — high lag 1, dip at lags 2-3, weekly revival at lag 7
# — which the free-AR dynamics of BOTH countries favour (the former bound
# theta2 <= 0 was binding in ~100% of Italian windows and zeroed the weekly
# lag).  RMIDAS_THETA2_MAX=0 reproduces the earlier log-concave restriction
# (flat / decaying / single hump).  Applied identically to both countries.
THETA2_MAX <- suppressWarnings(as.numeric(Sys.getenv("RMIDAS_THETA2_MAX", "8")))
if (!is.finite(THETA2_MAX)) THETA2_MAX <- 8
THETA_LOWER <- c(-8, -8)
THETA_UPPER <- c(8, THETA2_MAX)
BOUND_TOL   <- 1e-3      # |theta - bound| below this counts as "at the bound"
DELTA_SMALL <- 0.05      # |delta_y| below this: theta weakly identified
W_DEGENERATE <- 0.99     # max Almon weight above this: all mass on one lag

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_r_midas <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== %s (Almon on y lags {%s}, linear HF) — %s (n_cores=%d, n_starts=%d, theta2_max=%g, vst=%s) ===\n",
              RM_FAMILY, paste(RM_LAGS, collapse = ","), country, n_cores,
              N_STARTS_R_MIDAS, THETA2_MAX, VST_MODE))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  P       <- daily$elec_price              # RAW prices; normalised asinh per window

  # Almon AR block over the lag set RM_LAGS (raw lagged prices, transformed
  # per origin with that window's (a, b))
  ar_mat <- build_ar_lags(P, RM_LAGS)
  Z_full <- list(y_lags = ar_mat)
  AR_LAG_SETS <- list(y_lags = RM_LAGS)

  # Single-lag HF exo
  gas_lag1 <- c(NA_real_, daily$gas_log[1:(n - 1L)])
  brt_lag1 <- c(NA_real_, daily$brent_log[1:(n - 1L)])

  # LF (daily, release-aware: each monthly value becomes usable on its
  # release date and is carried forward until the next release — see
  # make_daily_lf in common_utils.R)
  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- LF_PUB_LAG
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  # The locf flag applies to the LEGACY calendar-lag mode only (surveys:
  # carry over missing monthly values; IPI: drop the row).  In the default
  # release-aware mode the latest released value is carried forward for
  # every variable by construction.
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
                         sprintf("%s_spec_dictionary_%s.csv", RM_PREFIX, cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS,
                          # Max lag actually used by this model: the Almon
                          # block runs on RM_LAGS and the fuel columns on lag
                          # 1, so max(RM_LAGS) pre-sample days suffice.
                          pre_lag_max = max(RM_LAGS),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))

  # Normalised-asinh location / scale per origin (identical to the values
  # written by ar_benchmark.R: same windows, same prices).
  VP <- vst_params_by_origin(P, origins, WINDOW_DAYS)
  write.csv(data.frame(origin_date = dates[origins], a = VP$a, b = VP$b),
            file.path(FORECAST_DIR, sprintf("vst_params_%s.csv", cc)), row.names = FALSE)

  # Per-window record: total AR loading, Almon shape parameters, the implied
  # AR coefficients phi_j = delta * b(j; theta) at every lag of RM_LAGS, the
  # stability diagnostic, optimizer / identification flags, and the LF
  # loadings.
  phi_cols  <- paste0("phi_", RM_LAGS)
  diag_cols <- c("theta1", "theta2", phi_cols, "min_root_mod",
                 "at_bound", "flag_delta_small", "flag_degenerate",
                 "convergence")
  coef_cols <- c("delta_y", diag_cols, lf_names)

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
    # Horizon-specific forecast-error variances (referee point 3): row =
    # origin, column = step 1..H_MAX; filled via iterated_error_variance().
    sigma2_mat <- matrix(NA_real_, n_orig, H_MAX)

    # ==== MAIN ROLLING-ORIGIN LOOP =========================================
    # One iteration = one forecast origin t0: estimate the model on the
    # WINDOW_DAYS days ending at t0, then forecast t0+1 ... t0+H_MAX.
    for (oi in seq_along(origins)) {
      t0 <- origins[oi]                     # row index of the forecast origin
      win_lo <- t0 - WINDOW_DAYS + 1L      # first row of the rolling window
      win_rows <- win_lo:t0                # estimation window (ends AT t0)
      p_o  <- list(a = VP$a[oi], b = VP$b[oi])   # this window's (a, b)
      ywin <- vst(P[win_rows], p_o)        # response: asinh((P - a)/b) in window

      # Almon block: the window rows of the lagged-price matrix (columns =
      # lags RM_LAGS), on the window's transformed scale.  The NLS weights
      # these columns with the exp-Almon weights b(j; theta).
      Z_list <- list(y_lags = vst(ar_mat[win_rows, , drop = FALSE], p_o))

      # Linear part of the design, assembled column-block by column-block:
      # intercept + 18 calendar dummies (always) ...
      lin_cols  <- list(X_base[win_rows, , drop = FALSE])
      lin_names <- base_cols
      # ... + one publication-aware daily column per LF variable in this spec
      if (length(lf_used) > 0L) {
        lin_cols[[length(lin_cols) + 1L]] <- lf_mat[win_rows, lf_used, drop = FALSE]
        lin_names <- c(lin_names, lf_used)
      }
      # ... + one-day lags of the fuel prices when the spec includes them
      # (linear, per model.tex: nothing to compress in a single lag)
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
      X_lin <- do.call(cbind, lin_cols)    # window x p linear design
      colnames(X_lin) <- lin_names

      # Keep only rows where response, every linear column and every Almon
      # lag are observed (drops e.g. days before the first LF release).
      keep_train <- stats::complete.cases(X_lin) &
                    rows_all_finite(Z_list$y_lags) &
                    is.finite(ywin)
      if (sum(keep_train) < ncol(X_lin) + 3L) next   # too few obs: skip origin

      Xlin_fit <- X_lin[keep_train, , drop = FALSE]
      yfit     <- ywin[keep_train]
      Z_fit    <- list(y_lags = Z_list$y_lags[keep_train, , drop = FALSE])
      # Keep UNweighted copies: the bias-correction variance must be the
      # innovation variance of the original transformed-price process, so it
      # is recomputed on these raw rows after the (possibly weighted)
      # estimation (referee point 4).
      Xlin_raw <- Xlin_fit
      yfit_raw <- yfit
      Zraw     <- Z_fit$y_lags

      # Discounted (weighted) LS via row scaling: WLS = OLS on sqrt-weighted
      # rows.  Coefficients are then used on UNscaled forecast rows as usual.
      if (DISCOUNT_LAMBDA < 1) {
        sw <- sqrt(DISCOUNT_LAMBDA ^ (t0 - win_rows))[keep_train]
        Xlin_fit <- Xlin_fit * sw
        yfit     <- yfit * sw
        Z_fit    <- lapply(Z_fit, function(M) M * sw)
      }

      # Warm start: theta from the PREVIOUS origin of this spec (windows
      # overlap by all but one day, so the optimum moves slowly).  Guard
      # against a stale cache of the wrong dimension.
      if (!is.null(theta_warm) && length(theta_warm) != 2L) theta_warm <- NULL

      # THE ESTIMATION.  Concentrated NLS (fit_r_midas_fast.R): for each
      # candidate theta the linear coefficients are profiled out by OLS and
      # the concentrated SSE is minimised over theta = (theta1, theta2) only,
      # by L-BFGS-B on the box [-8,8] x [-8, THETA2_MAX] (stability of the
      # implied AR dynamics is checked separately via min_root_mod below),
      # from the warm start + canonical shape restarts.
      # lag_sets tells the estimator the Almon exponent uses j in RM_LAGS.
      fit <- fit_r_midas(yfit, Z_fit, Xlin_fit,
                         theta_init  = theta_warm,
                         theta_upper_block = THETA_UPPER,
                         optim_maxit = 2000L,
                         optim_factr = 5e8,
                         n_starts    = N_STARTS_R_MIDAS,
                         lag_sets    = AR_LAG_SETS)
      if (is.null(fit)) next               # rank-deficient / degenerate window

      # ---- optional one-step Huber reweighting (RMIDAS_HUBER=1) -----------
      # NOTE: when RMIDAS_DISCOUNT is also active, the Huber residuals are
      # computed on the DISCOUNTED observations (robustness applied after
      # discounting) - state this explicitly in the methodology if both
      # options are enabled together.
      # Stage 1 = the fit above.  Compute its in-window residuals, build the
      # Huber weights (k = RMIDAS_HUBER_K * MAD/0.6745), and refit ONCE with
      # sqrt(w)-scaled rows (WLS = OLS on scaled data; same trick as the
      # discount).  Warm-started single-start refit, so cost ~ one extra fit.
      if (HUBER_ON) {
        w_alm  <- exp_almon_weights_j(RM_LAGS, fit$theta["y_lags", ])
        z_in   <- as.numeric(Z_fit$y_lags %*% w_alm)
        yhat_in <- as.numeric(cbind(Xlin_fit, z_in) %*%
                              c(fit$beta_lin, fit$delta["y_lags"]))
        hw <- huber_weights(yfit - yhat_in, HUBER_K)
        if (!is.null(hw)) {
          swh <- sqrt(hw)
          fit2 <- fit_r_midas(yfit * swh,
                              lapply(Z_fit, function(M) M * swh),
                              Xlin_fit * swh,
                              theta_init  = as.numeric(t(fit$theta)),
                              theta_upper_block = THETA_UPPER,
                              optim_maxit = 2000L,
                              optim_factr = 5e8,
                              n_starts    = 1L,
                              lag_sets    = AR_LAG_SETS)
          if (!is.null(fit2)) fit <- fit2
        }
      }
      theta_warm <- as.numeric(t(fit$theta))   # cache theta for next origin

      # ---- bias-correction variance (referee points 3 and 4) --------------
      # (a) Recompute the one-step residual variance on the ORIGINAL,
      #     unweighted observations at the final coefficients (identical to
      #     fit$sigma2 when discounting and Huber are off).
      w_bc   <- exp_almon_weights_j(RM_LAGS, fit$theta["y_lags", ])
      z_bc   <- as.numeric(Zraw %*% w_bc)
      res_bc <- yfit_raw - as.numeric(cbind(Xlin_raw, z_bc) %*%
                                      c(fit$beta_lin, fit$delta["y_lags"]))
      # (b) Propagate to each horizon through the fitted AR dynamics: the
      #     implied AR coefficients are phi_j = delta * b(j; theta_hat) at
      #     lags j in RM_LAGS; sigma2_h = sigma2 * sum_{s<h} psi_s^2.
      phi_bc <- as.numeric(fit$delta["y_lags"]) * w_bc
      sigma2_mat[oi, ] <- iterated_error_variance(RM_LAGS, phi_bc,
                                                  mean(res_bc^2), H_MAX)

      # Record the rolling coefficients: delta = total AR loading (the Almon
      # weights sum to 1, so delta is the sum of the four implied AR
      # coefficients) and the free loading of each LF variable.  Also record
      # the stability diagnostic: the minimum root modulus of the
      # characteristic polynomial of the implied AR coefficients phi_j =
      # delta_y * b_j(theta_hat) (> 1 = stationary window).
      th <- as.numeric(fit$theta["y_lags", ])
      coef_mat[oi, "delta_y"]      <- fit$delta["y_lags"]
      coef_mat[oi, c("theta1", "theta2")] <- th
      coef_mat[oi, phi_cols]       <- phi_bc
      coef_mat[oi, "min_root_mod"] <- ar_min_root_modulus(RM_LAGS, phi_bc)
      # Optimizer at the box boundary (a warm start can keep later windows in
      # a corner, so this is logged per window and summarised per country).
      coef_mat[oi, "at_bound"] <- as.integer(
        any(abs(th - THETA_LOWER) < BOUND_TOL | abs(th - THETA_UPPER) < BOUND_TOL))
      # Identification flags: theta is not identified when delta ~ 0, and the
      # SSE surface is flat when all the weight sits on a single lag.
      coef_mat[oi, "flag_delta_small"] <- as.integer(abs(fit$delta[["y_lags"]]) < DELTA_SMALL)
      coef_mat[oi, "flag_degenerate"]  <- as.integer(max(w_bc) > W_DEGENERATE)
      coef_mat[oi, "convergence"]      <- fit$convergence
      if (length(lf_used) > 0L) {
        for (v in lf_used) {
          b <- fit$beta_lin[v]
          if (!is.na(b)) coef_mat[oi, v] <- as.numeric(b)
        }
      }

      # ---- FORECAST-ROW CONSTRUCTION (information set at origin t0) -------
      # Build the single linear row used at every iteration step:
      #   calendar part taken at t0+1 (it is refreshed to the true target
      #   date inside the forecaster anyway);
      lin_row_parts <- list(X_base[t0 + 1L, , drop = TRUE])
      lin_row_names <- base_cols
      if (length(lf_used) > 0L) {
        # LF block FROZEN at row t0: the latest monthly value whose release
        # date is <= the origin (real-time information set; row t0+1 could
        # embed a release that appears only on the target day).
        lf_forecast <- lf_mat[t0, lf_used, drop = FALSE]
        lin_row_parts[[length(lin_row_parts) + 1L]] <- drop(lf_forecast)
        lin_row_names <- c(lin_row_names, lf_used)
      }
      # Fuel lag-1 columns at row t0+1 = the fuel price OBSERVED at t0
      # (the lag-1 column at row t is the value of day t-1), then frozen
      # across all iteration steps: future fuel prices are unknown at t0.
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
      if (anyNA(X_lin_row)) next           # regressor missing at the origin

      # ---- ITERATED MULTI-STEP FORECAST (option c.1) ----------------------
      # One-step model rolled forward h = 1..H_MAX: the AR-Almon block is
      # rebuilt at each step from realised prices where available and from
      # the model's own previous-step forecasts otherwise; calendar dummies
      # are set to the true target dates; LF and fuel terms stay frozen.
      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      yhat_mat[oi, ] <- iterated_forecast_r_midas(
        fit_obj         = fit,
        Z_list_full     = Z_full,
        X_lin_row       = X_lin_row,
        X_lin_colnames  = names(X_lin_row),
        history_y_aug   = vst(P[1:t0], p_o),   # past prices on the origin's scale
        origin_idx      = t0,
        target_dates    = target_dates,
        ar_block_name   = "y_lags",
        hf_block_names  = character(0),
        h_max           = H_MAX,
        holiday_dates   = hol_set,
        lag_sets        = AR_LAG_SETS
      )
    }
    list(yhat_mat = yhat_mat, coef_mat = coef_mat, sigma2 = sigma2_mat)
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec,
                                  n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.2fh) — wall-clock\n", el, el / 3600))

  for (h in HORIZONS) {
    yhat_h_mat <- vapply(seq_along(spec_ids),
                         function(si) vst_inv(spec_results[[si]]$yhat_mat[, h], VP),
                         numeric(n_orig))                   # median forecast, EUR/MWh
    colnames(yhat_h_mat) <- spec_ids
    # Feasible origin-by-origin inverse-transform correction using the
    # residual variance estimated inside that origin's rolling window.
    yhat_bc_mat <- vapply(seq_along(spec_ids), function(si) {
      mu <- spec_results[[si]]$yhat_mat[, h]
      s2 <- spec_results[[si]]$sigma2[, h]   # h-step error variance
      vst_inv(mu, VP, s2)                    # mean forecast, EUR/MWh
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
    keep_rows <- !is.na(out$y_actual)
    out <- out[keep_rows, , drop = FALSE]
    fp  <- file.path(FORECAST_DIR,
                     sprintf("%s_%s_h%02d.csv", RM_PREFIX, cc, h))
    write.csv(out, fp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", fp, nrow(out)))

    # Predictive standard deviation (transformed scale) of the h-step
    # Gaussian density N(mu_h, sigma2_h) per spec: input of the density
    # evaluation (log score / CRPS / ACPS on the price scale).
    sd_mat <- vapply(seq_along(spec_ids),
                     function(si) sqrt(spec_results[[si]]$sigma2[, h]),
                     numeric(n_orig))
    colnames(sd_mat) <- spec_ids
    sd_out <- cbind(data.frame(origin_date = dates[origins],
                               target_date = dates[origins + h],
                               stringsAsFactors = FALSE),
                    as.data.frame(sd_mat, stringsAsFactors = FALSE))
    sd_out <- sd_out[keep_rows, , drop = FALSE]
    fp_sd <- file.path(FORECAST_DIR,
                       sprintf("%s_%s_h%02d_sd.csv", RM_PREFIX, cc, h))
    write.csv(sd_out, fp_sd, row.names = FALSE)
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
  coef_path <- file.path(COEF_DIR, sprintf("%s_%s.csv", RM_PREFIX, cc))
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d, cols = %d)\n",
              coef_path, nrow(coef_df), ncol(coef_df)))

  # Almon-shape / optimizer diagnostics across all (origin, spec) windows.
  pct <- function(x) sprintf("%.1f%%", 100 * mean(x, na.rm = TRUE))
  cat(sprintf(paste0(
    "  theta diagnostics (%d windows): U-shape (theta2 > 0) %s | ",
    "at box bound %s | delta_small %s | degenerate weights %s | ",
    "optim code 1 (iteration limit) %s\n",
    "  median theta = (%.2f, %.2f); median implied phi(", paste(RM_LAGS, collapse = ","), ") = (%s); ",
    "median delta = %.3f; min root modulus < 1 in %s\n"),
    nrow(coef_df), pct(coef_df$theta2 > 0), pct(coef_df$at_bound == 1L),
    pct(coef_df$flag_delta_small == 1L), pct(coef_df$flag_degenerate == 1L),
    pct(coef_df$convergence == 1L),
    median(coef_df$theta1, na.rm = TRUE), median(coef_df$theta2, na.rm = TRUE),
    paste(sprintf("%.3f", vapply(phi_cols, function(cn) median(coef_df[[cn]], na.rm = TRUE),
                                 numeric(1))), collapse = ", "),
    median(coef_df$delta_y, na.rm = TRUE), pct(coef_df$min_root_mod < 1)))
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_r_midas(cc_)
cat(sprintf("\n%s (Almon on y lags {%s}, linear HF) done.\n", RM_FAMILY, paste(RM_LAGS, collapse = ",")))
