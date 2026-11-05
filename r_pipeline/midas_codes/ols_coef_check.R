#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — PLAIN-OLS COEFFICIENT CHECK OF EVERY REGRESSOR
#
# Run BEFORE any MIDAS estimation: are the regressors sensible when they
# enter a plain linear model?  For each country, the "kitchen-sink" ARX
#
#   y_t = a + calendar_t + sum_{j in {1,2,3,7}} phi_j y_{t-j}
#         + g gas_lag1 + b brent_lag1 + sum_k beta_k LF_k,t + e_t,
#   y = asinh((P - a) / b)   (normalised asinh; VST_MODE=raw -> plain asinh)
#
# (LF_k,t = release-aware daily value of monthly indicator k) is estimated
# by OLS and the coefficients are reported in four ways:
#
#   A. full sample and the three subperiods (pre <= 2020, crisis 2021-22,
#      post >= 2023), with Newey-West HAC standard errors (Bartlett kernel,
#      NW_LAG lags) — are signs / magnitudes stable?  (a, b) of the
#      transform are taken over the sample of each regression.
#   B. each fuel / LF variable ALONE, added to AR(1,2,3,7) + calendar
#      (full sample, HAC) — marginal relevance without collinearity among
#      the indicators;
#   C. rolling WINDOW_DAYS-day OLS of the kitchen sink, with the window's own
#      (a, b): median, 10th and 90th percentile of each coefficient across
#      windows and the share of windows with |t| > 1.96 (classical OLS t);
#   D. standardised coefficients (coef x sd(regressor) / sd(y)) for the full
#      sample — comparable economic magnitudes.
#
# Second block (for the R_MIDAS_L7 comparison): AR(1..7) + calendar with all
# seven lags free, full sample (HAC) and rolling medians, so the implied
# Almon weights of r_midas_l7.R can be compared with the free lag profile.
#
# OUTPUT: results/ols_coef_check_<cc>.csv (one row per regressor),
#         results/ols_ar7_check_<cc>.csv (free AR(1..7) profile) and console
#         tables.  Nothing downstream depends on this stage.
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

DATA_DIR   <- file.path(script_dir, "data")
RESULT_DIR <- file.path(script_dir, "results")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

NW_LAG <- suppressWarnings(as.integer(Sys.getenv("NW_LAG", "10")))
if (is.na(NW_LAG) || NW_LAG < 0L) NW_LAG <- 10L

PERIODS <- list(
  full   = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
  pre    = c(as.Date("1900-01-01"), as.Date("2020-12-31")),
  crisis = c(as.Date("2021-01-01"), as.Date("2022-12-31")),
  post   = c(as.Date("2023-01-01"), as.Date("2100-01-01"))
)

# Newey-West (Bartlett) HAC covariance of OLS coefficients.
nw_vcov <- function(X, e, L) {
  n  <- nrow(X)
  Xe <- X * e
  S  <- crossprod(Xe)
  if (L > 0L) for (l in seq_len(L)) {
    w <- 1 - l / (L + 1)
    G <- crossprod(Xe[(l + 1L):n, , drop = FALSE], Xe[1:(n - l), , drop = FALSE])
    S <- S + w * (G + t(G))
  }
  XtXi <- solve(crossprod(X))
  XtXi %*% S %*% XtXi
}

# OLS of the transformed price on X; the AR-lag columns (lag_idx) of X hold
# RAW lagged prices and are transformed with the same (a, b) as the response.
ols_hac <- function(X, P, lag_idx, L = NW_LAG) {
  p  <- vst_params(P)
  y  <- vst(P, p)
  X[, lag_idx] <- vst(X[, lag_idx, drop = FALSE], p)
  ok <- rows_all_finite(X) & is.finite(y)
  Xk <- X[ok, , drop = FALSE]; yk <- y[ok]
  qrf <- qr(Xk)
  if (qrf$rank < ncol(Xk)) return(NULL)
  b <- qr.coef(qrf, yk)
  e <- yk - drop(Xk %*% b)
  V <- nw_vcov(Xk, e, L)
  se <- sqrt(diag(V))
  list(coef = b, se = se, t = b / se, n = length(yk), a = p$a, b_scale = p$b,
       r2 = 1 - sum(e^2) / sum((yk - mean(yk))^2), sigma = sd(e),
       sdx = apply(Xk, 2, sd), sdy = sd(yk))
}

# Rolling OLS: per-window (a, b); returns coefficient and t matrices.
rolling_ols <- function(X, P, lag_idx, report_vars, origins) {
  p_ <- length(report_vars)
  B <- matrix(NA_real_, length(origins), p_, dimnames = list(NULL, report_vars))
  Tt <- B
  for (oi in seq_along(origins)) {
    t0 <- origins[oi]; wr <- (t0 - WINDOW_DAYS + 1L):t0
    pw <- vst_params(P[wr])
    Xw <- X[wr, , drop = FALSE]; Xw[, lag_idx] <- vst(Xw[, lag_idx, drop = FALSE], pw)
    yw <- vst(P[wr], pw)
    ok <- rows_all_finite(Xw) & is.finite(yw)
    if (sum(ok) < ncol(Xw) + 5L) next
    qrf <- qr(Xw[ok, , drop = FALSE])
    if (qrf$rank < ncol(Xw)) next
    b <- qr.coef(qrf, yw[ok])
    e <- yw[ok] - drop(Xw[ok, , drop = FALSE] %*% b)
    s2 <- sum(e^2) / (sum(ok) - ncol(Xw))
    Ri <- backsolve(qr.R(qrf), diag(ncol(Xw)))
    se <- sqrt(s2 * rowSums(Ri^2))[order(qrf$pivot)]
    names(se) <- colnames(Xw)
    B[oi, ]  <- b[report_vars]
    Tt[oi, ] <- (b / se)[report_vars]
  }
  list(B = B, Tt = Tt)
}

run_check <- function(country) {
  cc      <- if (country == "Germany") "de" else "it"
  cat(sprintf("\n=== OLS coefficient check — %s (HAC lag %d, vst=%s) ===\n", country, NW_LAG, VST_MODE))
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily; monthly <- dat$monthly
  dates   <- daily$date; n <- nrow(daily); P <- daily$elec_price

  Plag <- build_ar_lags(P, AR_LAGS)                 # raw lagged prices
  gas_lag1 <- c(NA_real_, daily$gas_log[1:(n - 1L)])
  brt_lag1 <- c(NA_real_, daily$brent_log[1:(n - 1L)])
  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_mat <- matrix(NA_real_, n, length(lf_names), dimnames = list(NULL, lf_names))
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(dates, monthly, v, LF_PUB_LAG[[v]],
                                                   locf = (LF_PUB_LAG[[v]] == PUB_LAG_SURVEY))
  cal <- build_calendar_dummies(dates, country)
  X_base <- cbind(intercept = 1, cal, Plag)
  X_exo  <- cbind(gas_lag1 = gas_lag1, brent_lag1 = brt_lag1, lf_mat)
  X_all  <- cbind(X_base, X_exo)
  lag_idx_all  <- match(colnames(Plag), colnames(X_all))
  lag_idx_base <- match(colnames(Plag), colnames(X_base))
  # weekday dummies reported too (Sunday = reference): the weekday profile
  # should be clearly positive vs Sunday once real weekends are in the data
  report_vars <- c(grep("^wd_", colnames(cal), value = TRUE), "holiday",
                   colnames(Plag), colnames(X_exo))
  res <- data.frame(variable = report_vars, stringsAsFactors = FALSE)

  # ---- A. full sample + subperiods, HAC ------------------------------------
  for (pn in names(PERIODS)) {
    sel <- dates >= PERIODS[[pn]][1] & dates <= PERIODS[[pn]][2]
    f <- ols_hac(X_all[sel, , drop = FALSE], P[sel], lag_idx_all)
    if (is.null(f)) next
    res[[paste0("coef_", pn)]] <- as.numeric(f$coef[report_vars])
    res[[paste0("t_",    pn)]] <- as.numeric(f$t[report_vars])
    if (pn == "full") {
      cat(sprintf("  kitchen-sink OLS, full sample: n = %d, R2 = %.4f, resid sd = %.4f, vst a = %.2f b = %.2f\n",
                  f$n, f$r2, f$sigma, f$a, f$b_scale))
      res$std_coef_full <- as.numeric(f$coef[report_vars]) * f$sdx[report_vars] / f$sdy   # D
    } else cat(sprintf("  subperiod %-6s: n = %d, R2 = %.4f\n", pn, f$n, f$r2))
  }

  # ---- B. each exogenous variable alone on AR + calendar -------------------
  res$coef_alone <- NA_real_; res$t_alone <- NA_real_
  for (v in colnames(X_exo)) {
    f <- ols_hac(cbind(X_base, X_exo[, v, drop = FALSE]), P, lag_idx_base)
    if (is.null(f)) next
    res$coef_alone[res$variable == v] <- as.numeric(f$coef[v])
    res$t_alone[res$variable == v]    <- as.numeric(f$t[v])
  }
  fb <- ols_hac(X_base, P, lag_idx_base)
  res$coef_alone[res$variable %in% colnames(Plag)] <- as.numeric(fb$coef[colnames(Plag)])
  res$t_alone[res$variable %in% colnames(Plag)]    <- as.numeric(fb$t[colnames(Plag)])

  # ---- C. rolling-window OLS of the kitchen sink ---------------------------
  origins <- make_origins(n, WINDOW_DAYS, pre_lag_max = max(AR_LAGS), h_max = H_MAX)
  ro <- rolling_ols(X_all, P, lag_idx_all, report_vars, origins)
  res$roll_median <- apply(ro$B, 2, median, na.rm = TRUE)
  res$roll_p10    <- apply(ro$B, 2, quantile, probs = 0.10, na.rm = TRUE)
  res$roll_p90    <- apply(ro$B, 2, quantile, probs = 0.90, na.rm = TRUE)
  res$roll_share_sig <- colMeans(abs(ro$Tt) > 1.96, na.rm = TRUE)
  res$roll_share_pos <- colMeans(ro$B > 0, na.rm = TRUE)

  cat(sprintf("\n  %-16s %9s %7s | %9s %7s | %8s %8s %8s | %8s %8s %8s %6s %6s\n",
              "variable", "coef_full", "t_HAC", "alone", "t_HAC",
              "pre", "crisis", "post", "roll_med", "p10", "p90", "%sig", "%pos"))
  for (i in seq_len(nrow(res))) {
    cat(sprintf("  %-16s %9.4f %7.2f | %9.4f %7.2f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f %5.0f%% %5.0f%%\n",
                res$variable[i], res$coef_full[i], res$t_full[i],
                res$coef_alone[i], res$t_alone[i],
                res$coef_pre[i], res$coef_crisis[i], res$coef_post[i],
                res$roll_median[i], res$roll_p10[i], res$roll_p90[i],
                100 * res$roll_share_sig[i], 100 * res$roll_share_pos[i]))
  }
  res$country <- cc
  out <- file.path(RESULT_DIR, sprintf("ols_coef_check_%s.csv", cc))
  write.csv(res, out, row.names = FALSE)
  cat(sprintf("  wrote %s\n", out))

  # ---- AR(1..7) + calendar, all lags free (reference for R_MIDAS_L7) -------
  Plag7 <- build_ar_lags(P, 1:7)
  X7 <- cbind(intercept = 1, cal, Plag7)
  li7 <- match(colnames(Plag7), colnames(X7))
  f7 <- ols_hac(X7, P, li7)
  ro7 <- rolling_ols(X7, P, li7, colnames(Plag7),
                     make_origins(n, WINDOW_DAYS, pre_lag_max = 7L, h_max = H_MAX))
  r7 <- data.frame(country = cc, lag = 1:7,
                   coef_full = as.numeric(f7$coef[colnames(Plag7)]),
                   t_full = as.numeric(f7$t[colnames(Plag7)]),
                   roll_median = apply(ro7$B, 2, median, na.rm = TRUE),
                   roll_p10 = apply(ro7$B, 2, quantile, probs = 0.10, na.rm = TRUE),
                   roll_p90 = apply(ro7$B, 2, quantile, probs = 0.90, na.rm = TRUE),
                   roll_share_sig = colMeans(abs(ro7$Tt) > 1.96, na.rm = TRUE),
                   stringsAsFactors = FALSE)
  cat(sprintf("\n  free AR(1..7) + calendar, full sample (n = %d, R2 = %.4f):\n", f7$n, f7$r2))
  cat("   lag:        ", paste(sprintf("%7d", 1:7), collapse = " "), "\n")
  cat("   coef (full):", paste(sprintf("%7.3f", r7$coef_full), collapse = " "), "  sum =", round(sum(r7$coef_full), 3), "\n")
  cat("   t (HAC):    ", paste(sprintf("%7.2f", r7$t_full), collapse = " "), "\n")
  cat("   rolling med:", paste(sprintf("%7.3f", r7$roll_median), collapse = " "), "  sum =", round(sum(r7$roll_median), 3), "\n")
  out7 <- file.path(RESULT_DIR, sprintf("ols_ar7_check_%s.csv", cc))
  write.csv(r7, out7, row.names = FALSE)
  cat(sprintf("  wrote %s\n", out7))
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL") c("Germany", "Italy") else country_arg
for (cc_ in countries) run_check(cc_)
cat("\nOLS coefficient check done.\n")
