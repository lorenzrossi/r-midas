#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — FORECAST COMBINATION, SUBPERIOD EVALUATION AND
#               GIACOMINI-ROSSI FLUCTUATION TEST
#
# Post-processing only: works on the forecast CSVs already written by
# r_midas.R / r_midas_extended.R and ar_benchmark.R.  No re-estimation.
# Applied IDENTICALLY to Germany and Italy (same combination rules, same
# subperiods, same test), so the methodology is symmetric across countries.
#
# What it does, per country, family and horizon:
#
# 1. FORECAST COMBINATION across the 256 specs (Bates-Granger 1969;
#    Bordignon et al. 2013; averaging literature surveyed in Weron 2014):
#      comb_eq      equal-weight mean of all spec forecasts (family only)
#      comb_trim    20%-trimmed mean (robust to bad specs; family only)
#      comb_invmse  pool = specs + AR candidates (COMB_INCLUDE_AR, default
#                   "ar_dum,ar_dum_gas"): adaptive fallback to the pure AR
#      comb_invmse  inverse-MSE weights, computed in REAL TIME: at origin o
#                   only forecast errors whose target date is <= o (i.e.
#                   resolved errors, origin <= o - h) enter the MSE, over a
#                   trailing window of TRAIL_DAYS days; specs need at least
#                   MIN_RESOLVED resolved errors to receive a weight,
#                   otherwise the combination falls back to the equal mean.
#    Combined forecasts are written to forecasts/comb_<family>_<cc>_h<h>.csv
#    so they can be scored like any other model.
#
# 2. SUBPERIOD EVALUATION against the AR benchmark (ar_dum):
#      full      whole evaluation sample
#      pre       target dates 2015-01-01 .. 2020-12-31
#      crisis    target dates 2021-01-01 .. 2022-12-31
#      post      target dates 2023-01-01 .. end
#    RMSE / MAE ratios in levels and on the asinh scale, plus full-sample
#    DM-HLN (negative statistic = model better, same convention as
#    summary_table.R).  Output: results/summary_subperiod_<cc>.csv
#
# 3. GIACOMINI-ROSSI (2010, J. Applied Econometrics) FLUCTUATION TEST of the
#    model against the AR benchmark: rolling sums of the squared-error loss
#    differential d_t = e_bench^2 - e_model^2 (positive = model better) over
#    a TRAILING window of a fraction mu = 0.3 of the sample (each statistic
#    is timestamped by the window END date), studentised with the same
#    rectangular-kernel HAC variance (truncation h-1) used in dm_hln.
#    Output: results/gr_fluctuation_<cc>.csv with the max and min of the
#    fluctuation statistic, the two-sided critical values, and the dates at
#    which the extrema occur.
#    NOTE: critical values below are transcribed from Giacomini & Rossi
#    (2010), Table 1 (two-sided).  Verify against the published table before
#    final use in the paper.
#
# Models evaluated: ar_dum_gas, ar_dum_brent, ar_dum_gas_brent (from the AR
# file), the best FULL-SAMPLE individual spec per family (flagged: selected
# ex post, look-ahead), and the three combinations (real-time, no look-ahead).
#
# Env knobs: COUNTRY_NAME ("Germany" | "Italy" | "ALL"), FAMILIES
# (comma-separated, default "r_midas,r_midas_extended"), TRAIL_DAYS (90),
# MIN_RESOLVED (30), TRIM (0.2), GR_MU (0.3).
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

FORECAST_DIR <- file.path(script_dir, "forecasts")
RESULT_DIR   <- file.path(script_dir, "results")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

FAMILIES     <- strsplit(Sys.getenv("FAMILIES", "r_midas,r_midas_extended"), ",")[[1]]
TRAIL_DAYS   <- as.integer(Sys.getenv("TRAIL_DAYS",   "90"))
TOP_N_SPECS  <- as.integer(Sys.getenv("TOP_N_SPECS",  "0"))  # 0 = ALL specs
MIN_RESOLVED <- as.integer(Sys.getenv("MIN_RESOLVED", "30"))
TRIM         <- as.numeric(Sys.getenv("TRIM",         "0.2"))
GR_MU        <- as.numeric(Sys.getenv("GR_MU",        "0.3"))
BENCH_COL    <- "ar_dum"

# AR candidates added to the inverse-MSE combination POOL (not to the
# benchmark, which stays ar_dum, and not to comb_eq / comb_trim, which remain
# equal/trimmed averages of the R-MIDAS family only).  With the AR forecasts
# in the pool, the real-time inverse-MSE weights migrate towards the pure AR
# whenever the macro block turns harmful (e.g. the 2021-22 crisis) and back
# when it pays off.  No look-ahead: weights use only resolved errors.
# Set COMB_INCLUDE_AR="" to disable.
COMB_INCLUDE_AR <- {
  raw <- Sys.getenv("COMB_INCLUDE_AR", unset = "ar_dum,ar_dum_gas")
  v <- trimws(strsplit(raw, ",")[[1]])
  v[nzchar(v)]
}

PERIODS <- list(
  full   = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
  pre    = c(as.Date("1900-01-01"), as.Date("2020-12-31")),
  crisis = c(as.Date("2021-01-01"), as.Date("2022-12-31")),
  post   = c(as.Date("2023-01-01"), as.Date("2100-01-01"))
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
rmse_ <- function(e) { e <- e[is.finite(e)]; if (!length(e)) NA_real_ else sqrt(mean(e^2)) }
mae_  <- function(e) { e <- e[is.finite(e)]; if (!length(e)) NA_real_ else mean(abs(e)) }

# Real-time inverse-MSE (Bates-Granger) combination.
# Fm: n x S forecast matrix (rows ordered by origin date), y: actuals,
# odates: origin dates.  Errors of the forecast made at origin j resolve at
# target j + h, so at origin i only rows with odates[j] <= odates[i] - h are
# usable; the MSE window is (odates[i] - h - trail, odates[i] - h].
combine_invmse <- function(Fm, y, odates, h,
                           trail = TRAIL_DAYS, min_res = MIN_RESOLVED,
                           return_weights = FALSE) {
  n  <- nrow(Fm)
  E2 <- (Fm - y)^2
  ok <- is.finite(E2)
  E2z <- E2; E2z[!ok] <- 0
  CS <- apply(E2z, 2, cumsum)              # cumulative SSE per spec
  CN <- apply(ok,  2, cumsum)              # cumulative count per spec
  od <- as.numeric(odates)
  out <- rep(NA_real_, n)
  eqm <- rowMeans(Fm, na.rm = TRUE)        # fallback
  W <- if (return_weights) matrix(NA_real_, n, ncol(Fm),
                                  dimnames = list(NULL, colnames(Fm))) else NULL
  for (i in seq_len(n)) {
    fi <- Fm[i, ]
    b <- findInterval(od[i] - h, od)
    a <- findInterval(od[i] - h - trail, od)
    if (b <= a || b < 1L) {                 # no resolved errors yet:
      out[i] <- eqm[i]                      # equal-weight fallback
      if (return_weights) {
        fin <- is.finite(fi)
        if (any(fin)) W[i, fin] <- 1 / sum(fin)
      }
      next
    }
    cnt <- CN[b, ] - if (a >= 1L) CN[a, ] else 0
    sse <- CS[b, ] - if (a >= 1L) CS[a, ] else 0
    mse <- ifelse(cnt >= min_res, sse / pmax(cnt, 1L), NA_real_)
    w   <- 1 / mse
    use <- is.finite(fi) & is.finite(w) & w > 0
    if (any(use)) {
      out[i] <- sum(w[use] * fi[use]) / sum(w[use])
      if (return_weights) W[i, use] <- w[use] / sum(w[use])
    } else {
      out[i] <- eqm[i]
      if (return_weights) {
        fin <- is.finite(fi)
        if (any(fin)) W[i, fin] <- 1 / sum(fin)
      }
    }
  }
  if (return_weights) list(out = out, W = W) else out
}

# Giacomini-Rossi (2010) two-sided fluctuation test, squared-error loss.
# d_t = e_bench^2 - e_model^2  (positive = model better).  Rolling sums over
# m = floor(mu * n) observations, studentised by the full-sample HAC variance
# (rectangular kernel, truncation h-1, as in dm_hln).
GR_CV <- data.frame(
  mu   = seq(0.1, 0.9, by = 0.1),
  cv05 = c(3.393, 3.179, 3.012, 2.890, 2.779, 2.634, 2.560, 2.433, 2.331),
  cv10 = c(3.170, 2.948, 2.766, 2.626, 2.500, 2.356, 2.252, 2.130, 2.020)
)

gr_fluctuation <- function(e_bench, e_model, tdates, h, mu = GR_MU) {
  ok <- is.finite(e_bench) & is.finite(e_model)
  d  <- e_bench[ok]^2 - e_model[ok]^2
  td <- tdates[ok]
  n  <- length(d)
  if (n < 100L) return(NULL)
  m  <- max(20L, floor(mu * n))
  dbar <- mean(d)
  s <- mean((d - dbar)^2)
  if (h > 1L) for (k in 1:(h - 1L))
    s <- s + 2 * sum((d[(k + 1L):n] - dbar) * (d[1:(n - k)] - dbar)) / n
  if (!is.finite(s) || s <= 0) return(NULL)
  cs <- c(0, cumsum(d))
  Ft <- (cs[(m + 1L):(n + 1L)] - cs[1:(n - m + 1L)]) / (sqrt(s) * sqrt(m))
  end_dates <- td[m:n]                     # window ends
  i_max <- which.max(Ft); i_min <- which.min(Ft)
  mu_eff <- m / n
  list(stat_max = Ft[i_max], date_max = end_dates[i_max],
       stat_min = Ft[i_min], date_min = end_dates[i_min],
       mu = mu_eff, m = m, n = n,
       cv05 = approx(GR_CV$mu, GR_CV$cv05, xout = mu_eff, rule = 2)$y,
       cv10 = approx(GR_CV$mu, GR_CV$cv10, xout = mu_eff, rule = 2)$y)
}

# Subperiod + full-sample evaluation of one forecast vector vs the benchmark.
# Per-model evaluation over the full sample and the three subperiods.
#   yhat_bc / ybench_bc : bias-corrected forecasts (optional; NA metrics if
#                         absent), so the bc block is CONSISTENT with
#                         summary_table across every table.
#   variables           : human-readable model composition, written into
#                         every output row (specs: their variable subset;
#                         AR models: their content; combinations: an
#                         explicit description of the weighting scheme).
#   a, b                : normalised-asinh location / scale of each row's
#                         origin window (transformed-scale errors are
#                         asinh((. - a)/b) differences; NULL = plain asinh).
#   Also reports the Giacomini-White (2006) conditional test (instruments:
#   constant and d_{t-h}) in every period, on squared and absolute errors,
#   EUR and transformed scale; d_{t-h} is built on the FULL origin sequence
#   and the period is applied afterwards (keep).
eval_model <- function(label, family, h, yhat, y, ybench, tdates,
                       yhat_bc = NULL, ybench_bc = NULL,
                       variables = "", a = NULL, b = NULL) {
  if (is.null(a)) a <- rep(0, length(y))
  if (is.null(b)) b <- rep(1, length(y))
  e_m <- yhat - y
  e_b <- ybench - y
  ea_m <- vst_z(yhat, a, b) - vst_z(y, a, b)        # model's own transformed scale
  ea_b <- vst_z(ybench, a, b) - vst_z(y, a, b)
  has_bc <- !is.null(yhat_bc) && !is.null(ybench_bc)
  if (has_bc) {
    eb_m <- yhat_bc - y
    eb_b <- ybench_bc - y
  }
  rows <- list()
  for (pn in names(PERIODS)) {
    sel <- tdates >= PERIODS[[pn]][1] & tdates <= PERIODS[[pn]][2]
    ok_lvl <- sel & is.finite(e_m) & is.finite(e_b)
    ok_as  <- sel & is.finite(ea_m) & is.finite(ea_b)
    ok_bc  <- if (has_bc) sel & is.finite(eb_m) & is.finite(eb_b) else
                rep(FALSE, length(sel))
    if (!any(ok_lvl)) next
    # DM-HLN computed within EVERY period (not only the full sample): on a
    # subsample it is the unconditional equal-accuracy test restricted to
    # those target dates, with the same h-1 HAC truncation.  The GR test
    # remains the instrument for the full time-varying picture.
    # Minimum sample length before computing the statistic (referee point
    # 10): with overlapping h-step errors, short subsamples give unreliable
    # HAC variances.  Require at least max(30, 10*h) observations.
    n_min_dm <- max(30L, 10L * h)
    dm_lvl <- if (sum(ok_lvl) >= n_min_dm)
      dm_hln(e_m[ok_lvl], e_b[ok_lvl], h) else
      list(stat = NA_real_, pvalue = NA_real_)
    dm_as  <- if (any(ok_as) && sum(ok_as) >= n_min_dm)
      dm_hln(ea_m[ok_as], ea_b[ok_as], h) else
      list(stat = NA_real_, pvalue = NA_real_)
    dm_bc  <- if (any(ok_bc) && sum(ok_bc) >= n_min_dm)
      dm_hln(eb_m[ok_bc], eb_b[ok_bc], h) else
      list(stat = NA_real_, pvalue = NA_real_)
    dm_mae <- if (sum(ok_lvl) >= n_min_dm)
      dm_hln(e_m[ok_lvl], e_b[ok_lvl], h, loss = "abserr") else
      list(stat = NA_real_, pvalue = NA_real_)
    gw_sq  <- gw_test(e_m, e_b, h, "sqerr",  keep = ok_lvl)
    gw_ab  <- gw_test(e_m, e_b, h, "abserr", keep = ok_lvl)
    gw_sqa <- gw_test(ea_m, ea_b, h, "sqerr",  keep = ok_as)
    gw_aba <- gw_test(ea_m, ea_b, h, "abserr", keep = ok_as)
    rows[[pn]] <- data.frame(
      family = family, model = label, variables = variables,
      h = h, period = pn,
      n = sum(ok_lvl),
      rmse        = rmse_(e_m[ok_lvl]),
      rmse_bench  = rmse_(e_b[ok_lvl]),
      rmse_ratio  = rmse_(e_m[ok_lvl]) / rmse_(e_b[ok_lvl]),
      mae_ratio   = mae_(e_m[ok_lvl])  / mae_(e_b[ok_lvl]),
      # Explicit guard (referee point 7): undefined when no valid pair of
      # transformed-scale observations remains in the subsample.
      rmse_ratio_asinh = if (any(ok_as))
        rmse_(ea_m[ok_as]) / rmse_(ea_b[ok_as]) else NA_real_,
      mae_ratio_asinh = if (any(ok_as))
        mae_(ea_m[ok_as]) / mae_(ea_b[ok_as]) else NA_real_,
      rmse_bc_ratio = if (any(ok_bc))
        rmse_(eb_m[ok_bc]) / rmse_(eb_b[ok_bc]) else NA_real_,
      dm_stat = dm_lvl$stat, dm_pval = dm_lvl$pvalue,
      dm_mae_stat = dm_mae$stat, dm_mae_pval = dm_mae$pvalue,
      dm_stat_asinh = dm_as$stat, dm_pval_asinh = dm_as$pvalue,
      dm_bc_stat = dm_bc$stat, dm_bc_pval = dm_bc$pvalue,
      gw_stat = gw_sq$stat, gw_pval = gw_sq$pvalue, gw_dbar = gw_sq$dbar,
      gw_mae_stat = gw_ab$stat, gw_mae_pval = gw_ab$pvalue, gw_mae_dbar = gw_ab$dbar,
      gw_stat_asinh = gw_sqa$stat, gw_pval_asinh = gw_sqa$pvalue, gw_dbar_asinh = gw_sqa$dbar,
      gw_mae_stat_asinh = gw_aba$stat, gw_mae_pval_asinh = gw_aba$pvalue, gw_mae_dbar_asinh = gw_aba$dbar,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# Weighted combination of the MEMBERS' bias-corrected forecasts, using the
# weight matrix estimated from the ORDINARY forecast errors (renormalised
# over members whose bc forecast is available at each origin).  This defines
# the bias-corrected combination consistently: same weights, corrected
# members.  Returns NA rows where no bc member is available.
apply_weights_bc <- function(W, Fbc) {
  num <- rowSums(ifelse(is.finite(Fbc), W * Fbc, NA), na.rm = TRUE)
  den <- rowSums(ifelse(is.finite(Fbc), W, NA), na.rm = TRUE)
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_country <- function(country) {
  cc <- if (country == "Germany") "de" else "it"
  cat(sprintf("\n=== Combination + subperiod + GR fluctuation — %s ===\n", country))

  sum_rows <- list()
  gr_rows  <- list()
  w1_rows  <- list()   # stage-1 average combination weights (diagnostics)

  # Human-readable model descriptions written into the `variables` column
  # of the subperiod and GR tables.
  AR_DESC <- c(
    ar_dum           = "AR(1,2,3,7) + calendar dummies [BENCHMARK]",
    ar_dum_gas       = "AR(1,2,3,7) + calendar + gas lag1",
    ar_dum_brent     = "AR(1,2,3,7) + calendar + brent lag1",
    ar_dum_gas_brent = "AR(1,2,3,7) + calendar + gas lag1 + brent lag1"
  )
  COMB_DESC <- c(
    comb_eq         = "equal-weight comb. of all R-MIDAS specs",
    comb_trim       = "20%-trimmed-mean comb. of all R-MIDAS specs",
    comb_invmse_fam = "inverse-MSE comb. of all R-MIDAS specs (90d trail)",
    comb_invmse     = "two-stage inverse-MSE comb.: {family comb., ar_dum, ar_dum_gas}"
  )
  # Spec composition lookup from the spec dictionary (per family, loaded in
  # the family loop below): spec_id -> "hf_vars + lf_vars".

  for (h in HORIZONS) {
    ar_fp <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    if (!file.exists(ar_fp)) { cat("  missing", ar_fp, "— skipped\n"); next }
    ar <- read.csv(ar_fp, stringsAsFactors = FALSE)
    ar$origin_date <- as.Date(ar$origin_date)
    ar$target_date <- as.Date(ar$target_date)
    if (!BENCH_COL %in% names(ar)) stop("Benchmark column not found: ", BENCH_COL)
    ar <- ar[order(ar$origin_date), , drop = FALSE]        # consecutive origins (GW lag)
    vp <- load_vst_params(FORECAST_DIR, cc, ar$origin_date)
    ar$vst_a <- vp$a; ar$vst_b <- vp$b

    # AR variants vs ar_dum (bc columns passed when present in the CSV)
    bench_bc <- if ("ar_dum__bc" %in% names(ar)) ar[["ar_dum__bc"]] else NULL
    for (mcol in setdiff(names(ar), c("origin_date", "target_date",
                                      "y_actual", BENCH_COL, "vst_a", "vst_b",
                                      grep("__bc$", names(ar), value = TRUE)))) {
      mdesc <- if (mcol %in% names(AR_DESC)) AR_DESC[[mcol]] else mcol
      mbc   <- if (paste0(mcol, "__bc") %in% names(ar))
        ar[[paste0(mcol, "__bc")]] else NULL
      sum_rows[[length(sum_rows) + 1L]] <-
        eval_model(mcol, "AR", h, ar[[mcol]], ar$y_actual,
                   ar[[BENCH_COL]], ar$target_date,
                   yhat_bc = mbc, ybench_bc = bench_bc,
                   variables = mdesc, a = ar$vst_a, b = ar$vst_b)
      gr <- gr_fluctuation(ar[[BENCH_COL]] - ar$y_actual,
                           ar[[mcol]] - ar$y_actual, ar$target_date, h)
      if (!is.null(gr))
        gr_rows[[length(gr_rows) + 1L]] <-
          data.frame(family = "AR", model = mcol, variables = mdesc,
                     h = h, as.data.frame(gr),
                     stringsAsFactors = FALSE)
    }

    for (fam in FAMILIES) {
      fp <- file.path(FORECAST_DIR, sprintf("%s_%s_h%02d.csv", fam, cc, h))
      if (!file.exists(fp)) { cat("  missing", fp, "— skipped\n"); next }
      fc <- read.csv(fp, stringsAsFactors = FALSE)
      fc$origin_date <- as.Date(fc$origin_date)
      fc$target_date <- as.Date(fc$target_date)

      ar_pool <- intersect(COMB_INCLUDE_AR, names(ar))
      ar_pool <- ar_pool[!grepl("__bc$", ar_pool)]
      # Also carry the bias-corrected twins of the benchmark and of the AR
      # pool candidates into the merge (needed for the bc combinations and
      # the bc metric block).
      ar_bc_cols <- intersect(paste0(c(BENCH_COL, ar_pool), "__bc"), names(ar))
      mg <- merge(ar[, unique(c("origin_date", "target_date", "y_actual",
                                BENCH_COL, ar_pool, ar_bc_cols, "vst_a", "vst_b"))],
                  fc[, setdiff(names(fc), "y_actual")],
                  by = c("origin_date", "target_date"))
      mg <- mg[order(mg$origin_date), , drop = FALSE]
      specs <- grep("^spec_", names(mg), value = TRUE)
      specs <- specs[!grepl("__bc$", specs)]   # exclude bias-corrected twins
      if (length(specs) == 0L) next
      Fm <- as.matrix(mg[, specs])
      y  <- mg$y_actual
      yb <- mg[[BENCH_COL]]
      td <- mg$target_date

      # ---- combinations (real time) ----
      # comb_eq / comb_trim: averages of the R-MIDAS family only.
      # comb_invmse_fam: inverse-MSE combination of the family only.
      # comb_invmse: TWO-STAGE.  Stage 1 = comb_invmse_fam.  Stage 2 =
      # inverse-MSE over the SMALL pool {stage-1 combo, AR candidates}.
      # With only 3 candidates the weights can concentrate on the AR when
      # the macro block hurts; putting the AR directly into the 256-strong
      # pool dilutes its weight to ~2/258 and has no effect.
      comb_eq     <- rowMeans(Fm, na.rm = TRUE)
      comb_trim   <- apply(Fm, 1, mean, trim = TRIM, na.rm = TRUE)
      cf1 <- combine_invmse(Fm, y, mg$origin_date, h, return_weights = TRUE)
      comb_invmse_fam <- cf1$out
      # Diagnostics: average stage-1 weight per spec across origins — shows
      # which specs de facto drive the family combination.
      w1_rows[[length(w1_rows) + 1L]] <- data.frame(
        family = fam, h = h, spec = colnames(Fm),
        avg_weight = colMeans(cf1$W, na.rm = TRUE),
        stringsAsFactors = FALSE)
      # ---- bias-corrected combinations ---------------------------------
      # Same weights (estimated on ORDINARY resolved errors), applied to the
      # members' bias-corrected forecasts (renormalised over available bc
      # members).  Defines comb_*__bc consistently: corrected members,
      # identical weighting scheme, no look-ahead.
      spec_bc_cols <- paste0(specs, "__bc")
      have_spec_bc <- all(spec_bc_cols %in% names(mg))
      Fbc <- if (have_spec_bc) as.matrix(mg[, spec_bc_cols]) else NULL
      comb_eq_bc   <- if (have_spec_bc) rowMeans(Fbc, na.rm = TRUE) else
                        rep(NA_real_, length(y))
      comb_trim_bc <- if (have_spec_bc)
        apply(Fbc, 1, mean, trim = TRIM, na.rm = TRUE) else
        rep(NA_real_, length(y))
      comb_fam_bc  <- if (have_spec_bc) apply_weights_bc(cf1$W, Fbc) else
                        rep(NA_real_, length(y))

      if (length(ar_pool) > 0L) {
        F_pool <- cbind(comb_fam = comb_invmse_fam,
                        as.matrix(mg[, ar_pool, drop = FALSE]))
        cf2 <- combine_invmse(F_pool, y, mg$origin_date, h,
                              return_weights = TRUE)
        comb_invmse <- cf2$out
        # Stage-2 bc: same stage-2 weights on the bc members {family bc
        # combo, ar bc candidates} where available.
        pool_bc <- intersect(paste0(ar_pool, "__bc"), names(mg))
        comb_invmse_bc <- if (have_spec_bc && length(pool_bc) == length(ar_pool)) {
          F_pool_bc <- cbind(comb_fam_bc,
                             as.matrix(mg[, pool_bc, drop = FALSE]))
          apply_weights_bc(cf2$W, F_pool_bc)
        } else rep(NA_real_, length(y))
        # Diagnostics: full stage-2 weight path (small: 3 candidates) — shows
        # WHEN the combination migrates to the AR and back.
        w2 <- data.frame(origin_date = format(mg$origin_date, "%Y-%m-%d"),
                         cf2$W, check.names = FALSE,
                         stringsAsFactors = FALSE)
        w2p <- file.path(FORECAST_DIR,
                         sprintf("comb_w2_%s_%s_h%02d.csv", fam, cc, h))
        write.csv(w2, w2p, row.names = FALSE)
      } else {
        comb_invmse    <- comb_invmse_fam
        comb_invmse_bc <- comb_fam_bc
      }

      combo <- data.frame(origin_date = mg$origin_date,
                          target_date = td, y_actual = y,
                          comb_eq = comb_eq, comb_trim = comb_trim,
                          comb_invmse_fam = comb_invmse_fam,
                          comb_invmse = comb_invmse,
                          comb_eq__bc = comb_eq_bc,
                          comb_trim__bc = comb_trim_bc,
                          comb_invmse_fam__bc = comb_fam_bc,
                          comb_invmse__bc = comb_invmse_bc,
                          check.names = FALSE,
                          stringsAsFactors = FALSE)
      cfp <- file.path(FORECAST_DIR, sprintf("comb_%s_%s_h%02d.csv", fam, cc, h))
      write.csv(combo, cfp, row.names = FALSE)

      # ---- individual specs to evaluate --------------------------------
      # TOP_N_SPECS <= 0 (default): evaluate ALL specs, plain labels (no
      # selection, hence no "_expost" flag).  TOP_N_SPECS = N > 0: only the
      # N best full-sample specs, labelled "_expost" because the ranking
      # uses the whole evaluation sample (look-ahead, disclose in the text).
      full_rmse <- apply(Fm, 2, function(f) rmse_(f - y))
      ord <- order(full_rmse)
      all_mode <- TOP_N_SPECS <= 0L || TOP_N_SPECS >= length(specs)
      sel_specs <- if (all_mode) specs else
        specs[ord[seq_len(TOP_N_SPECS)]]
      best_spec <- specs[ord[1L]]

      # Spec composition lookup ("gas|brent + pmi_serv|Energy" style) from
      # the spec dictionary, so the `variables` column of every table says
      # WHAT each model contains rather than only its id.
      dict_fp <- file.path(FORECAST_DIR,
                           sprintf("%s_spec_dictionary_%s.csv", fam, cc))
      spec_desc <- setNames(rep("", length(specs)), specs)
      if (file.exists(dict_fp)) {
        dd <- read.csv(dict_fp, stringsAsFactors = FALSE)
        dd$hf_vars[!nzchar(dd$hf_vars) | is.na(dd$hf_vars)] <- "none"
        dd$lf_vars[!nzchar(dd$lf_vars) | is.na(dd$lf_vars)] <- "none"
        spec_desc[dd$spec_id] <- paste0("HF: ", dd$hf_vars,
                                        " + LF: ", dd$lf_vars)
      }

      models    <- list()
      models_bc <- list()   # bias-corrected twin of each model (or NULL)
      desc      <- list()   # `variables` text of each model
      for (s in sel_specs) {
        lbl <- if (all_mode) s else paste0(s, "_expost")
        models[[lbl]] <- Fm[, s]
        bc_col <- paste0(s, "__bc")
        models_bc[[lbl]] <- if (bc_col %in% names(mg)) mg[[bc_col]] else NULL
        desc[[lbl]] <- spec_desc[[s]]
      }
      models[["comb_eq"]]          <- comb_eq
      models[["comb_trim"]]        <- comb_trim
      models[["comb_invmse_fam"]]  <- comb_invmse_fam
      models[["comb_invmse"]]      <- comb_invmse
      models_bc[["comb_eq"]]         <- comb_eq_bc
      models_bc[["comb_trim"]]       <- comb_trim_bc
      models_bc[["comb_invmse_fam"]] <- comb_fam_bc
      models_bc[["comb_invmse"]]     <- comb_invmse_bc
      for (nm in names(COMB_DESC)) desc[[nm]] <- COMB_DESC[[nm]]

      bench_bc2 <- if ("ar_dum__bc" %in% names(mg)) mg[["ar_dum__bc"]] else NULL
      for (lbl in names(models)) {
        sum_rows[[length(sum_rows) + 1L]] <-
          eval_model(lbl, fam, h, models[[lbl]], y, yb, td,
                     yhat_bc = models_bc[[lbl]], ybench_bc = bench_bc2,
                     variables = desc[[lbl]], a = mg$vst_a, b = mg$vst_b)
        gr <- gr_fluctuation(yb - y, models[[lbl]] - y, td, h)
        if (!is.null(gr))
          gr_rows[[length(gr_rows) + 1L]] <-
            data.frame(family = fam, model = lbl,
                       variables = desc[[lbl]], h = h, as.data.frame(gr),
                       stringsAsFactors = FALSE)
      }
      cat(sprintf("  %s h=%02d: best ex-post spec %s; combos written to %s\n",
                  fam, h, best_spec, basename(cfp)))
    }
  }

  sum_df <- do.call(rbind, sum_rows)
  sum_df$country <- cc
  sp <- file.path(RESULT_DIR, sprintf("summary_subperiod_%s.csv", cc))
  write.csv(sum_df, sp, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d)\n", sp, nrow(sum_df)))

  if (length(gr_rows) > 0L) {
    gr_df <- do.call(rbind, gr_rows)
    gr_df$country <- cc
    gr_df$reject_05 <- abs(pmax(gr_df$stat_max, -gr_df$stat_min)) > gr_df$cv05
    gp <- file.path(RESULT_DIR, sprintf("gr_fluctuation_%s.csv", cc))
    write.csv(gr_df, gp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", gp, nrow(gr_df)))
  }

  if (length(w1_rows) > 0L) {
    w1_df <- do.call(rbind, w1_rows)
    w1_df$country <- cc
    wp <- file.path(RESULT_DIR, sprintf("comb_weights_stage1_%s.csv", cc))
    write.csv(w1_df, wp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", wp, nrow(w1_df)))
  }
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_country(cc_)
cat("\nCombination + subperiod + GR fluctuation done.\n")
