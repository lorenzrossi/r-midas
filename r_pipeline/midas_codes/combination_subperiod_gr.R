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
#    a centered fraction mu = 0.3 of the sample, studentised with the same
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
# (comma-separated, default "r_midas,r_midas_extended"), TRAIL_DAYS (365),
# MIN_RESOLVED (60), TRIM (0.2), GR_MU (0.3).
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
TRAIL_DAYS   <- as.integer(Sys.getenv("TRAIL_DAYS",   "365"))
MIN_RESOLVED <- as.integer(Sys.getenv("MIN_RESOLVED", "60"))
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
                           trail = TRAIL_DAYS, min_res = MIN_RESOLVED) {
  n  <- nrow(Fm)
  E2 <- (Fm - y)^2
  ok <- is.finite(E2)
  E2z <- E2; E2z[!ok] <- 0
  CS <- apply(E2z, 2, cumsum)              # cumulative SSE per spec
  CN <- apply(ok,  2, cumsum)              # cumulative count per spec
  od <- as.numeric(odates)
  out <- rep(NA_real_, n)
  eqm <- rowMeans(Fm, na.rm = TRUE)        # fallback
  for (i in seq_len(n)) {
    b <- findInterval(od[i] - h, od)
    a <- findInterval(od[i] - h - trail, od)
    if (b <= a || b < 1L) { out[i] <- eqm[i]; next }
    cnt <- CN[b, ] - if (a >= 1L) CN[a, ] else 0
    sse <- CS[b, ] - if (a >= 1L) CS[a, ] else 0
    mse <- ifelse(cnt >= min_res, sse / pmax(cnt, 1L), NA_real_)
    w   <- 1 / mse
    fi  <- Fm[i, ]
    use <- is.finite(fi) & is.finite(w) & w > 0
    out[i] <- if (any(use)) sum(w[use] * fi[use]) / sum(w[use]) else eqm[i]
  }
  out
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
eval_model <- function(label, family, h, yhat, y, ybench, tdates) {
  e_m <- yhat - y
  e_b <- ybench - y
  ea_m <- asinh(yhat) - asinh(y)
  ea_b <- asinh(ybench) - asinh(y)
  rows <- list()
  for (pn in names(PERIODS)) {
    sel <- tdates >= PERIODS[[pn]][1] & tdates <= PERIODS[[pn]][2]
    ok_lvl <- sel & is.finite(e_m) & is.finite(e_b)
    ok_as  <- sel & is.finite(ea_m) & is.finite(ea_b)
    if (!any(ok_lvl)) next
    dm_lvl <- if (pn == "full") dm_hln(e_m[ok_lvl], e_b[ok_lvl], h) else
                list(stat = NA_real_, pvalue = NA_real_)
    dm_as  <- if (pn == "full") dm_hln(ea_m[ok_as], ea_b[ok_as], h) else
                list(stat = NA_real_, pvalue = NA_real_)
    rows[[pn]] <- data.frame(
      family = family, model = label, h = h, period = pn,
      n = sum(ok_lvl),
      rmse        = rmse_(e_m[ok_lvl]),
      rmse_bench  = rmse_(e_b[ok_lvl]),
      rmse_ratio  = rmse_(e_m[ok_lvl]) / rmse_(e_b[ok_lvl]),
      mae_ratio   = mae_(e_m[ok_lvl])  / mae_(e_b[ok_lvl]),
      rmse_ratio_asinh = rmse_(ea_m[ok_as]) / rmse_(ea_b[ok_as]),
      dm_stat = dm_lvl$stat, dm_pval = dm_lvl$pvalue,
      dm_stat_asinh = dm_as$stat, dm_pval_asinh = dm_as$pvalue,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_country <- function(country) {
  cc <- if (country == "Germany") "de" else "it"
  cat(sprintf("\n=== Combination + subperiod + GR fluctuation — %s ===\n", country))

  sum_rows <- list()
  gr_rows  <- list()

  for (h in HORIZONS) {
    ar_fp <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    if (!file.exists(ar_fp)) { cat("  missing", ar_fp, "— skipped\n"); next }
    ar <- read.csv(ar_fp, stringsAsFactors = FALSE)
    ar$origin_date <- as.Date(ar$origin_date)
    ar$target_date <- as.Date(ar$target_date)
    if (!BENCH_COL %in% names(ar)) stop("Benchmark column not found: ", BENCH_COL)

    # AR variants vs ar_dum
    for (mcol in setdiff(names(ar), c("origin_date", "target_date",
                                      "y_actual", BENCH_COL,
                                      grep("__bc$", names(ar), value = TRUE)))) {
      sum_rows[[length(sum_rows) + 1L]] <-
        eval_model(mcol, "AR", h, ar[[mcol]], ar$y_actual,
                   ar[[BENCH_COL]], ar$target_date)
      gr <- gr_fluctuation(ar[[BENCH_COL]] - ar$y_actual,
                           ar[[mcol]] - ar$y_actual, ar$target_date, h)
      if (!is.null(gr))
        gr_rows[[length(gr_rows) + 1L]] <-
          data.frame(family = "AR", model = mcol, h = h, as.data.frame(gr),
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
      mg <- merge(ar[, unique(c("origin_date", "target_date", "y_actual",
                                BENCH_COL, ar_pool))],
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
      # comb_invmse: pool = R-MIDAS specs + AR candidates, so the adaptive
      # weights can fall back to the pure AR when the macro block hurts.
      comb_eq     <- rowMeans(Fm, na.rm = TRUE)
      comb_trim   <- apply(Fm, 1, mean, trim = TRIM, na.rm = TRUE)
      F_pool <- if (length(ar_pool) > 0L)
        cbind(Fm, as.matrix(mg[, ar_pool, drop = FALSE])) else Fm
      comb_invmse <- combine_invmse(F_pool, y, mg$origin_date, h)

      combo <- data.frame(origin_date = mg$origin_date,
                          target_date = td, y_actual = y,
                          comb_eq = comb_eq, comb_trim = comb_trim,
                          comb_invmse = comb_invmse,
                          stringsAsFactors = FALSE)
      cfp <- file.path(FORECAST_DIR, sprintf("comb_%s_%s_h%02d.csv", fam, cc, h))
      write.csv(combo, cfp, row.names = FALSE)

      # ---- best full-sample individual spec (ex post; look-ahead flag) ----
      full_rmse <- apply(Fm, 2, function(f) rmse_(f - y))
      best_spec <- specs[which.min(full_rmse)]

      models <- list()
      models[[paste0(best_spec, "_expost")]] <- Fm[, best_spec]
      models[["comb_eq"]]     <- comb_eq
      models[["comb_trim"]]   <- comb_trim
      models[["comb_invmse"]] <- comb_invmse

      for (lbl in names(models)) {
        sum_rows[[length(sum_rows) + 1L]] <-
          eval_model(lbl, fam, h, models[[lbl]], y, yb, td)
        gr <- gr_fluctuation(yb - y, models[[lbl]] - y, td, h)
        if (!is.null(gr))
          gr_rows[[length(gr_rows) + 1L]] <-
            data.frame(family = fam, model = lbl, h = h, as.data.frame(gr),
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
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_country(cc_)
cat("\nCombination + subperiod + GR fluctuation done.\n")
