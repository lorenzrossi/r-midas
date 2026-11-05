#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — PREDICTIVE-ABILITY TESTS: EVERY MODEL vs THE AR BENCHMARK
#               (Table A) AND EVERY MACRO SPECIFICATION vs ITS DAILY-ONLY TWIN
#               (Table B), BY SUBPERIOD, WITH DM-HLN AND GIACOMINI-WHITE
#
# Families: AR (ar_dum_gas, ar_dum_brent, ar_dum_gas_brent), ARX, R_MIDAS,
# R_MIDAS_L7 (all 256 specifications each) and the real-time combinations
# of every family (COMB_<family>: comb_eq, comb_trim, comb_invmse_fam,
# comb_invmse).
#
# Table A — reference ar_dum (AR(1,2,3,7) + calendar).  Per
#   (country, h, period, model): n, rmse, mae, rmse_ratio, mae_ratio,
#   rmse_bc_ratio, rmse_ratio_asinh, mae_ratio_asinh, then for each of the
#   four loss definitions {squared, absolute} x {EUR/MWh, transformed scale}
#   the DM-HLN statistic / p-value and the GW statistic / p-value / dbar:
#     dm_stat, dm_pval, gw_stat, gw_pval, gw_dbar            squared, EUR
#     *_mae                                                   absolute, EUR
#     *_asinh                                                 squared, transformed
#     *_mae_asinh                                             absolute, transformed
#   plus the density metrics of the same period from density_eval.R
#   (logS_mean, CRPS_mean, share_in_90) when available.
#   -> results/tests_vs_benchmark_<cc>.csv
#      results/tests_l7_vs_benchmark_<cc>.csv  (R_MIDAS_L7 + its combinations)
#
# Table B — twin = specification of the SAME family with the same fuel set
#   and no monthly indicator (nested; isolates the monthly information).
#   Combinations are compared with the family's daily-only gas+brent spec.
#   Same test columns, ratios "_vs_twin".
#   -> results/tests_vs_twin_<cc>.csv, results/tests_l7_vs_twin_<cc>.csv
#      results/tests_vs_twin_summary.csv: per (country, h, period, family)
#      the median RMSE ratio vs twin and the share of the 252 macro
#      specifications beating their twin at 5% (DM, GW), overall and for the
#      specifications containing each indicator.
#
# Sign conventions: DM negative = model better; GW is a chi^2(2) statistic
# without sign — read the direction from gw_dbar (negative = model better).
# Transformed-scale errors use the normalised asinh with each origin's
# (a, b) from forecasts/vst_params_<cc>.csv.
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

BENCH <- "ar_dum"
FAM_FILES <- c(ARX = "arx", R_MIDAS = "r_midas", R_MIDAS_L7 = "r_midas_l7")
PERIODS <- list(
  full   = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
  pre    = c(as.Date("1900-01-01"), as.Date("2020-12-31")),
  crisis = c(as.Date("2021-01-01"), as.Date("2022-12-31")),
  post   = c(as.Date("2023-01-01"), as.Date("2100-01-01"))
)
rmse_ <- function(e) { e <- e[is.finite(e)]; if (!length(e)) NA_real_ else sqrt(mean(e^2)) }
mae_  <- function(e) { e <- e[is.finite(e)]; if (!length(e)) NA_real_ else mean(abs(e)) }

twin_map <- function(dict) {
  lf <- ifelse(is.na(dict$lf_vars), "", dict$lf_vars)
  hf <- ifelse(is.na(dict$hf_vars), "", dict$hf_vars)
  twins <- data.frame(hf_vars = hf[lf == ""], twin_id = dict$spec_id[lf == ""],
                      stringsAsFactors = FALSE)
  macro <- data.frame(spec_id = dict$spec_id[lf != ""], hf_vars = hf[lf != ""],
                      lf_vars = lf[lf != ""], stringsAsFactors = FALSE)
  merge(macro, twins, by = "hf_vars", all.x = TRUE)
}

# All tests of forecast f against reference r (both EUR/MWh), per period.
#   e_* : EUR errors, ea_* : transformed-scale errors (full sequences, origin order)
test_rows <- function(f, r, y, a, b, td, h, f_bc = NULL, r_bc = NULL, base, suffix = "") {
  e_m <- f - y; e_r <- r - y
  ya <- vst_z(y, a, b)                                   # model's own transformed scale
  ea_m <- vst_z(f, a, b) - ya; ea_r <- vst_z(r, a, b) - ya
  has_bc <- !is.null(f_bc) && !is.null(r_bc)
  out <- list()
  for (pn in names(PERIODS)) {
    sel <- td >= PERIODS[[pn]][1] & td <= PERIODS[[pn]][2]
    ok  <- sel & is.finite(e_m) & is.finite(e_r)
    oka <- ok & is.finite(ea_m) & is.finite(ea_r)
    n_ok <- sum(ok)
    if (n_ok < 5L) next
    n_min <- max(30L, 10L * h)
    row <- cbind(base, data.frame(period = pn, n = n_ok, stringsAsFactors = FALSE))
    row$rmse <- rmse_(e_m[ok]); row$mae <- mae_(e_m[ok])
    row[[paste0("rmse_ratio", suffix)]] <- rmse_(e_m[ok]) / rmse_(e_r[ok])
    row[[paste0("mae_ratio",  suffix)]] <- mae_(e_m[ok])  / mae_(e_r[ok])
    row[[paste0("rmse_bc_ratio", suffix)]] <- if (has_bc) {
      okb <- ok & is.finite(f_bc) & is.finite(r_bc)
      rmse_((f_bc - y)[okb]) / rmse_((r_bc - y)[okb])
    } else NA_real_
    row[[paste0("rmse_ratio_asinh", suffix)]] <- rmse_(ea_m[oka]) / rmse_(ea_r[oka])
    row[[paste0("mae_ratio_asinh",  suffix)]] <- mae_(ea_m[oka])  / mae_(ea_r[oka])
    tests <- list(
      sq        = list(sfx = "",           e1 = e_m,  e2 = e_r,  loss = "sqerr",  k = ok),
      mae       = list(sfx = "_mae",       e1 = e_m,  e2 = e_r,  loss = "abserr", k = ok),
      sq_asinh  = list(sfx = "_asinh",     e1 = ea_m, e2 = ea_r, loss = "sqerr",  k = oka),
      mae_asinh = list(sfx = "_mae_asinh", e1 = ea_m, e2 = ea_r, loss = "abserr", k = oka))
    for (tt in tests) {
      sfx <- tt$sfx
      if (sum(tt$k) >= n_min) {
        dm <- dm_hln(tt$e1[tt$k], tt$e2[tt$k], h, loss = tt$loss)
        gw <- gw_test(tt$e1, tt$e2, h, loss = tt$loss, keep = tt$k)
      } else dm <- gw <- list(stat = NA_real_, pvalue = NA_real_, dbar = NA_real_)
      row[[paste0("dm_stat", sfx)]] <- dm$stat
      row[[paste0("dm_pval", sfx)]] <- dm$pvalue
      row[[paste0("gw_stat", sfx)]] <- gw$stat
      row[[paste0("gw_pval", sfx)]] <- gw$pvalue
      row[[paste0("gw_dbar", sfx)]] <- gw$dbar
    }
    out[[pn]] <- row
  }
  do.call(rbind, out)
}

AR_DESC <- c(ar_dum_gas = "gas", ar_dum_brent = "brent", ar_dum_gas_brent = "gas|brent")

run_country <- function(cc) {
  cat(sprintf("\n=== DM / GW tests — %s ===\n", toupper(cc)))
  n_cores <- get_n_cores()
  dens_fp <- file.path(RESULT_DIR, sprintf("density_scores_%s.csv", cc))
  dens <- if (file.exists(dens_fp)) read.csv(dens_fp, stringsAsFactors = FALSE) else NULL
  A_rows <- list(); B_rows <- list()
  for (h in HORIZONS) {
    ar_fp <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    if (!file.exists(ar_fp)) { cat("  missing", ar_fp, "\n"); next }
    ar <- read.csv(ar_fp, stringsAsFactors = FALSE)
    ar$origin_date <- as.Date(ar$origin_date); ar$target_date <- as.Date(ar$target_date)
    ar <- ar[order(ar$origin_date), , drop = FALSE]
    key <- paste(ar$origin_date, ar$target_date)
    y <- ar$y_actual; td <- ar$target_date
    vp <- load_vst_params(FORECAST_DIR, cc, ar$origin_date); a <- vp$a; b <- vp$b
    bench <- ar[[BENCH]]; bench_bc <- ar[[paste0(BENCH, "__bc")]]
    t0 <- Sys.time()

    # ---- AR variants vs benchmark --------------------------------------------
    for (m in names(AR_DESC)) if (m %in% names(ar)) {
      A_rows[[length(A_rows) + 1L]] <- test_rows(
        ar[[m]], bench, y, a, b, td, h, ar[[paste0(m, "__bc")]], bench_bc,
        data.frame(country = cc, h = h, family = "AR", model = m,
                   hf_vars = AR_DESC[[m]], lf_vars = "", stringsAsFactors = FALSE))
    }

    # ---- spec families ----------------------------------------------------------
    for (famlab in names(FAM_FILES)) {
      fam <- FAM_FILES[[famlab]]
      fp <- file.path(FORECAST_DIR, sprintf("%s_%s_h%02d.csv", fam, cc, h))
      dp <- file.path(FORECAST_DIR, sprintf("%s_spec_dictionary_%s.csv", fam, cc))
      if (!file.exists(fp) || !file.exists(dp)) { cat("  missing", basename(fp), "— skipped\n"); next }
      fc <- read.csv(fp, stringsAsFactors = FALSE)
      fc <- fc[match(key, paste(fc$origin_date, fc$target_date)), , drop = FALSE]
      dict <- read.csv(dp, stringsAsFactors = FALSE)
      dict$hf_vars[is.na(dict$hf_vars)] <- ""; dict$lf_vars[is.na(dict$lf_vars)] <- ""
      specs <- dict$spec_id[dict$spec_id %in% names(fc)]
      tm <- twin_map(dict); tm <- tm[tm$spec_id %in% specs & tm$twin_id %in% specs, ]
      twin_of <- setNames(tm$twin_id, tm$spec_id)
      Fm  <- as.matrix(fc[, specs, drop = FALSE])
      Fbc <- as.matrix(fc[, paste0(specs, "__bc"), drop = FALSE])
      res <- parallel_lapply(seq_along(specs), function(j) {
        s <- specs[j]; i <- match(s, dict$spec_id)
        base <- data.frame(country = cc, h = h, family = famlab, model = s,
                           hf_vars = dict$hf_vars[i], lf_vars = dict$lf_vars[i],
                           stringsAsFactors = FALSE)
        A <- test_rows(Fm[, j], bench, y, a, b, td, h, Fbc[, j], bench_bc, base)
        B <- NULL
        if (s %in% names(twin_of)) {
          jt <- match(twin_of[[s]], specs)
          B <- test_rows(Fm[, j], Fm[, jt], y, a, b, td, h, Fbc[, j], Fbc[, jt],
                         cbind(base, data.frame(twin = twin_of[[s]], stringsAsFactors = FALSE)),
                         suffix = "_vs_twin")
        }
        list(A = A, B = B)
      }, n_cores = n_cores)
      A_rows[[length(A_rows) + 1L]] <- do.call(rbind, lapply(res, `[[`, "A"))
      Bl <- Filter(Negate(is.null), lapply(res, `[[`, "B"))
      if (length(Bl)) B_rows[[length(B_rows) + 1L]] <- do.call(rbind, Bl)

      # ---- the family's combinations: vs benchmark and vs the daily-only gas+brent spec
      cfp <- file.path(FORECAST_DIR, sprintf("comb_%s_%s_h%02d.csv", fam, cc, h))
      if (file.exists(cfp)) {
        cb <- read.csv(cfp, stringsAsFactors = FALSE)
        cb <- cb[match(key, paste(cb$origin_date, cb$target_date)), , drop = FALSE]
        twin_gb <- dict$spec_id[dict$hf_vars == "gas|brent" & dict$lf_vars == ""]
        jt <- if (length(twin_gb)) match(twin_gb[1], specs) else NA
        for (m in intersect(c("comb_eq", "comb_trim", "comb_invmse_fam", "comb_invmse"), names(cb))) {
          base <- data.frame(country = cc, h = h, family = paste0("COMB_", famlab), model = m,
                             hf_vars = "comb", lf_vars = "comb", stringsAsFactors = FALSE)
          mbc <- if (paste0(m, "__bc") %in% names(cb)) cb[[paste0(m, "__bc")]] else NULL
          A_rows[[length(A_rows) + 1L]] <- test_rows(cb[[m]], bench, y, a, b, td, h, mbc, bench_bc, base)
          if (!is.na(jt))
            B_rows[[length(B_rows) + 1L]] <- test_rows(
              cb[[m]], Fm[, jt], y, a, b, td, h, mbc, Fbc[, jt],
              cbind(base, data.frame(twin = twin_gb[1], stringsAsFactors = FALSE)),
              suffix = "_vs_twin")
        }
      }
    }
    cat(sprintf("  h=%02d done in %.0fs\n", h, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  A <- do.call(rbind, A_rows); B <- do.call(rbind, B_rows)
  # density metrics of the same (h, period, family, model)
  if (!is.null(dens)) {
    dk <- dens[, c("h", "period", "family", "model", "logS_mean", "CRPS_mean", "share_in_90")]
    A <- merge(A, dk, by = c("h", "period", "family", "model"), all.x = TRUE, sort = FALSE)
  }
  A <- A[order(A$h, A$family, A$model, match(A$period, names(PERIODS))), ]
  write.csv(A, file.path(RESULT_DIR, sprintf("tests_vs_benchmark_%s.csv", cc)), row.names = FALSE)
  write.csv(A[A$family %in% c("R_MIDAS_L7", "COMB_R_MIDAS_L7"), ],
            file.path(RESULT_DIR, sprintf("tests_l7_vs_benchmark_%s.csv", cc)), row.names = FALSE)
  if (!is.null(B)) {
    B <- B[order(B$h, B$family, B$model, match(B$period, names(PERIODS))), ]
    write.csv(B, file.path(RESULT_DIR, sprintf("tests_vs_twin_%s.csv", cc)), row.names = FALSE)
    write.csv(B[B$family %in% c("R_MIDAS_L7", "COMB_R_MIDAS_L7"), ],
              file.path(RESULT_DIR, sprintf("tests_l7_vs_twin_%s.csv", cc)), row.names = FALSE)
  }
  cat(sprintf("  wrote tests_vs_benchmark_%s.csv (%d rows), tests_vs_twin_%s.csv (%d rows)\n",
              cc, nrow(A), cc, if (is.null(B)) 0L else nrow(B)))
  list(A = A, B = B)
}

res <- lapply(eval_countries(), run_country)
B <- do.call(rbind, lapply(res, `[[`, "B"))
if (!is.null(B) && nrow(B)) {
  B <- B[grepl("^spec_", B$model), ]                  # macro specifications only
  indicators <- sort(unique(unlist(strsplit(B$lf_vars, "|", fixed = TRUE))))
  summ <- function(g) data.frame(
    n_specs = nrow(g),
    median_rmse_ratio_vs_twin = median(g$rmse_ratio_vs_twin, na.rm = TRUE),
    share_beat_twin           = mean(g$rmse_ratio_vs_twin < 1, na.rm = TRUE),
    share_dm_sig_better = mean(g$dm_pval < 0.05 & g$dm_stat < 0, na.rm = TRUE),
    share_gw_sig_better = mean(g$gw_pval < 0.05 & g$gw_dbar < 0, na.rm = TRUE),
    share_dm_sig_worse  = mean(g$dm_pval < 0.05 & g$dm_stat > 0, na.rm = TRUE),
    share_gw_sig_worse  = mean(g$gw_pval < 0.05 & g$gw_dbar > 0, na.rm = TRUE),
    stringsAsFactors = FALSE)
  keys <- unique(B[, c("country", "h", "period", "family")])
  out <- list()
  for (i in seq_len(nrow(keys))) {
    g <- B[B$country == keys$country[i] & B$h == keys$h[i] &
           B$period == keys$period[i] & B$family == keys$family[i], ]
    out[[length(out) + 1L]] <- cbind(keys[i, ], indicator = "all", summ(g))
    for (v in indicators) {
      gi <- g[vapply(strsplit(g$lf_vars, "|", fixed = TRUE), function(s) v %in% s, logical(1)), ]
      if (nrow(gi)) out[[length(out) + 1L]] <- cbind(keys[i, ], indicator = v, summ(gi))
    }
  }
  S <- do.call(rbind, out)
  sfile <- sprintf("tests_vs_twin_summary_%s.csv", paste(eval_countries(), collapse = "_"))
  write.csv(S, file.path(RESULT_DIR, sfile), row.names = FALSE)
  cat(sprintf("\nwrote %s (%d rows)\n", sfile, nrow(S)))
  cat("\n  macro specs vs twin, full sample, indicator = all (median RMSE ratio | share GW-sig better | DM-sig better):\n")
  sf <- S[S$period == "full" & S$indicator == "all", ]
  for (i in seq_len(nrow(sf)))
    cat(sprintf("    %s h=%02d %-10s  %.3f | %.2f | %.2f\n", sf$country[i], sf$h[i], sf$family[i],
                sf$median_rmse_ratio_vs_twin[i], sf$share_gw_sig_better[i], sf$share_dm_sig_better[i]))
}
cat("\nDM / GW tests done.\n")
