#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — DENSITY-FORECAST EVALUATION
#
# Every model in the pipeline is a Gaussian linear model on the transformed
# scale z = g(P) selected by VST_MODE (normalized asinh with the origin's
# window (a, b); plain asinh; log; level — see common_utils.R), so its h-step
# predictive distribution is  z ~ N(mu_h, sigma2_h)  with mu_h the iterated
# point forecast and sigma2_h the horizon-specific error variance.  On the
# PRICE scale the predictive CDF is  P(u) = Phi((g(u) - mu_h) / sigma_h)  in
# closed form and the density is  phi(.) / sigma_h * g'(u).  No
# retransformation bias arises: the density is mapped exactly (the point
# forecast g^-1(mu) is its median).  VST_MODE must be the one the forecasts
# were produced with.
#
# Scores (all on the price scale, euros/MWh; "+" = higher is better):
#   logS        log predictive likelihood  log p(y)                        (+)
#   CRPS        EXACT continuous ranked probability score                   (-)
#                 int (P(u) - 1{u >= y})^2 du = 2 int_0^1 (1{y<Q(t)} - t)(Q(t) - y) dt
#               closed form for level / log, K-point quantile rule otherwise
#               (crps_exact in common_utils.R)
#   CRPS_trunc  the same integral truncated to the common range [umin, umax]
#               (a threshold-weighted CRPS with indicator weight, Gneiting-
#               Ranjan 2011; kept because ACPS uses the same range)         (-)
#   ACPS_c      asymmetric continuous probability score of Iacopini,
#               Ravazzolo & Rossini (2023, JBES; arXiv 2006.11265), eq. (1): (+)
#                 ACPS(P,y;c) = int_{-inf}^{y} (c^2 - P(u)^2) w(u) du
#                             + int_{y}^{+inf} ((1-c)^2 - (1-P(u))^2) w(u) du,
#                 w(u) = 1/(1-c)^2 if P(u) > c,  1/c^2 if P(u) <= c.
#               c = 0.5 is symmetric; c < 0.5 penalises more an outcome in
#               the LEFT tail of the forecast, c > 0.5 in the RIGHT tail.
#               Truncated to [umin, umax] as in eq. (2): LEVELS depend on the
#               range, DIFFERENCES between models do not.
#   PIT         probability integral transform P(y): mean, sd (uniform: 0.5,
#               0.2887), Kolmogorov-Smirnov p-value, 90% interval coverage.
#
# Tests, per period (full / pre / crisis / post, >= max(30, 10h) obs):
#   vs the AR benchmark ar_dum — DM-HLN (Amisano-Giacomini) and GW (2006)
#     on the score differential, for every score;
#   vs the DAILY-ONLY TWIN (same family, same fuel set, no monthly
#     indicator) for every macro specification of ARX / R_MIDAS / R_MIDAS_L7
#     — on logS and CRPS.
#   NEGATIVE dm statistic = model better (loss = -score); GW has no sign:
#   gw_dbar < 0 = model better.
#
# INPUT : forecasts/<fam>_<cc>_h<h>.csv (+ _sd.csv), vst_params_<cc>.csv,
#         <fam>_spec_dictionary_<cc>.csv
# OUTPUT: results/density_scores_<cc>.csv          (h, family, model, period)
#         results/density_tests_vs_twin_<cc>.csv   (macro spec vs twin)
#         results/crps_summary_<cc>.csv            exact CRPS: ratio / skill vs
#                                                  ar_dum, DM + GW, twin columns
#         results/crps_twin_summary_<cc>.csv       median CRPS ratio vs twin and
#                                                  share of macro specs with a
#                                                  significantly lower CRPS,
#                                                  overall and by indicator
#         results/density_top5_per_country_h_<cc...>.csv  best 5 by logS
# Env   : ACPS_C (default "0.05,0.25,0.5,0.75,0.95"), CRPS_K (1000),
#         DENSITY_FAMILIES (default "ar,arx,r_midas,r_midas_l7"), COUNTRY_NAME.
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

# ---- unit tests of the CRPS code: stop if wrong --------------------------------
local({
  old <- VST_MODE
  VST_MODE <<- "level"
  stopifnot(abs(crps_quantile(71, 50, 12, 0, 1) - crps_norm(71, 50, 12)) < 1e-3)
  VST_MODE <<- "log"
  stopifnot(abs(crps_quantile(95, log(60), .35, 0, 1) - crps_lnorm(95, log(60), .35)) < 1e-3)
  VST_MODE <<- old
  cat(sprintf("CRPS unit tests passed (K = %d); VST_MODE = %s\n", CRPS_K, VST_MODE))
})

ACPS_C <- as.numeric(strsplit(Sys.getenv("ACPS_C", "0.05,0.25,0.5,0.75,0.95"), ",")[[1]])
ACPS_C <- ACPS_C[is.finite(ACPS_C) & ACPS_C > 0 & ACPS_C < 1]
FAMS   <- trimws(strsplit(Sys.getenv("DENSITY_FAMILIES", "ar,arx,r_midas,r_midas_l7"), ",")[[1]])
FAM_LABEL <- c(ar = "AR", arx = "ARX", r_midas = "R_MIDAS", r_midas_l7 = "R_MIDAS_L7")
BENCH  <- "ar_dum"
SCORE_COLS <- c("logS", "CRPS", "CRPS_trunc", paste0("ACPS_", ACPS_C))
LOSS_COLS  <- c("CRPS", "CRPS_trunc")                 # negatively oriented

PERIODS <- list(
  full   = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
  pre    = c(as.Date("1900-01-01"), as.Date("2020-12-31")),
  crisis = c(as.Date("2021-01-01"), as.Date("2022-12-31")),
  post   = c(as.Date("2023-01-01"), as.Date("2100-01-01"))
)

# Fixed z-grid; nodes at Phi^-1(c) so the ACPS weight switch sits on a node.
ZGRID <- sort(unique(c(seq(-8, 8, by = 0.05), qnorm(ACPS_C))))
PGRID <- pnorm(ZGRID)
G     <- length(ZGRID)

g_crps <- list(L = function(P) P^2, R = function(P) (1 - P)^2)
g_acps <- function(c) {
  w <- function(P) ifelse(P > c, 1 / (1 - c)^2, 1 / c^2)
  list(L = function(P) (c^2 - P^2) * w(P),
       R = function(P) ((1 - c)^2 - (1 - P)^2) * w(P))
}
GFUNS <- c(list(CRPS_trunc = g_crps), setNames(lapply(ACPS_C, g_acps), paste0("ACPS_", ACPS_C)))

# Integral over [umin, umax] of the side-dependent integrand for a block of
# forecasts (mu, sigma on the z scale; a, b the origin's VST parameters; y
# realised price), trapezoid on the z-grid with the segment containing y
# split at y and the constant tails added analytically.
integrate_scores <- function(mu, sigma, a, b, y, umin, umax) {
  n_ok <- length(y)
  U <- vst_inv0(outer(mu, rep(1, G)) + outer(sigma, ZGRID), a, b)
  U <- pmin(pmax(U, umin), umax)
  dU <- U[, -1, drop = FALSE] - U[, -G, drop = FALSE]
  z_y <- (vst_z(y, a, b) - mu) / sigma
  P_y <- pnorm(z_y)
  k_y <- findInterval(z_y, ZGRID)
  seg <- matrix(seq_len(G - 1L), n_ok, G - 1L, byrow = TRUE)
  left_seg  <- seg <  k_y
  right_seg <- seg >  k_y
  mid <- k_y >= 1L & k_y <= G - 1L
  im  <- which(mid)
  Uk  <- rep(NA_real_, n_ok); Uk1 <- Uk
  Uk[im]  <- U[cbind(im, k_y[im])]
  Uk1[im] <- U[cbind(im, k_y[im] + 1L)]
  yc <- pmin(pmax(y, umin), umax)
  U1 <- U[, 1]; UG <- U[, G]
  lo <- k_y == 0L; hi <- k_y == G
  out <- matrix(NA_real_, n_ok, length(GFUNS), dimnames = list(NULL, names(GFUNS)))
  for (s in names(GFUNS)) {
    gL <- GFUNS[[s]]$L(PGRID); gR <- GFUNS[[s]]$R(PGRID)
    trapL <- 0.5 * (gL[-G] + gL[-1]); trapR <- 0.5 * (gR[-G] + gR[-1])
    val <- rowSums(dU * left_seg  * rep(trapL, each = n_ok)) +
           rowSums(dU * right_seg * rep(trapR, each = n_ok))
    gLy <- GFUNS[[s]]$L(P_y); gRy <- GFUNS[[s]]$R(P_y)
    val[im] <- val[im] +
      0.5 * (gL[k_y[im]] + gLy[im]) * (yc[im] - Uk[im]) +
      0.5 * (gRy[im] + gR[k_y[im] + 1L]) * (Uk1[im] - yc[im])
    gL0 <- GFUNS[[s]]$L(0); gR0 <- GFUNS[[s]]$R(0)
    gL1 <- GFUNS[[s]]$L(1); gR1 <- GFUNS[[s]]$R(1)
    val <- val +
      ifelse(lo, gL0 * (yc - umin) + gR0 * (U1 - yc), gL0 * (U1 - umin)) +
      ifelse(hi, gL1 * (yc - UG) + gR1 * (umax - yc), gR1 * (umax - UG))
    out[, s] <- val
  }
  list(scores = out, pit = P_y, z_y = z_y)
}

# Per-observation scores of one model: n x (SCORE_COLS, PIT).
# point = MEDIAN forecast (the plain column, not __bc).
score_block <- function(point, sigma, a, b, y, umin, umax) {
  n <- length(y)
  S <- matrix(NA_real_, n, length(SCORE_COLS) + 1L,
              dimnames = list(NULL, c(SCORE_COLS, "PIT")))
  mu <- vst_z(point, a, b)
  ok <- is.finite(mu) & is.finite(sigma) & sigma > 0 & is.finite(y) & is.finite(a) & is.finite(b)
  if (!any(ok)) return(S)
  r <- integrate_scores(mu[ok], sigma[ok], a[ok], b[ok], y[ok], umin, umax)
  S[ok, "logS"] <- dnorm(r$z_y, log = TRUE) - log(sigma[ok]) + vst_logjac(y[ok], a[ok], b[ok])
  S[ok, colnames(r$scores)] <- r$scores
  S[ok, "CRPS"] <- crps_exact(y[ok], mu[ok], sigma[ok], a[ok], b[ok])
  S[ok, "PIT"]  <- r$pit
  S
}

# Period statistics of one model's scores against a reference's scores.
period_rows <- function(S, S_ref, td, h, base, with_pit = TRUE, scores = SCORE_COLS) {
  out <- list()
  for (pn in names(PERIODS)) {
    sel <- td >= PERIODS[[pn]][1] & td <= PERIODS[[pn]][2]
    ok  <- sel & is.finite(S[, "logS"]) & is.finite(S_ref[, "logS"])
    n_ok <- sum(ok)
    if (n_ok < 5L) next
    n_min <- max(30L, 10L * h)
    row <- cbind(base, data.frame(period = pn, n = n_ok, stringsAsFactors = FALSE))
    for (sc in scores) {
      sign_ <- if (sc %in% LOSS_COLS) 1 else -1              # loss = -score
      row[[paste0(sc, "_mean")]]  <- mean(S[ok, sc])
      row[[paste0(sc, "_ref")]]   <- mean(S_ref[ok, sc])
      row[[paste0(sc, "_diff")]]  <- mean(S[ok, sc]) - mean(S_ref[ok, sc])
      if (n_ok >= n_min) {
        dm <- dm_hln_loss(sign_ * S[ok, sc], sign_ * S_ref[ok, sc], h)
        gw <- gw_test_loss(sign_ * S[, sc], sign_ * S_ref[, sc], h, keep = ok)
      } else dm <- gw <- list(stat = NA_real_, pvalue = NA_real_, dbar = NA_real_)
      row[[paste0(sc, "_dm_stat")]] <- dm$stat
      row[[paste0(sc, "_dm_pval")]] <- dm$pvalue
      row[[paste0(sc, "_gw_stat")]] <- gw$stat
      row[[paste0(sc, "_gw_pval")]] <- gw$pvalue
      row[[paste0(sc, "_gw_dbar")]] <- gw$dbar
    }
    if (with_pit) {
      pit <- S[ok, "PIT"]
      row$pit_mean <- mean(pit); row$pit_sd <- sd(pit)
      row$pit_ks_pval <- tryCatch(
        suppressWarnings(ks.test(pit, "punif")$p.value), error = function(e) NA_real_)
      row$share_in_90 <- mean(pit > 0.05 & pit < 0.95)   # nominal 0.90
    }
    out[[pn]] <- row
  }
  do.call(rbind, out)
}

twin_map <- function(dict) {
  lf <- ifelse(is.na(dict$lf_vars), "", dict$lf_vars)
  hf <- ifelse(is.na(dict$hf_vars), "", dict$hf_vars)
  twins <- data.frame(hf_vars = hf[lf == ""], twin_id = dict$spec_id[lf == ""],
                      stringsAsFactors = FALSE)
  macro <- data.frame(spec_id = dict$spec_id[lf != ""], hf_vars = hf[lf != ""],
                      lf_vars = lf[lf != ""], stringsAsFactors = FALSE)
  merge(macro, twins, by = "hf_vars", all.x = TRUE)
}

AR_DESC <- c(ar_dum = "AR(1,2,3,7) + calendar dummies [BENCHMARK]",
             ar_dum_gas = "AR(1,2,3,7) + calendar + gas lag1",
             ar_dum_brent = "AR(1,2,3,7) + calendar + brent lag1",
             ar_dum_gas_brent = "AR(1,2,3,7) + calendar + gas lag1 + brent lag1")

# -----------------------------------------------------------------------------
run_country <- function(cc) {
  cat(sprintf("\n=== Density evaluation — %s (VST_MODE = %s) ===\n", toupper(cc), VST_MODE))
  y_all <- unlist(lapply(HORIZONS, function(h) {
    fp <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    if (file.exists(fp)) read.csv(fp)$y_actual else numeric(0)
  }))
  y_all <- y_all[is.finite(y_all)]
  if (!length(y_all)) { cat("  no AR forecasts found — skipped\n"); return(NULL) }
  umin <- floor(min(y_all)) - 100
  umax <- ceiling(max(y_all)) + 200
  cat(sprintf("  ACPS / CRPS_trunc range [%g, %g] euros/MWh; c = %s\n",
              umin, umax, paste(ACPS_C, collapse = ", ")))

  rows <- list(); twin_rows <- list(); n_viol <- 0L; n_tot <- 0L
  for (h in HORIZONS) {
    ar_fp <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d.csv", cc, h))
    ar_sd <- file.path(FORECAST_DIR, sprintf("ar_%s_h%02d_sd.csv", cc, h))
    if (!file.exists(ar_fp) || !file.exists(ar_sd)) { cat("  missing AR files for h =", h, "\n"); next }
    ar  <- read.csv(ar_fp, stringsAsFactors = FALSE)
    ars <- read.csv(ar_sd, stringsAsFactors = FALSE)
    ar$target_date <- as.Date(ar$target_date)
    key <- paste(ar$origin_date, ar$target_date)
    ars <- ars[match(key, paste(ars$origin_date, ars$target_date)), , drop = FALSE]
    y   <- ar$y_actual
    td  <- ar$target_date
    vp  <- load_vst_params(FORECAST_DIR, cc, ar$origin_date)
    S_bench <- score_block(ar[[BENCH]], ars[[BENCH]], vp$a, vp$b, y, umin, umax)

    for (fam in FAMS) {
      fp <- file.path(FORECAST_DIR, sprintf("%s_%s_h%02d.csv", fam, cc, h))
      sp <- file.path(FORECAST_DIR, sprintf("%s_%s_h%02d_sd.csv", fam, cc, h))
      if (!file.exists(fp) || !file.exists(sp)) { cat("  missing", basename(fp), "/ _sd — skipped\n"); next }
      fc  <- read.csv(fp, stringsAsFactors = FALSE)
      fcs <- read.csv(sp, stringsAsFactors = FALSE)
      fc  <- fc[match(key, paste(fc$origin_date, fc$target_date)), , drop = FALSE]
      fcs <- fcs[match(key, paste(fcs$origin_date, fcs$target_date)), , drop = FALSE]
      models <- setdiff(colnames(fc), c("origin_date", "target_date", "y_actual",
                                        grep("__bc$", colnames(fc), value = TRUE)))
      models <- models[models %in% colnames(fcs)]
      Fm <- as.matrix(fc[, models, drop = FALSE]); Sm <- as.matrix(fcs[, models, drop = FALSE])
      rm(fc, fcs)
      t_start <- Sys.time()
      S_all <- parallel_lapply(seq_along(models), function(j)
        score_block(Fm[, j], Sm[, j], vp$a, vp$b, y, umin, umax))
      names(S_all) <- models
      # the exact CRPS adds the tails, so it cannot fall below the truncated
      # one beyond the integration noise of the trapezoid rule (~2e-4
      # relative) and of the K-point quantile rule
      for (m in models) {
        d <- S_all[[m]][, "CRPS"] - S_all[[m]][, "CRPS_trunc"]
        tol <- 2e-3 * S_all[[m]][, "CRPS_trunc"] + 1e-3
        n_tot  <- n_tot + sum(is.finite(d)); n_viol <- n_viol + sum(d < -tol, na.rm = TRUE)
      }
      famlab <- FAM_LABEL[[fam]]
      desc <- setNames(rep(NA_character_, length(models)), models)
      dict <- NULL
      dp <- file.path(FORECAST_DIR, sprintf("%s_spec_dictionary_%s.csv", fam, cc))
      if (file.exists(dp)) {
        dict <- read.csv(dp, stringsAsFactors = FALSE)
        dict$hf_vars[is.na(dict$hf_vars)] <- ""; dict$lf_vars[is.na(dict$lf_vars)] <- ""
        i <- match(models, dict$spec_id)
        desc[] <- paste0("HF: ", ifelse(nzchar(dict$hf_vars[i]), dict$hf_vars[i], "none"),
                         " + LF: ", ifelse(nzchar(dict$lf_vars[i]), dict$lf_vars[i], "none"))
      }
      is_ar <- models %in% names(AR_DESC); desc[is_ar] <- AR_DESC[models[is_ar]]
      blocks <- lapply(models, function(m)
        period_rows(S_all[[m]], S_bench, td, h,
                    data.frame(country = cc, h = h, family = famlab, model = m,
                               variables = desc[[m]],
                               hf_vars = if (!is.null(dict)) dict$hf_vars[match(m, dict$spec_id)] else NA_character_,
                               lf_vars = if (!is.null(dict)) dict$lf_vars[match(m, dict$spec_id)] else NA_character_,
                               stringsAsFactors = FALSE)))
      rows[[length(rows) + 1L]] <- do.call(rbind, blocks)
      if (!is.null(dict)) {
        tm <- twin_map(dict)
        tm <- tm[tm$spec_id %in% models & tm$twin_id %in% models, , drop = FALSE]
        tb <- lapply(seq_len(nrow(tm)), function(i)
          period_rows(S_all[[tm$spec_id[i]]], S_all[[tm$twin_id[i]]], td, h,
                      data.frame(country = cc, h = h, family = famlab,
                                 model = tm$spec_id[i], twin = tm$twin_id[i],
                                 hf_vars = tm$hf_vars[i], lf_vars = tm$lf_vars[i],
                                 stringsAsFactors = FALSE),
                      with_pit = FALSE, scores = c("logS", "CRPS")))
        if (length(tb)) twin_rows[[length(twin_rows) + 1L]] <- do.call(rbind, tb)
      }
      cat(sprintf("  h=%02d %-10s %3d models scored in %.0fs\n", h, famlab,
                  length(models), as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
      rm(S_all, Fm, Sm); invisible(gc())
    }
  }
  cat(sprintf("  exact CRPS < truncated CRPS in %d of %d observations (expected 0)\n", n_viol, n_tot))
  df <- do.call(rbind, rows)
  names(df) <- sub("_ref$", "_bench", names(df))
  df$acps_umin <- umin; df$acps_umax <- umax; df$vst_mode <- VST_MODE
  out <- file.path(RESULT_DIR, sprintf("density_scores_%s.csv", cc))
  write.csv(df, out, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d)\n", out, nrow(df)))
  tw <- NULL
  if (length(twin_rows)) {
    tw <- do.call(rbind, twin_rows)
    names(tw) <- sub("_ref$", "_twin", names(tw))
    outt <- file.path(RESULT_DIR, sprintf("density_tests_vs_twin_%s.csv", cc))
    write.csv(tw, outt, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", outt, nrow(tw)))
  }

  # ---- CRPS tables for the paper ------------------------------------------------
  cs <- df[, c("country", "h", "period", "family", "model", "variables", "hf_vars", "lf_vars", "n",
               "CRPS_mean", "CRPS_bench", "CRPS_dm_stat", "CRPS_dm_pval",
               "CRPS_gw_stat", "CRPS_gw_pval", "CRPS_gw_dbar",
               "pit_sd", "share_in_90", "pit_ks_pval")]
  names(cs) <- c("country", "h", "period", "family", "model", "variables", "hf_vars", "lf_vars", "n",
                 "crps_mean", "crps_bench", "dm_stat", "dm_pval", "gw_stat", "gw_pval", "gw_dbar",
                 "pit_sd", "share_in_90", "pit_ks_pval")
  cs$crps_ratio <- cs$crps_mean / cs$crps_bench
  cs$crps_skill <- 1 - cs$crps_ratio
  if (!is.null(tw)) {
    tt <- tw[, c("h", "period", "family", "model", "twin", "CRPS_twin",
                 "CRPS_dm_stat", "CRPS_dm_pval", "CRPS_gw_stat", "CRPS_gw_pval", "CRPS_gw_dbar")]
    names(tt) <- c("h", "period", "family", "model", "twin", "crps_twin",
                   "dm_stat_vs_twin", "dm_pval_vs_twin", "gw_stat_vs_twin", "gw_pval_vs_twin", "gw_dbar_vs_twin")
    cs <- merge(cs, tt, by = c("h", "period", "family", "model"), all.x = TRUE, sort = FALSE)
    cs$crps_ratio_vs_twin <- cs$crps_mean / cs$crps_twin
  }
  cs <- cs[order(cs$h, match(cs$period, names(PERIODS)), cs$family, cs$crps_ratio), ]
  outc <- file.path(RESULT_DIR, sprintf("crps_summary_%s.csv", cc))
  write.csv(cs, outc, row.names = FALSE)
  cat(sprintf("  wrote %s  (rows = %d)\n", outc, nrow(cs)))

  if (!is.null(tw)) {
    # twin summary: median CRPS ratio and share of macro specs with a
    # significantly lower CRPS than their twin, overall and by indicator
    B <- cs[grepl("^spec_", cs$model) & !is.na(cs$twin), ]
    indicators <- sort(unique(unlist(strsplit(B$lf_vars, "|", fixed = TRUE))))
    summ <- function(g) data.frame(
      n_specs = nrow(g),
      median_crps_ratio_vs_twin = median(g$crps_ratio_vs_twin, na.rm = TRUE),
      share_beat_twin     = mean(g$crps_ratio_vs_twin < 1, na.rm = TRUE),
      share_dm_sig_better = mean(g$dm_pval_vs_twin < 0.05 & g$dm_stat_vs_twin < 0, na.rm = TRUE),
      share_gw_sig_better = mean(g$gw_pval_vs_twin < 0.05 & g$gw_dbar_vs_twin < 0, na.rm = TRUE),
      share_dm_sig_worse  = mean(g$dm_pval_vs_twin < 0.05 & g$dm_stat_vs_twin > 0, na.rm = TRUE),
      share_gw_sig_worse  = mean(g$gw_pval_vs_twin < 0.05 & g$gw_dbar_vs_twin > 0, na.rm = TRUE),
      stringsAsFactors = FALSE)
    keys <- unique(B[, c("country", "h", "period", "family")])
    out <- list()
    for (i in seq_len(nrow(keys))) {
      g <- B[B$h == keys$h[i] & B$period == keys$period[i] & B$family == keys$family[i], ]
      out[[length(out) + 1L]] <- cbind(keys[i, ], indicator = "all", summ(g))
      for (v in indicators) {
        gi <- g[vapply(strsplit(g$lf_vars, "|", fixed = TRUE), function(s) v %in% s, logical(1)), ]
        if (nrow(gi)) out[[length(out) + 1L]] <- cbind(keys[i, ], indicator = v, summ(gi))
      }
    }
    TS <- do.call(rbind, out)
    outs <- file.path(RESULT_DIR, sprintf("crps_twin_summary_%s.csv", cc))
    write.csv(TS, outs, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", outs, nrow(TS)))
  }

  cat("\n  full-sample digest (best logS per family; logS diff vs ar_dum, DM p, exact CRPS ratio, 90% coverage):\n")
  full <- df[df$period == "full", ]
  for (h in HORIZONS) for (fam in unique(full$family)) {
    g <- full[full$h == h & full$family == fam, ]
    g <- g[is.finite(g$logS_mean), ]
    if (!nrow(g)) next
    g <- g[order(-g$logS_mean), ][1, ]
    cat(sprintf("    h=%02d %-10s %-18s logS %+.4f (p=%.3f)  CRPS %.2f ratio %.3f  cov90 %.2f\n",
                h, fam, g$model, g$logS_diff, g$logS_dm_pval,
                g$CRPS_mean, g$CRPS_mean / g$CRPS_bench, g$share_in_90))
  }
  df
}

all_df <- do.call(rbind, Filter(Negate(is.null), lapply(eval_countries(), run_country)))
if (!is.null(all_df)) {
  full <- all_df[all_df$period == "full" & is.finite(all_df$logS_mean), ]
  best <- do.call(rbind, by(full, list(full$country, full$h), function(g)
    g[order(-g$logS_mean)[1:min(5, nrow(g))], , drop = FALSE]))
  bp <- file.path(RESULT_DIR, sprintf("density_top5_per_country_h_%s.csv",
                                      paste(eval_countries(), collapse = "_")))
  write.csv(best, bp, row.names = FALSE)
  cat(sprintf("\nwrote %s  (rows = %d)\n", bp, nrow(best)))
}
cat("\nDensity evaluation done.\n")
