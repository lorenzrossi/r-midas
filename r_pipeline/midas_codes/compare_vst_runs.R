#!/usr/bin/env Rscript
# =============================================================================
# COMPARE THE ITALIAN RUNS ACROSS TRANSFORMATIONS (README_v4b, section 4)
#
# Reads the Italian forecasts of several run folders, each produced under a
# different VST_MODE (normalized asinh in estimation_v4; level / log / asinh
# in estimation_v4_it_<mode>), and compares them ON THE EUR SCALE ONLY:
# RMSE, MAE, mean exact CRPS, mean log score (with the Jacobian), 90%
# coverage and PIT sd, for the full sample and the subperiods, for a fixed
# set of models:
#   ar_dum, ar_dum_gas, ar_dum_gas_brent
#   comb_eq and comb_invmse of every family (ARX, R_MIDAS, R_MIDAS_L7)
#   pre-specified macro models of every family: gas + all surveys,
#     gas + all IPI, gas + all monthly indicators
#   R_MIDAS_L7 spec_000 (lag-structure effect, no macro)
# and tests each transformation against the normalized baseline with DM-HLN
# and GW on squared errors, CRPS and log score (the targets are identical
# across runs, so the tests apply directly).
#
# Usage (from midas_codes/):  Rscript compare_vst_runs.R
# Env: VST_RUNS = "normalized=estimation_v4,level=estimation_v4_it_level,
#                  log=estimation_v4_it_log,asinh=estimation_v4_it_asinh"
#      VST_BASE = "normalized"; COMPARE_CC = "it"
# Output: results_compare/compare_vst_<cc>.csv (metrics, long by mode) and
#         results_compare/compare_vst_tests_<cc>.csv (tests vs the baseline).
# =============================================================================
args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
here <- if (length(farg) > 0) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE))) else getwd()
Sys.setenv(V2_SCRIPT_DIR = here)
source(file.path(here, "common_utils.R"))

runs_env <- Sys.getenv("VST_RUNS",
  "normalized=estimation_v4,level=estimation_v4_it_level,log=estimation_v4_it_log,asinh=estimation_v4_it_asinh")
RUNS <- do.call(rbind, lapply(strsplit(strsplit(runs_env, ",")[[1]], "="), function(x)
  data.frame(mode = trimws(x[1]), dir = file.path(here, trimws(x[2])), stringsAsFactors = FALSE)))
RUNS <- RUNS[dir.exists(file.path(RUNS$dir, "forecasts")), ]
BASE <- Sys.getenv("VST_BASE", "normalized")
cc   <- Sys.getenv("COMPARE_CC", "it")
OUT  <- file.path(here, "results_compare"); dir.create(OUT, showWarnings = FALSE)
cat("runs found:", paste(RUNS$mode, collapse = ", "), "\n")
stopifnot(BASE %in% RUNS$mode)

PERIODS <- list(full = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
                pre = c(as.Date("1900-01-01"), as.Date("2020-12-31")),
                crisis = c(as.Date("2021-01-01"), as.Date("2022-12-31")),
                post = c(as.Date("2023-01-01"), as.Date("2100-01-01")))
FAMS <- c(ARX = "arx", R_MIDAS = "r_midas", R_MIDAS_L7 = "r_midas_l7")
SURVEYS <- c("pmi_mfg_it", "pmi_serv_it", "iesi_it", "istat_bci_it", "pmi_mfg_de", "pmi_serv_de", "ifo_bci_de")
IPI     <- c("Consumer_Goods", "Manufacturing", "Energy")

read_al <- function(fp, key) { if (!file.exists(fp)) return(NULL)
  d <- read.csv(fp, stringsAsFactors = FALSE); d[match(key, paste(d$origin_date, d$target_date)), , drop = FALSE] }

# per run and horizon: list of model -> list(point, sd, a, b); y, td, key from the AR file
load_run <- function(dir, h) {
  fd <- file.path(dir, "forecasts")
  ar <- read.csv(file.path(fd, sprintf("ar_%s_h%02d.csv", cc, h)), stringsAsFactors = FALSE)
  ar$target_date <- as.Date(ar$target_date); ar <- ar[order(ar$origin_date), ]
  key <- paste(ar$origin_date, ar$target_date)
  ars <- read_al(file.path(fd, sprintf("ar_%s_h%02d_sd.csv", cc, h)), key)
  vp  <- load_vst_params(fd, cc, ar$origin_date)
  M <- list()
  for (m in c("ar_dum", "ar_dum_gas", "ar_dum_gas_brent"))
    M[[m]] <- list(point = ar[[m]], sd = ars[[m]])
  for (famlab in names(FAMS)) {
    fam <- FAMS[[famlab]]
    fc <- read_al(file.path(fd, sprintf("%s_%s_h%02d.csv", fam, cc, h)), key)
    fs <- read_al(file.path(fd, sprintf("%s_%s_h%02d_sd.csv", fam, cc, h)), key)
    cb <- read_al(file.path(fd, sprintf("comb_%s_%s_h%02d.csv", fam, cc, h)), key)
    dp <- file.path(fd, sprintf("%s_spec_dictionary_%s.csv", fam, cc))
    if (is.null(fc) || !file.exists(dp)) next
    dict <- read.csv(dp, stringsAsFactors = FALSE)
    dict$hf_vars[is.na(dict$hf_vars)] <- ""; dict$lf_vars[is.na(dict$lf_vars)] <- ""
    lfset <- strsplit(dict$lf_vars, "|", fixed = TRUE)
    lf_all <- sort(unique(unlist(lfset)))
    pick <- function(hf, lf) { i <- which(dict$hf_vars == hf & vapply(lfset, function(s) setequal(s, lf), logical(1)))
      if (length(i)) dict$spec_id[i[1]] else NA }
    sel <- c(gas_surveys = pick("gas", intersect(lf_all, SURVEYS)),
             gas_ipi     = pick("gas", intersect(lf_all, IPI)),
             gas_allmonthly = pick("gas", lf_all))
    if (famlab == "R_MIDAS_L7") sel <- c(sel, nomacro = "spec_000")
    for (nm in names(sel)) if (!is.na(sel[[nm]]) && sel[[nm]] %in% names(fc))
      M[[paste(famlab, nm, sel[[nm]], sep = ":")]] <- list(point = fc[[sel[[nm]]]], sd = if (!is.null(fs)) fs[[sel[[nm]]]] else NULL)
    if (!is.null(cb)) for (m in c("comb_eq", "comb_invmse"))
      M[[paste(famlab, m, sep = ":")]] <- list(point = cb[[m]], sd = NULL)   # combinations: point only
  }
  list(y = ar$y_actual, td = ar$target_date, key = key, a = vp$a, b = vp$b, M = M)
}

rows <- list(); trows <- list()
for (h in HORIZONS) {
  R <- list()
  for (i in seq_len(nrow(RUNS))) {
    VST_MODE <<- RUNS$mode[i]
    R[[RUNS$mode[i]]] <- tryCatch(load_run(RUNS$dir[i], h), error = function(e) { cat("  ", RUNS$mode[i], "h", h, ":", conditionMessage(e), "\n"); NULL })
  }
  R <- Filter(Negate(is.null), R)
  if (!length(R) || is.null(R[[BASE]])) next
  base <- R[[BASE]]
  # per-observation scores of every (mode, model), aligned to the baseline's keys
  SC <- list()
  for (mode in names(R)) {
    VST_MODE <<- mode
    r <- R[[mode]]; idx <- match(base$key, r$key)
    for (m in names(r$M)) {
      pt <- r$M[[m]]$point[idx]; s <- if (!is.null(r$M[[m]]$sd)) r$M[[m]]$sd[idx] else NULL
      a <- r$a[idx]; b <- r$b[idx]
      e <- pt - base$y
      sc <- list(e = e)
      if (!is.null(s)) {
        ok <- is.finite(pt) & is.finite(s) & s > 0 & is.finite(base$y)
        mu <- vst_z(pt, a, b)
        crps <- rep(NA_real_, length(e)); ls <- crps; pit <- crps
        crps[ok] <- crps_exact(base$y[ok], mu[ok], s[ok], a[ok], b[ok])
        ls[ok]   <- logscore_price(base$y[ok], pt[ok], s[ok], a[ok], b[ok])
        pit[ok]  <- pnorm((vst_z(base$y[ok], a[ok], b[ok]) - mu[ok]) / s[ok])
        sc$crps <- crps; sc$logS <- ls; sc$pit <- pit
      }
      SC[[mode]][[m]] <- sc
    }
  }
  VST_MODE <<- BASE
  for (pn in names(PERIODS)) {
    sel <- base$td >= PERIODS[[pn]][1] & base$td <= PERIODS[[pn]][2]
    for (mode in names(SC)) for (m in names(SC[[mode]])) {
      sc <- SC[[mode]][[m]]; ok <- sel & is.finite(sc$e)
      row <- data.frame(country = cc, h = h, period = pn, mode = mode, model = m, n = sum(ok),
                        rmse = sqrt(mean(sc$e[ok]^2)), mae = mean(abs(sc$e[ok])),
                        crps = if (!is.null(sc$crps)) mean(sc$crps[ok], na.rm = TRUE) else NA_real_,
                        logS = if (!is.null(sc$logS)) mean(sc$logS[ok], na.rm = TRUE) else NA_real_,
                        cov90 = if (!is.null(sc$pit)) mean(sc$pit[ok] > 0.05 & sc$pit[ok] < 0.95, na.rm = TRUE) else NA_real_,
                        pit_sd = if (!is.null(sc$pit)) sd(sc$pit[ok], na.rm = TRUE) else NA_real_,
                        stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- row
      if (mode != BASE && m %in% names(SC[[BASE]])) {
        b0 <- SC[[BASE]][[m]]; n_min <- max(30L, 10L * h)
        tr <- data.frame(country = cc, h = h, period = pn, mode = mode, base = BASE, model = m, n = sum(ok), stringsAsFactors = FALSE)
        add <- function(tag, l1, l2) {
          k <- ok & is.finite(l1) & is.finite(l2)
          dm <- if (sum(k) >= n_min) dm_hln_loss(l1[k], l2[k], h) else list(stat = NA, pvalue = NA)
          gw <- if (sum(k) >= n_min) gw_test_loss(l1, l2, h, keep = k) else list(stat = NA, pvalue = NA, dbar = NA)
          tr[[paste0(tag, "_dm_stat")]] <<- dm$stat; tr[[paste0(tag, "_dm_pval")]] <<- dm$pvalue
          tr[[paste0(tag, "_gw_stat")]] <<- gw$stat; tr[[paste0(tag, "_gw_pval")]] <<- gw$pvalue
          tr[[paste0(tag, "_gw_dbar")]] <<- gw$dbar
        }
        add("sqerr", sc$e^2, b0$e^2)
        if (!is.null(sc$crps) && !is.null(b0$crps)) { add("crps", sc$crps, b0$crps); add("logS", -sc$logS, -b0$logS) }
        trows[[length(trows) + 1L]] <- tr
      }
    }
  }
  cat(sprintf("h=%02d done (%d modes)\n", h, length(SC)))
}
D <- do.call(rbind, rows); write.csv(D, file.path(OUT, sprintf("compare_vst_%s.csv", cc)), row.names = FALSE)
if (length(trows)) { T_ <- do.call(rbind, trows); write.csv(T_, file.path(OUT, sprintf("compare_vst_tests_%s.csv", cc)), row.names = FALSE) }
cat(sprintf("wrote %s (%d rows)\n", file.path(OUT, sprintf("compare_vst_%s.csv", cc)), nrow(D)))

# console digest: full sample, ar_dum / ar_dum_gas / combos per mode
cat("\nfull sample, RMSE | exact CRPS | logS | cov90 by transformation (negative dm = mode better than", BASE, "):\n")
for (h in HORIZONS) for (m in c("ar_dum", "ar_dum_gas", "R_MIDAS:comb_invmse", "R_MIDAS_L7:comb_invmse", "ARX:comb_invmse")) {
  g <- D[D$h == h & D$period == "full" & D$model == m, ]
  if (!nrow(g)) next
  cat(sprintf("  h=%02d %-24s ", h, m))
  for (i in seq_len(nrow(g))) cat(sprintf("%s: %.2f | %.2f | %.3f | %.2f   ", g$mode[i], g$rmse[i], g$crps[i], g$logS[i], g$cov90[i]))
  cat("\n")
}
cat("\nCompare done.\n")
