#!/usr/bin/env Rscript
# =============================================================================
# v2/v4 PIPELINE — ORCHESTRATOR (stage-selectable)
#
# Every stage runs in its OWN Rscript process (not source()d into this one):
# the estimation stages leave hundreds of MB of objects behind, and forking
# mclapply workers from such a parent made the later stages crawl (the v3
# density stage took ~10 h instead of ~10 min).  A fresh process per stage
# also isolates the environment variables the family wrappers set.
#
# DEFAULT order:
#   ols_check    ols_coef_check.R    plain-OLS coefficients of every regressor
#                                    (full sample, subperiods, rolling, HAC)
#                                    + free AR(1..7) profile
#   ar           ar_benchmark.R      AR(1,2,3,7) + calendar; ar_dum = benchmark
#   arx          arx_free_lags.R     same 256 specs as r_midas with the four
#                                    AR coefficients FREE (OLS).  Family: "ARX".
#   r_midas      r_midas.R           STANDARD R-MIDAS — exp-Almon on the AR lag
#                                    set J = {1,2,3,7}; gas_lag1 / brent_lag1
#                                    linear; LF linear; 18 calendar dummies.
#   r_midas_l7   r_midas_l7.R        R-MIDAS with the Almon polynomial on the
#                                    consecutive lags 1..7.  Family: "R_MIDAS_L7".
#   combination  combination_subperiod_gr.R   combos + subperiods (DM, GW) +
#                                    GR fluctuation test, every family
#   summary      summary_table_ar_r_midas.R   point-forecast metrics, all families
#   density      density_eval.R      log score / CRPS / ACPS / PIT, DM + GW vs
#                                    benchmark and vs daily-only twin
#   tests        tests_dm_gw.R       Table A (vs benchmark) / Table B (vs twin),
#                                    DM-HLN + GW by subperiod, all families
#   cssed        cssed_plots.R       cumulative squared-error-difference plots
#
# Other stages, selectable via STAGES: ru_midas, r_midas_ext, summary_all.
#
# Usage:
#   Rscript run_all.R                                    # default order
#   STAGES="ar,r_midas,summary" Rscript run_all.R
#
# Knobs: COUNTRY_NAME ("Germany" | "Italy" | "ALL"), N_CORES, SPEC_LIMIT,
# WINDOW_DAYS, VST_MODE (normalized | raw), LF_IT_SURVEY (iesi_it |
# istat_bci_it), RMIDAS_THETA2_MAX (8 = theta2 free, 0 = old restriction),
# RMIDAS_DISCOUNT, RMIDAS_HUBER, FAMILIES (combination stage; default = the
# families whose estimation stages are in the run).
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
            "v2_paths.R")
} else file.path(getwd(), "v2_paths.R")
source(v2.boot.path)

STAGE_FILES <- c(
  ols_check   = "ols_coef_check.R",
  ar          = "ar_benchmark.R",
  arx         = "arx_free_lags.R",
  r_midas     = "r_midas.R",
  r_midas_l7  = "r_midas_l7.R",
  ru_midas    = "ru_midas.R",
  r_midas_ext = "r_midas_extended.R",
  combination = "combination_subperiod_gr.R",
  summary     = "summary_table_ar_r_midas.R",
  summary_all = "summary_table.R",
  density     = "density_eval.R",
  tests       = "tests_dm_gw.R",
  cssed       = "cssed_plots.R"
)
DEFAULT_ORDER <- c("ols_check", "ar", "arx", "r_midas", "r_midas_l7",
                   "combination", "summary", "density", "tests", "cssed")

stages_env <- Sys.getenv("STAGES", unset = "")
order <- if (nzchar(stages_env)) trimws(strsplit(stages_env, ",")[[1]]) else DEFAULT_ORDER

# Families handed to the combination stage (unless FAMILIES is set by the user)
fam_stage <- c(r_midas = "r_midas", arx = "arx", r_midas_l7 = "r_midas_l7")
families_default <- paste(fam_stage[names(fam_stage) %in% order], collapse = ",")
if (!nzchar(families_default)) families_default <- "r_midas"

rscript <- file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript)) rscript <- Sys.which("Rscript")

t0 <- Sys.time()
cat("=== v2 pipeline start: ", format(t0), " ===\n")
cat("script dir: ", script_dir, "\n")
cat("stages: ", paste(order, collapse = " -> "), "\n")
cat(sprintf("settings: VST_MODE=%s LF_IT_SURVEY=%s RMIDAS_THETA2_MAX=%s WINDOW_DAYS=%s N_CORES=%s FAMILIES=%s\n",
            Sys.getenv("VST_MODE", "normalized"), Sys.getenv("LF_IT_SURVEY", "iesi_it"),
            Sys.getenv("RMIDAS_THETA2_MAX", "8"), Sys.getenv("WINDOW_DAYS", "3650"),
            Sys.getenv("N_CORES", "auto"), Sys.getenv("FAMILIES", families_default)))

failed <- character(0)
for (st in order) {
  if (!st %in% names(STAGE_FILES)) {
    cat(sprintf("\nUnknown stage '%s' (known: %s) — skipped.\n",
                st, paste(names(STAGE_FILES), collapse = ", ")))
    next
  }
  fp <- file.path(script_dir, STAGE_FILES[[st]])
  cat(sprintf("\n[%s] %s\n", st, STAGE_FILES[[st]]))
  if (!file.exists(fp)) { cat("  file not found — stage skipped.\n"); next }
  env <- character(0)
  if (st == "combination" && !nzchar(Sys.getenv("FAMILIES", unset = "")))
    env <- c(env, paste0("FAMILIES=", families_default))
  ts <- Sys.time()
  status <- system2(rscript, args = shQuote(fp), env = env,
                    stdout = "", stderr = "")
  el <- format(round(difftime(Sys.time(), ts, units = "mins"), 1))
  if (!identical(as.integer(status), 0L)) {
    failed <- c(failed, st)
    cat(sprintf("  stage '%s' FAILED (exit status %s) after %s\n", st, status, el))
  } else cat(sprintf("  stage '%s' done in %s\n", st, el))
}

t1 <- Sys.time()
cat(sprintf("\n=== v2 pipeline done in %s ===\n",
            format(round(difftime(t1, t0, units = "mins"), 1))))
if (length(failed)) {
  cat("FAILED stages: ", paste(failed, collapse = ", "), "\n")
  quit(status = 1L)
}
