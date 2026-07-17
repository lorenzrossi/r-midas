#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — ORCHESTRATOR (stage-selectable)
#
# DEFAULT (current project focus): AR benchmark + updated standard R-MIDAS
# only, followed by the AR+R-MIDAS summary and the combination / subperiod /
# Giacomini-Rossi evaluation restricted to the r_midas family.
#
#   ar           ar_benchmark.R      AR(1,2,3,7) + calendar; ar_dum = benchmark
#   r_midas      r_midas.R           STANDARD R-MIDAS — exp-Almon on the AR lag
#                                    set J = {1,2,3,7} of y (the j values enter
#                                    the Almon exponent, see model.tex);
#                                    gas_lag1 / brent_lag1 linear; LF linear;
#                                    18 calendar dummies.  Family: "R_MIDAS".
#   summary      summary_table_ar_r_midas.R   AR + R_MIDAS metrics only
#   combination  combination_subperiod_gr.R   combos + subperiods + GR test
#                                    (FAMILIES defaults to "r_midas" here)
#
# Other stages, NOT run by default but selectable via the STAGES env var:
#   ru_midas     ru_midas.R            RU-MIDAS position dummies
#   r_midas_ext  r_midas_extended.R    extended R-MIDAS (Almon on y AND HF)
#   summary_all  summary_table.R       summary across all families
#
# Usage:
#   Rscript run_all.R                                    # default subset
#   STAGES="ar,r_midas,ru_midas,summary_all" Rscript run_all.R
#
# Knobs: COUNTRY_NAME ("Germany" | "Italy" | "ALL"), N_CORES, SPEC_LIMIT,
# WINDOW_DAYS, RMIDAS_DISCOUNT — see common_utils.R.
#
# Runtime estimate (one country, N_CORES = 9, 256 specs):
#   ar          : < 1 min
#   r_midas     : with the FWL-accelerated NLS expect tens of minutes rather
#                 than the 8-35 h of the pre-acceleration runs
#   summary     : < 30 sec
#   combination : ~1-2 min
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
  ar          = "ar_benchmark.R",
  r_midas     = "r_midas.R",
  ru_midas    = "ru_midas.R",
  r_midas_ext = "r_midas_extended.R",
  summary     = "summary_table_ar_r_midas.R",
  summary_all = "summary_table.R",
  combination = "combination_subperiod_gr.R"
)
DEFAULT_ORDER <- c("ar", "r_midas", "summary", "combination")

stages_env <- Sys.getenv("STAGES", unset = "")
order <- if (nzchar(stages_env)) {
  trimws(strsplit(stages_env, ",")[[1]])
} else DEFAULT_ORDER

t0 <- Sys.time()
cat("=== v2 pipeline start: ", format(t0), " ===\n")
cat("stages: ", paste(order, collapse = " -> "), "\n")

for (st in order) {
  if (!st %in% names(STAGE_FILES)) {
    cat(sprintf("\nUnknown stage '%s' (known: %s) — skipped.\n",
                st, paste(names(STAGE_FILES), collapse = ", ")))
    next
  }
  fp <- file.path(script_dir, STAGE_FILES[[st]])
  cat(sprintf("\n[%s] %s\n", st, STAGE_FILES[[st]]))
  if (!file.exists(fp)) {
    cat("  file not found — stage skipped.\n")
    next
  }
  # Restrict the combination stage to the r_midas family unless the user
  # already set FAMILIES explicitly.
  old_families <- Sys.getenv("FAMILIES", unset = NA_character_)
  set_families_here <- st == "combination" && is.na(old_families)
  if (set_families_here) Sys.setenv(FAMILIES = "r_midas")
  ts <- Sys.time()
  tryCatch(
    source(fp),
    finally = if (set_families_here) Sys.unsetenv("FAMILIES")
  )
  cat(sprintf("  stage '%s' done in %s\n", st,
              format(round(difftime(Sys.time(), ts, units = "mins"), 1))))
}

t1 <- Sys.time()
cat(sprintf("\n=== v2 pipeline done in %s ===\n",
            format(round(difftime(t1, t0, units = "mins"), 1))))
