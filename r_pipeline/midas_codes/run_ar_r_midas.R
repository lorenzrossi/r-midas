#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — AR + R-MIDAS STANDARD ONLY
#
# Subset orchestrator mirroring run_all.R, restricted to:
#   1) ar_benchmark.R   AR family (benchmark: ar_dum)
#   2) r_midas.R        Standard R-MIDAS (Almon on y, linear HF)
#   3) summary_table_ar_r_midas.R
#
# Skips RU-MIDAS and R-MIDAS extended.  Useful for faster reruns or while the
# full pipeline is busy on other steps.
#
# Knobs: same as run_all.R (common_utils.R, N_CORES, SPEC_LIMIT, COUNTRY_NAME).
#
# Country: COUNTRY_NAME env var ("Germany", "Italy", or "ALL" [default]).
#
# Runtime estimate (COUNTRY_NAME = ALL, N_CORES = 11, 256 specs, 5-start NLS):
#   AR family        : < 1 min
#   R-MIDAS standard : ~8–35 h per country (multi-start NLS; varies by machine)
#   summary          : < 30 sec
#
# Outputs:
#   forecasts/ar_<cc>_h<h>.csv
#   forecasts/r_midas_<cc>_h<h>.csv
#   forecasts/r_midas_spec_dictionary_<cc>.csv
#   results/summary_table_ar_r_midas.csv
#   results/summary_top5_ar_r_midas_per_country_h.csv
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
            "v2_paths.R")
} else file.path(getwd(), "v2_paths.R")
source(v2.boot.path)

t0 <- Sys.time()
cat("=== v2 AR + R-MIDAS pipeline start: ", format(t0), " ===\n")

cat("\n[1/3] AR benchmark family\n")
source(file.path(script_dir, "ar_benchmark.R"))

cat("\n[2/3] R-MIDAS standard (Almon on y, linear HF)\n")
source(file.path(script_dir, "r_midas.R"))

cat("\n[3/3] Summary table (AR + R-MIDAS only)\n")
source(file.path(script_dir, "summary_table_ar_r_midas.R"))

t1 <- Sys.time()
cat(sprintf("\n=== v2 AR + R-MIDAS pipeline done in %s mins ===\n",
            format(round(difftime(t1, t0, units = "mins"), 1))))
