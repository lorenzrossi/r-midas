#!/usr/bin/env Rscript
# =============================================================================
# R-MIDAS with the exponential Almon polynomial on the CONSECUTIVE daily lags
# 1..7 (family R_MIDAS_L7).  Same 256 specifications, theta box, starts,
# iterated forecasts, normalised asinh and calendar dummies as r_midas.R —
# this wrapper only sets the family parameters and sources the shared code,
# so the two families cannot drift apart.
#
# Why: the exp-Almon curve is designed for a run of consecutive lags; on the
# gappy set {1,2,3,7} the smoothness restriction has little meaning.  With
# lags 1..7 it compresses 7 coefficients into (delta, theta1, theta2) and
# uses the same information as the benchmark (prices up to one week back).
#
# Outputs carry the prefix r_midas_l7 (forecasts, _sd, spec dictionary,
# results/coefficients/r_midas_l7_<cc>.csv with phi_1 ... phi_7).
# =============================================================================
old <- Sys.getenv(c("RMIDAS_LAGS", "RMIDAS_PREFIX", "RMIDAS_FAMILY"), unset = NA)
Sys.setenv(RMIDAS_LAGS = "1:7", RMIDAS_PREFIX = "r_midas_l7", RMIDAS_FAMILY = "R_MIDAS_L7")
args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
here <- if (length(farg) > 0) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE),
                        mustWork = FALSE))
} else getwd()
env_dir <- Sys.getenv("V2_SCRIPT_DIR", unset = "")
if (nzchar(env_dir) && file.exists(file.path(env_dir, "r_midas.R"))) here <- env_dir
source(file.path(here, "r_midas.R"), chdir = TRUE)
# restore the environment so that a later stage in the same session is not
# affected (stages normally run in separate Rscript processes anyway)
for (nm in names(old)) if (is.na(old[[nm]])) Sys.unsetenv(nm) else
  do.call(Sys.setenv, as.list(setNames(old[[nm]], nm)))
