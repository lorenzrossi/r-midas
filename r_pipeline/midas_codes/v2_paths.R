# Resolve v2_new_dataset root (iCloud-safe: ~+~ in --file= paths).
# Sets script_dir in the caller's environment; honours V2_SCRIPT_DIR.
resolve_v2_script_dir <- function() {
  env <- Sys.getenv("V2_SCRIPT_DIR", unset = "")
  if (nzchar(env)) {
    cand <- suppressWarnings(normalizePath(env, winslash = "/", mustWork = FALSE))
    if (!is.na(cand) && nzchar(cand) && dir.exists(cand)) return(cand)
  }
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grep("^--file=", args)]
  if (length(farg) > 0) {
    raw <- sub("^--file=", "", farg[1])
    raw <- gsub("~+~", " ", raw, fixed = TRUE)
    cand <- suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))
    if (!is.na(cand) && nzchar(cand) && file.exists(cand)) return(dirname(cand))
  }
  getwd()
}

script_dir <- resolve_v2_script_dir()
if (dir.exists(script_dir)) setwd(script_dir)
Sys.setenv(V2_SCRIPT_DIR = script_dir)
