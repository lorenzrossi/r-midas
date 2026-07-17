#!/usr/bin/env Rscript

# =============================================================================
# Plot Germany LF coefficient paths:
#   1) LF variable alone (no HF, no other LF)
#   2) LF variable with gas only (HF = gas, single LF)
#
# Default inputs:
#   - results/coefficients/r_midas_de.csv
#   - forecasts/r_midas_spec_dictionary_de.csv
#
# Override via env vars:
#   COEF_FILE, SPEC_DICT_FILE, OUT_DIR
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(
    dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
    "v2_paths.R"
  )
} else {
  file.path(getwd(), "v2_paths.R")
}
source(v2.boot.path)

script_root <- script_dir

coef_file <- Sys.getenv(
  "COEF_FILE",
  unset = file.path(script_root, "results", "coefficients", "r_midas_de.csv")
)
dict_file <- Sys.getenv(
  "SPEC_DICT_FILE",
  unset = file.path(script_root, "forecasts", "r_midas_spec_dictionary_de.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(script_root, "results", "plots", "lf_alone_vs_gas_de")
)

if (!file.exists(coef_file)) stop("Missing coefficient file: ", coef_file)
if (!file.exists(dict_file)) stop("Missing spec dictionary file: ", dict_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coefs <- read.csv(coef_file, stringsAsFactors = FALSE, check.names = FALSE)
specs <- read.csv(dict_file, stringsAsFactors = FALSE, check.names = FALSE)

required_coef <- c("origin_date", "spec_id")
required_dict <- c("spec_id", "hf_vars", "lf_vars")
if (!all(required_coef %in% names(coefs))) {
  stop("Coefficients file must contain: ", paste(required_coef, collapse = ", "))
}
if (!all(required_dict %in% names(specs))) {
  stop("Dictionary file must contain: ", paste(required_dict, collapse = ", "))
}

coefs$origin_date <- as.Date(coefs$origin_date)

sanitize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

pretty_name <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

find_spec_id <- function(hf_vars, lf_vars) {
  idx <- which(specs$hf_vars == hf_vars & specs$lf_vars == lf_vars)
  if (length(idx) == 0L) return(NA_character_)
  sort(specs$spec_id[idx])[1]
}

lf_vars_all <- sort(unique(specs$lf_vars[nchar(specs$lf_vars) > 0L & !grepl("\\|", specs$lf_vars)]))
lf_vars_all <- lf_vars_all[lf_vars_all %in% names(coefs)]
if (length(lf_vars_all) == 0L) {
  stop("No single LF variables found in dictionary that match coefficient columns.")
}

map_rows <- list()

for (i in seq_along(lf_vars_all)) {
  lf_var <- lf_vars_all[i]
  sid_alone <- find_spec_id("", lf_var)
  sid_gas <- find_spec_id("gas", lf_var)

  if (is.na(sid_alone) || is.na(sid_gas)) {
    warning("Skipping ", lf_var, ": missing spec for alone and/or gas setup.")
    next
  }

  b0 <- coefs[coefs$spec_id == sid_alone, c("origin_date", lf_var), drop = FALSE]
  b1 <- coefs[coefs$spec_id == sid_gas, c("origin_date", lf_var), drop = FALSE]
  names(b0)[2] <- "coef_alone"
  names(b1)[2] <- "coef_gas"

  merged <- merge(b0, b1, by = "origin_date", all = TRUE, sort = TRUE)
  if (nrow(merged) == 0L) {
    warning("No rows for ", lf_var, " after merging alone vs gas paths.")
    next
  }

  out_file <- file.path(out_dir, sprintf("%02d_%s_alone_vs_gas.png", i, sanitize_name(lf_var)))
  png(out_file, width = 1400, height = 800, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(5.0, 5.0, 4.8, 1.5))

  year_ticks <- as.Date(sprintf("%d-01-01", 2015:2025))
  x_rng <- as.Date(c("2015-01-01", "2025-12-31"))
  y_min <- min(c(merged$coef_alone, merged$coef_gas), na.rm = TRUE)
  y_max <- max(c(merged$coef_alone, merged$coef_gas), na.rm = TRUE)
  if (!is.finite(y_min) || !is.finite(y_max)) {
    y_min <- -1
    y_max <- 1
  }
  if (y_min == y_max) {
    y_min <- y_min - 0.1
    y_max <- y_max + 0.1
  }

  plot(
    merged$origin_date, merged$coef_alone,
    type = "l",
    lwd = 2,
    col = "#1f77b4",
    xaxt = "n",
    xlim = x_rng,
    ylim = c(y_min, y_max),
    xlab = "Forecast origin date",
    ylab = "Estimated coefficient",
    main = sprintf("LF coefficient path: %s", pretty_name(lf_var))
  )
  lines(merged$origin_date, merged$coef_gas, lwd = 2, col = "#d62728")
  axis.Date(1, at = year_ticks, format = "%Y")
  abline(h = 0, col = "gray50", lty = 2, lwd = 1)
  legend(
    "topright",
    legend = c("LF alone", "LF + gas"),
    col = c("#1f77b4", "#d62728"),
    lty = 1,
    lwd = 2,
    bty = "n"
  )
  dev.off()

  map_rows[[length(map_rows) + 1L]] <- data.frame(
    lf_variable = lf_var,
    spec_id_alone = sid_alone,
    spec_id_gas = sid_gas,
    output_file = out_file,
    stringsAsFactors = FALSE
  )
}

if (length(map_rows) == 0L) {
  stop("No plots produced. Check dictionary/spec coverage.")
}

map_df <- do.call(rbind, map_rows)
map_path <- file.path(out_dir, "lf_alone_vs_gas_map_de.csv")
write.csv(map_df, map_path, row.names = FALSE)

cat("Saved plots in: ", out_dir, "\n", sep = "")
cat("Saved map: ", map_path, "\n", sep = "")
