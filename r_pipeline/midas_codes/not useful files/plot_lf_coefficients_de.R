#!/usr/bin/env Rscript

# =============================================================================
# Plot low-frequency coefficients over time (Germany)
#
# Creates one plot per screenshot model using Germany coefficients.
# Default inputs:
#   - results/coefficients/r_midas_de.csv
#   - forecasts/r_midas_spec_dictionary_de.csv
#
# Override with env vars:
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
  unset = file.path(script_root, "results", "plots", "lf_coefficients_de")
)
summary_file <- Sys.getenv(
  "SUMMARY_FILE",
  unset = file.path(script_root, "results", "summary_table_ar_r_midas.csv")
)

if (!file.exists(coef_file)) stop("Missing coefficient file: ", coef_file)
if (!file.exists(dict_file)) stop("Missing spec dictionary file: ", dict_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coefs <- read.csv(coef_file, stringsAsFactors = FALSE, check.names = FALSE)
specs <- read.csv(dict_file, stringsAsFactors = FALSE, check.names = FALSE)
summary_df <- if (file.exists(summary_file)) {
  read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}

if (!("origin_date" %in% names(coefs))) stop("Column 'origin_date' not found in coefficients file.")
if (!("spec_id" %in% names(coefs))) stop("Column 'spec_id' not found in coefficients file.")
if (!all(c("spec_id", "hf_vars", "lf_vars") %in% names(specs))) {
  stop("Dictionary must contain columns: spec_id, hf_vars, lf_vars.")
}

coefs$origin_date <- as.Date(coefs$origin_date)

model_defs <- list(
  list(
    label = "R-MIDAS + gas + brent + PMI MFG + ifo BCI + IPI (C.Goods, Mfg, Energy)",
    hf_vars = "gas|brent",
    lf_vars = "pmi_mfg_de|ifo_bci_de|Consumer_Goods|Manufacturing|Energy"
  ),
  list(
    label = "R-MIDAS + gas + PMI MFG + ifo BCI + IPI (C.Goods, Mfg)",
    hf_vars = "gas",
    lf_vars = "pmi_mfg_de|ifo_bci_de|Consumer_Goods|Manufacturing"
  ),
  list(
    label = "R-MIDAS + gas + PMI MFG + ifo BCI + IPI (Mfg)",
    hf_vars = "gas",
    lf_vars = "pmi_mfg_de|ifo_bci_de|Manufacturing"
  ),
  list(
    label = "R-MIDAS + gas + PMI MFG + ifo BCI",
    hf_vars = "gas",
    lf_vars = "pmi_mfg_de|ifo_bci_de"
  )
)

sanitize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

format_var_label <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

resolve_spec <- function(hf_vars, lf_vars) {
  idx <- which(specs$hf_vars == hf_vars & specs$lf_vars == lf_vars)
  if (length(idx) == 0L) return(NA_character_)
  if (length(idx) > 1L) {
    # Defensive tie-break: pick lexicographically first spec_id.
    sid <- sort(specs$spec_id[idx])[1]
    return(sid)
  }
  specs$spec_id[idx]
}

palette_cols <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf")

summary_rows <- list()
resolved_specs <- rep(NA_character_, length(model_defs))

for (k in seq_along(model_defs)) {
  md <- model_defs[[k]]
  resolved_specs[k] <- resolve_spec(md$hf_vars, md$lf_vars)
}

winning_h_by_spec <- setNames(vector("list", length(resolved_specs)), resolved_specs)
for (sid in resolved_specs) winning_h_by_spec[[sid]] <- integer(0)

if (!is.null(summary_df) && all(c("country", "family", "model", "h", "rmse") %in% names(summary_df))) {
  horizon_set <- c(1, 7, 14, 21, 28)
  cand <- summary_df[
    summary_df$country == "de" &
      summary_df$family == "R_MIDAS" &
      summary_df$model %in% resolved_specs &
      summary_df$h %in% horizon_set &
      is.finite(summary_df$rmse),
    c("model", "h", "rmse"),
    drop = FALSE
  ]
  if (nrow(cand) > 0L) {
    for (hh in horizon_set) {
      ch <- cand[cand$h == hh, , drop = FALSE]
      if (nrow(ch) == 0L) next
      ch <- ch[order(ch$rmse, ch$model), , drop = FALSE]
      best_model <- as.character(ch$model[1])
      winning_h_by_spec[[best_model]] <- c(winning_h_by_spec[[best_model]], hh)
    }
  }
}

for (k in seq_along(model_defs)) {
  md <- model_defs[[k]]
  sid <- resolved_specs[k]
  if (!is.finite(match(sid, coefs$spec_id)) && is.na(sid)) {
    warning("Spec not found for model: ", md$label)
    next
  }

  block <- coefs[coefs$spec_id == sid, , drop = FALSE]
  if (nrow(block) == 0L) {
    warning("No coefficient rows found for spec_id ", sid, " (", md$label, ").")
    next
  }

  lf_vars <- strsplit(md$lf_vars, "\\|", fixed = FALSE)[[1]]
  lf_vars <- lf_vars[lf_vars %in% names(block)]
  if (length(lf_vars) == 0L) {
    warning("No low-frequency coefficient columns found for spec_id ", sid, ".")
    next
  }

  block <- block[order(block$origin_date), , drop = FALSE]
  ymat <- as.matrix(block[, lf_vars, drop = FALSE])
  mode(ymat) <- "numeric"
  lf_vars_pretty <- vapply(lf_vars, format_var_label, character(1))
  winning_h <- winning_h_by_spec[[sid]]
  subtitle <- if (length(winning_h) == 0L) {
    sprintf("%s (lowest RMSE at h = n/a)", paste(lf_vars_pretty, collapse = ", "))
  } else {
    sprintf(
      "%s (lowest RMSE at h = %s)",
      paste(lf_vars_pretty, collapse = ", "),
      paste(winning_h, collapse = ", ")
    )
  }
  year_ticks <- as.Date(sprintf("%d-01-01", 2015:2025))

  out_file <- file.path(out_dir, paste0(sprintf("%02d", k), "_", sanitize_name(md$label), ".png"))
  png(out_file, width = 1400, height = 800, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(5.0, 5.0, 5.5, 11.0))

  matplot(
    x = block$origin_date,
    y = ymat,
    type = "l",
    lty = 1,
    lwd = 2,
    col = palette_cols[seq_len(ncol(ymat))],
    xaxt = "n",
    xlim = as.Date(c("2015-01-01", "2025-12-31")),
    xlab = "Forecast origin date",
    ylab = "Estimated coefficient",
    main = "R-MIDAS LF coefficients"
  )
  axis.Date(1, at = year_ticks, format = "%Y")
  abline(h = 0, col = "gray50", lty = 2, lwd = 1)
  mtext(subtitle, side = 3, line = 0.8, cex = 0.85)
  legend(
    "topright",
    inset = c(-0.34, 0),
    xpd = NA,
    legend = lf_vars_pretty,
    col = palette_cols[seq_len(ncol(ymat))],
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.9
  )
  dev.off()

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    model_label = md$label,
    spec_id = sid,
    hf_vars = md$hf_vars,
    lf_vars = md$lf_vars,
    lowest_rmse_h = if (length(winning_h) == 0L) NA_character_ else paste(winning_h, collapse = ","),
    output_file = out_file,
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows) == 0L) {
  stop("No plots produced. Check model definitions and dictionary coverage.")
}

map_df <- do.call(rbind, summary_rows)
map_path <- file.path(out_dir, "model_to_spec_map_de.csv")
write.csv(map_df, map_path, row.names = FALSE)

cat("Saved plots in: ", out_dir, "\n", sep = "")
cat("Saved model map: ", map_path, "\n", sep = "")
