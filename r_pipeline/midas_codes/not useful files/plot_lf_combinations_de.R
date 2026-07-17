#!/usr/bin/env Rscript

# =============================================================================
# Germany LF coefficient combinations — mix of "alone" and "with gas" series
#
# Eight plots, all reading from the same r_midas_de.csv coefficient file
# and the matching dictionary:
#
#  1) "pmi_mfg_de + ifo_bci_de (mixed scope)":
#       ifo_bci_de    → coefficient from the "alone" spec   (no HF, only LF)
#       pmi_mfg_de    → coefficient from the "with gas" spec (HF=gas, single LF)
#     Title/subtitle/colors match the user-provided PNG.
#
#  2-5) Singular plots, one per OTHER LF variable, each from the "with gas"
#       single-LF spec:
#          pmi_serv_de, Consumer_Goods, Manufacturing, Energy.
#
#  6-8) Group plots — coefficients pulled from a SINGLE spec where the listed
#       LF variables and gas are all in the same regression:
#          6) pmi_mfg_de + ifo_bci_de + Manufacturing
#          7) pmi_mfg_de + ifo_bci_de + Manufacturing + Consumer_Goods
#          8) pmi_mfg_de + ifo_bci_de + Manufacturing + Consumer_Goods + Energy
#
# Same per-variable colour assignment across every plot so the same indicator
# is always the same colour.
#
# Path resolution priority (any one of these works):
#   1) V2_DATA_DIR env var set to the v2_new_dataset folder
#   2) v2_paths.R sitting next to this script (the original convention)
#   3) Run from inside v2_new_dataset/ (getwd() pickup)
#   4) Explicit COEF_FILE / SPEC_DICT_FILE / OUT_DIR env vars
#
# So either:
#   cd "/path/to/v2_new_dataset"
#   Rscript plot_lf_combinations_de.R
# or:
#   V2_DATA_DIR="/path/to/v2_new_dataset" Rscript plot_lf_combinations_de.R
# =============================================================================

# -----------------------------------------------------------------------------
# Bootstrap — robust path resolution
# -----------------------------------------------------------------------------
resolve_script_root <- function() {
  # 1) Explicit env var wins
  env_root <- Sys.getenv("V2_DATA_DIR", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) return(normalizePath(env_root, mustWork = FALSE))
  
  # 2) Look for v2_paths.R sitting next to this script (original convention)
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grep("^--file=", args)]
  if (length(farg) > 0L) {
    raw <- sub("^--file=", "", farg[1L])
    raw <- gsub("~+~", " ", raw, fixed = TRUE)
    script_dir_cand <- tryCatch(
      dirname(normalizePath(raw, winslash = "/", mustWork = FALSE)),
      error = function(e) NA_character_
    )
    if (!is.na(script_dir_cand) && nzchar(script_dir_cand) && dir.exists(script_dir_cand)) {
      boot_path <- file.path(script_dir_cand, "v2_paths.R")
      if (file.exists(boot_path)) {
        # Source it in the parent frame so any `script_dir` it sets is visible
        eval(parse(file = boot_path), envir = globalenv())
        if (exists("script_dir", envir = globalenv())) {
          sd <- get("script_dir", envir = globalenv())
          if (nzchar(sd) && dir.exists(sd)) return(normalizePath(sd, mustWork = FALSE))
        }
      }
      # No v2_paths.R but we have the script directory — use it directly
      return(script_dir_cand)
    }
  }
  
  # 3) Fall back to the current working directory
  wd <- getwd()
  if (dir.exists(wd)) return(wd)
  
  stop("Cannot resolve project root.  Set V2_DATA_DIR to the v2_new_dataset folder.")
}

script_root <- resolve_script_root()

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
  unset = file.path(script_root, "results", "plots", "lf_combinations_de")
)

# Sanity diagnostics so a misconfigured run is easy to debug
cat("Project root :", script_root, "\n")
cat("Coef file    :", coef_file, "\n")
cat("Dict file    :", dict_file, "\n")
cat("Output dir   :", out_dir,   "\n\n")

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

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
sanitize_name <- function(x) {
  x <- tolower(x); x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}
pretty_name <- function(x) gsub("_", " ", x, fixed = TRUE)

# Resolve a spec_id by exact (hf_vars, lf_vars) match in the dictionary.
# hf_vars / lf_vars are pipe-separated strings; "" means "no HF" / "no LF".
find_spec_id <- function(hf_vars, lf_vars) {
  idx <- which(specs$hf_vars == hf_vars & specs$lf_vars == lf_vars)
  if (length(idx) == 0L) return(NA_character_)
  sort(specs$spec_id[idx])[1]
}

# Pull a single LF coefficient time series from a given spec_id.
pull_series <- function(sid, lf_var) {
  if (is.na(sid)) return(NULL)
  if (!(lf_var %in% names(coefs))) return(NULL)
  blk <- coefs[coefs$spec_id == sid, c("origin_date", lf_var), drop = FALSE]
  if (nrow(blk) == 0L) return(NULL)
  blk <- blk[order(blk$origin_date), , drop = FALSE]
  names(blk)[2] <- "value"
  blk
}

# Fixed colour per LF variable.
COL_BY_LF <- c(
  pmi_mfg_de      = "#1f77b4",
  ifo_bci_de      = "#d62728",
  Consumer_Goods  = "#2ca02c",
  Manufacturing   = "#9467bd",
  Energy          = "#ff7f0e",
  pmi_serv_de     = "#17becf"
)

YEAR_TICKS <- as.Date(sprintf("%d-01-01", 2015:2025))
X_RANGE    <- as.Date(c("2015-01-01", "2025-12-31"))

plot_panel <- function(series, subtitle, legend_lab, out_file,
                       title = "R-MIDAS LF coefficients",
                       coef_plot_scale = 1) {
  cols <- COL_BY_LF[names(series)]
  x_all <- sort(unique(do.call(c, lapply(series, function(s) s$origin_date))))
  ymat  <- vapply(
    names(series),
    function(nm) {
      s <- series[[nm]]
      s$value[match(x_all, s$origin_date)]
    },
    numeric(length(x_all))
  )
  if (!is.matrix(ymat)) ymat <- matrix(ymat, ncol = length(series))
  colnames(ymat) <- names(series)
  if (is.finite(coef_plot_scale) && coef_plot_scale != 1) {
    ymat <- ymat * coef_plot_scale
  }
  
  y_min <- suppressWarnings(min(ymat, na.rm = TRUE))
  y_max <- suppressWarnings(max(ymat, na.rm = TRUE))
  if (!is.finite(y_min) || !is.finite(y_max)) { y_min <- -1; y_max <- 1 }
  if (y_min == y_max) { y_min <- y_min - 0.1; y_max <- y_max + 0.1 }
  
  png(out_file, width = 1400, height = 800, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(5.0, 5.0, 5.5, 11.0))
  
  matplot(
    x = x_all, y = ymat,
    type = "l", lty = 1, lwd = 2, col = cols,
    xaxt = "n",
    xlim = X_RANGE, ylim = c(y_min, y_max),
    xlab = "Time",
    ylab = "Estimated coefficient",
    main = title
  )
  axis.Date(1, at = YEAR_TICKS, format = "%Y")
  abline(h = 0, col = "gray50", lty = 2, lwd = 1)
  sub_line <- subtitle
  if (is.finite(coef_plot_scale) && coef_plot_scale != 1) {
    sub_line <- paste0(subtitle, sprintf(" (coefficients ×%.0f)", coef_plot_scale))
  }
  mtext(sub_line, side = 3, line = 0.8, cex = 0.85)
  legend(
    "topright",
    inset  = c(-0.34, 0),
    xpd    = NA,
    legend = unname(legend_lab),
    col    = cols,
    lty    = 1, lwd = 2,
    bty    = "n", cex = 0.9
  )
  dev.off()
}

# -----------------------------------------------------------------------------
# Plot 1.  ifo BCI (alone) + PMI MFG (with gas), mirroring the user-PNG.
# -----------------------------------------------------------------------------
sid_ifo_alone <- find_spec_id("",     "ifo_bci_de")
sid_pmi_gas   <- find_spec_id("gas",  "pmi_mfg_de")
if (is.na(sid_ifo_alone) || is.na(sid_pmi_gas)) {
  stop("Could not resolve the mixed-scope specs (ifo alone / pmi+gas).")
}
plot_panel(
  series = list(
    pmi_mfg_de = pull_series(sid_pmi_gas,   "pmi_mfg_de"),
    ifo_bci_de = pull_series(sid_ifo_alone, "ifo_bci_de")
  ),
  subtitle   = "PMI Manufacturing DE, ifo BCI DE",
  legend_lab = c(pmi_mfg_de = "pmi mfg de", ifo_bci_de = "ifo bci de"),
  out_file   = file.path(out_dir, "01_pmi_mfg_gas_plus_ifo_alone.png"),
  coef_plot_scale = 10
)

# -----------------------------------------------------------------------------
# Plots 2-5.  Each OTHER LF variable, alone with gas (single LF + HF=gas).
# -----------------------------------------------------------------------------
other_lf <- c("pmi_serv_de", "Consumer_Goods", "Manufacturing", "Energy")
for (i in seq_along(other_lf)) {
  v <- other_lf[i]
  sid <- find_spec_id("gas", v)
  if (is.na(sid)) {
    warning("Skipping singular plot for ", v, ": no 'gas + ", v, "' spec found.")
    next
  }
  s <- pull_series(sid, v)
  if (is.null(s)) {
    warning("Skipping singular plot for ", v, ": empty coefficient block.")
    next
  }
  plot_panel(
    series = setNames(list(s), v),
    subtitle = sprintf(pretty_name(v)),
    legend_lab = setNames(pretty_name(v), v),
    out_file = file.path(out_dir, sprintf("%02d_%s_with_gas.png", 1 + i, sanitize_name(v)))
  )
}

# -----------------------------------------------------------------------------
# Plots 6-8.  Group plots — all LF coefs from ONE spec (multi-LF + gas).
# -----------------------------------------------------------------------------
group_plot <- function(lf_vars_vec, out_basename, plot_idx) {
  # Dictionary canonical order:
  canonical <- c("pmi_mfg_de", "pmi_serv_de", "ifo_bci_de",
                 "Consumer_Goods", "Manufacturing", "Energy")
  lf_dict <- paste(intersect(canonical, lf_vars_vec), collapse = "|")
  sid <- find_spec_id("gas", lf_dict)
  if (is.na(sid)) {
    warning("Skipping group plot ", out_basename,
            ": no 'gas + ", lf_dict, "' spec found.")
    return(invisible(NULL))
  }
  series <- setNames(
    lapply(lf_vars_vec, function(v) pull_series(sid, v)),
    lf_vars_vec
  )
  ok <- !vapply(series, is.null, logical(1))
  if (!any(ok)) {
    warning("Skipping group plot ", out_basename, ": no series resolved.")
    return(invisible(NULL))
  }
  series <- series[ok]
  plot_panel(
    series     = series,
    subtitle   = sprintf("%s all LFs in same spec)",
                         paste(vapply(names(series), pretty_name, character(1)),
                               collapse = ", ")),
    legend_lab = setNames(vapply(names(series), pretty_name, character(1)),
                          names(series)),
    out_file   = file.path(out_dir, sprintf("%02d_%s.png", plot_idx, out_basename))
  )
}

group_plot(
  lf_vars_vec  = c("pmi_mfg_de", "ifo_bci_de", "Manufacturing"),
  out_basename = "group_pmi_ifo_mfg_with_gas",
  plot_idx     = 6
)
group_plot(
  lf_vars_vec  = c("pmi_mfg_de", "ifo_bci_de", "Manufacturing", "Consumer_Goods"),
  out_basename = "group_pmi_ifo_mfg_cgs_with_gas",
  plot_idx     = 7
)
group_plot(
  lf_vars_vec  = c("pmi_mfg_de", "ifo_bci_de", "Manufacturing", "Consumer_Goods", "Energy"),
  out_basename = "group_pmi_ifo_mfg_cgs_energy_with_gas",
  plot_idx     = 8
)

cat("\nSaved plots in: ", out_dir, "\n", sep = "")