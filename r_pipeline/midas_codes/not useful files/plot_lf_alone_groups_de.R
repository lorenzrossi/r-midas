#!/usr/bin/env Rscript

# =============================================================================
# Germany LF coefficient plots (legacy coefficients)
#
# What this script creates:
# 1) One plot per LF variable using LF-alone specs only (no HF, single LF).
# 2) Group plots using LF-alone specs only:
#    - pmi_mfg_de + ifo_bci_de
#    - pmi_mfg_de + ifo_bci_de + Manufacturing
#    - pmi_mfg_de + ifo_bci_de + Manufacturing + Consumer_Goods
#    - pmi_mfg_de + ifo_bci_de + Manufacturing + Consumer_Goods + Energy
# 3) One final plot with the same 5 LF variables where each line is the
#    average coefficient across all model specs (by origin date).
#
# Defaults:
#   - results/coefficients/r_midas_de.csv
#   - forecasts/r_midas_spec_dictionary_de.csv
#   - results/plots/lf_alone_groups_de/
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
  unset = file.path(script_root, "results", "plots", "lf_alone_groups_de")
)

if (!file.exists(coef_file)) stop("Missing coefficient file: ", coef_file)
if (!file.exists(dict_file)) stop("Missing dictionary file: ", dict_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coefs <- read.csv(coef_file, stringsAsFactors = FALSE, check.names = FALSE)
specs <- read.csv(dict_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("origin_date", "spec_id") %in% names(coefs))) {
  stop("Coefficients file must contain: origin_date, spec_id")
}
if (!all(c("spec_id", "hf_vars", "lf_vars") %in% names(specs))) {
  stop("Dictionary must contain: spec_id, hf_vars, lf_vars")
}

coefs$origin_date <- as.Date(coefs$origin_date)

split_lf <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  strsplit(x, "\\|", fixed = FALSE)[[1]]
}

pretty_var <- function(x) gsub("_", " ", x, fixed = TRUE)
slug <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

find_spec <- function(lf_set, hf_vars = "") {
  lf_set <- as.character(lf_set)
  lf_set <- lf_set[nzchar(lf_set)]
  idx <- which(specs$hf_vars == hf_vars)
  if (length(idx) == 0L) return(NA_character_)

  ok <- logical(length(idx))
  for (i in seq_along(idx)) {
    this <- split_lf(specs$lf_vars[idx[i]])
    ok[i] <- length(this) == length(lf_set) && setequal(this, lf_set)
  }
  hits <- specs$spec_id[idx[ok]]
  if (length(hits) == 0L) return(NA_character_)
  sort(hits)[1]
}

plot_block <- function(df, vars, out_file, main_title) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) == 0L) return(FALSE)

  d <- df[, c("origin_date", vars), drop = FALSE]
  d <- d[order(d$origin_date), , drop = FALSE]

  ymat <- as.matrix(d[, vars, drop = FALSE])
  mode(ymat) <- "numeric"

  y_min <- suppressWarnings(min(ymat, na.rm = TRUE))
  y_max <- suppressWarnings(max(ymat, na.rm = TRUE))
  if (!is.finite(y_min) || !is.finite(y_max)) {
    y_min <- -1
    y_max <- 1
  }
  if (y_min == y_max) {
    y_min <- y_min - 0.1
    y_max <- y_max + 0.1
  }

  palette_cols <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf")
  cols <- palette_cols[seq_len(ncol(ymat))]
  x_ticks <- as.Date(sprintf("%d-01-01", 2015:2025))

  png(out_file, width = 1400, height = 800, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(5, 5, 5, 2))

  matplot(
    x = d$origin_date,
    y = ymat,
    type = "l",
    lty = 1,
    lwd = 2,
    col = cols,
    xaxt = "n",
    xlim = as.Date(c("2015-01-01", "2025-12-31")),
    xlab = "Forecast origin date",
    ylab = "Estimated coefficient",
    main = main_title
  )
  axis.Date(1, at = x_ticks, format = "%Y")
  abline(h = 0, col = "gray50", lty = 2, lwd = 1)
  legend(
    "topright",
    legend = vapply(vars, pretty_var, character(1)),
    col = cols,
    lty = 1,
    lwd = 2,
    bty = "n"
  )
  dev.off()
  TRUE
}

map_rows <- list()
row_i <- 1L

# ---------------------------------------------------------------------------
# 1) Single-LF plots (LF alone only)
# ---------------------------------------------------------------------------
single_specs <- specs[
  specs$hf_vars == "" &
    nzchar(specs$lf_vars) &
    !grepl("\\|", specs$lf_vars),
  c("spec_id", "lf_vars"),
  drop = FALSE
]
single_vars <- sort(unique(single_specs$lf_vars))
single_vars <- single_vars[single_vars %in% names(coefs)]

for (v in single_vars) {
  sid <- find_spec(v, hf_vars = "")
  if (is.na(sid)) next
  block <- coefs[coefs$spec_id == sid, c("origin_date", v), drop = FALSE]
  if (nrow(block) == 0L) next
  out_file <- file.path(out_dir, paste0("single_", slug(v), "_alone.png"))
  ok <- plot_block(
    df = block,
    vars = v,
    out_file = out_file,
    main_title = sprintf("LF coefficient alone: %s", pretty_var(v))
  )
  if (!ok) next

  map_rows[[row_i]] <- data.frame(
    plot_type = "single_alone",
    spec_id = sid,
    lf_vars = v,
    output_file = out_file,
    stringsAsFactors = FALSE
  )
  row_i <- row_i + 1L
}

# ---------------------------------------------------------------------------
# 1b) Combined plot: ifo_bci_de + pmi_mfg_de (each from its own LF-alone spec)
# ---------------------------------------------------------------------------
combo_vars <- c("ifo_bci_de", "pmi_mfg_de")
combo_vars <- combo_vars[combo_vars %in% names(coefs)]
if (length(combo_vars) == 2L) {
  sid_ifo <- find_spec("ifo_bci_de", hf_vars = "")
  sid_pmi <- find_spec("pmi_mfg_de", hf_vars = "")
  if (!is.na(sid_ifo) && !is.na(sid_pmi)) {
    b_ifo <- coefs[coefs$spec_id == sid_ifo, c("origin_date", "ifo_bci_de"), drop = FALSE]
    b_pmi <- coefs[coefs$spec_id == sid_pmi, c("origin_date", "pmi_mfg_de"), drop = FALSE]
    combo_df <- merge(b_ifo, b_pmi, by = "origin_date", all = TRUE, sort = TRUE)

    combo_out <- file.path(out_dir, "single_ifo_bci_de_and_pmi_mfg_de_alone.png")
    ok <- plot_block(
      df = combo_df,
      vars = c("ifo_bci_de", "pmi_mfg_de"),
      out_file = combo_out,
      main_title = "LF coefficients alone: ifo bci de and pmi mfg de"
    )
    if (ok) {
      map_rows[[row_i]] <- data.frame(
        plot_type = "single_pair_alone",
        spec_id = paste(sid_ifo, sid_pmi, sep = "|"),
        lf_vars = "ifo_bci_de|pmi_mfg_de",
        output_file = combo_out,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
}

# ---------------------------------------------------------------------------
# 2) Group plots from single-LF-alone specs (hf_vars == "")
# ---------------------------------------------------------------------------
build_group_from_single_alone <- function(vars) {
  vars <- vars[vars %in% names(coefs)]
  if (length(vars) == 0L) return(NULL)

  merged <- data.frame(origin_date = sort(unique(coefs$origin_date)))
  sid_map <- character(0)
  kept_vars <- character(0)

  for (v in vars) {
    sid <- find_spec(v, hf_vars = "")
    if (is.na(sid)) next
    block <- coefs[coefs$spec_id == sid, c("origin_date", v), drop = FALSE]
    if (nrow(block) == 0L) next
    merged <- merge(merged, block, by = "origin_date", all.x = TRUE, sort = TRUE)
    sid_map <- c(sid_map, sid)
    kept_vars <- c(kept_vars, v)
  }

  if (length(kept_vars) == 0L) return(NULL)
  list(df = merged, vars = kept_vars, spec_ids = sid_map)
}

group_defs <- list(
  c("pmi_mfg_de", "ifo_bci_de"),
  c("pmi_mfg_de", "ifo_bci_de", "Manufacturing"),
  c("pmi_mfg_de", "ifo_bci_de", "Manufacturing", "Consumer_Goods"),
  c("pmi_mfg_de", "ifo_bci_de", "Manufacturing", "Consumer_Goods", "Energy")
)

for (g in group_defs) {
  grp <- build_group_from_single_alone(g)
  if (is.null(grp)) {
    warning("No single-LF-alone data found for group: ", paste(g, collapse = "|"))
    next
  }

  vars <- grp$vars
  block <- grp$df
  out_file <- file.path(out_dir, paste0("group_", length(vars), "_", slug(paste(vars, collapse = "_")), "_alone.png"))
  ok <- plot_block(
    df = block,
    vars = vars,
    out_file = out_file,
    main_title = sprintf("LF coefficients alone (from single-LF specs, %d vars)", length(vars))
  )
  if (!ok) next

  map_rows[[row_i]] <- data.frame(
    plot_type = "group_alone",
    spec_id = paste(paste(vars, grp$spec_ids, sep = "="), collapse = "|"),
    lf_vars = paste(vars, collapse = "|"),
    output_file = out_file,
    stringsAsFactors = FALSE
  )
  row_i <- row_i + 1L
}

# ---------------------------------------------------------------------------
# 3) Average across all specs for the 5 requested LF variables
# ---------------------------------------------------------------------------
avg_vars <- c("pmi_mfg_de", "ifo_bci_de", "Manufacturing", "Consumer_Goods", "Energy")
avg_vars <- avg_vars[avg_vars %in% names(coefs)]
if (length(avg_vars) > 0L) {
  mean_finite <- function(x) {
    z <- x[is.finite(x)]
    if (length(z) == 0L) NA_real_ else mean(z)
  }

  avg_df <- data.frame(origin_date = sort(unique(coefs$origin_date)))
  for (v in avg_vars) {
    tmp <- aggregate(coefs[[v]], by = list(origin_date = coefs$origin_date), FUN = mean_finite)
    names(tmp)[2] <- v
    avg_df <- merge(avg_df, tmp, by = "origin_date", all.x = TRUE, sort = TRUE)
  }

  avg_out <- file.path(out_dir, "group_5_average_across_all_specs.png")
  plot_block(
    df = avg_df,
    vars = avg_vars,
    out_file = avg_out,
    main_title = "LF coefficients average across all specs"
  )

  map_rows[[row_i]] <- data.frame(
    plot_type = "average_all_specs",
    spec_id = NA_character_,
    lf_vars = paste(avg_vars, collapse = "|"),
    output_file = avg_out,
    stringsAsFactors = FALSE
  )
}

if (length(map_rows) == 0L) {
  stop("No plots produced. Check inputs and dictionary coverage.")
}

map_df <- do.call(rbind, map_rows)
map_path <- file.path(out_dir, "lf_alone_groups_map_de.csv")
write.csv(map_df, map_path, row.names = FALSE)

cat("Saved plots in: ", out_dir, "\n", sep = "")
cat("Saved map: ", map_path, "\n", sep = "")
