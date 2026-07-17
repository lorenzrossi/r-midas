#!/usr/bin/env Rscript

# Average LF coefficients by day-of-month position:
# for day d, mean coef across all calendar months where that day exists
# (e.g. average over ~130 values for every "1st of month" in the sample).
#
# Uses standard R-MIDAS coefficients: r_midas_de.csv
#
# Series:
#   - ifo_bci_de     : spec_008 (LF alone)
#   - pmi_mfg_de     : spec_160 (gas + pmi_mfg_de only)
#   - Manufacturing  : spec_130 (gas + Manufacturing only)
#   - Consumer_Goods : spec_132 (gas + Consumer_Goods only)
#   - Energy         : spec_129 (gas + Energy only)

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
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(script_root, "results", "plots", "lf_monthly_avg_de")
)
max_dom <- as.integer(Sys.getenv("MAX_DAY_OF_MONTH", unset = "28"))

if (!file.exists(coef_file)) stop("Missing coefficient file: ", coef_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coefs <- read.csv(coef_file, stringsAsFactors = FALSE, check.names = FALSE)
coefs$origin_date <- as.Date(coefs$origin_date)

series_defs <- list(
  list(spec_id = "spec_008", col = "ifo_bci_de",     label = "ifo bci de (LF alone)"),
  list(spec_id = "spec_160", col = "pmi_mfg_de",     label = "pmi mfg de (gas)"),
  list(spec_id = "spec_130", col = "Manufacturing",  label = "IPI manufacturing (gas)"),
  list(spec_id = "spec_132", col = "Consumer_Goods", label = "IPI consumer goods (gas)"),
  list(spec_id = "spec_129", col = "Energy",         label = "IPI energy (gas)")
)

pretty_var <- function(x) gsub("_", " ", x, fixed = TRUE)

dom_mean <- function(x) {
  ok <- is.finite(x)
  if (!any(ok)) return(NA_real_)
  mean(x[ok])
}

dom_rows <- list()
map_rows <- list()

for (i in seq_along(series_defs)) {
  sd <- series_defs[[i]]
  if (!sd$col %in% names(coefs)) stop("Missing column in coef file: ", sd$col)

  block <- coefs[coefs$spec_id == sd$spec_id, c("origin_date", sd$col), drop = FALSE]
  if (nrow(block) == 0L) {
    warning("No rows for ", sd$spec_id)
    next
  }
  block <- block[is.finite(block[[sd$col]]), , drop = FALSE]
  block$day_of_month <- as.integer(format(block$origin_date, "%d"))
  block <- block[block$day_of_month <= max_dom, , drop = FALSE]

  by_dom <- split(block[[sd$col]], block$day_of_month)
  agg <- data.frame(
    day_of_month = as.integer(names(by_dom)),
    coef_dom_avg = vapply(by_dom, dom_mean, numeric(1)),
    n_obs = vapply(by_dom, function(x) sum(is.finite(x)), integer(1)),
    stringsAsFactors = FALSE
  )
  agg$series <- sd$label
  agg$spec_id <- sd$spec_id
  agg$variable <- sd$col

  dom_rows[[i]] <- agg
  map_rows[[i]] <- data.frame(
    spec_id = sd$spec_id,
    variable = sd$col,
    label = sd$label,
    stringsAsFactors = FALSE
  )
}

if (length(dom_rows) == 0L) stop("No day-of-month series built.")

long_df <- do.call(rbind, dom_rows)
map_df <- do.call(rbind, map_rows)

wide <- long_df[, c("day_of_month", "variable", "coef_dom_avg")]
wide <- reshape(
  wide,
  direction = "wide",
  idvar = "day_of_month",
  timevar = "variable",
  v.names = "coef_dom_avg"
)
coef_cols <- grep("^coef_dom_avg", names(wide), value = TRUE)
wide <- wide[order(wide$day_of_month), , drop = FALSE]

palette_cols <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e")

png_out <- file.path(out_dir, "lf_dom_avg_five_indicators.png")
png(png_out, width = 1400, height = 800, res = 150)
op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)
par(mar = c(5, 5, 5, 2))

ymat <- as.matrix(wide[, coef_cols, drop = FALSE])
mode(ymat) <- "numeric"
leg_labels <- vapply(map_df$variable, pretty_var, character(1))

matplot(
  x = wide$day_of_month,
  y = ymat,
  type = "l",
  lty = 1,
  lwd = 2,
  col = palette_cols[seq_len(ncol(ymat))],
  xlab = sprintf("Day of month (1-%d)", max_dom),
  ylab = sprintf("Average coefficient (across months, day position)"),
  main = "R-MIDAS LF coefficients (average by day-of-month position)",
  xaxt = "n",
  xlim = c(1, max_dom)
)
axis(1, at = seq(1, max_dom, by = 1))
abline(h = 0, col = "gray50", lty = 2, lwd = 1)
legend(
  "topright",
  legend = leg_labels,
  col = palette_cols[seq_len(ncol(ymat))],
  lty = 1,
  lwd = 2,
  bty = "n",
  cex = 0.85
)
dev.off()

csv_out <- file.path(out_dir, "lf_dom_avg_five_indicators.csv")
write.csv(long_df, csv_out, row.names = FALSE)
map_path <- file.path(out_dir, "lf_dom_avg_series_map.csv")
write.csv(map_df, map_path, row.names = FALSE)

cat("Saved plot: ", png_out, "\n", sep = "")
cat("Saved data: ", csv_out, "\n", sep = "")
cat("Saved map: ", map_path, "\n", sep = "")
