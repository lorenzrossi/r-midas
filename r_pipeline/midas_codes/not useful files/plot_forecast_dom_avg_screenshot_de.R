#!/usr/bin/env Rscript

# Average Germany R-MIDAS screenshot forecasts by day-of-month position
# (excludes AR benchmark models). One PNG per spec × horizon; legend outside.

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(
    dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
    "v2_paths.R"
  )
} else file.path(getwd(), "v2_paths.R")
source(v2.boot.path)

script_root <- script_dir
forecast_dir <- file.path(script_root, "forecasts")
out_dir <- file.path(script_root, "results", "plots", "screenshot_models_de")
plot_dir <- file.path(out_dir, "forecast_dom_avg_by_model")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

HORIZONS <- c(1, 7, 14, 21, 28)
max_dom <- as.integer(Sys.getenv("MAX_DAY_OF_MONTH", unset = "28"))

model_defs <- list(
  list(
    label = "R-MIDAS + gas + brent + PMI + ifo + IPI (all)",
    family = "R_MIDAS",
    spec_id = "spec_239",
    plot_col = "#1f77b4"
  ),
  list(
    label = "R-MIDAS + gas + PMI + ifo + IPI (C.Goods, Mfg)",
    family = "R_MIDAS",
    spec_id = "spec_174",
    plot_col = "#ff7f0e"
  ),
  list(
    label = "R-MIDAS + gas + PMI + ifo + IPI (Mfg)",
    family = "R_MIDAS",
    spec_id = "spec_170",
    plot_col = "#2ca02c"
  ),
  list(
    label = "R-MIDAS + gas + PMI + ifo",
    family = "R_MIDAS",
    spec_id = "spec_168",
    plot_col = "#d62728"
  )
)

dom_mean <- function(x) {
  ok <- is.finite(x)
  if (!any(ok)) return(NA_real_)
  mean(x[ok])
}

dom_aggregate <- function(dates, values, max_day) {
  df <- data.frame(
    target_date = as.Date(dates),
    value = as.numeric(values),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$value), , drop = FALSE]
  df$day_of_month <- as.integer(format(df$target_date, "%d"))
  df <- df[df$day_of_month <= max_day, , drop = FALSE]
  by_dom <- split(df$value, df$day_of_month)
  data.frame(
    day_of_month = as.integer(names(by_dom)),
    dom_avg = vapply(by_dom, dom_mean, numeric(1)),
    n_obs = vapply(by_dom, function(x) sum(is.finite(x)), integer(1)),
    stringsAsFactors = FALSE
  )
}

load_forecasts <- function(h) {
  ar_path <- file.path(forecast_dir, sprintf("ar_de_h%02d.csv", h))
  rm_path <- file.path(forecast_dir, sprintf("r_midas_de_h%02d.csv", h))
  if (!file.exists(ar_path) || !file.exists(rm_path)) {
    return(NULL)
  }
  ar <- read.csv(ar_path, stringsAsFactors = FALSE)
  rm <- read.csv(rm_path, stringsAsFactors = FALSE)
  ar$target_date <- as.Date(ar$target_date)
  rm$target_date <- as.Date(rm$target_date)
  base <- ar[order(ar$target_date), c("target_date", "y_actual"), drop = FALSE]

  forecasts <- list()
  for (m in model_defs) {
    if (!m$spec_id %in% names(rm)) next
    v <- rm[[m$spec_id]]
    names(v) <- as.character(rm$target_date)
    forecasts[[m$spec_id]] <- v[as.character(base$target_date)]
  }
  list(
    dates = base$target_date,
    actual = base$y_actual,
    forecasts = forecasts
  )
}

plot_dom_pair <- function(dom_actual, dom_fc, m, h, out_path) {
  x <- dom_actual$day_of_month
  y_act <- dom_actual$dom_avg
  y_fc <- dom_fc$dom_avg[match(x, dom_fc$day_of_month)]
  y_rng <- range(c(y_act, y_fc), na.rm = TRUE)

  png(out_path, width = 1500, height = 750, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mar = c(5, 5, 4, 14), xpd = FALSE)

  plot(
    x, y_act,
    type = "l", lwd = 2.2, col = "black",
    xlab = sprintf("Day of month (1-%d)", max_dom),
    ylab = "Average EUR/MWh (across months at day position)",
    main = sprintf("%s — h = %d", m$label, h),
    xlim = c(1, max_dom),
    ylim = y_rng,
    xaxt = "n"
  )
  axis(1, at = seq(1, max_dom, by = 1))
  lines(x, y_fc, lwd = 2.2, col = m$plot_col)

  par(xpd = TRUE)
  legend(
    x = max_dom + 0.55,
    y = mean(y_rng),
    legend = c("Actual", "Forecast"),
    col = c("black", m$plot_col),
    lty = 1,
    lwd = 2.2,
    bty = "n",
    xjust = 0,
    yjust = 0.5,
    cex = 0.9
  )
}

all_long <- list()

for (h in HORIZONS) {
  loaded <- load_forecasts(h)
  if (is.null(loaded)) {
    warning("Missing forecast files for h=", h)
    next
  }

  dom_actual <- dom_aggregate(loaded$dates, loaded$actual, max_dom)
  dom_actual$series <- "Actual"
  dom_actual$spec_id <- NA_character_
  dom_actual$horizon <- h

  rows_h <- list(dom_actual)

  for (m in model_defs) {
    if (!m$spec_id %in% names(loaded$forecasts)) next
    dom_fc <- dom_aggregate(loaded$dates, loaded$forecasts[[m$spec_id]], max_dom)
    dom_fc$series <- m$label
    dom_fc$spec_id <- m$spec_id
    dom_fc$horizon <- h
    rows_h[[length(rows_h) + 1L]] <- dom_fc

    plot_dom_pair(
      dom_actual,
      dom_fc,
      m,
      h,
      file.path(plot_dir, sprintf("forecast_dom_avg_%s_h%02d.png", m$spec_id, h))
    )
  }

  long_h <- do.call(rbind, rows_h)
  long_h <- long_h[order(long_h$series, long_h$day_of_month), , drop = FALSE]
  all_long[[as.character(h)]] <- long_h

  write.csv(
    long_h,
    file.path(out_dir, sprintf("forecast_dom_avg_rmidas_de_h%02d.csv", h)),
    row.names = FALSE
  )
}

if (length(all_long) == 0L) stop("No horizons processed.")

combined <- do.call(rbind, all_long)
write.csv(
  combined,
  file.path(out_dir, "forecast_dom_avg_rmidas_de_all_horizons.csv"),
  row.names = FALSE
)

old_combined <- file.path(out_dir, "forecast_dom_avg_de_h01.png")
if (file.exists(old_combined)) {
  for (pat in c("forecast_dom_avg_de_h*.png", "forecast_dom_avg_de_h*.csv")) {
    old_files <- Sys.glob(file.path(out_dir, pat))
    invisible(file.remove(old_files))
  }
  invisible(file.remove(file.path(out_dir, "forecast_dom_avg_de_all_horizons.csv")))
}

cat("Saved per-model plots in: ", plot_dir, "\n", sep = "")
