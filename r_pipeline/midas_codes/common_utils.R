# =============================================================================
# v2 PIPELINE — COMMON UTILITIES
#
# Loads daily prices (y = elec, X_HF = gas, brent) and monthly LF blocks for
# Germany or Italy, applies the agreed transformations, builds publication-
# aware daily LF series via LOCF on the (pub_lag-shifted) monthly index, and
# exposes helpers used by ar_benchmark.R, ru_midas.R (Cimadomo principal),
# r_midas.R (Almon-on-y-only), r_midas_extended.R (Almon on y AND HF) and
# r_midas_legacy.R (Almon-on-y-only, legacy-style NLS).
#
# Transformations (Weron 2014, Marcellino-Stock-Watson 2006 style):
#   elec_price   : asinh   (handles negative DE prices, near-linear at low |y|)
#   gas_price    : log
#   brent_price  : log
#   IPI series   : Δlog    (Consumer_Goods, Manufacturing, Energy)
#   PMI/IFO/ISTAT: scale() (full-sample z-score; safe — OLS forecasts are
#                           invariant to affine rescaling of regressors)
#
# Forecast scoring is done on the elec_price scale by back-transforming with
# sinh(.). RMSE / DM-HLN are then computed on euros/MWh.
#
# Calendar controls (built by build_calendar_dummies(dates, country)):
#   6 weekday dummies (wd_mon..wd_sat, Sunday = reference) — standard EPF
#                                                            controls for the
#                                                            weekday cycle in
#                                                            electricity demand
#  11 month-of-year dummies (mon02..mon12, January = reference) — long-run
#                                                                  seasonality
#   1 holiday dummy (national holidays per country, fixed + Easter-derived
#                    movable feasts)
#
# Publication lags:
#   IPI                       : pub_lag = 2 months (Eurostat ~T-2 release)
#   PMI / IFO BCI / ISTAT BCI : pub_lag = 1 month  (released end of ref. month)
#
# AR structure used everywhere (Weron 2014, JoF):
#   AR lags = {1, 2, 3, 7}  — free coefficients, no Almon restriction.
#
# NLS for R-MIDAS:
#   Concentrated NLS with multi-start L-BFGS-B (Ghysels-Wright 2009 protocol).
#   exp-Almon weights computed with log-sum-exp normalisation for stability.
# =============================================================================

suppressWarnings(suppressPackageStartupMessages(library(stats)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))

# -----------------------------------------------------------------------------
# 0. PARALLEL HELPER
# -----------------------------------------------------------------------------
get_n_cores <- function() {
  raw <- Sys.getenv("N_CORES", unset = "")
  if (nchar(raw) > 0L) {
    n <- suppressWarnings(as.integer(raw))
    if (!is.na(n) && n >= 1L) return(n)
  }
  max(1L, parallel::detectCores(logical = TRUE) - 1L)
}

parallel_lapply <- function(X, FUN, ..., n_cores = NULL,
                            mc_preschedule = TRUE) {
  if (is.null(n_cores)) n_cores <- get_n_cores()
  if (n_cores <= 1L || .Platform$OS.type == "windows") {
    return(lapply(X, FUN, ...))
  }
  parallel::mclapply(X, FUN, ..., mc.cores = n_cores,
                     mc.preschedule = mc_preschedule,
                     mc.allow.recursive = FALSE)
}


# -----------------------------------------------------------------------------
# 1. CONSTANTS
# -----------------------------------------------------------------------------
AR_LAGS        <- c(1L, 2L, 3L, 7L)
# Rolling-window length in days.  Overridable via the WINDOW_DAYS env var so
# that alternative window lengths (e.g. 730 for faster regime adaptation) can
# be run with the SAME methodology for both countries.
WINDOW_DAYS    <- local({
  raw <- suppressWarnings(as.integer(Sys.getenv("WINDOW_DAYS", "")))
  if (!is.na(raw) && raw >= 365L) raw else 3650L
})
HORIZONS       <- c(1L, 7L, 14L, 21L, 28L)
H_MAX          <- max(HORIZONS)
Q_HF           <- 14L
PUB_LAG_IPI    <- 2L
PUB_LAG_SURVEY <- 1L

HF_CANDIDATES  <- c("gas", "brent")
LF_CANDIDATES_DE <- c("pmi_mfg_de", "pmi_serv_de", "ifo_bci_de",
                      "Consumer_Goods", "Manufacturing", "Energy")
LF_CANDIDATES_IT <- c("pmi_mfg_it", "pmi_serv_it", "istat_bci_it",
                      "Consumer_Goods", "Manufacturing", "Energy")

# -----------------------------------------------------------------------------
# 2. DATA LOADING + TRANSFORMATION
# -----------------------------------------------------------------------------
load_country_data <- function(country, data_dir) {
  stopifnot(country %in% c("Germany", "Italy"))
  cc <- if (country == "Germany") "de" else "it"
  daily_path   <- file.path(data_dir, paste0("prices_", cc, "_daily.csv"))
  monthly_path <- file.path(data_dir, paste0("dataset_", cc, "_m.csv"))
  if (!file.exists(daily_path))   stop("Daily file not found: ", daily_path)
  if (!file.exists(monthly_path)) stop("Monthly file not found: ", monthly_path)

  daily   <- read.csv(daily_path,   stringsAsFactors = FALSE)
  monthly <- read.csv(monthly_path, stringsAsFactors = FALSE)
  daily$date   <- as.Date(daily$date)
  monthly$date <- as.Date(monthly$date)
  daily   <- daily[order(daily$date), , drop = FALSE]
  monthly <- monthly[order(monthly$date), , drop = FALSE]

  daily$y_asinh   <- asinh(daily$elec_price)
  daily$gas_log   <- log(daily$gas_price)
  daily$brent_log <- log(daily$brent_price)

  names(monthly)[names(monthly) == "Electricity_Gas_Steam"] <- "Energy"
  for (col in c("Consumer_Goods", "Manufacturing", "Energy")) {
    if (!col %in% names(monthly)) stop("Missing monthly column: ", col)
    monthly[[col]] <- c(0, diff(log(monthly[[col]])))
  }
  survey_cols <- if (country == "Germany")
    c("pmi_mfg_de", "pmi_serv_de", "ifo_bci_de") else
    c("pmi_mfg_it", "pmi_serv_it", "istat_bci_it")
  for (col in survey_cols) {
    if (!col %in% names(monthly)) stop("Missing monthly column: ", col)
    monthly[[col]] <- as.numeric(scale(monthly[[col]]))
  }
  monthly$ref_month <- as.Date(format(monthly$date, "%Y-%m-01"))
  list(daily = daily, monthly = monthly, country = country)
}

# -----------------------------------------------------------------------------
# 3. PUBLICATION-AWARE DAILY LF SERIES
# -----------------------------------------------------------------------------
# The preferred mode maps each monthly observation to the date on which it
# would normally have become public, then carries the latest released value
# forward day by day.  Exact historical release dates take precedence when a
# column named <variable>_release_date is present in monthly_df.
#
# Default conventional schedules (all overridable through environment vars):
#   Manufacturing PMI : first working day of the following month
#   Services PMI      : third working day of the following month
#   ifo BCI           : 25th calendar day of the reference month
#   ISTAT BCI         : 27th calendar day of the reference month
#   IPI components    : 45 calendar days after reference-month end
# Weekends are rolled forward to Monday.  These rules avoid anticipatory use
# of a value but remain a pseudo-real-time approximation when vintage-specific
# release dates are unavailable.
#
# Set LF_RELEASE_MODE=calendar_lag to reproduce the former whole-month shift.

roll_to_weekday <- function(d, direction = c("forward", "backward")) {
  direction <- match.arg(direction)
  d <- as.Date(d)
  wd <- as.POSIXlt(d)$wday  # Sunday=0, Saturday=6
  if (direction == "forward") {
    d[wd == 6L] <- d[wd == 6L] + 2L
    d[wd == 0L] <- d[wd == 0L] + 1L
  } else {
    d[wd == 6L] <- d[wd == 6L] - 1L
    d[wd == 0L] <- d[wd == 0L] - 2L
  }
  d
}

nth_workday_of_month <- function(month_first, n) {
  month_first <- as.Date(month_first)
  stopifnot(length(n) == 1L, n >= 1L)
  out <- rep(as.Date(NA), length(month_first))
  for (i in seq_along(month_first)) {
    cand <- seq(month_first[i], shift_month(month_first[i], 1L) - 1L, by = "day")
    cand <- cand[!as.POSIXlt(cand)$wday %in% c(0L, 6L)]
    if (length(cand) >= n) out[i] <- cand[n]
  }
  out
}

usual_lf_release_date <- function(ref_month, lf_col) {
  ref_month <- as.Date(format(as.Date(ref_month), "%Y-%m-01"))
  next_month <- shift_month(ref_month, 1L)
  month_end  <- next_month - 1L

  env_int <- function(name, default) {
    z <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
    if (is.na(z)) default else z
  }

  if (grepl("^pmi_mfg_", lf_col)) {
    n <- env_int("LF_PMI_MFG_WORKDAY", 1L)
    return(nth_workday_of_month(next_month, n))
  }
  if (grepl("^pmi_serv_", lf_col)) {
    n <- env_int("LF_PMI_SERV_WORKDAY", 3L)
    return(nth_workday_of_month(next_month, n))
  }
  if (identical(lf_col, "ifo_bci_de")) {
    day <- env_int("LF_IFO_DAY", 25L)
    return(roll_to_weekday(ref_month + (day - 1L), "forward"))
  }
  if (identical(lf_col, "istat_bci_it")) {
    day <- env_int("LF_ISTAT_BCI_DAY", 27L)
    return(roll_to_weekday(ref_month + (day - 1L), "forward"))
  }
  if (lf_col %in% c("Consumer_Goods", "Manufacturing", "Energy")) {
    delay <- env_int("LF_IPI_DELAY_DAYS", 45L)
    return(roll_to_weekday(month_end + delay, "forward"))
  }
  stop("No usual release-date rule defined for LF variable: ", lf_col)
}

make_daily_lf <- function(daily_dates, monthly_df, lf_col, pub_lag,
                          locf = FALSE) {
  stopifnot(lf_col %in% names(monthly_df))
  mode <- tolower(Sys.getenv("LF_RELEASE_MODE", unset = "usual_release"))

  # Backward-compatible calendar-month mapping.
  if (mode %in% c("calendar_lag", "legacy")) {
    ref_first <- monthly_df$ref_month
    vals <- monthly_df[[lf_col]]
    if (locf) {
      last_val <- NA_real_
      for (ii in seq_along(vals)) {
        if (is.na(vals[ii])) {
          if (!is.na(last_val)) vals[ii] <- last_val
        } else last_val <- vals[ii]
      }
    }
    d_month_first <- as.Date(format(daily_dates, "%Y-%m-01"))
    target_month  <- shift_month(d_month_first, -pub_lag)
    idx <- match(target_month, ref_first)
    return(vals[idx])
  }

  vals <- monthly_df[[lf_col]]
  exact_col <- paste0(lf_col, "_release_date")
  release_dates <- if (exact_col %in% names(monthly_df)) {
    as.Date(monthly_df[[exact_col]])
  } else {
    usual_lf_release_date(monthly_df$ref_month, lf_col)
  }

  releases <- data.frame(
    release_date = as.Date(release_dates),
    value = as.numeric(vals),
    row_id = seq_along(vals)
  )
  releases <- releases[!is.na(releases$release_date) & is.finite(releases$value), , drop = FALSE]
  releases <- releases[order(releases$release_date, releases$row_id), , drop = FALSE]
  if (!nrow(releases)) return(rep(NA_real_, length(daily_dates)))

  # If two observations share a release date, retain the latest row.
  releases <- releases[!duplicated(releases$release_date, fromLast = TRUE), , drop = FALSE]
  idx <- findInterval(as.numeric(as.Date(daily_dates)),
                      as.numeric(releases$release_date))
  out <- rep(NA_real_, length(daily_dates))
  available <- idx > 0L
  out[available] <- releases$value[idx[available]]
  out
}

# Vectorised row filter: TRUE where every entry of the row is finite.
# Replaces apply(M, 1, function(r) all(is.finite(r))), which is slow on
# (3650 x 14) blocks evaluated at every rolling origin.
rows_all_finite <- function(M) rowSums(!is.finite(M)) == 0L

# One-step Huber reweighting (used when RMIDAS_HUBER=1 in the runners).
# Given the residuals of a first-stage fit, returns observation weights
#   w_t = 1                 if |e_t| <= k
#   w_t = k / |e_t|         if |e_t| >  k
# with k = k_mult * MAD(resid)/0.6745 (robust scale; k_mult = 1.345 keeps
# 95% efficiency under Gaussian errors).  The refit with these weights is
# the standard one-step reweighted M-estimator: crisis-day residuals stop
# dominating the coefficient estimates quadratically.  Returns NULL when
# the robust scale is degenerate (all residuals ~ equal).
huber_weights <- function(resid, k_mult = 1.345) {
  s <- stats::median(abs(resid - stats::median(resid))) / 0.6745
  if (!is.finite(s) || s <= 0) return(NULL)
  k <- k_mult * s
  a <- abs(resid)
  w <- rep(1, length(resid))
  big <- a > k
  w[big] <- k / a[big]
  w
}

shift_month <- function(d, k) {
  yrs <- as.integer(format(d, "%Y"))
  mos <- as.integer(format(d, "%m"))
  tot <- yrs * 12L + (mos - 1L) + k
  ny  <- tot %/% 12L
  nm  <- (tot %% 12L) + 1L
  as.Date(sprintf("%04d-%02d-01", ny, nm))
}

# -----------------------------------------------------------------------------
# 4. DESIGN-MATRIX BUILDING BLOCKS
# -----------------------------------------------------------------------------
build_ar_lags <- function(y, lags = AR_LAGS) {
  n <- length(y)
  out <- matrix(NA_real_, n, length(lags))
  for (j in seq_along(lags)) {
    L <- lags[j]
    if (L < n) out[(L + 1L):n, j] <- y[1:(n - L)]
  }
  colnames(out) <- paste0("y_lag", lags)
  out
}

build_hf_exo_lags <- function(x, prefix, Q = Q_HF) {
  n <- length(x)
  out <- matrix(NA_real_, n, Q)
  for (j in 1:Q) {
    out[(j + 1L):n, j] <- x[1:(n - j)]
  }
  colnames(out) <- paste0(prefix, "_lag", 1:Q)
  out
}

# Month-of-year dummies (legacy / back-compat).  11 columns, January = reference.
build_month_dummies <- function(dates) {
  m <- as.integer(format(dates, "%m"))
  M <- matrix(0L, length(m), 11L)
  for (k in 2:12) M[, k - 1L] <- as.integer(m == k)
  colnames(M) <- paste0("mon", sprintf("%02d", 2:12))
  M
}

# -----------------------------------------------------------------------------
# 4b. CALENDAR CONTROLS — weekday + month + holiday
# -----------------------------------------------------------------------------
# Easter Sunday by the Anonymous Gregorian (Butcher/Meeus) algorithm.
# Returns a Date.  Easter-derived movable feasts:
#     Good Friday    = Easter Sunday - 2
#     Easter Monday  = Easter Sunday + 1
#     Ascension Day  = Easter Sunday + 39
#     Whit Monday    = Easter Sunday + 50
easter_sunday <- function(year) {
  a <- year %% 19
  b <- year %/% 100
  c <- year %% 100
  d <- b %/% 4
  e <- b %% 4
  f <- (b + 8) %/% 25
  g <- (b - f + 1) %/% 3
  h <- (19 * a + b - d - g + 15) %% 30
  i <- c %/% 4
  k <- c %% 4
  l <- (32 + 2 * e + 2 * i - h - k) %% 7
  m <- (a + 11 * h + 22 * l) %/% 451
  month <- (h + l - 7 * m + 114) %/% 31
  day   <- ((h + l - 7 * m + 114) %% 31) + 1
  as.Date(sprintf("%04d-%02d-%02d", year, month, day))
}

# National holiday set per (country, year).  Returns a Date vector.
# Italy:  Jan 1, Jan 6 (Epiphany), Apr 25 (Liberation Day), Easter Monday,
#         May 1 (Labour), Jun 2 (Republic), Aug 15 (Ferragosto),
#         Nov 1 (All Saints), Dec 8 (Immaculate Conception), Dec 25, Dec 26.
# Germany: Jan 1, Good Friday, Easter Monday, May 1, Ascension Day, Whit Monday,
#          Oct 3 (German Unity Day), Dec 25, Dec 26.
country_holidays <- function(years, country) {
  stopifnot(country %in% c("Germany", "Italy"))
  dts <- as.Date(character(0))
  for (y in years) {
    e <- easter_sunday(y)
    if (country == "Germany") {
      dts <- c(dts, as.Date(c(
        sprintf("%04d-01-01", y),
        format(e - 2, "%Y-%m-%d"),
        format(e + 1, "%Y-%m-%d"),
        sprintf("%04d-05-01", y),
        format(e + 39, "%Y-%m-%d"),
        format(e + 50, "%Y-%m-%d"),
        sprintf("%04d-10-03", y),
        sprintf("%04d-12-25", y),
        sprintf("%04d-12-26", y)
      )))
    } else {
      dts <- c(dts, as.Date(c(
        sprintf("%04d-01-01", y),
        sprintf("%04d-01-06", y),
        sprintf("%04d-04-25", y),
        format(e + 1, "%Y-%m-%d"),
        sprintf("%04d-05-01", y),
        sprintf("%04d-06-02", y),
        sprintf("%04d-08-15", y),
        sprintf("%04d-11-01", y),
        sprintf("%04d-12-08", y),
        sprintf("%04d-12-25", y),
        sprintf("%04d-12-26", y)
      )))
    }
  }
  sort(unique(dts))
}

# Calendar-dummy design matrix:  6 weekday (Sun = reference)
#                              + 11 month-of-year (Jan = reference)
#                              +  1 holiday  (1 if date in country_holidays)
# Returns an n x 18 integer matrix with named columns.
build_calendar_dummies <- function(dates, country) {
  stopifnot(country %in% c("Germany", "Italy"))
  n <- length(dates)
  # Weekday: format(., "%u") gives 1=Mon..7=Sun (ISO 8601)
  wd <- as.integer(format(dates, "%u"))
  W <- matrix(0L, n, 6L)
  for (k in 1:6) W[, k] <- as.integer(wd == k)
  colnames(W) <- c("wd_mon", "wd_tue", "wd_wed", "wd_thu", "wd_fri", "wd_sat")
  # Month-of-year
  M <- build_month_dummies(dates)
  # Holiday
  yrs   <- unique(as.integer(format(dates, "%Y")))
  hd    <- country_holidays(yrs, country)
  H     <- matrix(as.integer(dates %in% hd), n, 1L,
                  dimnames = list(NULL, "holiday"))
  cbind(W, M, H)
}

# -----------------------------------------------------------------------------
# 5. ROLLING-ORIGIN INDEX GENERATOR
# -----------------------------------------------------------------------------
make_origins <- function(n_rows, window = WINDOW_DAYS,
                         pre_lag_max = max(c(AR_LAGS, Q_HF)),
                         h_max = H_MAX) {
  min_o <- window + pre_lag_max
  max_o <- n_rows - h_max
  if (min_o > max_o) stop("Insufficient data: min_origin=", min_o,
                          " max_origin=", max_o)
  min_o:max_o
}

# -----------------------------------------------------------------------------
# 6. DIEBOLD–MARIANO HARVEY–LEYBOURNE–NEWBOLD TEST
# -----------------------------------------------------------------------------
dm_hln <- function(e1, e2, h, loss = c("sqerr", "abserr")) {
  loss <- match.arg(loss)
  ok <- is.finite(e1) & is.finite(e2)
  e1 <- e1[ok]; e2 <- e2[ok]
  T_ <- length(e1)
  if (T_ < 1L) return(list(stat = NA_real_, pvalue = NA_real_, T = T_))
  d  <- if (loss == "sqerr") e1^2 - e2^2 else abs(e1) - abs(e2)
  dbar <- mean(d)
  g0   <- mean((d - dbar)^2)
  acov_sum <- 0
  if (h > 1L) {
    for (k in 1:(h - 1L)) {
      gk <- sum((d[(k + 1L):T_] - dbar) * (d[1:(T_ - k)] - dbar)) / T_
      acov_sum <- acov_sum + 2 * gk
    }
  }
  v_d <- (g0 + acov_sum) / T_
  if (!is.finite(v_d) || v_d <= 0) return(list(stat = NA_real_, pvalue = NA_real_, T = T_))
  dm <- dbar / sqrt(v_d)
  corr <- sqrt((T_ + 1 - 2 * h + h * (h - 1) / T_) / T_)
  dm_hln <- dm * corr
  pval <- 2 * pt(-abs(dm_hln), df = T_ - 1L)
  list(stat = dm_hln, pvalue = pval, T = T_, dbar = dbar)
}

# -----------------------------------------------------------------------------
# 7. CALENDAR-COLUMN REFRESH AT TARGET ROW (used in iterated multi-step)
# -----------------------------------------------------------------------------
# In iterated multi-step option (c.1), the calendar dummies (weekday, month-
# of-year, holiday) at the target date are different from those at row
# (origin + 1).  This helper updates the appropriate columns of a frozen row
# in place at each iteration step.
#
#   row_lin        : named numeric vector — the linear-design row to update
#   colnames_      : character vector of column names of the design
#   target_date    : Date — the date being forecast at this step
#   holiday_dates  : Date vector — pre-computed national-holiday set; if NULL
#                    the holiday column (if present) is left untouched
#
# Refreshes any of {mon02..mon12, wd_mon..wd_sat, holiday} that are present
# in colnames_, leaving everything else alone.
refresh_calendar_row <- function(row_lin, colnames_, target_date,
                                 holiday_dates = NULL) {
  # Month-of-year dummies
  mon_idx <- grep("^mon[0-9]{2}$", colnames_)
  if (length(mon_idx) > 0L) {
    row_lin[mon_idx] <- 0
    m <- as.integer(format(target_date, "%m"))
    if (m >= 2L) {
      col_name <- sprintf("mon%02d", m)
      ci <- match(col_name, colnames_)
      if (!is.na(ci)) row_lin[ci] <- 1
    }
  }
  # Weekday dummies
  wd_idx <- grep("^wd_(mon|tue|wed|thu|fri|sat)$", colnames_)
  if (length(wd_idx) > 0L) {
    row_lin[wd_idx] <- 0
    wd_iso <- as.integer(format(target_date, "%u"))   # 1=Mon..7=Sun
    if (wd_iso >= 1L && wd_iso <= 6L) {
      wd_name <- c("wd_mon", "wd_tue", "wd_wed", "wd_thu", "wd_fri", "wd_sat")[wd_iso]
      ci <- match(wd_name, colnames_)
      if (!is.na(ci)) row_lin[ci] <- 1
    }
  }
  # Holiday dummy
  hol_idx <- match("holiday", colnames_)
  if (!is.na(hol_idx) && !is.null(holiday_dates)) {
    row_lin[hol_idx] <- as.integer(target_date %in% holiday_dates)
  }
  row_lin
}

# -----------------------------------------------------------------------------
# 8. EXPONENTIAL-ALMON HELPERS FOR R-MIDAS (NLS)
# -----------------------------------------------------------------------------
# Foroni-Marcellino-Schumacher (2015, JoE) Eq. (11), with log-sum-exp
# normalisation for numerical stability at the box boundary.
#
# exp_almon_weights_j takes the ACTUAL lag numbers (j values) entering the
# exponent, so the weights can be computed over a non-consecutive lag set
# such as J = {1, 2, 3, 7} (the set used by the paper's R-MIDAS AR block,
# see model.tex eq. (almon)).  exp_almon_weights keeps the old consecutive
# 1..Q interface for back-compatibility (extended HF blocks, legacy code).
exp_almon_weights_j <- function(jvec, theta) {
  log_w <- theta[1] * jvec + theta[2] * jvec^2
  log_w <- log_w - max(log_w)
  w <- exp(log_w)
  w / sum(w)
}

exp_almon_weights <- function(Q, theta) exp_almon_weights_j(1:Q, theta)

# Concentrated NLS for R-MIDAS, with optional multi-start L-BFGS-B
# (Ghysels-Wright 2009 protocol).  See r_midas.R / r_midas_extended.R for
# the per-spec callers that set n_starts and seed the random draws.
fit_r_midas <- function(y, Z_list, X_lin,
                        theta_init = NULL,
                        theta_lower_block = c(-8, -8),
                        theta_upper_block = c( 8,  0),
                        optim_maxit = 2000L,
                        optim_factr = 5e8,
                        n_starts = 1L) {
  K <- length(Z_list)
  if (K < 1L) stop("Need at least one Z block (the AR block).")
  Q <- ncol(Z_list[[1]])
  n <- length(y)
  if (is.null(X_lin))
    X_lin <- matrix(1, n, 1L, dimnames = list(NULL, "intercept"))

  theta_lower <- rep(theta_lower_block, K)
  theta_upper <- rep(theta_upper_block, K)
  dim_theta   <- 2L * K

  if (n_starts < 1L) n_starts <- 1L
  starts <- matrix(NA_real_, nrow = n_starts, ncol = dim_theta)
  if (!is.null(theta_init) && length(theta_init) == dim_theta &&
      all(is.finite(theta_init))) {
    starts[1, ] <- as.numeric(theta_init)
  } else {
    starts[1, ] <- rep(0, dim_theta)
  }
  if (n_starts > 1L) {
    for (s in 2:n_starts) {
      starts[s, ] <- runif(dim_theta, min = theta_lower, max = theta_upper)
    }
  }

  build_full_design <- function(theta_v) {
    zmat <- matrix(NA_real_, n, K)
    for (k in 1:K) {
      th_k <- theta_v[(2L * k - 1L):(2L * k)]
      w_k  <- exp_almon_weights(Q, th_k)
      zmat[, k] <- as.numeric(Z_list[[k]] %*% w_k)
    }
    cbind(X_lin, zmat)
  }
  obj <- function(theta_v) {
    Xfull <- build_full_design(theta_v)
    keep  <- stats::complete.cases(Xfull) & is.finite(y)
    if (sum(keep) < ncol(Xfull) + 5L) return(1e12)
    qrfit <- qr(Xfull[keep, , drop = FALSE])
    if (qrfit$rank < ncol(Xfull)) return(1e12)
    beta  <- qr.coef(qrfit, y[keep])
    res   <- y[keep] - drop(Xfull[keep, , drop = FALSE] %*% beta)
    sum(res^2)
  }

  best <- NULL
  best_val <- Inf
  for (s in seq_len(n_starts)) {
    cand <- tryCatch(
      optim(par = as.numeric(starts[s, ]), fn = obj,
            method = "L-BFGS-B",
            lower = theta_lower, upper = theta_upper,
            control = list(maxit = optim_maxit, factr = optim_factr)),
      error = function(e) NULL
    )
    if (is.null(cand)) next
    if (!is.finite(cand$value) || cand$value >= 1e11) next
    if (!(cand$convergence %in% c(0L, 1L))) next
    if (cand$value < best_val) {
      best_val <- cand$value
      best     <- cand
    }
  }
  if (is.null(best)) return(NULL)
  theta_hat <- best$par

  Xfull <- build_full_design(theta_hat)
  keep  <- stats::complete.cases(Xfull) & is.finite(y)
  qrfit <- qr(Xfull[keep, , drop = FALSE])
  if (qrfit$rank < ncol(Xfull)) return(NULL)
  beta_full <- qr.coef(qrfit, y[keep])
  res       <- y[keep] - drop(Xfull[keep, , drop = FALSE] %*% beta_full)
  p_lin     <- ncol(X_lin)
  beta_lin  <- beta_full[1:p_lin]
  names(beta_lin) <- colnames(X_lin)
  delta     <- beta_full[(p_lin + 1L):(p_lin + K)]
  names(delta) <- names(Z_list)
  theta_mat <- matrix(theta_hat, nrow = K, ncol = 2L, byrow = TRUE,
                      dimnames = list(names(Z_list), c("theta1", "theta2")))

  list(theta    = theta_mat,
       beta_lin = beta_lin,
       delta    = delta,
       sse      = best$value,
       sigma2   = mean(res^2),
       converged = best$convergence %in% c(0L, 1L),
       n_eff    = sum(keep),
       n_starts = n_starts)
}

# Iterated multi-step forecast for an R-MIDAS fit.  Refreshes calendar columns
# (mon, weekday, holiday) at each target step; HF Almon blocks and LF columns
# are frozen at row (origin + 1).
# lag_sets: optional named list giving, per Almon block, the ACTUAL lag
# numbers of its columns (e.g. list(y_lags = c(1,2,3,7))).  Defaults to
# consecutive 1..ncol per block (the legacy convention).
iterated_forecast_r_midas <- function(fit_obj, Z_list_full, X_lin_row,
                                      X_lin_colnames, history_y_aug,
                                      origin_idx, target_dates,
                                      ar_block_name = "y_lags",
                                      hf_block_names = character(0),
                                      h_max = H_MAX,
                                      holiday_dates = NULL,
                                      lag_sets = NULL) {
  theta_mat <- fit_obj$theta
  delta_v   <- fit_obj$delta
  beta_lin  <- fit_obj$beta_lin
  block_names <- rownames(theta_mat)
  if (is.null(lag_sets)) lag_sets <- list()
  for (nm in block_names) {
    if (is.null(lag_sets[[nm]]))
      lag_sets[[nm]] <- seq_len(ncol(Z_list_full[[nm]]))
  }

  W <- lapply(block_names, function(nm) {
    exp_almon_weights_j(lag_sets[[nm]], theta_mat[nm, ])
  })
  names(W) <- block_names
  ar_lag_set <- lag_sets[[ar_block_name]]

  frozen_z <- numeric(length(hf_block_names))
  names(frozen_z) <- hf_block_names
  for (nm in hf_block_names) {
    frozen_z[nm] <- sum(W[[nm]] * Z_list_full[[nm]][origin_idx + 1L, ])
  }

  y_hat_h <- numeric(h_max)
  for (s in 1:h_max) {
    row_lin <- refresh_calendar_row(X_lin_row, X_lin_colnames,
                                    target_dates[s], holiday_dates)
    lin_part <- sum(beta_lin * row_lin)

    n_hist <- length(history_y_aug)
    ar_idx_in_hist <- n_hist - ar_lag_set + 1L
    if (any(ar_idx_in_hist < 1L)) {
      y_hat_h[s] <- NA_real_
      history_y_aug <- c(history_y_aug, NA_real_)
      next
    }
    z_y <- sum(W[[ar_block_name]] * history_y_aug[ar_idx_in_hist])

    hf_part <- 0
    for (nm in hf_block_names) {
      hf_part <- hf_part + delta_v[nm] * frozen_z[nm]
    }

    yh <- lin_part + delta_v[ar_block_name] * z_y + hf_part
    y_hat_h[s] <- yh
    history_y_aug <- c(history_y_aug, yh)
  }
  y_hat_h
}

# -----------------------------------------------------------------------------
# 9. LINEAR ITERATED MULTI-STEP — kept for AR / RU-MIDAS families
# -----------------------------------------------------------------------------
iterated_forecast <- function(coef, colnames_, history_y_aug, frozen_row,
                              ar_lag_idx, ar_lags, mon_idx, target_dates,
                              h_max = H_MAX,
                              holiday_dates = NULL) {
  y_hat_h <- numeric(h_max)
  for (s in 1:h_max) {
    row <- frozen_row
    n_hist <- length(history_y_aug)
    for (jj in seq_along(ar_lags)) {
      L <- ar_lags[jj]
      idx <- n_hist - L + 1L
      row[ar_lag_idx[jj]] <- history_y_aug[idx]
    }
    # Refresh ALL calendar columns (mon, weekday, holiday) at target_dates[s].
    # The legacy `mon_idx` arg is now ignored — refresh_calendar_row uses
    # name-pattern matching on colnames_.
    row <- refresh_calendar_row(row, colnames_, target_dates[s], holiday_dates)
    yh <- sum(coef * row)
    y_hat_h[s] <- yh
    history_y_aug <- c(history_y_aug, yh)
  }
  y_hat_h
}
