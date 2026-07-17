#!/usr/bin/env Rscript
# =============================================================================
# v2 PIPELINE — RU-MIDAS WITH 28 WITHIN-MONTH POSITION DUMMIES (CIMADOMO-STYLE)
#
# Replicates the design structure of Foroni-Rossini-et-al. (ECB WP 2250 / 2019,
# "Forecasting daily electricity prices with monthly macroeconomic variables"),
# Eq. (6): every lag of every regressor is interacted with within-month
# position dummies D_1, ..., D_28.  Days whose strict within-month position
# is > 28 are dropped from estimation (paper convention).
#
# Single-equation OLS form:
#
#   x_t  =  Σ_{j} ψ_j  C_{j,t}                          calendar controls (linear)
#          + Σ_{i=1..28} α_i  D_i  z^{LF}_{v,t}        per LF v in spec
#          + Σ_{i=1..28} β_{1,i} D_i  x_{t-1}          AR lag 1
#          + Σ_{i=1..28} β_{2,i} D_i  x_{t-2}          AR lag 2
#          + Σ_{i=1..28} β_{3,i} D_i  x_{t-7}          AR lag 7 (Weron weekly lag)
#          + Σ_{i=1..28} γ_{1,i} D_i  gas_{t-1}        per HF lag of gas in spec
#          + Σ_{i=1..28} γ_{2,i} D_i  brent_{t-1}      per HF lag of brent in spec
#          + v_t
#
# Two distinct dummy roles:
#   D_1, ..., D_28   ← within-month POSITION dummies, the MIDAS-structural
#                       element; interact multiplicatively with every lag block.
#   C_j              ← CALENDAR controls (Foroni-Rossini-style daily base):
#                        6 weekday dummies (Sun = reference),
#                       11 month-of-year dummies (Jan = reference),
#                        1 holiday dummy (country-specific national holidays).
#                       These are plain linear regressors; they do NOT interact
#                       with position dummies.
#
# AR-y lag set is paper-faithful {1, 2, 7}.  HF exogenous default to a single
# lag (gas_{t-1}, brent_{t-1}), matching the paper's "one lag of daily oil"
# specification.  Position dummies are SATURATED — there is no separate
# global intercept, since the 28 dummies span the constant.
#
# Spec grid: 2^8 = 256 subsets of {gas, brent} x {pmi_mfg, pmi_serv, bci,
# Consumer_Goods, Manufacturing, Energy}.  SPEC_LIMIT env var to truncate.
#
# Iterated multi-step (option c.1):
#   - y AR feedback via history_y_aug, weighted by the position dummies at
#     each target step (i.e., the AR coefficients used at step s are those
#     of position i_s = within_month_position(origin + s))
#   - HF exo (gas, brent) FROZEN at row (origin + 1)
#   - LF FROZEN at row (origin + 1)
#   - Position dummies REFRESHED at the target row
#   - Calendar controls (weekday, MOY, holiday) REFRESHED at target row
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
farg <- args[grep("^--file=", args)]
v2.boot.path <- if (length(farg) > 0) {
  raw <- gsub("~+~", " ", sub("^--file=", "", farg[1]), fixed = TRUE)
  file.path(dirname(suppressWarnings(normalizePath(raw, winslash = "/", mustWork = FALSE))),
            "v2_paths.R")
} else file.path(getwd(), "v2_paths.R")
source(v2.boot.path)

source(file.path(script_dir, "common_utils.R"))

DATA_DIR     <- file.path(script_dir, "data")
FORECAST_DIR <- file.path(script_dir, "forecasts")
dir.create(FORECAST_DIR, showWarnings = FALSE, recursive = TRUE)

AR_LAGS_POS <- c(1L, 2L, 7L)
Q_HF_POS    <- 1L
K_POS       <- 28L
SPEC_LIMIT  <- suppressWarnings(as.integer(Sys.getenv("SPEC_LIMIT", "256")))
if (is.na(SPEC_LIMIT) || SPEC_LIMIT <= 0L) SPEC_LIMIT <- 256L

within_month_position <- function(dates, k = K_POS) {
  pos <- as.integer(format(dates, "%d"))
  ifelse(pos > k, NA_integer_, pos)
}

build_position_dummies <- function(dates, k = K_POS) {
  pos <- within_month_position(dates, k)
  D <- matrix(0L, length(dates), k)
  for (i in 1:k) D[, i] <- as.integer(!is.na(pos) & pos == i)
  colnames(D) <- paste0("D", sprintf("%02d", 1:k))
  D
}

interact_with_positions <- function(src, D, src_prefix) {
  n <- nrow(src); p <- ncol(src); k <- ncol(D)
  out <- matrix(NA_real_, n, p * k)
  cn  <- character(p * k)
  idx <- 1L
  for (j in seq_len(p)) {
    src_col <- src[, j]
    for (i in seq_len(k)) {
      out[, idx] <- src_col * D[, i]
      cn[idx]    <- sprintf("%s_%s_x_%s", src_prefix, colnames(src)[j], colnames(D)[i])
      idx <- idx + 1L
    }
  }
  colnames(out) <- cn
  out
}

# -----------------------------------------------------------------------------
# RUNNER
# -----------------------------------------------------------------------------
run_pos_country <- function(country) {
  n_cores <- get_n_cores()
  cat(sprintf("\n=== RU-MIDAS POS-DUMMIES (Cimadomo-style) — %s (n_cores=%d) ===\n",
              country, n_cores))
  cc      <- if (country == "Germany") "de" else "it"
  dat     <- load_country_data(country, DATA_DIR)
  daily   <- dat$daily
  monthly <- dat$monthly
  dates   <- daily$date
  n       <- nrow(daily)
  y       <- daily$y_asinh
  pos     <- within_month_position(dates, K_POS)

  ar_lag_mat <- matrix(NA_real_, n, length(AR_LAGS_POS),
                       dimnames = list(NULL, paste0("y_lag", AR_LAGS_POS)))
  for (j in seq_along(AR_LAGS_POS)) {
    L <- AR_LAGS_POS[j]
    if (L < n) ar_lag_mat[(L + 1L):n, j] <- y[1:(n - L)]
  }
  gas_lag <- c(NA_real_, daily$gas_log[1:(n - 1L)])
  brt_lag <- c(NA_real_, daily$brent_log[1:(n - 1L)])
  hf_lag_mat <- cbind(gas = gas_lag, brent = brt_lag)
  hf_lookup_lag1 <- list(gas = gas_lag, brent = brt_lag)

  lf_names <- if (country == "Germany") LF_CANDIDATES_DE else LF_CANDIDATES_IT
  lf_pub   <- c(
    pmi_mfg_de = PUB_LAG_SURVEY, pmi_serv_de = PUB_LAG_SURVEY, ifo_bci_de = PUB_LAG_SURVEY,
    pmi_mfg_it = PUB_LAG_SURVEY, pmi_serv_it = PUB_LAG_SURVEY, istat_bci_it = PUB_LAG_SURVEY,
    Consumer_Goods = PUB_LAG_IPI, Manufacturing = PUB_LAG_IPI, Energy = PUB_LAG_IPI
  )
  lf_mat <- matrix(NA_real_, n, length(lf_names),
                   dimnames = list(NULL, lf_names))
  for (v in lf_names) lf_mat[, v] <- make_daily_lf(dates, monthly, v, lf_pub[[v]])

  # ---- Position dummies + calendar (weekday/month/holiday) controls ------
  # Position dummies are the MIDAS-structural element (interact with every lag
  # block).  Calendar dummies are plain linear seasonal controls — 6 weekday
  # (Sun = ref) + 11 month-of-year (Jan = ref) + 1 holiday.
  D_mat   <- build_position_dummies(dates, K_POS)
  X_CAL   <- build_calendar_dummies(dates, country)
  cal_cols <- colnames(X_CAL)
  hol_set <- country_holidays(unique(as.integer(format(dates, "%Y"))), country)

  # AR block: 3 lags × 28 positions = 84 cols
  X_AR  <- interact_with_positions(ar_lag_mat,  D_mat, "ar")
  X_HF  <- list(
    gas   = interact_with_positions(matrix(hf_lag_mat[, "gas"],   ncol = 1, dimnames = list(NULL, "gas")),   D_mat, "hf"),
    brent = interact_with_positions(matrix(hf_lag_mat[, "brent"], ncol = 1, dimnames = list(NULL, "brent")), D_mat, "hf")
  )
  X_LF  <- list()
  for (v in lf_names) {
    X_LF[[v]] <- interact_with_positions(
      matrix(lf_mat[, v], ncol = 1, dimnames = list(NULL, v)),
      D_mat, "lf"
    )
  }

  hf_names   <- c("gas", "brent")
  flag_names <- c(hf_names, lf_names)

  masks <- expand.grid(replicate(length(flag_names), 0:1, simplify = FALSE),
                       KEEP.OUT.ATTRS = FALSE)
  masks <- as.matrix(masks); colnames(masks) <- flag_names
  masks_int <- as.integer(masks %*% (2^(rev(seq_len(ncol(masks))) - 1L)))
  ord <- order(masks_int)
  masks <- masks[ord, , drop = FALSE]
  spec_ids <- sprintf("spec_%03d", 0:(nrow(masks) - 1L))
  rownames(masks) <- spec_ids
  if (SPEC_LIMIT < nrow(masks)) {
    masks <- masks[1:SPEC_LIMIT, , drop = FALSE]
    spec_ids <- spec_ids[1:SPEC_LIMIT]
    cat(sprintf("  SPEC_LIMIT=%d (truncating from 256)\n", SPEC_LIMIT))
  }

  dict <- data.frame(
    spec_id = spec_ids,
    n_hf    = rowSums(masks[, hf_names, drop = FALSE]),
    n_lf    = rowSums(masks[, lf_names, drop = FALSE]),
    hf_vars = vapply(seq_len(nrow(masks)), function(i) {
      v <- hf_names[masks[i, hf_names] == 1L]
      if (length(v) == 0L) "" else paste(v, collapse = "|")
    }, character(1)),
    lf_vars = vapply(seq_len(nrow(masks)), function(i) {
      v <- lf_names[masks[i, lf_names] == 1L]
      if (length(v) == 0L) "" else paste(v, collapse = "|")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  dict_path <- file.path(FORECAST_DIR,
                         sprintf("ru_midas_spec_dictionary_%s.csv", cc))
  write.csv(dict, dict_path, row.names = FALSE)
  cat(sprintf("  wrote %s  (specs = %d)\n", dict_path, nrow(dict)))

  origins <- make_origins(n, WINDOW_DAYS,
                          pre_lag_max = max(c(AR_LAGS_POS, Q_HF_POS)),
                          h_max = H_MAX)
  n_orig <- length(origins)
  cat(sprintf("  rolling origins: %d (from %s to %s)\n",
              n_orig, format(dates[origins[1]]),
              format(dates[origins[n_orig]])))
  cat(sprintf("  base AR block: %d cols (3 lags × 28 pos)\n", ncol(X_AR)))
  cat(sprintf("  calendar block: %d cols (6 wd + 11 mon + 1 hol)\n", ncol(X_CAL)))

  fit_one_spec <- function(si) {
    sid <- spec_ids[si]
    mk  <- masks[si, ]
    hf_used <- character(0)
    if (mk["gas"]   == 1L) hf_used <- c(hf_used, "gas")
    if (mk["brent"] == 1L) hf_used <- c(hf_used, "brent")
    lf_used <- lf_names[mk[lf_names] == 1L]

    # Assemble full design: AR + calendar (linear) + HF/LF blocks per mask.
    X_blocks <- list(X_AR, X_CAL)
    for (nm in hf_used) X_blocks[[length(X_blocks) + 1L]] <- X_HF[[nm]]
    for (v in lf_used)  X_blocks[[length(X_blocks) + 1L]] <- X_LF[[v]]
    X <- do.call(cbind, X_blocks)

    yhat_mat <- matrix(NA_real_, n_orig, H_MAX)

    for (oi in seq_along(origins)) {
      t0 <- origins[oi]
      win_lo <- t0 - WINDOW_DAYS + 1L
      win_rows <- win_lo:t0
      keep_pos <- !is.na(pos[win_rows])
      Xtr <- X[win_rows, , drop = FALSE]
      ytr <- y[win_rows]
      ok  <- keep_pos & stats::complete.cases(Xtr) & is.finite(ytr)
      if (sum(ok) < ncol(Xtr) + 5L) next
      Xfit <- Xtr[ok, , drop = FALSE]
      yfit <- ytr[ok]
      qrfit <- qr(Xfit)
      if (qrfit$rank < ncol(Xfit)) next
      beta <- qr.coef(qrfit, yfit)
      if (anyNA(beta)) next
      names(beta) <- colnames(X)

      target_dates <- dates[(t0 + 1L):(t0 + H_MAX)]
      target_pos   <- within_month_position(target_dates, K_POS)
      frozen_hf <- list()
      for (nm in hf_used) frozen_hf[[nm]] <- hf_lookup_lag1[[nm]][t0 + 1L]
      frozen_lf <- numeric(length(lf_used)); names(frozen_lf) <- lf_used
      if (length(lf_used) > 0L) {
        for (v in lf_used) frozen_lf[v] <- lf_mat[t0 + 1L, v]
      }
      if (any(!is.finite(unlist(frozen_hf))) || any(!is.finite(frozen_lf))) next

      history_y_aug <- y[1:t0]
      h_path <- rep(NA_real_, H_MAX)

      for (s in 1:H_MAX) {
        i_s <- target_pos[s]
        if (is.na(i_s)) {
          h_path[s] <- if (s > 1L) h_path[s - 1L] else NA_real_
          history_y_aug <- c(history_y_aug,
                             if (is.finite(h_path[s])) h_path[s] else 0)
          next
        }
        n_hist <- length(history_y_aug)
        ar_vals <- numeric(length(AR_LAGS_POS))
        for (jj in seq_along(AR_LAGS_POS)) {
          ar_vals[jj] <- history_y_aug[n_hist - AR_LAGS_POS[jj] + 1L]
        }
        D_tag <- sprintf("D%02d", i_s)
        # AR coefficient lookup
        ar_part <- 0
        for (jj in seq_along(AR_LAGS_POS)) {
          col <- sprintf("ar_y_lag%d_x_%s", AR_LAGS_POS[jj], D_tag)
          b   <- beta[col]
          if (!is.na(b)) ar_part <- ar_part + b * ar_vals[jj]
        }
        # HF coefficient lookup
        hf_part <- 0
        for (nm in hf_used) {
          col <- sprintf("hf_%s_x_%s", nm, D_tag)
          b   <- beta[col]
          if (!is.na(b)) hf_part <- hf_part + b * frozen_hf[[nm]]
        }
        # LF coefficient lookup
        lf_part <- 0
        for (v in lf_used) {
          col <- sprintf("lf_%s_x_%s", v, D_tag)
          b   <- beta[col]
          if (!is.na(b)) lf_part <- lf_part + b * frozen_lf[v]
        }
        # Calendar contributions (plain linear, refreshed at target_date):
        #   - 11 MOY dummies (Jan = ref, so no contribution if target month = 1)
        #   - 6 weekday dummies (Sun = ref, no contribution on Sundays)
        #   - 1 holiday dummy (0 / 1 depending on holiday set membership)
        target_date <- dates[t0 + s]
        # MOY
        moy_part <- 0
        target_month <- as.integer(format(target_date, "%m"))
        if (target_month >= 2L) {
          col <- sprintf("mon%02d", target_month)
          b   <- beta[col]
          if (!is.na(b)) moy_part <- b
        }
        # Weekday
        wd_part <- 0
        wd_iso <- as.integer(format(target_date, "%u"))   # 1=Mon..7=Sun
        if (wd_iso >= 1L && wd_iso <= 6L) {
          wd_name <- c("wd_mon","wd_tue","wd_wed","wd_thu","wd_fri","wd_sat")[wd_iso]
          b <- beta[wd_name]
          if (!is.na(b)) wd_part <- b
        }
        # Holiday
        hol_part <- 0
        if (target_date %in% hol_set) {
          b <- beta["holiday"]
          if (!is.na(b)) hol_part <- b
        }

        yh <- ar_part + hf_part + lf_part + moy_part + wd_part + hol_part
        h_path[s] <- yh
        history_y_aug <- c(history_y_aug, yh)
      }

      for (h in HORIZONS) {
        if (h <= length(h_path) && is.finite(h_path[h]))
          yhat_mat[oi, h] <- sinh(h_path[h])
      }
    }
    yhat_mat
  }

  t_start <- Sys.time()
  spec_results <- parallel_lapply(seq_along(spec_ids), fit_one_spec,
                                  n_cores = n_cores)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  fitting done in %.1fs (%.2fh) — wall-clock\n", el, el / 3600))

  target_pos_full <- within_month_position(dates, K_POS)
  for (h in HORIZONS) {
    yhat_h_mat <- vapply(seq_along(spec_ids),
                         function(si) spec_results[[si]][, h],
                         numeric(n_orig))
    colnames(yhat_h_mat) <- spec_ids
    out <- data.frame(
      origin_date = dates[origins],
      target_date = dates[origins + h],
      target_pos  = target_pos_full[origins + h],
      y_actual    = daily$elec_price[origins + h],
      stringsAsFactors = FALSE
    )
    out <- cbind(out, as.data.frame(yhat_h_mat, stringsAsFactors = FALSE))
    out <- out[!is.na(out$y_actual), , drop = FALSE]
    fp  <- file.path(FORECAST_DIR,
                     sprintf("ru_midas_%s_h%02d.csv", cc, h))
    write.csv(out, fp, row.names = FALSE)
    cat(sprintf("  wrote %s  (rows = %d)\n", fp, nrow(out)))
  }
}

country_arg <- Sys.getenv("COUNTRY_NAME", unset = "ALL")
countries   <- if (toupper(country_arg) == "ALL")
  c("Germany", "Italy") else country_arg
for (cc_ in countries) run_pos_country(cc_)
cat("\nRU-MIDAS POS-DUMMIES (Cimadomo-style) done.\n")
