# =============================================================================
# PYTHON PIPELINE - COMMON UTILITIES
#
# Line-by-line port of r_pipeline/midas_codes/common_utils.R (including
# fit_r_midas_fast.R).  Same constants, same transformations, same
# publication-aware LF alignment, same calendar controls, same DM-HLN test,
# same concentrated NLS (FWL-accelerated, deterministic shape restarts) and
# the same iterated multi-step forecasting conventions (option c.1).
#
# Transformations (Weron 2014 style):
#   elec_price   : asinh     gas/brent : log
#   IPI series   : dlog      PMI/IFO/ISTAT : full-sample z-score (ddof=1)
#
# Env knobs honoured (same names as the R pipeline):
#   WINDOW_DAYS      rolling window length in days     (default 3650)
#   N_CORES          worker processes                  (default cpu-1)
#   SPEC_LIMIT       truncate the 256-spec grid        (default 256)
#   RMIDAS_DISCOUNT  discounted WLS lambda             (default 1 = off)
#   COUNTRY_NAME     Germany | Italy | ALL
#   ORIGIN_LIMIT     python-only: use only the first N rolling origins
#                    (quick tests / validation against the R outputs)
# =============================================================================

import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize

# -----------------------------------------------------------------------------
# 0. PATHS AND PARALLEL HELPERS
# -----------------------------------------------------------------------------
PIPE_DIR = os.path.dirname(os.path.abspath(__file__))


def default_data_dir():
    env = os.environ.get("V2_DATA_DIR", "")
    if env and os.path.isdir(env):
        return env
    cand = os.path.normpath(os.path.join(
        PIPE_DIR, "..", "r_pipeline", "midas_codes", "data"))
    return cand


FORECAST_DIR = os.path.join(PIPE_DIR, "forecasts")
RESULT_DIR   = os.path.join(PIPE_DIR, "results")
COEF_DIR     = os.path.join(RESULT_DIR, "coefficients")
for _d in (FORECAST_DIR, RESULT_DIR, COEF_DIR):
    os.makedirs(_d, exist_ok=True)


def get_n_cores():
    raw = os.environ.get("N_CORES", "")
    if raw:
        try:
            n = int(raw)
            if n >= 1:
                return n
        except ValueError:
            pass
    return max(1, (os.cpu_count() or 2) - 1)


# -----------------------------------------------------------------------------
# 1. CONSTANTS  (identical to common_utils.R)
# -----------------------------------------------------------------------------
AR_LAGS  = [1, 2, 3, 7]
HORIZONS = [1, 7, 14, 21, 28]
H_MAX    = max(HORIZONS)
Q_HF     = 14
PUB_LAG_IPI    = 2
PUB_LAG_SURVEY = 1

HF_CANDIDATES    = ["gas", "brent"]
LF_CANDIDATES_DE = ["pmi_mfg_de", "pmi_serv_de", "ifo_bci_de",
                    "Consumer_Goods", "Manufacturing", "Energy"]
LF_CANDIDATES_IT = ["pmi_mfg_it", "pmi_serv_it", "istat_bci_it",
                    "Consumer_Goods", "Manufacturing", "Energy"]
LF_PUB = {
    "pmi_mfg_de": PUB_LAG_SURVEY, "pmi_serv_de": PUB_LAG_SURVEY,
    "ifo_bci_de": PUB_LAG_SURVEY,
    "pmi_mfg_it": PUB_LAG_SURVEY, "pmi_serv_it": PUB_LAG_SURVEY,
    "istat_bci_it": PUB_LAG_SURVEY,
    "Consumer_Goods": PUB_LAG_IPI, "Manufacturing": PUB_LAG_IPI,
    "Energy": PUB_LAG_IPI,
}


def window_days():
    raw = os.environ.get("WINDOW_DAYS", "")
    try:
        v = int(raw)
        if v >= 365:
            return v
    except ValueError:
        pass
    return 3650


WINDOW_DAYS = window_days()


def discount_lambda():
    try:
        v = float(os.environ.get("RMIDAS_DISCOUNT", "1"))
        if np.isfinite(v) and 0 < v <= 1:
            return v
    except ValueError:
        pass
    return 1.0


def spec_limit():
    try:
        v = int(os.environ.get("SPEC_LIMIT", "256"))
        if v > 0:
            return v
    except ValueError:
        pass
    return 256


# -----------------------------------------------------------------------------
# 2. DATA LOADING + TRANSFORMATION  (= load_country_data in R)
# -----------------------------------------------------------------------------
def load_country_data(country, data_dir=None):
    assert country in ("Germany", "Italy")
    if data_dir is None:
        data_dir = default_data_dir()
    cc = "de" if country == "Germany" else "it"
    daily_path   = os.path.join(data_dir, f"prices_{cc}_daily.csv")
    monthly_path = os.path.join(data_dir, f"dataset_{cc}_m.csv")
    for p in (daily_path, monthly_path):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    daily   = pd.read_csv(daily_path)
    monthly = pd.read_csv(monthly_path)
    daily["date"]   = pd.to_datetime(daily["date"])
    monthly["date"] = pd.to_datetime(monthly["date"])
    daily   = daily.sort_values("date").reset_index(drop=True)
    monthly = monthly.sort_values("date").reset_index(drop=True)

    daily["y_asinh"]   = np.arcsinh(daily["elec_price"].to_numpy(float))
    daily["gas_log"]   = np.log(daily["gas_price"].to_numpy(float))
    daily["brent_log"] = np.log(daily["brent_price"].to_numpy(float))

    monthly = monthly.rename(columns={"Electricity_Gas_Steam": "Energy"})
    for col in ("Consumer_Goods", "Manufacturing", "Energy"):
        if col not in monthly.columns:
            raise KeyError(f"Missing monthly column: {col}")
        x = np.log(monthly[col].to_numpy(float))
        monthly[col] = np.concatenate([[0.0], np.diff(x)])
    survey_cols = (["pmi_mfg_de", "pmi_serv_de", "ifo_bci_de"]
                   if country == "Germany"
                   else ["pmi_mfg_it", "pmi_serv_it", "istat_bci_it"])
    for col in survey_cols:
        if col not in monthly.columns:
            raise KeyError(f"Missing monthly column: {col}")
        x = monthly[col].to_numpy(float)
        monthly[col] = (x - np.nanmean(x)) / np.nanstd(x, ddof=1)  # R scale()
    monthly["ref_month"] = monthly["date"].dt.to_period("M").dt.to_timestamp()
    return {"daily": daily, "monthly": monthly, "country": country}


# -----------------------------------------------------------------------------
# 3. PUBLICATION-AWARE DAILY LF SERIES
# -----------------------------------------------------------------------------
# Preferred mode: assign each monthly observation its usual public release
# date and carry the latest released value forward.  An exact historical
# column named <variable>_release_date takes precedence when supplied.
# Set LF_RELEASE_MODE=calendar_lag to reproduce the former whole-month shift.


def _roll_weekend_forward(dates):
    x = pd.DatetimeIndex(pd.to_datetime(dates))
    add = np.where(x.weekday == 5, 2, np.where(x.weekday == 6, 1, 0))
    return pd.DatetimeIndex(x + pd.to_timedelta(add, unit="D"))


def _nth_workday(month_first, n):
    out = []
    for d in pd.DatetimeIndex(pd.to_datetime(month_first)):
        days = pd.date_range(d, d + pd.offsets.MonthEnd(0), freq="B")
        out.append(days[n - 1] if len(days) >= n else pd.NaT)
    return pd.DatetimeIndex(out)


def usual_lf_release_date(ref_month, lf_col):
    ref = pd.DatetimeIndex(pd.to_datetime(ref_month)).to_period("M").to_timestamp()
    next_month = ref + pd.offsets.MonthBegin(1)
    month_end = ref + pd.offsets.MonthEnd(0)

    def env_int(name, default):
        try:
            return int(os.environ.get(name, str(default)))
        except ValueError:
            return default

    if lf_col.startswith("pmi_mfg_"):
        return _nth_workday(next_month, env_int("LF_PMI_MFG_WORKDAY", 1))
    if lf_col.startswith("pmi_serv_"):
        return _nth_workday(next_month, env_int("LF_PMI_SERV_WORKDAY", 3))
    if lf_col == "ifo_bci_de":
        day = env_int("LF_IFO_DAY", 25)
        return _roll_weekend_forward(ref + pd.to_timedelta(day - 1, unit="D"))
    if lf_col == "istat_bci_it":
        day = env_int("LF_ISTAT_BCI_DAY", 27)
        return _roll_weekend_forward(ref + pd.to_timedelta(day - 1, unit="D"))
    if lf_col in ("Consumer_Goods", "Manufacturing", "Energy"):
        delay = env_int("LF_IPI_DELAY_DAYS", 45)
        return _roll_weekend_forward(month_end + pd.to_timedelta(delay, unit="D"))
    raise KeyError(f"No usual release-date rule defined for LF variable: {lf_col}")


def make_daily_lf(daily_dates, monthly_df, lf_col, pub_lag, locf=False):
    mode = os.environ.get("LF_RELEASE_MODE", "usual_release").lower()
    vals = monthly_df[lf_col].to_numpy(float).copy()

    if mode in ("calendar_lag", "legacy"):
        if locf:
            last = np.nan
            for i in range(len(vals)):
                if np.isnan(vals[i]):
                    if not np.isnan(last):
                        vals[i] = last
                else:
                    last = vals[i]
        ref = monthly_df["ref_month"].dt.to_period("M")
        ref_index = pd.Series(np.arange(len(ref)), index=ref)
        target = pd.Series(pd.to_datetime(daily_dates)).dt.to_period("M") - pub_lag
        idx = ref_index.reindex(target).to_numpy()
        out = np.full(len(daily_dates), np.nan)
        ok = ~pd.isna(idx)
        out[ok] = vals[idx[ok].astype(int)]
        return out

    exact_col = f"{lf_col}_release_date"
    if exact_col in monthly_df.columns:
        release = pd.to_datetime(monthly_df[exact_col], errors="coerce")
    else:
        release = usual_lf_release_date(monthly_df["ref_month"], lf_col)

    rel = pd.DataFrame({
        "release_date": pd.to_datetime(release),
        "value": vals,
        "row_id": np.arange(len(vals)),
    })
    rel = rel[np.isfinite(rel["value"]) & rel["release_date"].notna()]
    rel = rel.sort_values(["release_date", "row_id"]).drop_duplicates(
        "release_date", keep="last")
    if rel.empty:
        return np.full(len(daily_dates), np.nan)

    left = pd.DataFrame({"date": pd.to_datetime(daily_dates),
                         "row_id": np.arange(len(daily_dates))})
    joined = pd.merge_asof(left.sort_values("date"),
                           rel[["release_date", "value"]],
                           left_on="date", right_on="release_date",
                           direction="backward")
    return joined.sort_values("row_id")["value"].to_numpy(float)


# -----------------------------------------------------------------------------
# 4. DESIGN-MATRIX BUILDING BLOCKS
# -----------------------------------------------------------------------------
def build_ar_lags(y, lags=AR_LAGS):
    n = len(y)
    out = np.full((n, len(lags)), np.nan)
    for j, L in enumerate(lags):
        if L < n:
            out[L:, j] = y[:n - L]
    cols = [f"y_lag{L}" for L in lags]
    return out, cols


def build_hf_exo_lags(x, prefix, Q=Q_HF):
    n = len(x)
    out = np.full((n, Q), np.nan)
    for j in range(1, Q + 1):
        out[j:, j - 1] = x[:n - j]
    cols = [f"{prefix}_lag{j}" for j in range(1, Q + 1)]
    return out, cols


def build_month_dummies(dates):
    m = dates.dt.month.to_numpy()
    M = np.zeros((len(m), 11), dtype=float)
    for k in range(2, 13):
        M[:, k - 2] = (m == k).astype(float)
    cols = [f"mon{k:02d}" for k in range(2, 13)]
    return M, cols


# Easter Sunday, Anonymous Gregorian (Butcher/Meeus) algorithm.
def easter_sunday(year):
    a = year % 19
    b = year // 100
    c = year % 100
    d = b // 4
    e = b % 4
    f = (b + 8) // 25
    g = (b - f + 1) // 3
    h = (19 * a + b - d - g + 15) % 30
    i = c // 4
    k = c % 4
    l = (32 + 2 * e + 2 * i - h - k) % 7
    m = (a + 11 * h + 22 * l) // 451
    month = (h + l - 7 * m + 114) // 31
    day   = ((h + l - 7 * m + 114) % 31) + 1
    return pd.Timestamp(year=year, month=month, day=day)


def country_holidays(years, country):
    assert country in ("Germany", "Italy")
    out = set()
    for y in years:
        e = easter_sunday(int(y))
        if country == "Germany":
            ds = [pd.Timestamp(y, 1, 1), e - pd.Timedelta(days=2),
                  e + pd.Timedelta(days=1), pd.Timestamp(y, 5, 1),
                  e + pd.Timedelta(days=39), e + pd.Timedelta(days=50),
                  pd.Timestamp(y, 10, 3), pd.Timestamp(y, 12, 25),
                  pd.Timestamp(y, 12, 26)]
        else:
            ds = [pd.Timestamp(y, 1, 1), pd.Timestamp(y, 1, 6),
                  pd.Timestamp(y, 4, 25), e + pd.Timedelta(days=1),
                  pd.Timestamp(y, 5, 1), pd.Timestamp(y, 6, 2),
                  pd.Timestamp(y, 8, 15), pd.Timestamp(y, 11, 1),
                  pd.Timestamp(y, 12, 8), pd.Timestamp(y, 12, 25),
                  pd.Timestamp(y, 12, 26)]
        out.update(ds)
    return out


def build_calendar_dummies(dates, country):
    n = len(dates)
    wd = dates.dt.dayofweek.to_numpy() + 1          # ISO 1=Mon..7=Sun
    W = np.zeros((n, 6))
    for k in range(1, 7):
        W[:, k - 1] = (wd == k).astype(float)
    wd_cols = ["wd_mon", "wd_tue", "wd_wed", "wd_thu", "wd_fri", "wd_sat"]
    M, mon_cols = build_month_dummies(dates)
    years = dates.dt.year.unique()
    hol = country_holidays(years, country)
    H = dates.isin(hol).to_numpy(float).reshape(-1, 1)
    X = np.hstack([W, M, H])
    cols = wd_cols + mon_cols + ["holiday"]
    return X, cols


# -----------------------------------------------------------------------------
# 5. ROLLING-ORIGIN INDEX GENERATOR  (0-based, = make_origins in R minus 1)
# -----------------------------------------------------------------------------
def make_origins(n_rows, window=None, pre_lag_max=None, h_max=H_MAX):
    if window is None:
        window = WINDOW_DAYS
    if pre_lag_max is None:
        pre_lag_max = max(AR_LAGS + [Q_HF])
    # R: min_o = window + pre_lag_max (1-based)  ->  0-based subtract 1
    min_o = window + pre_lag_max - 1
    max_o = n_rows - h_max - 1
    if min_o > max_o:
        raise ValueError(f"Insufficient data: min_origin={min_o} max={max_o}")
    origins = np.arange(min_o, max_o + 1)
    lim = os.environ.get("ORIGIN_LIMIT", "")
    if lim:
        try:
            L = int(lim)
            if L > 0:
                origins = origins[:L]
        except ValueError:
            pass
    return origins


# -----------------------------------------------------------------------------
# 6. DIEBOLD-MARIANO HARVEY-LEYBOURNE-NEWBOLD TEST  (= dm_hln in R)
# -----------------------------------------------------------------------------
def dm_hln(e1, e2, h, loss="sqerr"):
    from scipy import stats as st
    e1 = np.asarray(e1, float)
    e2 = np.asarray(e2, float)
    ok = np.isfinite(e1) & np.isfinite(e2)
    e1, e2 = e1[ok], e2[ok]
    T = len(e1)
    if T < 1:
        return {"stat": np.nan, "pvalue": np.nan, "T": T}
    d = e1**2 - e2**2 if loss == "sqerr" else np.abs(e1) - np.abs(e2)
    dbar = d.mean()
    g0 = ((d - dbar)**2).mean()
    acov = 0.0
    if h > 1:
        for k in range(1, h):
            gk = np.sum((d[k:] - dbar) * (d[:-k] - dbar)) / T
            acov += 2 * gk
    v = (g0 + acov) / T
    if not np.isfinite(v) or v <= 0:
        return {"stat": np.nan, "pvalue": np.nan, "T": T}
    dm = dbar / np.sqrt(v)
    corr = np.sqrt((T + 1 - 2 * h + h * (h - 1) / T) / T)
    stat = dm * corr
    pval = 2 * st.t.sf(abs(stat), df=T - 1)
    return {"stat": stat, "pvalue": pval, "T": T, "dbar": dbar}


# -----------------------------------------------------------------------------
# 7. CALENDAR REFRESH AT TARGET ROW  (= refresh_calendar_row in R)
# -----------------------------------------------------------------------------
WD_NAMES = ["wd_mon", "wd_tue", "wd_wed", "wd_thu", "wd_fri", "wd_sat"]


def refresh_calendar_row(row, colnames, target_date, holiday_set=None):
    row = row.copy()
    cidx = {c: i for i, c in enumerate(colnames)}
    mon_cols = [c for c in colnames if c.startswith("mon") and len(c) == 5
                and c[3:].isdigit()]
    if mon_cols:
        for c in mon_cols:
            row[cidx[c]] = 0.0
        m = target_date.month
        if m >= 2:
            c = f"mon{m:02d}"
            if c in cidx:
                row[cidx[c]] = 1.0
    wd_present = [c for c in WD_NAMES if c in cidx]
    if wd_present:
        for c in wd_present:
            row[cidx[c]] = 0.0
        wd_iso = target_date.dayofweek + 1
        if 1 <= wd_iso <= 6:
            c = WD_NAMES[wd_iso - 1]
            if c in cidx:
                row[cidx[c]] = 1.0
    if "holiday" in cidx and holiday_set is not None:
        row[cidx["holiday"]] = 1.0 if target_date in holiday_set else 0.0
    return row


# -----------------------------------------------------------------------------
# 8. EXP-ALMON + CONCENTRATED NLS  (= fit_r_midas_fast.R)
# -----------------------------------------------------------------------------
def exp_almon_weights_j(jvec, theta):
    """Exp-Almon weights over an arbitrary lag set (the actual j values enter
    the exponent), e.g. J = {1, 2, 3, 7} for the paper's AR block."""
    j = np.asarray(jvec, dtype=float)
    lw = theta[0] * j + theta[1] * j * j
    lw -= lw.max()
    w = np.exp(lw)
    return w / w.sum()


def exp_almon_weights(Q, theta):
    return exp_almon_weights_j(np.arange(1, Q + 1), theta)


CANONICAL_SHAPES = np.array([
    [0.0,  0.00],
    [-1.0, 0.00],
    [-0.2, 0.00],
    [0.5, -0.10],
    [0.1, -0.02],
])

_EPS = np.finfo(float).eps


def _qr_rank(R):
    d = np.abs(np.diag(R))
    if d.size == 0:
        return 0
    tol = d.max() * max(R.shape) * _EPS
    return int((d > tol).sum())


def fit_r_midas_fast(y, Z_list, X_lin, lin_cols,
                     theta_init=None,
                     theta_lower_block=(-8.0, -8.0),
                     theta_upper_block=(8.0, 0.0),
                     optim_maxit=2000, optim_factr=5e8,
                     n_starts=1, deterministic_starts=True,
                     lag_sets=None):
    """Concentrated NLS for R-MIDAS; identical estimator to the R version.

    y       : (n,) response (already filtered if desired; refiltered here)
    Z_list  : dict name -> (n, Q_k) lag matrix (Almon blocks; first = AR)
    X_lin   : (n, p) linear design (intercept first)
    lin_cols: list of p column names
    lag_sets: optional dict name -> lag numbers of the block's columns
              (e.g. {"y_lags": [1, 2, 3, 7]}); defaults to consecutive 1..Q_k
    """
    names = list(Z_list.keys())
    K = len(names)
    if lag_sets is None:
        lag_sets = {}
    lag_sets = {nm: np.asarray(lag_sets.get(
        nm, np.arange(1, Z_list[nm].shape[1] + 1)), float) for nm in names}
    Qs = [Z_list[nm].shape[1] for nm in names]
    y = np.asarray(y, float)
    n = len(y)
    keep = np.isfinite(y) & np.isfinite(X_lin).all(axis=1)
    for nm in names:
        keep &= np.isfinite(Z_list[nm]).all(axis=1)
    if keep.sum() < X_lin.shape[1] + K + 5:
        return None
    Xk = X_lin[keep]
    yk = y[keep]
    Zk = [Z_list[nm][keep] for nm in names]

    # ---- one-off FWL pre-computation ----
    Qx, Rx = np.linalg.qr(Xk)
    if _qr_rank(Rx) < Xk.shape[1]:
        return None
    Zcat = np.hstack(Zk)                       # n x KQ
    Zt = Zcat - Qx @ (Qx.T @ Zcat)
    yt = yk - Qx @ (Qx.T @ yk)
    G  = Zt.T @ Zt
    cv = Zt.T @ yt
    yy = float(yt @ yt)
    Qtot = sum(Qs)
    offs = np.concatenate([[0], np.cumsum(Qs)])

    def obj(theta_v):
        W = np.zeros((Qtot, K))
        for k in range(K):
            W[offs[k]:offs[k + 1], k] = exp_almon_weights_j(
                lag_sets[names[k]], theta_v[2 * k:2 * k + 2])
        A = W.T @ G @ W
        b = W.T @ cv
        try:
            d = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return 1e12
        if not np.all(np.isfinite(d)):
            return 1e12
        sse = yy - float(b @ d)
        if not np.isfinite(sse) or sse < -1e-6 * yy:
            return 1e12
        return max(sse, 0.0)

    lower = list(theta_lower_block) * K
    upper = list(theta_upper_block) * K
    bounds = list(zip(lower, upper))
    dim = 2 * K
    n_starts = max(1, int(n_starts))
    if deterministic_starts:
        n_starts = min(n_starts, 1 + len(CANONICAL_SHAPES))
    starts = []
    if theta_init is not None and len(theta_init) == dim and \
            np.all(np.isfinite(theta_init)):
        starts.append(np.asarray(theta_init, float))
    else:
        starts.append(np.zeros(dim))
    for s in range(1, n_starts):
        if deterministic_starts:
            starts.append(np.tile(CANONICAL_SHAPES[s - 1], K))
        else:
            starts.append(np.random.uniform(lower, upper))

    ftol = optim_factr * _EPS                  # R optim factr convention
    best = None
    best_val = np.inf
    for s0 in starts:
        try:
            res = minimize(obj, s0, method="L-BFGS-B", bounds=bounds,
                           options={"maxiter": optim_maxit, "ftol": ftol})
        except Exception:
            continue
        if not np.isfinite(res.fun) or res.fun >= 1e11:
            continue
        if res.status not in (0, 1):
            continue
        if res.fun < best_val:
            best_val = res.fun
            best = res
    if best is None:
        return None
    theta_hat = best.x

    # ---- final fit at theta_hat (same code path as R) ----
    zmat = np.column_stack([
        Zk[k] @ exp_almon_weights_j(lag_sets[names[k]],
                                    theta_hat[2 * k:2 * k + 2])
        for k in range(K)])
    Xfull = np.hstack([Xk, zmat])
    beta, _, rank, _ = np.linalg.lstsq(Xfull, yk, rcond=None)
    if rank < Xfull.shape[1]:
        return None
    res_v = yk - Xfull @ beta
    p_lin = Xk.shape[1]
    beta_lin = dict(zip(lin_cols, beta[:p_lin]))
    delta = dict(zip(names, beta[p_lin:p_lin + K]))
    theta_mat = {nm: theta_hat[2 * k:2 * k + 2].copy()
                 for k, nm in enumerate(names)}
    return {"theta": theta_mat, "beta_lin": beta_lin, "delta": delta,
            "sse": float(best.fun), "sigma2": float(np.mean(res_v**2)),
            "converged": best.status in (0, 1), "n_eff": int(keep.sum()),
            "n_starts": n_starts}


# -----------------------------------------------------------------------------
# 9. ITERATED MULTI-STEP FORECASTS  (option c.1, = R versions)
# -----------------------------------------------------------------------------
def iterated_forecast_r_midas(fit, Z_full, X_lin_row, lin_cols,
                              history_y, origin_idx, target_dates,
                              ar_block_name="y_lags", hf_block_names=(),
                              h_max=H_MAX, holiday_set=None, lag_sets=None):
    names = list(fit["theta"].keys())
    if lag_sets is None:
        lag_sets = {}
    lag_sets = {nm: np.asarray(lag_sets.get(
        nm, np.arange(1, Z_full[nm].shape[1] + 1)), int) for nm in names}
    W = {nm: exp_almon_weights_j(lag_sets[nm], fit["theta"][nm])
         for nm in names}
    delta = fit["delta"]
    beta_lin = np.array([fit["beta_lin"][c] for c in lin_cols])
    ar_lag_set = lag_sets[ar_block_name]
    max_ar_lag = int(ar_lag_set.max())

    frozen_z = {}
    for nm in hf_block_names:
        frozen_z[nm] = float(W[nm] @ Z_full[nm][origin_idx + 1])

    hist = list(np.asarray(history_y, float))
    out = np.full(h_max, np.nan)
    row0 = np.asarray(X_lin_row, float)
    for s in range(h_max):
        row = refresh_calendar_row(row0, lin_cols, target_dates[s],
                                   holiday_set)
        lin_part = float(beta_lin @ row)
        if len(hist) < max_ar_lag:
            hist.append(np.nan)
            continue
        lags = np.array([hist[-L] for L in ar_lag_set])
        z_y = float(W[ar_block_name] @ lags)
        hf_part = sum(delta[nm] * frozen_z[nm] for nm in hf_block_names)
        yh = lin_part + delta[ar_block_name] * z_y + hf_part
        out[s] = yh
        hist.append(yh)
    return out


def iterated_forecast_linear(coef, colnames, history_y, frozen_row,
                             ar_lag_idx, ar_lags, target_dates,
                             h_max=H_MAX, holiday_set=None):
    hist = list(np.asarray(history_y, float))
    out = np.full(h_max, np.nan)
    coef = np.asarray(coef, float)
    for s in range(h_max):
        row = np.asarray(frozen_row, float).copy()
        for jj, L in enumerate(ar_lags):
            row[ar_lag_idx[jj]] = hist[-L]
        row = refresh_calendar_row(row, colnames, target_dates[s],
                                   holiday_set)
        yh = float(coef @ row)
        out[s] = yh
        hist.append(yh)
    return out


# -----------------------------------------------------------------------------
# 10. SPEC GRID  (= masks construction in the R runners; spec_id = bit value)
# -----------------------------------------------------------------------------
def build_spec_grid(hf_names, lf_names, limit=None):
    """256 subsets of {gas, brent} x LF set, ordered exactly like the R
    scripts: integer value of the flag vector with the FIRST flag as the most
    significant bit, ascending; spec_id index equals that integer."""
    flags = list(hf_names) + list(lf_names)
    nf = len(flags)
    rows = []
    for v in range(2 ** nf):
        bits = [(v >> (nf - 1 - j)) & 1 for j in range(nf)]
        rows.append(bits)
    masks = pd.DataFrame(rows, columns=flags)
    masks.index = [f"spec_{i:03d}" for i in range(len(masks))]
    if limit is not None and limit < len(masks):
        masks = masks.iloc[:limit]
    return masks


def spec_dictionary(masks, hf_names, lf_names):
    recs = []
    for sid, row in masks.iterrows():
        hf = [v for v in hf_names if row[v] == 1]
        lf = [v for v in lf_names if row[v] == 1]
        recs.append({"spec_id": sid, "n_hf": len(hf), "n_lf": len(lf),
                     "hf_vars": "|".join(hf), "lf_vars": "|".join(lf)})
    return pd.DataFrame(recs)


# -----------------------------------------------------------------------------
# 11. OUTPUT HELPER  (forecast CSV with the same layout as the R scripts)
# -----------------------------------------------------------------------------
def write_forecast_csv(path, dates, origins, h, y_actual_level, yhat_cols):
    """yhat_cols: dict column_name -> (n_orig,) array (already in EUR/MWh)."""
    out = pd.DataFrame({
        "origin_date": dates.iloc[origins].dt.strftime("%Y-%m-%d").to_numpy(),
        "target_date": dates.iloc[origins + h].dt.strftime("%Y-%m-%d").to_numpy(),
        "y_actual":    y_actual_level[origins + h],
    })
    for cname, vals in yhat_cols.items():
        out[cname] = vals
    out = out[np.isfinite(out["y_actual"].to_numpy(float))]
    out.to_csv(path, index=False)
    return len(out)
