#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - RU-MIDAS WITH 28 WITHIN-MONTH POSITION DUMMIES
# Port of r_pipeline/midas_codes/ru_midas.R (Cimadomo / Foroni-Rossini
# style, ECB WP 2250 Eq. 6).
#
# Every lag of every regressor is interacted with within-month position
# dummies D1..D28; days with day-of-month > 28 are dropped from estimation.
# AR lags {1, 2, 7}; HF gas/brent at lag 1; LF publication-aware daily
# columns; calendar controls (6 weekday + 11 month + 1 holiday) enter as
# plain linear regressors.  Position dummies are saturated (no intercept).
#
# OLS per rolling window; iterated multi-step with position-specific
# coefficients at each target step; HF and LF frozen at (origin + 1).
#
# Outputs: forecasts/ru_midas_<cc>_h<h>.csv (+ spec dictionary), same layout
# as the R script (including the target_pos column).
# =============================================================================

import os
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

import common_utils as cu

AR_LAGS_POS = [1, 2, 7]
K_POS = 28

_CTX = {}


def within_month_position(dates, k=K_POS):
    pos = dates.dt.day.to_numpy().astype(float)
    pos[pos > k] = np.nan
    return pos


def _build_context(country):
    cc = "de" if country == "Germany" else "it"
    dat = cu.load_country_data(country)
    daily = dat["daily"]
    monthly = dat["monthly"]
    dates = daily["date"]
    n = len(daily)
    y = daily["y_asinh"].to_numpy()
    pos = within_month_position(dates)

    ar_lag_mat = np.full((n, len(AR_LAGS_POS)), np.nan)
    for j, L in enumerate(AR_LAGS_POS):
        ar_lag_mat[L:, j] = y[:n - L]
    gas_lag = np.concatenate([[np.nan], daily["gas_log"].to_numpy()[:-1]])
    brt_lag = np.concatenate([[np.nan], daily["brent_log"].to_numpy()[:-1]])

    lf_names = (cu.LF_CANDIDATES_DE if country == "Germany"
                else cu.LF_CANDIDATES_IT)
    # NOTE: ru_midas.R calls make_daily_lf without LOCF; kept identical.
    lf_mat = np.column_stack([
        cu.make_daily_lf(dates, monthly, v, cu.LF_PUB[v], locf=False)
        for v in lf_names])

    # Position dummies (n x 28)
    D = np.zeros((n, K_POS))
    for i in range(1, K_POS + 1):
        D[:, i - 1] = (pos == i).astype(float)

    cal_mat, cal_cols = cu.build_calendar_dummies(dates, country)
    hol_set = cu.country_holidays(dates.dt.year.unique(), country)

    # Interacted blocks.  Column order matches interact_with_positions in R:
    # for each source column j, all 28 positions.
    def interact(src):
        p = src.shape[1]
        out = np.empty((n, p * K_POS))
        for j in range(p):
            out[:, j * K_POS:(j + 1) * K_POS] = src[:, [j]] * D
        return out

    X_AR = interact(ar_lag_mat)                       # 84 cols
    X_HF = {"gas":   interact(gas_lag.reshape(-1, 1)),
            "brent": interact(brt_lag.reshape(-1, 1))}
    X_LF = {v: interact(lf_mat[:, [j]])
            for j, v in enumerate(lf_names)}

    origins = cu.make_origins(n, cu.WINDOW_DAYS,
                              pre_lag_max=max(AR_LAGS_POS + [1]),
                              h_max=cu.H_MAX)
    masks = cu.build_spec_grid(cu.HF_CANDIDATES, lf_names, cu.spec_limit())

    return {"cc": cc, "daily": daily, "dates": dates, "n": n, "y": y,
            "pos": pos, "lf_names": lf_names, "lf_mat": lf_mat,
            "gas_lag": gas_lag, "brt_lag": brt_lag,
            "X_AR": X_AR, "X_HF": X_HF, "X_LF": X_LF,
            "cal_mat": cal_mat, "cal_cols": cal_cols, "hol_set": hol_set,
            "origins": origins, "masks": masks}


def _init_worker(country):
    global _CTX
    _CTX = _build_context(country)


def fit_one_spec(si):
    c = _CTX
    mk = c["masks"].iloc[si]
    lf_names = c["lf_names"]
    hf_used = [v for v in ("gas", "brent") if mk[v] == 1]
    lf_used = [v for v in lf_names if mk[v] == 1]

    cal_cols = c["cal_cols"]
    n_ar = len(AR_LAGS_POS)
    blocks = [c["X_AR"], c["cal_mat"]]
    for nm in hf_used:
        blocks.append(c["X_HF"][nm])
    for v in lf_used:
        blocks.append(c["X_LF"][v])
    X = np.hstack(blocks)

    # Column index bookkeeping (replaces R's name-based lookup):
    #   AR block: lag jj, position i  ->  jj*K_POS + (i-1)
    #   calendar: offset n_ar*K_POS, order = cal_cols
    #   HF nm   : offset_hf[nm] + (i-1)
    #   LF v    : offset_lf[v] + (i-1)
    off_cal = n_ar * K_POS
    off = off_cal + len(cal_cols)
    off_hf = {}
    for nm in hf_used:
        off_hf[nm] = off
        off += K_POS
    off_lf = {}
    for v in lf_used:
        off_lf[v] = off
        off += K_POS
    cal_pos = {cname: off_cal + j for j, cname in enumerate(cal_cols)}

    y = c["y"]
    pos = c["pos"]
    dates = c["dates"]
    origins = c["origins"]
    yhat_mat = np.full((len(origins), cu.H_MAX), np.nan)

    for oi, t0 in enumerate(origins):
        lo = t0 - cu.WINDOW_DAYS + 1
        sl = slice(lo, t0 + 1)
        Xtr = X[sl]
        ytr = y[sl]
        ok = (~np.isnan(pos[sl])) & np.isfinite(Xtr).all(axis=1) & \
            np.isfinite(ytr)
        if ok.sum() < Xtr.shape[1] + 5:
            continue
        beta, _, rank, _ = np.linalg.lstsq(Xtr[ok], ytr[ok], rcond=None)
        if rank < Xtr.shape[1] or not np.all(np.isfinite(beta)):
            continue

        frozen_hf = {nm: (c["gas_lag"] if nm == "gas" else
                          c["brt_lag"])[t0 + 1] for nm in hf_used}
        frozen_lf = {v: c["lf_mat"][t0 + 1, lf_names.index(v)]
                     for v in lf_used}
        if (any(not np.isfinite(v) for v in frozen_hf.values()) or
                any(not np.isfinite(v) for v in frozen_lf.values())):
            continue

        hist = list(y[:t0 + 1])
        h_path = np.full(cu.H_MAX, np.nan)
        for s in range(1, cu.H_MAX + 1):
            tgt = t0 + s
            tdate = dates.iloc[tgt]
            i_s = tdate.day
            if i_s > K_POS:
                h_path[s - 1] = h_path[s - 2] if s > 1 else np.nan
                hist.append(h_path[s - 1] if np.isfinite(h_path[s - 1])
                            else 0.0)
                continue
            ar_part = 0.0
            for jj, L in enumerate(AR_LAGS_POS):
                ar_part += beta[jj * K_POS + (i_s - 1)] * hist[-L]
            hf_part = sum(beta[off_hf[nm] + (i_s - 1)] * frozen_hf[nm]
                          for nm in hf_used)
            lf_part = sum(beta[off_lf[v] + (i_s - 1)] * frozen_lf[v]
                          for v in lf_used)
            moy_part = 0.0
            m = tdate.month
            if m >= 2:
                moy_part = beta[cal_pos[f"mon{m:02d}"]]
            wd_part = 0.0
            wd_iso = tdate.dayofweek + 1
            if 1 <= wd_iso <= 6:
                wd_part = beta[cal_pos[cu.WD_NAMES[wd_iso - 1]]]
            hol_part = (beta[cal_pos["holiday"]]
                        if tdate in c["hol_set"] else 0.0)
            yh = ar_part + hf_part + lf_part + moy_part + wd_part + hol_part
            h_path[s - 1] = yh
            hist.append(yh)

        for h in cu.HORIZONS:
            if np.isfinite(h_path[h - 1]):
                yhat_mat[oi, h - 1] = np.sinh(h_path[h - 1])
    return si, yhat_mat


def run_country(country):
    global _CTX
    _CTX = _build_context(country)
    c = _CTX
    cc = c["cc"]
    spec_ids = list(c["masks"].index)
    n_cores = cu.get_n_cores()
    print(f"\n=== RU-MIDAS POS-DUMMIES - {country} (n_cores={n_cores}) ===")

    dict_df = cu.spec_dictionary(c["masks"], cu.HF_CANDIDATES, c["lf_names"])
    dp = os.path.join(cu.FORECAST_DIR, f"ru_midas_spec_dictionary_{cc}.csv")
    dict_df.to_csv(dp, index=False)
    origins = c["origins"]
    print(f"  rolling origins: {len(origins)}")

    t0c = time.time()
    results = [None] * len(spec_ids)
    if n_cores > 1 and len(spec_ids) > 1:
        with ProcessPoolExecutor(max_workers=n_cores,
                                 initializer=_init_worker,
                                 initargs=(country,)) as ex:
            for si, ym in ex.map(fit_one_spec, range(len(spec_ids))):
                results[si] = ym
    else:
        for si in range(len(spec_ids)):
            _, ym = fit_one_spec(si)
            results[si] = ym
    print(f"  fitting done in {time.time() - t0c:.1f}s")

    dates = c["dates"]
    y_lvl = c["daily"]["elec_price"].to_numpy(float)
    pos_full = within_month_position(dates)
    for h in cu.HORIZONS:
        out = pd.DataFrame({
            "origin_date": dates.iloc[origins].dt.strftime("%Y-%m-%d").to_numpy(),
            "target_date": dates.iloc[origins + h].dt.strftime("%Y-%m-%d").to_numpy(),
            "target_pos":  pos_full[origins + h],
            "y_actual":    y_lvl[origins + h],
        })
        for si, sid in enumerate(spec_ids):
            out[sid] = results[si][:, h - 1]     # already sinh-transformed
        out = out[np.isfinite(out["y_actual"].to_numpy(float))]
        fp = os.path.join(cu.FORECAST_DIR, f"ru_midas_{cc}_h{h:02d}.csv")
        out.to_csv(fp, index=False)
        print(f"  wrote {fp}  (rows = {len(out)})")


def main():
    country_arg = os.environ.get("COUNTRY_NAME", "ALL")
    countries = (["Germany", "Italy"] if country_arg.upper() == "ALL"
                 else [country_arg])
    for cn in countries:
        run_country(cn)
    print("\nRU-MIDAS POS-DUMMIES done.")


if __name__ == "__main__":
    main()
